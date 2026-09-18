#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "euroc.h"
#include "quat.h"
#include "eskf.h"
#include "image.h"
#include "frontend.h"
#include "msckf.h"

static size_t gt_nearest(const gt_sample_t *gt, size_t n, double t) {
        for (size_t i = 0; i < n; i++)
                if (gt[i].timestamp >= t)
                        return i;
        return n - 1;
}

int main(int argc, char *argv[]) {
        if (argc < 4) {
                fprintf(stderr, "usage: %s <imu.csv> <gt.csv> <cam0 dir>\n", argv[0]);
                return 1;
        }

        imu_sample_t *imu;
        size_t n = euroc_load_imu(argv[1], &imu);
        if (n == 0) return 1;

        gt_sample_t *gt;
        size_t m = euroc_load_gt(argv[2], &gt);
        if (m == 0) return 1;

        const char *cam_dir = argv[3];
        char path[512];
        snprintf(path, sizeof path, "%s/data.csv", cam_dir);
        cam_frame_t *cam;
        size_t ncam = euroc_load_cam(path, &cam);
        if (ncam == 0) return 1;

        size_t i0 = 0;
        while (i0 < n && imu[i0].timestamp < gt[0].timestamp)
                i0++;

        size_t j0 = gt_nearest(gt, m, imu[i0].timestamp);

        eskf_t f;
        eskf_init(&f, gt[j0].q, gt[j0].pos, gt[j0].vel, gt[0].accel_bias, gt[0].gyro_bias);

        frontend_t fe;
        frontend_init(&fe);
        size_t ic = 0;
        while (ic < ncam && cam[ic].timestamp < imu[i0].timestamp)
                ic++;

        int updates_ok = 0, updates_rej = 0;

        /* ---- diagnostics (passive: reads state, never changes it) ----
         * Optional 4th argument selects the CSV path. Default "diag.csv".
         * Every value written below is already computed by the run; nothing
         * here feeds back into the estimator.
         */
        const char *diag_path = (argc >= 5) ? argv[4] : "diag.csv";
        /* Optional 5th argument: stop after this many seconds of sequence.
         * Diagnostic-only. Window 1 lives at t = 20-45 s, so a 50 s run answers
         * most questions in a quarter of the time. The estimator path is
         * untouched; the same frames produce the same numbers. */
        double t_end = (argc >= 6) ? atof(argv[5]) : 0.0;
        int stop = 0;
        FILE *diag = fopen(diag_path, "w");
        if (diag)
                fprintf(diag, "t,frame,pos_err,sigma_1d,sigma_3d,ratio_1d,ratio_3d,"
                              "att_err_deg,n_live,n_dead,n_sub,n_trunc,nobs_sum,nobs_max,"
                              "nobs_all_max,trunc_min,trunc_max,n_clones,ok,rej,ok_frac,n_wtd,w_min,"
                              "r_k,r_tri,r_jac,r_null,r_chy,r_chi,iobs,vobs\n");
        else
                fprintf(stderr, "warning: cannot open %s — diagnostics disabled\n", diag_path);

        printf("%-8s %-12s %-12s %-10s\n", "t [s]", "pos err [m]", "pred +- [m]", "att err [deg]");

        for (size_t k = i0; k < n-1 && !stop; k++) {
                double dt = imu[k+1].timestamp - imu[k].timestamp;
                eskf_predict(&f, imu[k], dt);

                if (ic < ncam && cam[ic].timestamp <= imu[k+1].timestamp) {
                        if (t_end > 0.0 && cam[ic].timestamp - imu[i0].timestamp > t_end) {
                                stop = 1;
                                break;
                        }
                        snprintf(path, sizeof path, "%s/data/%s", cam_dir, cam[ic].filename);
                        image_t img;
                        if (image_load(path, &img) == 0) {
                                int f_ok = 0, f_rej = 0, f_sub = 0, f_trunc = 0;
                                int n_wtd = 0;      /* tracks down-weighted this frame */
                                double w_min = 1.0; /* worst weight applied this frame */
                                int rej_before[MSCKF_REJ_COUNT];
                                for (int r = 0; r < MSCKF_REJ_COUNT; r++)
                                        rej_before[r] = msckf_rej_count[r];
                                int iobs_before = msckf_invalid_obs;
                                int vobs_before = msckf_valid_obs;
                                int nobs_sum = 0, nobs_max = 0;
                                int all_max = 0;         /* longest track offered this frame */
                                int trunc_min = 0, trunc_max = 0;  /* length range of truncated tracks */

                                eskf_augment(&f, cam[ic].timestamp);
                                frontend_process(&fe, &img);

                                for (int d = 0; d < fe.n_dead; d++) {
                                        dead_track_t *tk = &fe.dead[d];
                                        int kk = tk->nobs;
                                        if (kk > all_max) all_max = kk;

                                        /* Truncate, never discard, and pair by frame.
                                         *
                                         * A dead track's newest observation is ALWAYS from
                                         * frame F-1, because harvest() copies hist[0..n-1]
                                         * before the current frame's observation is
                                         * appended. The newest clone, index n_clones-1, is
                                         * frame F. So the newest observation pairs with
                                         * clone n_clones-2, not n_clones-1:
                                         *
                                         *   clone index for frame f = n_clones - 1 - (F - f)
                                         *   obs[j] is frame F - ks + j  ->  ci[j] = n_clones - 1 - ks + j
                                         *
                                         * The window is therefore one effective clone
                                         * smaller than the array: kmax = n_clones - 1.
                                         * The oldest `off` observations are dropped.
                                         */
                                        int kmax = f.n_clones - 1;
                                        int ks   = (kk > kmax) ? kmax : kk;
                                        int off  = kk - ks;
                                        if (off > 0) {
                                                if (f_trunc == 0 || kk < trunc_min) trunc_min = kk;
                                                if (kk > trunc_max) trunc_max = kk;
                                                f_trunc++;
                                        }
                                        int ci[FE_HIST];
                                        for (int j = 0; j < ks; j++)
                                                ci[j] = f.n_clones - 1 - ks + j;
                                        f_sub++;
                                        nobs_sum += ks;
                                        if (ks > nobs_max) nobs_max = ks;
                                        double w = 1.0;
                                        if (msckf_update_track(&f, ci, tk->obs + off, ks, 3.0/458.0, &w)) {
                                                updates_ok++;
                                                f_ok++;
                                        } else {
                                                updates_rej++;
                                                f_rej++;
                                        }
                                        if (w < 1.0) {
                                                n_wtd++;
                                                if (w < w_min) w_min = w;
                                        }
                                }

                                if (diag) {
                                        gt_sample_t *gd = &gt[gt_nearest(gt, m, cam[ic].timestamp)];
                                        double dx = f.pos.x - gd->pos.x;
                                        double dy = f.pos.y - gd->pos.y;
                                        double dz = f.pos.z - gd->pos.z;
                                        double err = sqrt(dx*dx + dy*dy + dz*dz);
                                        double s1 = sqrt(mat_get(f.P, 0, 0));
                                        double s3 = sqrt(mat_get(f.P, 0, 0) +
                                                         mat_get(f.P, 1, 1) +
                                                         mat_get(f.P, 2, 2));
                                        double dot = f.q.w*gd->q.w + f.q.x*gd->q.x +
                                                     f.q.y*gd->q.y + f.q.z*gd->q.z;
                                        int nsub = f_ok + f_rej;
                                        fprintf(diag,
                                                "%.6f,%zu,%.4f,%.6f,%.6f,%.4f,%.4f,%.4f,"
                                                "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.4f,%d,%.4f,%d,%d,%d,%d,%d,%d,%d,%d\n",
                                                cam[ic].timestamp - imu[i0].timestamp, ic,
                                                err, s1, s3,
                                                (s1 > 0.0) ? err / s1 : 0.0,
                                                (s3 > 0.0) ? err / s3 : 0.0,
                                                2.0 * acos(fabs(dot)) * 180.0 / M_PI,
                                                fe.n, fe.n_dead, f_sub, f_trunc,
                                                nobs_sum, nobs_max,
                                                all_max, trunc_min, trunc_max, f.n_clones,
                                                f_ok, f_rej,
                                                (nsub > 0) ? (double)f_ok / (double)nsub : 0.0,
                                                n_wtd, w_min,
                                                msckf_rej_count[0] - rej_before[0],
                                                msckf_rej_count[1] - rej_before[1],
                                                msckf_rej_count[2] - rej_before[2],
                                                msckf_rej_count[3] - rej_before[3],
                                                msckf_rej_count[4] - rej_before[4],
                                                msckf_rej_count[5] - rej_before[5],
                                                msckf_invalid_obs - iobs_before,
                                                msckf_valid_obs - vobs_before);
                                }
                                image_free(&img);
                        }
                        ic++;
                }

                if ((k - i0) % 4000 == 0) {
                        gt_sample_t *g = &gt[gt_nearest(gt, m, imu[k].timestamp)];
                        double dx = f.pos.x - g->pos.x, dy = f.pos.y - g->pos.y, dz = f.pos.z - g->pos.z;
                        double dot = f.q.w*g->q.w + f.q.x*g->q.x + f.q.y*g->q.y + f.q.z*g->q.z;
                        printf("%-8.1f %-12.2f %-12.2f %-10.2f\n",
                               imu[k].timestamp - imu[i0].timestamp,
                               sqrt(dx*dx + dy*dy + dz*dz),
                               sqrt(mat_get(f.P, 0, 0)),
                               2.0 * acos(fabs(dot)) * 180.0 / M_PI);
                }
        }

        gt_sample_t *g = &gt[gt_nearest(gt, m, imu[n-1].timestamp)];
        double dx = f.pos.x - g->pos.x, dy = f.pos.y - g->pos.y, dz = f.pos.z - g->pos.z;
        printf("final: measured %.1f m, predicted +-%.1f m\n",
               sqrt(dx*dx + dy*dy + dz*dz),
               sqrt(mat_get(f.P, 0, 0)));
        printf("att 1-sigma: %.3f deg\n", sqrt(mat_get(f.P, 6, 6)) * 180.0 / M_PI);
        printf("updates: %d ok, %d rejected\n", updates_ok, updates_rej);
        printf("gyro bias est %.5f %.5f %.5f | true %.5f %.5f %.5f\n",
               f.bg.x, f.bg.y, f.bg.z,
               gt[m-1].gyro_bias.x, gt[m-1].gyro_bias.y, gt[m-1].gyro_bias.z);

        if (diag) fclose(diag);
        frontend_free(&fe);
        free(imu);
        free(gt);
        free(cam);
        return 0;
}

