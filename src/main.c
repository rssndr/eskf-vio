#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "euroc.h"
#include "quat.h"
#include "eskf.h"
#include "image.h"
#include "frontend.h"
#include "msckf.h"

/* Monotonic wall clock for the stage timing report. Diagnostic only: nothing
 * downstream reads it, so the estimator's numbers stay bit-identical. */
static double now_s(void) {
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static size_t gt_nearest(const gt_sample_t *gt, size_t n, double t) {
        for (size_t i = 0; i < n; i++)
                if (gt[i].timestamp >= t)
                        return i;
        return n - 1;
}

/* ---- observability diagnostic (passive: reads f, never writes it) -------
 *
 * A monocular VIO system cannot observe four directions of its own state: the
 * three global translations, and rotation about the gravity axis. Along those
 * directions the covariance must therefore NOT shrink. Prediction adds no
 * process noise there, and no measurement can carry information there, so a
 * correct filter simply cannot learn anything about them. If the covariance
 * does shrink, the filter is inventing information about a direction it cannot
 * see -- that is what "overconfident" means, and this turns it into a number.
 *
 * The directions are exact, not approximations. With the state laid out as in
 * eskf.c -- [0..2] pos, [3..5] vel, [6..8] attitude, [9..11] ba, [12..14] bg,
 * then 6 per clone (position, attitude) from index 15:
 *
 *   translation along world axis i:
 *       dp = e_i on the body position AND on every clone position; nothing
 *       else moves. F n = n exactly (the POS<-VEL coupling needs a velocity
 *       component, and there is none), and n'Qn = 0 exactly (Q is non-zero
 *       only on VEL, TH, BA, BG). So the variance along it is a CONSTANT, and
 *       any decrease is spurious information, full stop.
 *
 *   rotation about world z:
 *       dp = e_z x p, dv = e_z x v, dtheta = R' e_z on the body and every
 *       clone. F n = n to first order in dt. Here Q does inject noise, so the
 *       variance may grow -- but it must never DECREASE.
 *
 * R' e_z is the third row of the body-to-world rotation: the attitude error is
 * right-multiplied (q_true = q * exp(dtheta), see eskf_inject), so a world-frame
 * yaw corresponds to a body-frame axis R' e_z. The convention is checked once
 * against eskf.c's quat_to_R() at startup rather than assumed.
 *
 * n'Pn is read straight out of f.P.d[] because mat_get() copies all 51,200
 * bytes of mat_t by value -- a 75x75 loop through it would move a gigabyte per
 * frame and wreck the timing report. It runs outside the t_front and t_upd
 * spans; at ~22,500 mul-adds per frame against a 223 ms frame it is 0.003%.
 */

/* Third row of the body-to-world rotation, i.e. R' e_z. Same formula and
 * convention as quat_to_R() in eskf.c, repeated here because that function
 * returns a 51,200-byte mat_t and this runs once per clone per frame. */
static void R_row2(const quaternion_t *q, double out[3]) {
        out[0] = 2.0 * (q->x*q->z - q->y*q->w);
        out[1] = 2.0 * (q->y*q->z + q->x*q->w);
        out[2] = 1.0 - 2.0 * (q->x*q->x + q->y*q->y);
}

mat_t quat_to_R(quaternion_t q);   /* defined in eskf.c, not in eskf.h */

static size_t state_dim(const eskf_t *f) {
        return 15 + 6 * (size_t)f->n_clones;
}

static void dir_translation(const eskf_t *f, int axis, double *n) {
        size_t d = state_dim(f);
        for (size_t i = 0; i < d; i++) n[i] = 0.0;
        double s = 1.0 / sqrt((double)(1 + f->n_clones));
        n[(size_t)axis] = s;
        for (int c = 0; c < f->n_clones; c++)
                n[15 + 6*(size_t)c + (size_t)axis] = s;
}

static void dir_yaw(const eskf_t *f, double *n) {
        size_t d = state_dim(f);
        for (size_t i = 0; i < d; i++) n[i] = 0.0;

        double r[3];
        n[0] = -f->pos.y;  n[1] = f->pos.x;          /* e_z x p */
        n[3] = -f->vel.y;  n[4] = f->vel.x;          /* e_z x v */
        R_row2(&f->q, r);
        n[6] = r[0];  n[7] = r[1];  n[8] = r[2];     /* R' e_z */

        for (int c = 0; c < f->n_clones; c++) {
                size_t q = 15 + 6*(size_t)c;
                n[q+0] = -f->clones[c].pos.y;
                n[q+1] =  f->clones[c].pos.x;
                R_row2(&f->clones[c].q, r);
                n[q+3] = r[0];  n[q+4] = r[1];  n[q+5] = r[2];
        }

        double s = 0.0;
        for (size_t i = 0; i < d; i++) s += n[i]*n[i];
        s = sqrt(s);
        if (s > 0.0)
                for (size_t i = 0; i < d; i++) n[i] /= s;
}

static double dir_var(const eskf_t *f, const double *n) {
        size_t d = state_dim(f);
        double v = 0.0;
        for (size_t i = 0; i < d; i++) {
                if (n[i] == 0.0) continue;
                const double *row = f->P.d + i * MAT_MAX;
                double acc = 0.0;
                for (size_t j = 0; j < d; j++)
                        acc += row[j] * n[j];
                v += n[i] * acc;
        }
        return v;
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
        /* Optional 6th argument: measurement noise in pixels (default 3.0).
         * Expressed against the focal length to get normalized units, which is
         * what the filter wants. Exposed so the noise model can be swept
         * against the overconfidence ratio instead of being assumed. */
        double sigma_px = (argc >= 7) ? atof(argv[6]) : 3.0;
        int stop = 0;
        FILE *diag = fopen(diag_path, "w");
        if (diag)
                fprintf(diag, "t,frame,pos_err,sigma_1d,sigma_3d,ratio_1d,ratio_3d,"
                              "att_err_deg,n_live,n_dead,n_sub,n_trunc,nobs_sum,nobs_max,"
                              "nobs_all_max,trunc_min,trunc_max,n_clones,ok,rej,ok_frac,n_wtd,w_min,"
                              "r_k,r_tri,r_jac,r_null,r_chy,r_chi,iobs,vobs,nis_sum,nis_dof,"
                              "var_tx,var_ty,var_tz,var_yaw\n");
        else
                fprintf(stderr, "warning: cannot open %s — diagnostics disabled\n", diag_path);

        printf("%-8s %-12s %-12s %-10s\n", "t [s]", "pos err [m]", "pred +- [m]", "att err [deg]");

        /* Stage timing. Answers "what do I optimise" with a measurement instead
         * of arithmetic: how much of the run is the MSCKF update, and how much
         * is the front-end. RT frac is that stage's total divided by the flight
         * duration, so a value above 1.00 means the stage alone is slower than
         * real time. */
        double t_prop = 0.0, t_front = 0.0, t_upd = 0.0;
        size_t k_last = i0;   /* last IMU sample processed, for the true span */
        double t_wall0 = now_s();

        /* Observability diagnostic state. Index 0-2 = translation along world
         * x/y/z, index 3 = rotation about world z.
         *
         * The direction vector is an average over the body position and every
         * clone position, so its weights depend on the clone count. The raw value
         * therefore moves whenever the window changes size -- bookkeeping, not
         * information. Only two frames with the SAME clone count have identical
         * weights, so only those pairs are compared. Once the window is full
         * (MAX_CLONES, reached after the first 10 camera frames) the count stops
         * changing, so every later frame is comparable with its predecessor. */
        double ndir[MAT_MAX];
        double ob[4];
        double ob_ref[4];       /* value at the first full-window frame */
        double ob_min[4], ob_max[4];
        double ob_worst[4];     /* most negative single-frame change */
        double ob_negsum[4];    /* sum of every negative change */
        int    ob_ndrop[4];
        double ob_prev[4];
        int have_full = 0, have_prev = 0, prev_c = -1;
        for (int a = 0; a < 4; a++) {
                ob[a] = ob_ref[a] = ob_prev[a] = 0.0;
                ob_min[a] = 1e300;
                ob_max[a] = -1e300;
                ob_worst[a] = 0.0;
                ob_negsum[a] = 0.0;
                ob_ndrop[a] = 0;
        }
        {
                /* Check R_row2() against eskf.c's quat_to_R() once, so the yaw
                 * direction rests on a verified convention and not on a guess
                 * about how the attitude error is applied. */
                double r[3];
                mat_t R = quat_to_R(f.q);
                R_row2(&f.q, r);
                double dm = 0.0;
                for (int k = 0; k < 3; k++) {
                        double e = fabs(r[k] - mat_get(R, 2, (size_t)k));
                        if (e > dm) dm = e;
                }
                printf("convention check: |R_row2 - quat_to_R row 2| = %.3e\n\n", dm);
        }

        for (size_t k = i0; k < n-1 && !stop; k++) {
                k_last = k;
                double dt = imu[k+1].timestamp - imu[k].timestamp;
                double tp0 = now_s();
                eskf_predict(&f, imu[k], dt);
                t_prop += now_s() - tp0;

                if (ic < ncam && cam[ic].timestamp <= imu[k+1].timestamp) {
                        if (t_end > 0.0 && cam[ic].timestamp - imu[i0].timestamp > t_end) {
                                stop = 1;
                                break;
                        }
                        double tf0 = now_s();
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
                                double nis_before = msckf_nis_sum;
                                int dof_before = msckf_nis_dof;
                                int nobs_sum = 0, nobs_max = 0;
                                int all_max = 0;         /* longest track offered this frame */
                                int trunc_min = 0, trunc_max = 0;  /* length range of truncated tracks */

                                eskf_augment(&f, cam[ic].timestamp);
                                frontend_process(&fe, &img);
                                double tf1 = now_s();
                                t_front += tf1 - tf0;

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
                                        if (msckf_update_track(&f, ci, tk->obs + off, ks, sigma_px/458.0, &w)) {
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
                                t_upd += now_s() - tf1;

                                /* Unobservable-direction variances. Outside
                                 * both stage spans, reads f only. */
                                {
                                        for (int a = 0; a < 3; a++) {
                                                dir_translation(&f, a, ndir);
                                                ob[a] = dir_var(&f, ndir);
                                        }
                                        dir_yaw(&f, ndir);
                                        ob[3] = dir_var(&f, ndir);

                                        int c = f.n_clones;
                                        if (have_prev && c == prev_c)
                                                for (int a = 0; a < 4; a++) {
                                                        double d = ob[a] - ob_prev[a];
                                                        if (d < 0.0) {
                                                                ob_negsum[a] += d;
                                                                ob_ndrop[a]++;
                                                                if (d < ob_worst[a]) ob_worst[a] = d;
                                                        }
                                                }
                                        for (int a = 0; a < 4; a++) ob_prev[a] = ob[a];
                                        prev_c = c;
                                        have_prev = 1;
                                        if (c == MAX_CLONES) {
                                                if (!have_full)
                                                        for (int a = 0; a < 4; a++) ob_ref[a] = ob[a];
                                                for (int a = 0; a < 4; a++) {
                                                        if (ob[a] < ob_min[a]) ob_min[a] = ob[a];
                                                        if (ob[a] > ob_max[a]) ob_max[a] = ob[a];
                                                }
                                                have_full = 1;
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
                                                "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.4f,%d,%.4f,%d,%d,%d,%d,%d,%d,%d,%d,%.4f,%d,"
                                                "%.6e,%.6e,%.6e,%.6e\n",
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
                                                msckf_valid_obs - vobs_before,
                                                msckf_nis_sum - nis_before,
                                                msckf_nis_dof - dof_before,
                                                ob[0], ob[1], ob[2], ob[3]);
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
        /* gamma is chi-squared with m degrees of freedom when R is right, so a
         * mean of 1 means the noise model is calibrated and 4 means sigma^2 is
         * four times too small. */
        if (msckf_nis_dof)
                printf("NIS: mean gamma/dof = %.3f   mean gamma = %.1f   (%d dof, %d tracks)\n",
                       msckf_nis_sum / msckf_nis_dof,
                       msckf_nis_sum / (updates_ok ? updates_ok : 1),
                       msckf_nis_dof, updates_ok);

        /* The headline. A correct filter gains NO information along these
         * directions: prediction adds none (n'Qn = 0 for translation, exactly)
         * and no measurement can. So the value is either exactly constant or
         * growing through Q, and every fall is information the filter cannot
         * have. 'ref' is the first frame with a full window; only frames with the
         * same clone count are compared, which after those first MAX_CLONES
         * frames is all of them. */
        if (have_full) {
                static const char *nm[4] = { "translation x", "translation y",
                                             "translation z", "yaw about z" };
                printf("\nunobservable directions -- P along a direction the camera cannot see\n");
                printf("compared only between frames with the same clone count;\n"
                       "'worst' and 'total' are the falls, as %% of the first full-window value\n\n");
                printf("%-15s %12s %12s %12s %12s %12s %6s\n",
                       "direction", "ref", "min", "max", "worst fall", "total fall", "falls");
                for (int a = 0; a < 4; a++)
                        printf("%-15s %12.4e %12.4e %12.4e %11.4f%% %11.4f%% %6d\n",
                               nm[a], ob_ref[a], ob_min[a], ob_max[a],
                               (ob_ref[a] > 0.0) ? 100.0*ob_worst[a]  / ob_ref[a] : 0.0,
                               (ob_ref[a] > 0.0) ? 100.0*ob_negsum[a] / ob_ref[a] : 0.0,
                               ob_ndrop[a]);
        }

        double t_wall = now_s() - t_wall0;
        /* Span actually processed, not the span of the dataset: with the stop
         * option these differ, and dividing by the dataset length silently
         * halves every RT frac. */
        double t_flight = imu[k_last].timestamp - imu[i0].timestamp;
        double t_acc = t_prop + t_front + t_upd;
        printf("\nstage timing (one thread, this machine — not the P4)\n");
        printf("%-10s %11s %8s %9s\n", "stage", "total [s]", "% run", "RT frac");
        printf("%-10s %11.1f %7.1f%% %9.2f\n", "predict",   t_prop,  100*t_prop/t_wall,  t_prop/t_flight);
        printf("%-10s %11.1f %7.1f%% %9.2f\n", "front-end", t_front, 100*t_front/t_wall, t_front/t_flight);
        printf("%-10s %11.1f %7.1f%% %9.2f\n", "update",    t_upd,   100*t_upd/t_wall,   t_upd/t_flight);
        printf("%-10s %11.1f %7.1f%% %9.2f\n", "accounted", t_acc,   100*t_acc/t_wall,   t_acc/t_flight);
        printf("%-10s %11.1f %7.1f%% %9.2f\n", "wall",      t_wall,  100.0,              t_wall/t_flight);
        printf("flight %.1f s over %zu frames -> %.0f ms/frame, %.0f ms of it update\n",
               t_flight, ic, 1000.0*t_wall/(ic ? ic : 1), 1000.0*t_upd/(ic ? ic : 1));
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

