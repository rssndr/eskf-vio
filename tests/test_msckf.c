#include <stdio.h>
#include <math.h>
#include "msckf.h"
#include "ins.h"
#include "tri.h"

static int nfail = 0;
static void check(const char *name, int cond) {
        printf("%-44s %s\n", name, cond ? "ok" : "FAIL");
        if (!cond) nfail++;
}

int main(void) {
        const double s2 = 0.7071067811865476;
        quaternion_t qc = q_norm((quaternion_t){ s2, 0.1, -0.05, s2 });
        clone_t cl = { qc, { 0.4, -0.2, 0.15 }, 0.0 };
        vector_3d_t pf = { 1.5, 0.8, 4.0 };
        const double eps = 1e-7;

        obs_jac_t o = obs_jacobian(&cl, pf);
        check("jac: valid", o.valid);

        double maxerr = 0;
        for (int k = 0; k < 3; k++) {
                clone_t p = cl;
                if (k == 0) p.pos.x += eps;
                if (k == 1) p.pos.y += eps;
                if (k == 2) p.pos.z += eps;
                obs_jac_t op = obs_jacobian(&p, pf);
                double n0 = (op.z.x - o.z.x) / eps, n1 = (op.z.y - o.z.y) / eps;
                double e = fmax(fabs(n0 - o.Hp[k]), fabs(n1 - o.Hp[3+k]));
                if (e > maxerr) maxerr = e;
        }
        printf("  Hp  max err %.2e\n", maxerr);
        check("jac: Hp matches finite differences", maxerr < 1e-5);

        maxerr = 0;
        for (int k = 0; k < 3; k++) {
                vector_3d_t dth = { 0, 0, 0 };
                if (k == 0) dth.x = eps;
                if (k == 1) dth.y = eps;
                if (k == 2) dth.z = eps;
                clone_t p = cl;
                p.q = q_norm(q_mul_q(cl.q, gyro_to_q(dth, 1.0)));
                obs_jac_t op = obs_jacobian(&p, pf);
                double n0 = (op.z.x - o.z.x) / eps, n1 = (op.z.y - o.z.y) / eps;
                double e = fmax(fabs(n0 - o.Hth[k]), fabs(n1 - o.Hth[3+k]));
                if (e > maxerr) maxerr = e;
        }
        printf("  Hth max err %.2e\n", maxerr);
        check("jac: Hth matches finite differences", maxerr < 1e-5);

        maxerr = 0;
        for (int k = 0; k < 3; k++) {
                vector_3d_t p = pf;
                if (k == 0) p.x += eps;
                if (k == 1) p.y += eps;
                if (k == 2) p.z += eps;
                obs_jac_t op = obs_jacobian(&cl, p);
                double n0 = (op.z.x - o.z.x) / eps, n1 = (op.z.y - o.z.y) / eps;
                double e = fmax(fabs(n0 - o.Hf[k]), fabs(n1 - o.Hf[3+k]));
                if (e > maxerr) maxerr = e;
        }
        printf("  Hf  max err %.2e\n", maxerr);
        check("jac: Hf matches finite differences", maxerr < 1e-5);

        check("jac: Hp = -Hf", fabs(o.Hp[0] + o.Hf[0]) < 1e-15 &&
                               fabs(o.Hp[5] + o.Hf[5]) < 1e-15);

        quaternion_t id = {1,0,0,0};
        vector_3d_t zero = {0,0,0};
        eskf_t f;
        eskf_init(&f, id, zero, (vector_3d_t){0.6, 0.1, 0.05}, zero, zero);
        imu_sample_t hover = { 0, {0,0,0}, {0,0,9.81} };
        for (int c = 0; c < 6; c++) {
                for (int s = 0; s < 40; s++) eskf_predict(&f, hover, 0.005);
                eskf_augment(&f, c * 0.2);
        }

        vector_3d_t pft = { 0.3, 0.2, 5.0 };
        pt2_t obs[6];
        int ci[6] = { 0, 1, 2, 3, 4, 5 };
        int allv = 1;
        for (int i = 0; i < 6; i++) {
                obs_jac_t oi = obs_jacobian(&f.clones[i], pft);
                allv &= oi.valid;
                obs[i] = oi.z;
        }
        check("upd: all views valid", allv);

        {
                vector_3d_t tp;
                clone_t tmp[6];
                for (int i = 0; i < 6; i++) tmp[i] = f.clones[i];
                triangulate(tmp, obs, 6, &tp);
                mat_t Hf = mat_zero(12, 3);
                for (int i = 0; i < 6; i++) {
                        obs_jac_t oi = obs_jacobian(&f.clones[i], tp);
                        for (int rr = 0; rr < 2; rr++)
                                for (int cc = 0; cc < 3; cc++)
                                        mat_set(&Hf, 2*i+rr, cc, oi.Hf[rr*3+cc]);
                }
                mat_t A = msckf_nullspace(Hf);
                check("null: dims 9x12", A.rows == 9 && A.cols == 12);
                mat_t AH = mat_mul(A, Hf);
                double mx = 0;
                for (size_t i = 0; i < AH.rows * 3; i++)
                        if (fabs(AH.d[i]) > mx) mx = fabs(AH.d[i]);
                printf("  |A Hf| max %.2e\n", mx);
                check("null: A Hf = 0", mx < 1e-10);
                mat_t AAt = mat_mul(A, mat_transpose(A));
                int orth = 1;
                for (size_t i = 0; i < 9; i++)
                        for (size_t j = 0; j < 9; j++)
                                orth &= fabs(mat_get(AAt,i,j) - (i==j?1.0:0.0)) < 1e-10;
                check("null: rows orthonormal", orth);
        }

        double sigma = 0.5 / 458.0;

        double pv0 = mat_get(f.P, 15, 15);
        vector_3d_t before[6];
        for (int i = 0; i < 6; i++) before[i] = f.clones[i].pos;
        int ret = msckf_update_track(&f, ci, obs, 6, sigma, NULL);
        check("upd: perfect obs accepted", ret == 1);
        double mv = 0;
        for (int i = 0; i < 6; i++) {
                mv = fmax(mv, fabs(f.clones[i].pos.x - before[i].x));
                mv = fmax(mv, fabs(f.clones[i].pos.y - before[i].y));
                mv = fmax(mv, fabs(f.clones[i].pos.z - before[i].z));
        }
        printf("  perfect obs: max clone move %.2e m\n", mv);
        check("upd: near-zero correction", mv < 1e-6);
        check("upd: clone pos variance shrinks", mat_get(f.P, 15, 15) < pv0);
        int usym = 1;
        for (size_t i = 0; i < f.P.rows; i++)
                for (size_t j = 0; j < f.P.rows; j++)
                        usym &= fabs(mat_get(f.P,i,j) - mat_get(f.P,j,i)) < 1e-9;
        check("upd: P symmetric", usym);

        vector_3d_t truth[6];
        for (int i = 0; i < 6; i++) truth[i] = f.clones[i].pos;
        f.clones[2].pos.x += 0.004;  f.clones[2].pos.y -= 0.003;  f.clones[2].pos.z += 0.002;
        f.clones[4].pos.x -= 0.003;  f.clones[4].pos.y += 0.005;  f.clones[4].pos.z += 0.001;
        for (int c = 0; c < 6; c++)
                for (int j = 0; j < 3; j++) {
                        size_t d = 15 + 6*c + j;
                        mat_set(&f.P, d, d, mat_get(f.P, d, d) + 1e-4);
                }
        double e0 = 0;
        for (int i = 0; i < 6; i++)
                e0 += fabs(f.clones[i].pos.x - truth[i].x)
                    + fabs(f.clones[i].pos.y - truth[i].y)
                    + fabs(f.clones[i].pos.z - truth[i].z);
        ret = msckf_update_track(&f, ci, obs, 6, sigma, NULL);
        double e1 = 0;
        for (int i = 0; i < 6; i++)
                e1 += fabs(f.clones[i].pos.x - truth[i].x)
                    + fabs(f.clones[i].pos.y - truth[i].y)
                    + fabs(f.clones[i].pos.z - truth[i].z);
        printf("  corrupted clones: err %.4f -> %.4f m (ret %d)\n", e0, e1, ret);
        check("upd: corruption accepted", ret == 1);
        check("upd: error reduced > 40%", e1 < 0.6 * e0);

        double p00 = mat_get(f.P, 0, 0), p20 = mat_get(f.P, 20, 20);
        pt2_t bad[6];
        for (int i = 0; i < 6; i++) bad[i] = obs[i];
        bad[2].x += 0.1;

        /* One outlier in six observations. The Huber reweight went in with the
         * chi-squared gate replacement: the track is accepted and down-weighted
         * rather than thrown away. So the assertion is no longer "rejected and P
         * untouched" -- it is the stronger pair: the weight must actually bite,
         * and the outlier must move the covariance strictly less than the same
         * update would with clean observations. Two copies of the filter,
         * identical before the update, make that a direct comparison. */
        eskf_t fc = f;
        double w = 1.0;
        ret = msckf_update_track(&f, ci, bad, 6, sigma, &w);
        printf("  outlier: ret %d, Huber weight %.4f\n", ret, w);
        check("huber: outlier accepted, not rejected", ret == 1);
        check("huber: weight bites", w > 0.0 && w < 0.2);

        (void)msckf_update_track(&fc, ci, obs, 6, sigma, NULL);
        double drop_bad   = (p00 - mat_get(f.P, 0, 0))   + (p20 - mat_get(f.P, 20, 20));
        double drop_clean = (p00 - mat_get(fc.P, 0, 0))  + (p20 - mat_get(fc.P, 20, 20));
        printf("  P drop: outlier %.3e, clean obs %.3e\n", drop_bad, drop_clean);
        check("huber: P still moved by the outlier", drop_bad > 0.0);
        check("huber: outlier shrinks P less than clean obs", drop_bad < drop_clean);

        printf(nfail ? "FAIL (%d)\n" : "PASS\n", nfail);
        return nfail != 0;
}

