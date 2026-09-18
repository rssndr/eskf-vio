#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "msckf.h"
#include "cam.h"
#include "tri.h"

/* Upper-tail chi-squared quantiles: entry [d-1] is exceeded with probability
 * 0.01 by a chi-squared variable with d degrees of freedom. Named CHI2_95
 * until 2026-09-18 — the name was wrong, the values were always the 99%
 * quantiles, so the designed false-rejection rate was 1%, not 5%. Measured on
 * MH_01: 0.4-2.1% of healthy tracks were rejected, which matches 1% and
 * confirms the values.
 *
 * Now used as the Huber knee rather than as a hard rejection threshold. */
static const double CHI2_99[17] = {
        6.635, 9.210, 11.345, 13.277, 15.086, 16.812, 18.475, 20.090, 21.666,
        23.209, 24.725, 26.217, 27.688, 29.141, 30.578, 32.000, 33.409,
};

int msckf_rej_count[MSCKF_REJ_COUNT];
int msckf_invalid_obs;
int msckf_valid_obs;
double msckf_nis_sum;
int    msckf_nis_dof;

extern mat_t quat_to_R(quaternion_t q);
extern mat_t mat_skew(vector_3d_t v);

obs_jac_t obs_jacobian(const clone_t *cl, vector_3d_t pf) {
        obs_jac_t o = { .valid = 0 };

        mat_t Rwb = quat_to_R(cl->q);
        mat_t Rbs = mat_zero(3, 3);
        for (int k = 0; k < 9; k++) Rbs.d[k] = EUROC_T_BS.R[k];
        mat_t Rwc = mat_mul(Rwb, Rbs);
        mat_t Rcw = mat_transpose(Rwc);

        vector_3d_t t = EUROC_T_BS.t;
        vector_3d_t c = {
                cl->pos.x + Rwb.d[0]*t.x + Rwb.d[1]*t.y + Rwb.d[2]*t.z,
                cl->pos.y + Rwb.d[3]*t.x + Rwb.d[4]*t.y + Rwb.d[5]*t.z,
                cl->pos.z + Rwb.d[6]*t.x + Rwb.d[7]*t.y + Rwb.d[8]*t.z,
        };

        vector_3d_t d = { pf.x - c.x, pf.y - c.y, pf.z - c.z };
        double x = Rcw.d[0]*d.x + Rcw.d[1]*d.y + Rcw.d[2]*d.z;
        double y = Rcw.d[3]*d.x + Rcw.d[4]*d.y + Rcw.d[5]*d.z;
        double z = Rcw.d[6]*d.x + Rcw.d[7]*d.y + Rcw.d[8]*d.z;
        o.depth = z;
        if (z <= 0.1) return o;

        o.z.x = x / z;
        o.z.y = y / z;

        // Jpi = (1/z) [1 0 -u; 0 1 -v]  (2x3)
        double Jpi[6] = { 1/z, 0, -o.z.x/z,  0, 1/z, -o.z.y/z };

        /* Hf = Jpi * Rcw ; Hp = -Hf */
        for (int r = 0; r < 2; r++)
                for (int cc = 0; cc < 3; cc++) {
                        double s = 0;
                        for (int k = 0; k < 3; k++)
                                s += Jpi[r*3 + k] * Rcw.d[k*3 + cc];
                        o.Hf[r*3 + cc] = s;
                        o.Hp[r*3 + cc] = -s;
                }

        // Hth = Jpi * Rbs^T * [Rwb^T (pf - p)]x
        vector_3d_t pb = { pf.x - cl->pos.x, pf.y - cl->pos.y, pf.z - cl->pos.z };
        vector_3d_t lb = {
                Rwb.d[0]*pb.x + Rwb.d[3]*pb.y + Rwb.d[6]*pb.z,
                Rwb.d[1]*pb.x + Rwb.d[4]*pb.y + Rwb.d[7]*pb.z,
                Rwb.d[2]*pb.x + Rwb.d[5]*pb.y + Rwb.d[8]*pb.z,
        };
        mat_t A = mat_mul(mat_transpose(Rbs), mat_skew(lb));
        for (int r = 0; r < 2; r++)
                for (int cc = 0; cc < 3; cc++) {
                        double s = 0;
                        for (int k = 0; k < 3; k++)
                                s += Jpi[r*3 + k] * A.d[k*3 + cc];
                        o.Hth[r*3 + cc] = s;
                }

        o.valid = 1;
        return o;
}

mat_t msckf_nullspace(mat_t Hf) {
        int rows = (int)Hf.rows;
        double q[3][40];
        int nq = 0;

        for (int c = 0; c < 3; c++) {
                double v[40];
                for (int i = 0; i < rows; i++) v[i] = mat_get(Hf, i, c);
                for (int p = 0; p < nq; p++) {
                        double d = 0;
                        for (int i = 0; i < rows; i++) d += v[i] * q[p][i];
                        for (int i = 0; i < rows; i++) v[i] -= d * q[p][i];
                }
                double nn = 0;
                for (int i = 0; i < rows; i++) nn += v[i] * v[i];
                if (nn > 1e-12) {
                        nn = sqrt(nn);
                        for (int i = 0; i < rows; i++) q[nq][i] = v[i] / nn;
                        nq++;
                }
        }

        mat_t A = mat_zero(rows - nq, rows);
        double acc[40][40];
        int m = 0;
        for (int e = 0; e < rows && m < rows - nq; e++) {
                double v[40] = {0};
                v[e] = 1.0;
                for (int p = 0; p < nq; p++) {
                        double d = q[p][e];
                        for (int i = 0; i < rows; i++) v[i] -= d * q[p][i];
                }
                for (int p = 0; p < m; p++) {
                        double d = 0;
                        for (int i = 0; i < rows; i++) d += v[i] * acc[p][i];
                        for (int i = 0; i < rows; i++) v[i] -= d * acc[p][i];
                }
                double nn = 0;
                for (int i = 0; i < rows; i++) nn += v[i] * v[i];
                if (nn > 1e-8) {
                        nn = sqrt(nn);
                        for (int i = 0; i < rows; i++) acc[m][i] = v[i] / nn;
                        m++;
                }
        }
        for (int r = 0; r < m; r++)
                for (int i = 0; i < rows; i++)
                        mat_set(&A, r, i, acc[r][i]);
        return A;
}

int msckf_update_track(eskf_t *f, const int *ci, const pt2_t *obs, int k, double sigma, double *w_out) {
        /* Weight actually applied to this track: 1.0 when the residual is
         * inside the knee. Reported as 1.0 on every early return so the caller
         * never reads a stale value. */
        if (w_out) *w_out = 1.0;
        if (k < 2 || 2*k > 40) { msckf_rej_count[MSCKF_REJ_K_RANGE]++; return 0; }
        size_t n = 15 + 6 * (size_t)f->n_clones;

        clone_t cl[MAX_CLONES];
        for (int i = 0; i < k; i++) cl[i] = f->clones[ci[i]];

        vector_3d_t pf;
        if (!triangulate(cl, obs, k, &pf)) { msckf_rej_count[MSCKF_REJ_TRIANG]++; return 0; }

        mat_t r  = mat_zero(2*k, 1);
        mat_t Hx = mat_zero(2*k, n);
        mat_t Hf = mat_zero(2*k, 3);
        int n_invalid = 0;
        for (int i = 0; i < k; i++) {
                obs_jac_t o = obs_jacobian(&f->clones[ci[i]], pf);
                if (!o.valid) {
                        /* Diagnostic only: keep scanning so the counters can
                         * say how many observations of this track were
                         * unusable. The track is still discarded exactly as
                         * before — row i is simply left at zero. */
                        n_invalid++;
                        continue;
                }
                r.d[2*i]   = obs[i].x - o.z.x;
                r.d[2*i+1] = obs[i].y - o.z.y;
                size_t d = 15 + 6 * (size_t)ci[i];
                for (int rr = 0; rr < 2; rr++)
                        for (int cc = 0; cc < 3; cc++) {
                                mat_set(&Hx, 2*i+rr, d+cc,   o.Hp[rr*3+cc]);
                                mat_set(&Hx, 2*i+rr, d+3+cc, o.Hth[rr*3+cc]);
                                mat_set(&Hf, 2*i+rr, cc,     o.Hf[rr*3+cc]);
                        }
        }
        if (n_invalid) {
                msckf_rej_count[MSCKF_REJ_JACOBIAN]++;
                msckf_invalid_obs += n_invalid;
                msckf_valid_obs   += k - n_invalid;
                /* MSCKF_DUMP=<n> prints the observations and clone positions of
                 * the first n tracks discarded this way, then stops. Diagnostic
                 * only; unset means no cost beyond one getenv. */
                static int dump_left = -1;
                if (dump_left < 0) {
                        const char *e = getenv("MSCKF_DUMP");
                        dump_left = e ? atoi(e) : 0;
                }
                if (dump_left > 0) {
                        dump_left--;
                        fprintf(stderr, "REJ k=%d pf=(%.5f %.5f %.5f) invalid=%d\n",
                                k, pf.x, pf.y, pf.z, n_invalid);
                        for (int i = 0; i < k; i++) {
                                obs_jac_t od = obs_jacobian(&f->clones[ci[i]], pf);
                                fprintf(stderr, "   obs[%d]=(%.6f %.6f) z=%+.6f clone%d=(%.5f %.5f %.5f)\n",
                                        i, obs[i].x, obs[i].y, od.depth, ci[i],
                                        f->clones[ci[i]].pos.x, f->clones[ci[i]].pos.y,
                                        f->clones[ci[i]].pos.z);
                        }
                }
                return 0;
        }

        mat_t A = msckf_nullspace(Hf);
        size_t m = A.rows;
        if (m == 0 || m > 17) { msckf_rej_count[MSCKF_REJ_NULLSPACE]++; return 0; }
        mat_t rp = mat_mul(A, r);
        mat_t Hp = mat_mul(A, Hx);

        mat_t S = mat_mul(mat_mul(Hp, f->P), mat_transpose(Hp));
        double s2 = sigma * sigma;
        for (size_t i = 0; i < m; i++)
                mat_set(&S, i, i, mat_get(S, i, i) + s2);

        mat_t y;
        if (!mat_chol_solve(S, rp, &y)) { msckf_rej_count[MSCKF_REJ_CHOL_Y]++; return 0; }
        double gamma = 0;
        for (size_t i = 0; i < m; i++) gamma += rp.d[i] * y.d[i];
        /* Raw innovation statistic, before any weighting: mean gamma/m should be
         * 1 when R is right, and is the estimator for the sigma correction. */
        msckf_nis_sum += gamma;
        msckf_nis_dof += (int)m;

        /* Robust (Huber) weighting replaces the hard chi-squared gate.
         *
         * A binary gate assumes the model is right. When it is not — the pose
         * has drifted, so the reprojection residual grows — the gate rejects
         * the very measurements that would correct the drift, and the error
         * grows further. Measured on MH_01: 85-90% of submitted tracks were
         * rejected for 25 s while the track population stayed identical to the
         * healthy phases.
         *
         * So do not reject on residual size; down-weight instead, scaling the
         * measurement noise up by 1/w with
         *
         *     w = sqrt(knee / gamma)   for gamma > knee,   else 1
         *
         * the Huber influence function generalised to the Mahalanobis distance
         * sqrt(gamma). The knee is the quantile the old gate used, so every
         * track that passed before is bit-identical.
         *
         * Influence is bounded: as gamma grows, w -> 0, the s2/w diagonal
         * dominates S, S^-1 -> 0 and the correction dx -> 0. No second hard
         * threshold is needed.
         *
         * Only the diagonal changes. Entry i of S is s0_i + s2, so
         * s0_i = S_ii - s2 and the weighted entry is s0_i + s2/w. */
        double w = 1.0;
        double knee = CHI2_99[m-1];
        if (gamma > knee) {
                w = sqrt(knee / gamma);
                for (size_t i = 0; i < m; i++)
                        mat_set(&S, i, i, mat_get(S, i, i) - s2 + s2 / w);
        }
        if (w_out) *w_out = w;

        mat_t Sinv;
        if (!mat_chol_solve(S, mat_eye(m), &Sinv)) { msckf_rej_count[MSCKF_REJ_CHOL_INV]++; return 0; }
        mat_t PHt = mat_mul(f->P, mat_transpose(Hp));
        mat_t K   = mat_mul(PHt, Sinv);
        mat_t dx  = mat_mul(K, rp);
        mat_t KH  = mat_mul(K, Hp);
        f->P = mat_mul(mat_add(mat_eye(n), mat_scale(KH, -1.0)), f->P);
        eskf_inject(f, &dx);
        return 1;
}

