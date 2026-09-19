#include <math.h>
#include "eskf.h"
#include "ins.h"

enum { POS = 0, VEL = 3, TH = 6, BA = 9, BG = 12 };

mat_t mat_skew(vector_3d_t a) {
        mat_t M = mat_zero(3, 3);
        mat_set(&M, 0, 1, -a.z);
        mat_set(&M, 0, 2, a.y);
        mat_set(&M, 1, 0, a.z);
        mat_set(&M, 1, 2, -a.x);
        mat_set(&M, 2, 0, -a.y);
        mat_set(&M, 2, 1, a.x);
        return M;
}

mat_t quat_to_R(quaternion_t q) {
        mat_t R = mat_zero(3, 3);
        mat_set(&R, 0, 0, 1 - 2*(q.y*q.y + q.z*q.z));
        mat_set(&R, 0, 1, 2*(q.x*q.y - q.z*q.w));
        mat_set(&R, 0, 2, 2*(q.x*q.z + q.y*q.w));
        mat_set(&R, 1, 0, 2*(q.x*q.y + q.z*q.w));
        mat_set(&R, 1, 1, 1 - 2*(q.x*q.x + q.z*q.z));
        mat_set(&R, 1, 2, 2*(q.y*q.z - q.x*q.w));
        mat_set(&R, 2, 0, 2*(q.x*q.z - q.y*q.w));
        mat_set(&R, 2, 1, 2*(q.y*q.z + q.x*q.w));
        mat_set(&R, 2, 2, 1 - 2*(q.x*q.x + q.y*q.y));
        return R;
}

mat_t build_F(quaternion_t q, vector_3d_t a, vector_3d_t w, double dt) {
        mat_t F = mat_eye(15);
        mat_t R = quat_to_R(q);

        mat_set(&F, POS+0, VEL+0, dt);
        mat_set(&F, POS+1, VEL+1, dt);
        mat_set(&F, POS+2, VEL+2, dt);

        mat_t RA = mat_mul(R, mat_skew(a));
        for (size_t i = 0; i < 3; i++)
                for (size_t j = 0; j < 3; j++)
                        mat_set(&F, VEL+i, TH+j, -dt * mat_get(RA, i, j));

        for (size_t i = 0; i < 3; i++)
                for (size_t j = 0; j < 3; j++)
                        mat_set(&F, VEL+i, BA+j, -dt * mat_get(R, i, j));

        mat_t Rw = mat_transpose(quat_to_R(gyro_to_q(w, dt)));
        for (size_t i = 0; i < 3; i++)
                for (size_t j = 0; j < 3; j++)
                        mat_set(&F, TH+i, TH+j, mat_get(Rw, i, j));

        mat_set(&F, TH+0, BG+0, -dt);
        mat_set(&F, TH+1, BG+1, -dt);
        mat_set(&F, TH+2, BG+2, -dt);

        return F;
}

mat_t build_Q(double dt) {
        mat_t Q = mat_zero(15, 15);

        double qv  = SIG_A * SIG_A * dt;
        double qth = SIG_G * SIG_G * dt;
        double qba = SIG_BA * SIG_BA * dt;
        double qbg = SIG_BG * SIG_BG * dt;

        for (size_t i = 0; i < 3; i++) {
                mat_set(&Q, VEL+i, VEL+i, qv);
                mat_set(&Q, TH+i,  TH+i,  qth);
                mat_set(&Q, BA+i,  BA+i,  qba);
                mat_set(&Q, BG+i,  BG+i,  qbg);
        }
        return Q;
}

void eskf_init(eskf_t *f, quaternion_t q, vector_3d_t pos, vector_3d_t vel,
               vector_3d_t ba, vector_3d_t bg) {
        f->q = q;  f->pos = pos;  f->vel = vel;  f->ba = ba;  f->bg = bg;
        f->n_clones = 0;
        f->max_span = 0.0;
        f->P = mat_zero(15, 15);
        for (size_t i = POS; i < TH+3; i++)  mat_set(&f->P, i, i, 1e-5);
        for (size_t i = BA;  i < BA+3; i++)  mat_set(&f->P, i, i, 1e-6);
        for (size_t i = BG;  i < 15;   i++)  mat_set(&f->P, i, i, 1e-8);
}

void eskf_predict(eskf_t *f, imu_sample_t s, double dt) {
        vector_3d_t w = s.gyro;
        w.x -= f->bg.x;  w.y -= f->bg.y;  w.z -= f->bg.z;
        vector_3d_t a = s.accel;
        a.x -= f->ba.x;  a.y -= f->ba.y;  a.z -= f->ba.z;

        f->q = q_norm(q_mul_q(f->q, gyro_to_q(w, dt)));

        size_t n = 15 + 6 * (size_t)f->n_clones;
        mat_t F15 = build_F(f->q, a, w, dt);
        mat_t Q15 = build_Q(dt);
        mat_t F = mat_eye(n);
        mat_t Q = mat_zero(n, n);
        for (size_t i = 0; i < 15; i++)
                for (size_t j = 0; j < 15; j++) {
                        mat_set(&F, i, j, mat_get(F15, i, j));
                        mat_set(&Q, i, j, mat_get(Q15, i, j));
                }

        mat_t Ft = mat_transpose(F);
        f->P = mat_add(mat_mul(mat_mul(F, f->P), Ft), Q);

        rotate_vector(&a, f->q);
        a.z -= 9.81;
        f->vel.x += a.x*dt;  f->vel.y += a.y*dt;  f->vel.z += a.z*dt;
        f->pos.x += f->vel.x*dt;  f->pos.y += f->vel.y*dt;  f->pos.z += f->vel.z*dt;
}

void eskf_inject(eskf_t *f, const mat_t *dx) {
        f->pos.x += dx->d[POS+0];  f->pos.y += dx->d[POS+1];  f->pos.z += dx->d[POS+2];
        f->vel.x += dx->d[VEL+0];  f->vel.y += dx->d[VEL+1];  f->vel.z += dx->d[VEL+2];
        vector_3d_t dth = { dx->d[TH+0], dx->d[TH+1], dx->d[TH+2] };
        f->q = q_norm(q_mul_q(f->q, gyro_to_q(dth, 1.0)));
        f->ba.x += dx->d[BA+0];  f->ba.y += dx->d[BA+1];  f->ba.z += dx->d[BA+2];
        f->bg.x += dx->d[BG+0];  f->bg.y += dx->d[BG+1];  f->bg.z += dx->d[BG+2];

        for (int c = 0; c < f->n_clones; c++) {
                size_t d = 15 + 6 * (size_t)c;
                f->clones[c].pos.x += dx->d[d+0];
                f->clones[c].pos.y += dx->d[d+1];
                f->clones[c].pos.z += dx->d[d+2];
                vector_3d_t cth = { dx->d[d+3], dx->d[d+4], dx->d[d+5] };
                f->clones[c].q = q_norm(q_mul_q(f->clones[c].q, gyro_to_q(cth, 1.0)));
        }
}

/* Accelerometer as a gravity measurement: 9.81 * R' e_z, gated on | |a| - g |. */
int eskf_update_gravity(eskf_t *f, vector_3d_t a_m, double gate, double sigma_a) {
        double am = sqrt(a_m.x*a_m.x + a_m.y*a_m.y + a_m.z*a_m.z);
        if (gate <= 0.0 || fabs(am - 9.81) > gate)
                return 0;

        size_t n = 15 + 6 * (size_t)f->n_clones, st = f->P.cols;
        mat_t R = quat_to_R(f->q);
        double v[3] = { 9.81 * mat_get(R, 2, 0),
                        9.81 * mat_get(R, 2, 1),
                        9.81 * mat_get(R, 2, 2) };
        double sk[9] = {  0.0, -v[2],  v[1],
                        v[2],   0.0, -v[0],
                       -v[1],  v[0],   0.0 };

        static double H[3*MAT_MAX], PHt[MAT_MAX*3], tmp[MAT_MAX*MAT_MAX];
        for (size_t j = 0; j < n; j++) {
                H[0*n+j] = H[1*n+j] = H[2*n+j] = 0.0;
        }
        for (size_t i = 0; i < 3; i++) {
                for (size_t j = 0; j < 3; j++)
                        H[i*n + TH + j] = -sk[3*i+j];
                H[i*n + BA + i] = -1.0;
        }

        double y[3] = { a_m.x - f->ba.x - v[0],
                        a_m.y - f->ba.y - v[1],
                        a_m.z - f->ba.z - v[2] };

        for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < 3; j++) {
                        double s = 0.0;
                        for (size_t k = 0; k < n; k++)
                                s += f->P.d[i*st+k] * H[j*n+k];
                        PHt[i*3+j] = s;
                }

        mat_t S = mat_zero(3, 3);
        for (size_t i = 0; i < 3; i++)
                for (size_t j = 0; j < 3; j++) {
                        double s = 0.0;
                        for (size_t k = 0; k < n; k++)
                                s += H[i*n+k] * PHt[k*3+j];
                        S.d[i*3+j] = s + (i == j ? sigma_a*sigma_a : 0.0);
                }
        mat_t Si = mat3_inv(S);

        for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < n; j++) {
                        double s = 0.0;
                        for (size_t k = 0; k < 3; k++) {
                                double Kik = 0.0;
                                for (size_t m = 0; m < 3; m++)
                                        Kik += PHt[i*3+m] * Si.d[m*3+k];
                                s += Kik * PHt[j*3+k];
                        }
                        tmp[i*n+j] = f->P.d[i*st+j] - s;
                }
        for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < n; j++)
                        f->P.d[i*st+j] = tmp[i*n+j];

        mat_t dx = mat_zero(n, 1);
        for (size_t i = 0; i < n; i++) {
                double s = 0.0;
                for (size_t j = 0; j < 3; j++) {
                        double Kij = 0.0;
                        for (size_t m = 0; m < 3; m++)
                                Kij += PHt[i*3+m] * Si.d[m*3+j];
                        s += Kij * y[j];
                }
                dx.d[i] = s;
        }
        eskf_inject(f, &dx);
        return 1;
}

/* Velocity against a reference, H on the velocity block only; v_ref = 0 is a ZUPT. */
int eskf_update_vel(eskf_t *f, vector_3d_t v_ref, double sigma_v) {
        size_t n = 15 + 6 * (size_t)f->n_clones, st = f->P.cols;
        if (sigma_v <= 0.0)
                return 0;

        static double H[3*MAT_MAX], PHt[MAT_MAX*3], tmp[MAT_MAX*MAT_MAX];
        for (size_t j = 0; j < n; j++)
                H[0*n+j] = H[1*n+j] = H[2*n+j] = 0.0;
        for (size_t i = 0; i < 3; i++)
                H[i*n + VEL + i] = 1.0;

        double y[3] = { v_ref.x - f->vel.x, v_ref.y - f->vel.y, v_ref.z - f->vel.z };

        for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < 3; j++) {
                        double s = 0.0;
                        for (size_t k = 0; k < n; k++)
                                s += f->P.d[i*st+k] * H[j*n+k];
                        PHt[i*3+j] = s;
                }

        mat_t S = mat_zero(3, 3);
        for (size_t i = 0; i < 3; i++)
                for (size_t j = 0; j < 3; j++) {
                        double s = 0.0;
                        for (size_t k = 0; k < n; k++)
                                s += H[i*n+k] * PHt[k*3+j];
                        S.d[i*3+j] = s + (i == j ? sigma_v*sigma_v : 0.0);
                }
        mat_t Si = mat3_inv(S);

        for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < n; j++) {
                        double s = 0.0;
                        for (size_t k = 0; k < 3; k++) {
                                double Kik = 0.0;
                                for (size_t m = 0; m < 3; m++)
                                        Kik += PHt[i*3+m] * Si.d[m*3+k];
                                s += Kik * PHt[j*3+k];
                        }
                        tmp[i*n+j] = f->P.d[i*st+j] - s;
                }
        for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < n; j++)
                        f->P.d[i*st+j] = tmp[i*n+j];

        mat_t dx = mat_zero(n, 1);
        for (size_t i = 0; i < n; i++) {
                double s = 0.0;
                for (size_t j = 0; j < 3; j++) {
                        double Kij = 0.0;
                        for (size_t m = 0; m < 3; m++)
                                Kij += PHt[i*3+m] * Si.d[m*3+j];
                        s += Kij * y[j];
                }
                dx.d[i] = s;
        }
        eskf_inject(f, &dx);
        return 1;
}

int eskf_update_zupt(eskf_t *f, double sigma_v) {
        return eskf_update_vel(f, (vector_3d_t){ 0.0, 0.0, 0.0 }, sigma_v);
}

/* No parallax across the clone window leaves translation and metric scale unobserved.
   Add the acceleration uncertainty no measurement constrains, so P stops claiming a
   precision the camera cannot supply. Inert at or above base_ref. */
int eskf_inflate_noparallax(eskf_t *f, double base_ref, double sigma_a, double dt) {
        if (sigma_a <= 0.0 || base_ref <= 0.0 || f->n_clones < 2)
                return 0;

        double bx = f->clones[f->n_clones-1].pos.x - f->clones[0].pos.x;
        double by = f->clones[f->n_clones-1].pos.y - f->clones[0].pos.y;
        double bz = f->clones[f->n_clones-1].pos.z - f->clones[0].pos.z;
        double g = 1.0 - sqrt(bx*bx + by*by + bz*bz) / base_ref;
        if (g <= 0.0)
                return 0;

        double s = sigma_a * g, dv = s * dt, dp = 0.5 * s * dt * dt;
        size_t st = f->P.cols;
        for (size_t i = 0; i < 3; i++) {
                size_t p = POS + i, v = VEL + i;
                f->P.d[p*st+p] += dp*dp;
                f->P.d[p*st+v] += dp*dv;  f->P.d[v*st+p] += dp*dv;
                f->P.d[v*st+v] += dv*dv;
        }
        return 1;
}

void eskf_update_pos(eskf_t *f, vector_3d_t z, double sigma_z) {
        size_t n = 15 + 6 * (size_t)f->n_clones;

        mat_t y = mat_zero(3, 1);
        y.d[0] = z.x - f->pos.x;
        y.d[1] = z.y - f->pos.y;
        y.d[2] = z.z - f->pos.z;

        mat_t PHt = mat_zero(n, 3);
        mat_t S   = mat_zero(3, 3);
        for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < 3; j++)
                        mat_set(&PHt, i, j, mat_get(f->P, i, j));
        for (size_t i = 0; i < 3; i++)
                for (size_t j = 0; j < 3; j++)
                        mat_set(&S, i, j, mat_get(f->P, i, j) + (i == j ? sigma_z*sigma_z : 0.0));

        mat_t K  = mat_mul(PHt, mat3_inv(S));
        mat_t dx = mat_mul(K, y);

        mat_t KH = mat_zero(n, n);
        for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < 3; j++)
                        mat_set(&KH, i, j, mat_get(K, i, j));
        f->P = mat_mul(mat_add(mat_eye(n), mat_scale(KH, -1.0)), f->P);

        eskf_inject(f, &dx);
}

static void marginalize_oldest(eskf_t *f) {
        size_t n = 15 + 6 * (size_t)f->n_clones;
        mat_t Pn = mat_zero(n - 6, n - 6);
        size_t ii = 0;
        for (size_t i = 0; i < n; i++) {
                if (i >= 15 && i < 21) continue;
                size_t jj = 0;
                for (size_t j = 0; j < n; j++) {
                        if (j >= 15 && j < 21) continue;
                        mat_set(&Pn, ii, jj, mat_get(f->P, i, j));
                        jj++;
                }
                ii++;
        }
        f->P = Pn;

        for (int c = 0; c + 1 < f->n_clones; c++)
                f->clones[c] = f->clones[c + 1];
        f->n_clones--;
}

void eskf_augment(eskf_t *f, double timestamp) {
        /* Optional duration cap: without it the window spans more time at lower rates. */
        while (f->n_clones > 0 &&
               (f->n_clones >= MAX_CLONES ||
                (f->max_span > 0.0 && f->n_clones > 2 &&
                 timestamp - f->clones[0].timestamp > f->max_span)))
                marginalize_oldest(f);

        size_t d = 15 + 6 * (size_t)f->n_clones;
        size_t src[6] = { POS, POS+1, POS+2, TH, TH+1, TH+2 };

        mat_t Pn = mat_zero(d + 6, d + 6);
        for (size_t i = 0; i < d; i++)
                for (size_t j = 0; j < d; j++)
                        mat_set(&Pn, i, j, mat_get(f->P, i, j));
        for (size_t i = 0; i < 6; i++)
                for (size_t k = 0; k < d; k++) {
                        mat_set(&Pn, d+i, k, mat_get(f->P, src[i], k));
                        mat_set(&Pn, k, d+i, mat_get(f->P, k, src[i]));
                }
        for (size_t i = 0; i < 6; i++)
                for (size_t j = 0; j < 6; j++)
                        mat_set(&Pn, d+i, d+j, mat_get(f->P, src[i], src[j]));
        f->P = Pn;

        f->clones[f->n_clones] = (clone_t){ f->q, f->pos, timestamp };
        f->n_clones++;
}

