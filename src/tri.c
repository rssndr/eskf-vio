#include <math.h>
#include "tri.h"
#include "cam.h"

extern mat_t quat_to_R(quaternion_t q);

/* Views where x is behind the camera, and the mean camera centre. */
static int count_bad_depths(const clone_t *cl, int n, vector_3d_t x, vector_3d_t *cmean) {
        double mx = 0, my = 0, mz = 0;
        int bad = 0;

        for (int i = 0; i < n; i++) {
                mat_t Rwb = quat_to_R(cl[i].q);
                mat_t Rbs = mat_zero(3, 3);
                for (int k = 0; k < 9; k++) Rbs.d[k] = EUROC_T_BS.R[k];
                mat_t Rwc = mat_mul(Rwb, Rbs);

                vector_3d_t t = EUROC_T_BS.t;
                vector_3d_t c = {
                        cl[i].pos.x + Rwb.d[0]*t.x + Rwb.d[1]*t.y + Rwb.d[2]*t.z,
                        cl[i].pos.y + Rwb.d[3]*t.x + Rwb.d[4]*t.y + Rwb.d[5]*t.z,
                        cl[i].pos.z + Rwb.d[6]*t.x + Rwb.d[7]*t.y + Rwb.d[8]*t.z,
                };
                mx += c.x; my += c.y; mz += c.z;

                double dx = x.x - c.x, dy = x.y - c.y, dz = x.z - c.z;
                double z = Rwc.d[2]*dx + Rwc.d[5]*dy + Rwc.d[8]*dz;
                if (z <= 0.1) bad++;
        }

        cmean->x = mx / n; cmean->y = my / n; cmean->z = mz / n;
        return bad;
}

int triangulate(const clone_t *cl, const pt2_t *obs, int n, vector_3d_t *out) {
        mat_t AtA = mat_zero(3, 3);
        mat_t Atb = mat_zero(3, 1);

        for (int i = 0; i < n; i++) {
                mat_t Rwb = quat_to_R(cl[i].q);

                // R_WC = R_WB * R_BS
                mat_t Rbs = mat_zero(3, 3);
                for (int k = 0; k < 9; k++) Rbs.d[k] = EUROC_T_BS.R[k];
                mat_t Rwc = mat_mul(Rwb, Rbs);

                // camera centre c = p + R_WB * t_BS
                vector_3d_t t = EUROC_T_BS.t;
                vector_3d_t c = {
                        cl[i].pos.x + Rwb.d[0]*t.x + Rwb.d[1]*t.y + Rwb.d[2]*t.z,
                        cl[i].pos.y + Rwb.d[3]*t.x + Rwb.d[4]*t.y + Rwb.d[5]*t.z,
                        cl[i].pos.z + Rwb.d[6]*t.x + Rwb.d[7]*t.y + Rwb.d[8]*t.z,
                };

                // M = R_WC^T; rows of M are columns of R_WC
                double m1[3] = { Rwc.d[0], Rwc.d[3], Rwc.d[6] };
                double m2[3] = { Rwc.d[1], Rwc.d[4], Rwc.d[7] };
                double m3[3] = { Rwc.d[2], Rwc.d[5], Rwc.d[8] };

                double rows[2][3], rhs[2];
                for (int k = 0; k < 3; k++) {
                        rows[0][k] = obs[i].x * m3[k] - m1[k];
                        rows[1][k] = obs[i].y * m3[k] - m2[k];
                }
                rhs[0] = rows[0][0]*c.x + rows[0][1]*c.y + rows[0][2]*c.z;
                rhs[1] = rows[1][0]*c.x + rows[1][1]*c.y + rows[1][2]*c.z;

                for (int r = 0; r < 2; r++)
                        for (int a = 0; a < 3; a++) {
                                for (int b = 0; b < 3; b++)
                                        mat_set(&AtA, a, b, mat_get(AtA, a, b) + rows[r][a]*rows[r][b]);
                                Atb.d[a] += rows[r][a] * rhs[r];
                        }
        }

        double det = AtA.d[0]*(AtA.d[4]*AtA.d[8] - AtA.d[5]*AtA.d[7])
                   - AtA.d[1]*(AtA.d[3]*AtA.d[8] - AtA.d[5]*AtA.d[6])
                   + AtA.d[2]*(AtA.d[3]*AtA.d[7] - AtA.d[4]*AtA.d[6]);
        if (fabs(det) < 1e-12) return 0;

        mat_t x = mat_mul(mat3_inv(AtA), Atb);
        *out = (vector_3d_t){ x.d[0], x.d[1], x.d[2] };

        /* Homogeneous in (x - c_i): take the in-front branch if one exists, else return as-is. */
        vector_3d_t cmean;
        if (count_bad_depths(cl, n, *out, &cmean) > 0) {
                vector_3d_t flip = { 2.0*cmean.x - out->x,
                                     2.0*cmean.y - out->y,
                                     2.0*cmean.z - out->z };
                if (count_bad_depths(cl, n, flip, &cmean) == 0)
                        *out = flip;
        }
        return 1;
}

