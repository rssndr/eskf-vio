#include <stdio.h>
#include <math.h>
#include "eskf.h"
#include "ins.h"

extern mat_t mat_skew(vector_3d_t v);
extern mat_t quat_to_R(quaternion_t q);
extern mat_t build_F(quaternion_t q, vector_3d_t a, vector_3d_t w, double dt);
extern mat_t build_Q(double dt);

static int nfail = 0;

static void check(const char *name, int cond) {
        printf("%-34s %s\n", name, cond ? "ok" : "FAIL");
        if (!cond) nfail++;
}

static vector_3d_t mat_apply(mat_t R, vector_3d_t v) {
        mat_t u = mat_zero(3, 1);
        u.d[0] = v.x; u.d[1] = v.y; u.d[2] = v.z;
        mat_t r = mat_mul(R, u);
        return (vector_3d_t){ r.d[0], r.d[1], r.d[2] };
}

static int veq(vector_3d_t a, vector_3d_t b) {
        return fabs(a.x-b.x) < 1e-12 && fabs(a.y-b.y) < 1e-12 && fabs(a.z-b.z) < 1e-12;
}

int main(void) {
        const double s = 0.7071067811865476;

        vector_3d_t v = {1,2,3}, u = {4,5,6};
        vector_3d_t cross = {-3, 6, -3};
        check("skew: [v]x u = v x u", veq(mat_apply(mat_skew(v), u), cross));

        quaternion_t qs[3] = {{s,s,0,0},{s,0,s,0},{s,0,0,s}};
        vector_3d_t  vs[2] = {{1,0,0},{0.3,-1.2,2.5}};
        int weld = 1;
        for (int i = 0; i < 3; i++)
                for (int j = 0; j < 2; j++) {
                        vector_3d_t a = mat_apply(quat_to_R(qs[i]), vs[j]);
                        vector_3d_t b = vs[j];
                        rotate_vector(&b, qs[i]);
                        weld &= veq(a, b);
                }
        check("weld: R(q)v = rotate_vector", weld);

        mat_t R = quat_to_R(qs[0]);
        mat_t RtR = mat_mul(mat_transpose(R), R);
        int orth = 1;
        for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                        orth &= fabs(mat_get(RtR,i,j) - (i==j ? 1.0 : 0.0)) < 1e-12;
        check("R: orthogonal (RtR = I)", orth);

        quaternion_t id = {1,0,0,0};
        vector_3d_t a0 = {1,2,3}, w0 = {0,0,0};
        double dt = 0.01;
        mat_t F = build_F(id, a0, w0, dt);
        check("F: (0,1) block = dt",      fabs(mat_get(F,0,3) - dt)  < 1e-15);
        check("F: (2,4) block = -dt",     fabs(mat_get(F,6,12) + dt) < 1e-15);
        check("F: (1,3) block = -dt",     fabs(mat_get(F,3,9)  + dt) < 1e-15);
        check("F: (1,2) sample entry",    fabs(mat_get(F,3,7) - 3*dt) < 1e-15);
        check("F: diag stays 1",          fabs(mat_get(F,0,0) - 1.0) < 1e-15 &&
                                          fabs(mat_get(F,6,6) - 1.0) < 1e-15);
        check("F: untouched block is 0",  mat_get(F,0,6) == 0.0);

        vector_3d_t wr = {0.3, -0.2, 0.5};
        mat_t Fr = build_F(id, a0, wr, 0.005);
        mat_t B = mat_zero(3,3);
        for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                        mat_set(&B, i, j, mat_get(Fr, 6+i, 6+j));
        mat_t BtB = mat_mul(mat_transpose(B), B);
        int rot = 1;
        for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                        rot &= fabs(mat_get(BtB,i,j) - (i==j?1.0:0.0)) < 1e-9;
        check("F: (2,2) is a rotation at w!=0", rot);

        mat_t Q  = build_Q(0.005);
        mat_t Q2 = build_Q(0.010);
        check("Q: vel variance is sigma_a^2*dt", fabs(mat_get(Q,3,3) - SIG_A*SIG_A*0.005) < 1e-20);
        check("Q: vel variance scales with dt",  fabs(mat_get(Q2,3,3) - 2.0*mat_get(Q,3,3)) < 1e-20);
        check("Q: pos block empty",              mat_get(Q,0,0) == 0.0);

        eskf_t f;
        vector_3d_t zero = {0,0,0};
        eskf_init(&f, id, zero, zero, zero, zero);
        double p0v = mat_get(f.P, 3, 3);
        double p0p = mat_get(f.P, 0, 0);
        imu_sample_t hover = { .timestamp = 0, .gyro = {0,0,0}, .accel = {0,0,9.81} };
        for (int k = 0; k < 200; k++)
                eskf_predict(&f, hover, 0.005);
        int sym = 1;
        for (int i = 0; i < 15; i++)
                for (int j = 0; j < 15; j++)
                        sym &= fabs(mat_get(f.P,i,j) - mat_get(f.P,j,i)) < 1e-9;
        check("P: stays symmetric",       sym);
        check("P: vel variance grows",    mat_get(f.P, 3, 3) > p0v);
        check("P: pos variance grows",    mat_get(f.P, 0, 0) > p0p);

        eskf_t g;
        eskf_init(&g, id, zero, zero, zero, zero);
        for (int k = 0; k < 100; k++) eskf_predict(&g, hover, 0.005);
        eskf_augment(&g, 0.5);
        check("aug: P is 21x21", g.P.rows == 21 && g.P.cols == 21);

        int csym = 1, corner = 1, xrow = 1;
        size_t src[6] = { 0, 1, 2, 6, 7, 8 };
        for (size_t i = 0; i < 21; i++)
                for (size_t j = 0; j < 21; j++)
                        csym &= fabs(mat_get(g.P,i,j) - mat_get(g.P,j,i)) < 1e-12;
        for (size_t i = 0; i < 6; i++)
                for (size_t j = 0; j < 6; j++)
                        corner &= mat_get(g.P, 15+i, 15+j) == mat_get(g.P, src[i], src[j]);
        for (size_t i = 0; i < 6; i++)
                for (size_t k = 0; k < 15; k++)
                        xrow &= mat_get(g.P, 15+i, k) == mat_get(g.P, src[i], k);
        check("aug: symmetric", csym);
        check("aug: clone corner copied", corner);
        check("aug: cross rows copied", xrow);

        double clone_var = mat_get(g.P, 15, 15);
        double cross_iv  = mat_get(g.P, 15, 3);
        for (int k = 0; k < 400; k++) eskf_predict(&g, hover, 0.005);
        int psym = 1;
        for (size_t i = 0; i < 21; i++)
                for (size_t j = 0; j < 21; j++)
                        psym &= fabs(mat_get(g.P,i,j) - mat_get(g.P,j,i)) < 1e-9;
        check("aug: symmetric after predict", psym);
        check("aug: clone block frozen", mat_get(g.P, 15, 15) == clone_var);
        check("aug: IMU-clone cross evolves", mat_get(g.P, 15, 3) != cross_iv);

        for (int c = 0; c < 12; c++) eskf_augment(&g, 1.0 + c);
        check("marg: capped at 10 clones", g.n_clones == 10 && g.P.rows == 75);
        check("marg: oldest gone", g.clones[0].timestamp > 1.5);
        int msym = 1;
        for (size_t i = 0; i < 75; i++)
                for (size_t j = 0; j < 75; j++)
                        msym &= fabs(mat_get(g.P,i,j) - mat_get(g.P,j,i)) < 1e-9;
        check("marg: symmetric at full window", msym);

        printf(nfail ? "FAIL (%d)\n" : "PASS\n", nfail);
        return nfail != 0;
}

