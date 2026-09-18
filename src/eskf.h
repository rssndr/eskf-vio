#ifndef ESKF_H
#define ESKF_H
#include "imu.h"
#include "mat.h"
#include "quat.h"

#define MAX_CLONES 10

/* SIG_A is not the datasheet value: inflated to cover unmodelled error sources.
 * The others are the ADIS16448 datasheet densities. */
#define SIG_A   1.0000e-2
#define SIG_G   1.6968e-4
#define SIG_BA  3.0000e-3
#define SIG_BG  1.9393e-5

typedef struct {
        quaternion_t q;
        vector_3d_t  pos;
        double       timestamp;
} clone_t;

typedef struct {
        quaternion_t q;
        vector_3d_t pos, vel;
        vector_3d_t ba, bg;
        clone_t clones[MAX_CLONES];
        int n_clones;
        mat_t P;
} eskf_t;

void eskf_init(eskf_t *f, quaternion_t q, vector_3d_t pos, vector_3d_t vel,
               vector_3d_t ba, vector_3d_t bg);
void eskf_predict(eskf_t *f, imu_sample_t s, double dt);
void eskf_update_pos(eskf_t *f, vector_3d_t z, double sigma_z);
int  eskf_update_gravity(eskf_t *f, vector_3d_t a_m, double gate, double sigma_a);
int  eskf_update_zupt(eskf_t *f, double sigma_v);
int  eskf_update_vel(eskf_t *f, vector_3d_t v_ref, double sigma_v);
void eskf_augment(eskf_t *f, double timestamp);
void eskf_inject(eskf_t *f, const mat_t *dx);

#endif

