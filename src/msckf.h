#ifndef MSCKF_H
#define MSCKF_H
#include "eskf.h"
#include "klt.h"

typedef struct {
        double Hp[6];
        double Hth[6];
        double Hf[6];
        pt2_t  z;
        int    valid;
        /* Camera-frame depth of the point in this clone's frame. Always set,
         * valid or not, so a caller can tell "behind the camera" (negative)
         * from "on the lens" (small positive) from "off the top of the
         * frustum". */
        double depth;
} obs_jac_t;

obs_jac_t obs_jacobian(const clone_t *cl, vector_3d_t pf);
mat_t msckf_nullspace(mat_t Hf);
/* Rejection-reason counters. Diagnostics only: no estimator path reads these.
 * msckf_update_track() has six exit paths and the caller cannot tell them
 * apart, so a "the gate rejects too much" reading of a rejection count is a
 * guess. That guess was wrong once, at a cost of one 13-minute run. */
enum {
        MSCKF_REJ_K_RANGE = 0,   /* k < 2 or 2k > 40 */
        MSCKF_REJ_TRIANG,        /* triangulate() failed */
        MSCKF_REJ_JACOBIAN,      /* obs_jacobian() invalid (feature depth <= 0.1) */
        MSCKF_REJ_NULLSPACE,     /* null-space dimension out of range */
        MSCKF_REJ_CHOL_Y,        /* Cholesky solve S*y = rp failed */
        MSCKF_REJ_CHOL_INV,      /* Cholesky solve for S^-1 failed */
        MSCKF_REJ_COUNT
};
extern int msckf_rej_count[MSCKF_REJ_COUNT];
/* Diagnostic: observations whose Jacobian was invalid, and observations that
 * were valid, summed over the tracks discarded for that reason. The ratio
 * decides the fix: near 1 invalid of 9 means one bad observation kills an
 * otherwise usable track; near 9 of 9 means the triangulation itself is
 * garbage. */
extern int msckf_invalid_obs;
extern int msckf_valid_obs;

/* w_out (may be NULL) receives the Huber weight applied to this track: 1.0
 * when the residual is inside the knee, sqrt(knee/gamma) below it. */
int msckf_update_track(eskf_t *f, const int *ci, const pt2_t *obs, int k, double sigma,
                       double *w_out);

#endif

