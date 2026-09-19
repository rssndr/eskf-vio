#ifndef FRONTEND_H
#define FRONTEND_H
#include "image.h"
#include "klt.h"

#define FE_MAX 256
#define FE_HIST 12

typedef struct {
        int id;
        pt2_t pt;
        int age;
        pt2_t hist[FE_HIST];
        int nhist;
} feature_t;

typedef struct {
        pt2_t obs[FE_HIST];
        int nobs;
} dead_track_t;

typedef struct {
        image_t prev;
        int has_prev;
        feature_t f[FE_MAX];
        int n;
        int next_id;
        dead_track_t dead[FE_MAX];
        int n_dead;
        int n_in;        /* features offered to the tracker this frame */
        double fb_rms;   /* RMS forward-backward error, px */
        double disp;     /* mean pixel motion of the features that survived */
} frontend_t;

void frontend_init(frontend_t *fe);
void frontend_process(frontend_t *fe, const image_t *frame);
void frontend_free(frontend_t *fe);

#endif
