#include <math.h>
#include "frontend.h"
#include "fast.h"
#include "cam.h"

#define TARGET     200
#define FAST_T     30
#define WIN        7
#define LEVELS     4
#define FB_MAX     1.0
#define MIN_DIST   10

void frontend_init(frontend_t *fe) {
        fe->has_prev = 0;
        fe->n = 0;
        fe->next_id = 0;
        fe->n_dead = 0;
}

void frontend_free(frontend_t *fe) {
        if (fe->has_prev) image_free(&fe->prev);
        fe->has_prev = 0;
}

static void replenish(frontend_t *fe, const image_t *frame) {
        static corner_t raw[100000];
        int nr = fast_detect(frame, FAST_T, raw, 100000);

        int nk = 0;
        for (int i = 0; i < nr; i++) {
                int clear = 1;
                for (int j = 0; j < fe->n; j++) {
                        double dx = raw[i].x - fe->f[j].pt.x;
                        double dy = raw[i].y - fe->f[j].pt.y;
                        if (dx*dx + dy*dy < MIN_DIST*MIN_DIST) { clear = 0; break; }
                }
                if (clear) raw[nk++] = raw[i];
        }

        corner_t sel[FE_MAX];
        int ns = fast_select(raw, nk, frame->w, frame->h, 5, 8, 6, 5, sel, FE_MAX);

        for (int i = 0; i < ns && fe->n < TARGET; i++) {
                vector_3d_t ray = cam_unproject(&EUROC_CAM0, sel[i].x, sel[i].y);
                feature_t *nf = &fe->f[fe->n++];
                nf->id = fe->next_id++;
                nf->pt = (pt2_t){ sel[i].x, sel[i].y };
                nf->age = 0;
                nf->hist[0] = (pt2_t){ ray.x, ray.y };
                nf->nhist = 1;
        }
}

static void harvest(frontend_t *fe, const feature_t *ft) {
        if (ft->nhist < 3) return;
        dead_track_t *dtk = &fe->dead[fe->n_dead++];
        for (int h = 0; h < ft->nhist; h++)
                dtk->obs[h] = ft->hist[h];
        dtk->nobs = ft->nhist;
}

void frontend_process(frontend_t *fe, const image_t *frame) {
        fe->n_dead = 0;

        if (fe->has_prev && fe->n > 0) {
                pt2_t p0[FE_MAX], p1[FE_MAX], pb[FE_MAX];
                unsigned char st[FE_MAX], stb[FE_MAX];

                for (int i = 0; i < fe->n; i++) { p0[i] = fe->f[i].pt; p1[i] = p0[i]; }
                klt_track(&fe->prev, frame, p0, p1, st, fe->n, WIN, LEVELS);

                for (int i = 0; i < fe->n; i++) pb[i] = p1[i];
                klt_track(frame, &fe->prev, p1, pb, stb, fe->n, WIN, LEVELS);

                int m = 0;
                fe->n_in = fe->n;
                double fbsum = 0.0, dsum = 0.0;
                int nfb = 0, ndsp = 0;
                for (int i = 0; i < fe->n; i++) {
                        double ex = pb[i].x - p0[i].x, ey = pb[i].y - p0[i].y;
                        if (st[i] && stb[i] && ex*ex + ey*ey < FB_MAX*FB_MAX) {
                                double gx = p1[i].x - p0[i].x, gy = p1[i].y - p0[i].y;
                                fbsum += ex*ex + ey*ey;  nfb++;
                                dsum  += sqrt(gx*gx + gy*gy);  ndsp++;
                                fe->f[m] = fe->f[i];
                                fe->f[m].pt = p1[i];
                                fe->f[m].age++;

                                if (fe->f[m].nhist == FE_HIST) {
                                        harvest(fe, &fe->f[m]);
                                        fe->f[m].nhist = 0;
                                }
                                vector_3d_t ray = cam_unproject(&EUROC_CAM0, p1[i].x, p1[i].y);
                                fe->f[m].hist[fe->f[m].nhist++] = (pt2_t){ ray.x, ray.y };
                                m++;
                        } else {
                                harvest(fe, &fe->f[i]);
                        }
                }
                fe->n = m;
                fe->fb_rms = (nfb  > 0) ? sqrt(fbsum / nfb) : 0.0;
                fe->disp   = (ndsp > 0) ? dsum / ndsp : 0.0;
        }

        if (fe->n < TARGET)
                replenish(fe, frame);

        if (fe->has_prev) image_free(&fe->prev);
        image_copy(&fe->prev, frame);
        fe->has_prev = 1;
}

