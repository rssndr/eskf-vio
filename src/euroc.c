#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "euroc.h"
#include "imu.h"

size_t euroc_load_imu(const char *path, imu_sample_t **out) {
        FILE *f = fopen(path, "r");
        if (!f) { perror(path); return 0; }

        char line[256];
        size_t n = 0;
        while (fgets(line, sizeof line, f))
                if (line[0] != '#') n++;
        rewind(f);

        imu_sample_t *s = malloc(n * sizeof *s);
        if (!s) { perror("malloc"); fclose(f); return 0; }

        size_t i = 0;
        while (fgets(line, sizeof line, f)) {
                if (line[0] == '#') continue;

                long long ts;
                imu_sample_t *m = &s[i];
                if (sscanf(line, "%lld,%lf,%lf,%lf,%lf,%lf,%lf",
                                &ts,
                                &m->gyro.x,  &m->gyro.y,  &m->gyro.z,
                                &m->accel.x, &m->accel.y, &m->accel.z) != 7)
                        continue;

                m->timestamp = ts / 1e9;
                i++;
        }
        fclose(f);

        *out = s;

        return i;
}

size_t euroc_load_gt(const char *path, gt_sample_t **out) {
        FILE *f = fopen(path, "r");
        if (!f) { perror(path); return 0; }

        char line[512];
        size_t n = 0;
        while (fgets(line, sizeof line, f))
                if (line[0] != '#') n++;
        rewind(f);

        gt_sample_t *s = malloc(n * sizeof *s);
        if (!s) { perror("malloc"); fclose(f); return 0; }

        size_t i = 0;
        while (fgets(line, sizeof line, f)) {
                if (line[0] == '#') continue;

                long long ts;
                gt_sample_t *g = &s[i];
                if (sscanf(line,
                        "%lld,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                        &ts,
                        &g->pos.x, &g->pos.y, &g->pos.z,
                        &g->q.w, &g->q.x, &g->q.y, &g->q.z,
                        &g->vel.x, &g->vel.y, &g->vel.z,
                        &g->gyro_bias.x,  &g->gyro_bias.y,  &g->gyro_bias.z,
                        &g->accel_bias.x, &g->accel_bias.y, &g->accel_bias.z) != 17)
                        continue;

                /* EuRoC stores these up to 6.8e-5 off unit norm, 80% of samples above 1.
                   An attitude must be unit, or a dot product passes 1 and the angle
                   comparison clamps to zero for every small error. */
                double qn = sqrt(g->q.w*g->q.w + g->q.x*g->q.x +
                                 g->q.y*g->q.y + g->q.z*g->q.z);
                if (qn > 0.0) {
                        g->q.w /= qn;  g->q.x /= qn;
                        g->q.y /= qn;  g->q.z /= qn;
                }

                g->timestamp = ts / 1e9;
                i++;
        }
        fclose(f);

        *out = s;

        return i;
}

size_t euroc_load_cam(const char *path, cam_frame_t **out) {
        FILE *f = fopen(path, "r");
        if (!f) { perror(path); return 0; }

        char line[256];
        size_t n = 0;
        while (fgets(line, sizeof line, f))
                if (line[0] != '#') n++;
        rewind(f);

        cam_frame_t *s = malloc(n * sizeof *s);
        if (!s) { perror("malloc"); fclose(f); return 0; }

        size_t i = 0;
        while (fgets(line, sizeof line, f)) {
                if (line[0] == '#') continue;

                long long ts;
                cam_frame_t *m = &s[i];
                if (sscanf(line, "%lld,%63s", &ts, m->filename) != 2)
                        continue;

                m->timestamp = ts / 1e9;
                i++;
        }
        fclose(f);

        *out = s;

        return i;

}

