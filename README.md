# eskf-vio: MSCKF visual-inertial odometry in C

15-state ESKF prediction at 200 Hz fused with MSCKF visual corrections, from hand-rolled FAST-9
+ Lucas-Kanade tracking and a pinhole + radial-tangential camera model. IMU and one camera
only; no external input after the initial state.

EuRoC Machine Hall, shipped configuration: **ATE 0.376 m on MH_01**, and **0.22–0.96 m across
MH_01–05 with no per-sequence tuning**. Dead reckoning alone drifts 3808 m on the same sequence.

Built from scratch as a learning vehicle. Not competitive with mature systems: OpenVINS-class
filters are more accurate on these sequences — though that is the literature's claim and is not
measured here, since no reference implementation has been run.

## Results

MH_01, initialised from ground truth at t0 only. Ground truth ends 181.9 s into the run while the
sensors continue for another 1.07 s, so every figure below is scored over **0–181.9 s** and the 22
camera frames past that are dropped — scoring them compares the estimate against a frozen
reference. `ATE` is the RMS position error over that window. `ratio` is ‖error‖ / √tr(P_pos):
1.0 means the covariance tells the truth.

| t [s] | pos err [m] | predicted 1σ [m] | ratio |
|---|---|---|---|
| 0–20 | 0.04 | 0.009 | 4.8 |
| 20–40 | 0.06 | 0.021 | 3.1 |
| 40–60 | 0.08 | 0.034 | 2.4 |
| 60–80 | 0.17 | 0.067 | 2.6 |
| 80–100 | 0.22 | 0.078 | 2.9 |
| 100–120 | 0.25 | 0.090 | 2.8 |
| 120–140 | 0.19 | 0.102 | 1.9 |
| 140–160 | 0.43 | 0.120 | 3.5 |
| 160–180 | 0.87 | 0.144 | 6.0 |
| 180–181.9 | 1.13 | 0.155 | 7.3 |

ATE **0.376 m**, peak 1.181 m, final 1.2 m at t = 181.9 s. 62 407 tracks accepted, 233 rejected
(0.4%).

**Breadth.** All five Machine Hall sequences, full rate, no parameter change and no tuning:

| sequence | ATE [m] | peak [m] | rejected |
|---|---|---|---|
| MH_01_easy | 0.376 | 1.181 | 233 |
| MH_02_easy | **0.220** | 0.526 | 281 |
| MH_03_medium | 0.378 | 0.736 | 55 |
| MH_04_difficult | 0.559 | 1.062 | 519 |
| MH_05_difficult | 0.964 | 1.530 | 615 |

The error tracks the sequence difficulty rating, which is the shape a working estimator should
have. It is one scene family, though: the Vicon Room sequences are not in this set.

**Consistency.** The filter is **1.9–7.3× overconfident** on position, best in the middle of the
run. The measurement noise model, however, is calibrated: the innovation statistic
`gamma/dof` = 0.79 against the 1.0 that a correct `R` predicts, so `R` is right to within 21%
and the defect is in how `P` contracts. Two things localise it further. The ratio falls when the
vision rate falls — 4.46 → 2.61 halving 20 → 10 Hz at constant observations per frame — which is
the signature of **time-correlated measurement error**: `R` is diagonal while the tracker's
error is not white, so the covariance shrinks with the measurement count as if the errors were
independent. And `P` projected onto the unobservable directions (the three global translations
and rotation about gravity) **grows, it never collapses**, so the classic first-estimate-Jacobian
failure is absent here and FEJ is not what the ratio is made of. The correlated-error mechanism
is not fixable with a scalar and has no fix in the tree yet.

**Vision rate.** Decimating the camera has a **minimum at 10 Hz**, not at full rate:

| rate | frames | ATE [m] | peak [m] | CPU vs real time (one session) |
|---|---|---|---|---|
| 20 Hz | 3660 | 0.376 | 1.181 | 4.42× |
| **10 Hz** | 1830 | **0.236** | **0.615** | **2.57×** |
| 5 Hz | 915 | 0.706 | 1.395 | 1.25× |
| 3.33 Hz | 610 | 1.184 | 2.986 | — |
| 2 Hz | 366 | 11.41 | 30.67 | 0.51× |

Halving the rate buys 1.6× accuracy *and* 1.7× less CPU. The obvious confound — a fixed
clone-window spans more time at lower rates, which enlarges the baseline and improves depth
observability on its own — was controlled with an optional window duration cap, and the
advantage survives at matched baseline. Below 5 Hz the failure is the **tracker**: mean feature
displacement reaches 27 px while the forward-backward error stays at 0.086–0.210 px, far inside
its 1.0 px gate, even where the filter is 11 m wrong. A consistent mismatch closes the
round trip on itself, so the FB check is a *self-consistency* test and cannot see the collapse.
Displacement is the signal that can, and unlike the error it needs no ground truth.

## What's there

- **INS**, strapdown dead-reckoning: `gyro_to_q` (small-angle guard), accel to world frame, double integration
- **ESKF**, 15-state error propagation (pos, vel, att, ba, bg). Discrete error-state transition per Solà, exact per-tick rotation on the δθ diagonal, Q from the ADIS16448 noise densities published in `eskf.h`. Clone augmentation (6-dim P extension per pose), marginalization of the oldest clone
- **MSCKF**, Mourikis & Roumeliotis: linear multi-view triangulation with a cheirality reflection, per-observation Jacobians, null-space projection of the feature state (Gram-Schmidt), Cholesky solve, Huber robust reweight in place of a binary χ² gate
- **Degeneracy responses**: a zero-velocity update gated on sustained absence of linear acceleration (shipped), and a no-parallax covariance inflation that is measured, falsified and off by default — see the scale degeneracy below
- **VIO front-end**, FAST-9 + NMS + grid selection (~10k raw → 231 spread corners), pyramidal KLT, forward-backward gating, persistent feature IDs, dead-track management for MSCKF (max 12 obs per track), track-length and displacement diagnostics
- **Camera model**, pinhole + radial-tangential distortion, EuRoC cam0 calibration as compile-time constants, T_BS extrinsics

## Architecture

Two-rate: ESKF prediction at IMU rate (200 Hz), MSCKF correction at camera rate. Error state
injected into nominal via Solà right-multiplication, `q_true = q_nom ⊗ q(δθ)`. The covariance
reset Jacobian is deferred.

Clone poses are augmented into `P` each camera frame and marginalized when the window fills
(`MAX_CLONES` = 10, optionally also capped by duration). The null-space projection eliminates
the 3D feature state from the measurement, leaving only pose residuals — the feature is never
estimated and never stored.

**The scale degeneracy.** For about 24 s of the MH_01 run the platform is translation-free
while still rotating, and the camera baseline over the clone window collapses by a factor of
1000 (0.220 m → 0.00023 m). With no parallax there is no metric scale and no depth to observe,
and nothing warns of it — no residual can, since the filter is legitimately consistent
throughout, `gamma/dof` per frame ≈ 0.001. Under the shipped noise set, with no response
enabled, the episode costs up to **1.23 m** of position error and a 0.259 m ATE over the first
60 s. A zero-velocity update
detects the stationary stretch from the accelerometer (`| ‖a‖ − g | < 0.2 m/s²` sustained 1 s,
a detector that never fires before t = 18 s on MH_01) and holds velocity to zero: the window's
peak error falls **1.23 → 0.072 m** and the 60 s ATE **0.259 → 0.065 m**. See the limits below —
the premise is a scenario property.

**Inflating the covariance does not work.** Below a set parallax threshold the translation block
takes the acceleration uncertainty no measurement constrains (`base_ref`, `sigma_d`); the window's
consistency ratio is **3.00** without it and **2.77 / 3.12 / 3.80** at σ_d = 0.02 / 0.05 / 0.20,
so no strength improves the number it exists to fix. Every camera update re-shrinks `P`, and a
larger `P` only raises the gain, so the injected covariance is spent on the next frame and the
estimate takes on more of the measurement's correlated error. σ_d = 0.20 lowers the window peak
1.230 → 0.485 m and the 60 s ATE 0.259 → 0.175 m, and degrades what follows: 0.170 → 0.222 m and
1.47 → 3.34 in ratio. With the shipped zero-velocity update on, both error and ratio get worse
(peak 0.072 → 0.080 m, window ratio 3.07 → 4.28). The function stays in the tree with `base_ref`
and `sigma_d` at 0, and the six configurations are recorded beside it.

## Building and running

```bash
make test     # compile + run all unit tests
make run      # full EuRoC pipeline at the shipped configuration
make clean    # remove build artifacts
```

GCC and `libm` only. The dataset is expected at `~/datasets/euroc/MH_01_easy/` (override
`DATA_IMU`, `DATA_GT`, `DATA_CAM` in the Makefile).

`make run` uses the compiled-in shipped configuration. Everything is overridable on the
command line, positionally:

```bash
./build/main <imu.csv> <gt.csv> <cam0/> <out.csv> [t_end] [sigma_px] [grav_gate]
             [grav_sigma] [zupt_sigma] [hold_sigma] [decim] [max_span]
```

| arg | default | meaning |
|---|---|---|
| `out.csv` | `diag.csv` | one row per camera frame; `sigma_1d`, `ratio_1d`, `ratio_3d`, `att_err_deg`, `tilt_deg`, `yaw_deg`, `n_clones`, `ok`, `rej`, `base`, `disp` |
| `t_end` | 0 (all) | stop after this many seconds of flight; clamped to the reference extent, so the default scores everything ground truth covers |
| `sigma_px` | **0.65** | measurement noise in pixels; the shipped value is the fixed point of the NIS calibration rule |
| `grav_gate`, `grav_sigma` | 0, 0.5 | accelerometer-as-gravity update. Measured to degrade the estimate in six configurations — kept, off by default |
| `zupt_sigma` | **0.05** | zero-velocity update, m/s. On by default |
| `hold_sigma` | 0 | velocity-hold variant; the assumption-free alternative to the ZUPT |
| `decim` | 1 | keep every n-th camera frame |
| `max_span` | 0 | clone-window duration cap [s]; 0 caps by count only |
| `base_ref`, `sigma_d` | 0, 0 | no-parallax covariance inflation, m and m/s². Measured not to improve consistency at any strength — kept, off by default |

## Verification

- **11 test binaries, 99 checks**, `make test`. Jacobians against finite differences to 1e-8, null-space annihilation and orthonormality, perfect-observation zero-correction, corrupted-clone recovery, outlier rejection with `P` untouched, triangulation to 8.9e-16 m on exact observations.
- **`tools/xcheck.py`** reimplements the null-space projection and the update independently in numpy, from a dump the C emits when `MSCKF_XCHECK=<n>` is set. On 120 real tracks the two agree at machine precision: `m = 2k−3` exactly at every track length, the innovation statistic to 6.3e-13, the Huber weight to 6.7e-15, the state correction to 2.3e-9 and `P` to 2.3e-12. It compares only basis-independent quantities, so a mismatch means the subspaces differ rather than the bases. It does **not** cover triangulation or the observation Jacobian.
- **Instrumentation is passive.** Every diagnostic added to the driver was checked bit-identical against a frozen baseline binary before its numbers were trusted.

## Honest limits

**The zero-velocity update assumes the vehicle stops, and a drone in cruise does not.** The
accelerometer measures Δv = 0 — velocity *constant*, not zero, and the filter's velocity at that
moment is not zero either. The variant that avoids the assumption, holding velocity at the value
it had when the quiet stretch began, is **worse than doing nothing** under the shipped noise set:
window peak 1.230 → **1.638 m** and window ratio 3.00 → **12.06**, because it freezes a velocity
that is itself wrong and then shrinks the covariance around it. It improves this benchmark by
more than it would improve a real flight profile, and on a platform that never stops it will not
fire.

**One scene family.** All five sequences are the same machine hall, camera, lighting and motion
style. The Vicon Room is the genuinely different regime and is not measured.

**ATE only, and nothing to compare it against.** Relative pose error, the field's usual drift
metric, is not computed. No reference implementation has been run on these sequences either, so
the distance to OpenVINS or VINS-Mono is not measured. Running one is the only way an accuracy
claim means anything against the field.

**The late attitude error is yaw about gravity.** An earlier version of this page reported a
12.9° non-yaw tilt worth 2.2 m/s² of spurious horizontal acceleration. The number was wrong twice.
The diagnostic behind it measured the angle between
the filter's and the truth's **body-z** axes, which is a tilt only if body z is the vertical one:
in EuRoC the accelerometer at rest puts the near-vertical body axis in **x** (`+0.963, −0.054,
−0.264`), leaving body z 105–114° from gravity. A near-horizontal axis is swept almost entirely by
rotation *about* gravity, so the metric was tracking yaw. And the reference stream ends 1.07 s
before the sensors, so the last 22 frames were scored against a frozen sample — which is where the
largest values came from. A third defect sat under both: EuRoC's reference quaternions are not unit
norm (80% exceed 1, by up to 6.8e-5), so a quaternion dot product passed 1, a `nan` guard caught
it, and every error below 0.4–1.3° read as exactly 0.0000. The loader now normalises them.

Split at the gravity axis, over the scored window:

| | mean | max |
|---|---|---|
| tilt (perpendicular to gravity) | 0.28° | **0.64°** |
| yaw (about gravity) | 1.94° | **7.59°** |
| combined | 2.00° | **7.59°** |

The attitude error is 97% yaw about gravity, and the spurious horizontal acceleration from tilt is
**0.11 m/s²**, not 2.2. About two thirds of the endpoint yaw is the integrated residual gyro-bias
error; the rest sits in the weakly observable direction, so it is the same defect as the
consistency ratio above rather than a second one.

**Compute is the port blocker, and a separate axis from accuracy.** One thread, `make run`,
nothing else on the machine: wall 839.7 s for 181.9 s of flight, **4.62× real time**, split
predict 28.1 s (3.3%), front-end 201.7 s (24.0%), **update 608.4 s (72.5%)** — 229 ms per frame
with 166 ms of it the update. Two mechanical defects dominate it: `mat_t` is 51 200 bytes and
every operation passes it by value, and one update step pays 57% of its cost for an algebraic
form that is identical at 5× less. Never move a cost change and an accuracy change together.

## Road map

Done: the strapdown INS, the 15-state ESKF, the VIO front-end and the MSCKF update, the
update-path fixes, the Huber kernel, the zero-velocity update, the decimation curve, breadth
across MH_01–05, and the numpy cross-check.

Open, in rough priority order: implement first-estimate Jacobians, then re-run the noise-parameter
grid with and without them — closing the gap between the physically-derived settings and the
hand-fit ones is the acceptance test; test a per-track common-mode term in `R`
against the correlated-error hypothesis; compute RPE and run a reference implementation on the
same sequences; extend the cross-check to triangulation; and gyro-warm-started KLT to reach below 5 Hz.

Then the closed-source work: ESP32-P4 logger, a learned inertial model fused as an ordinary ESKF
measurement, the MAVLink `ODOMETRY` wire to ArduPilot or PX4, and int8 deployment.

## References

- **Solà** — *Quaternion kinematics for the error-state Kalman filter* — ESKF propagation and update
- **Mourikis & Roumeliotis** — *A multi-state constraint Kalman filter for vision-aided inertial navigation* — MSCKF, null-space projection, clone poses
