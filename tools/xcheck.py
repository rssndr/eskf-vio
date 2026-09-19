#!/usr/bin/env python3
"""Independent numpy reimplementation of the MSCKF null-space projection and update.

Reads the fixture dumped by MSCKF_XCHECK=<n> and recomputes the projection, the
innovation statistic, the Huber weight, the Kalman gain and the covariance update
from the same inputs.

What this covers: src/msckf.c from the null-space construction onward -- the
Gram-Schmidt basis, the two thresholds, the chi-squared knee and its `m` index,
the Huber reweight, the gain, the (I - K H) P update, and the state correction.

What it does not cover: triangulation (src/tri.c) and the observation Jacobian.
Inputs to those are taken from the C dump as given, so an error inside them is
invisible here.

Comparison is on quantities that do not depend on the choice of null-space basis:
`m`, gamma, w, dx and P after the update. Two different bases of the same subspace
give bit-comparable dx and P, so a mismatch means the subspaces differ, not that
the bases do.

Usage: python3 tools/xcheck.py <dump.txt>
"""
import sys
import numpy as np

# Must match CHI2_99 in src/msckf.c: entry [d-1] is the 99% quantile for d dof.
CHI2_99 = [6.635, 9.210, 11.345, 13.277, 15.086, 16.812, 18.475, 20.090,
           21.666, 23.209, 24.725, 26.217, 27.688, 29.141, 30.578, 32.000, 33.409]


def parse(path):
    tracks, cur = [], None
    for line in open(path):
        p = line.split()
        if not p:
            continue
        if p[0] == "XC":
            cur = {key: float(val) for key, val in
                   (kv.split("=") for kv in p[1:])}
            tracks.append(cur)
        elif p[0] in ("R", "HX", "HF", "P0", "DX", "P1") and cur is not None:
            rows, cols = int(p[1]), int(p[2])
            cur[p[0]] = np.array([float(x) for x in p[3:]]).reshape(rows, cols)
    return tracks


def recompute(t):
    """Recompute the update from the C inputs. Returns (m, gamma, w, dx, P1)."""
    sigma = t["sigma"]
    Hf, Hx, r, P0 = t["HF"], t["HX"], t["R"], t["P0"]
    rows = Hf.shape[0]
    n = P0.shape[0]

    # Left null space of Hf. The basis is not unique; the update is.
    U, s, _ = np.linalg.svd(Hf, full_matrices=True)
    rank = int((s > s[0] * 1e-12).sum()) if s.size else 0
    A = U[:, rank:].T
    m = rows - rank

    rp = A @ r
    Hp = A @ Hx
    S = Hp @ P0 @ Hp.T + sigma * sigma * np.eye(m)

    gamma = float((rp.T @ np.linalg.solve(S, rp))[0, 0])

    w = 1.0
    knee = CHI2_99[m - 1]
    if gamma > knee:
        w = float(np.sqrt(knee / gamma))
        S = S + (sigma * sigma / w - sigma * sigma) * np.eye(m)

    Sinv = np.linalg.inv(S)
    K = P0 @ Hp.T @ Sinv
    dx = K @ rp
    P1 = (np.eye(n) - K @ Hp) @ P0
    return m, gamma, w, dx, P1


def rel(a, b):
    """Relative difference, floored so a near-zero magnitude does not inflate it."""
    return np.abs(a - b).max() / max(np.abs(b).max(), 1e-12)


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    tracks = parse(sys.argv[1])
    if not tracks:
        print("no XC records found in %s" % sys.argv[1])
        return 2

    tol_m, tol_num = 0, 1e-6
    worst, fails = {}, 0
    print("%-4s %-4s %-6s %-6s %-12s %-12s %-12s %-12s" %
          ("trk", "k", "m(C)", "m(py)", "gamma rel", "w rel", "dx rel", "P1 rel"))
    for i, t in enumerate(tracks):
        m, gamma, w, dx, P1 = recompute(t)
        dg, dw = rel(gamma, t["gamma"]), rel(w, t["w"])
        dd, dp = rel(dx, t["DX"]), rel(P1, t["P1"])
        print("%-4d %-4d %-6d %-6d %-12.3e %-12.3e %-12.3e %-12.3e" %
              (i, int(t["k"]), int(t["m"]), m, dg, dw, dd, dp))
        for name, val in (("gamma", dg), ("w", dw), ("dx", dd), ("P1", dp)):
            worst[name] = max(worst.get(name, 0.0), val)
        if int(t["m"]) != m:
            print("  FAIL: null-space dimension m differs — the C projection spans "
                  "a different subspace than the exact one")
            fails += 1
        elif max(dg, dw, dd, dp) > tol_num:
            print("  FAIL: m agrees but a recomputed quantity exceeds %g" % tol_num)
            fails += 1

    print()
    print("tracks: %d   m mismatches: %d" % (len(tracks), fails))
    print("worst relative differences: " +
          "  ".join("%s %.3e" % (k, v) for k, v in sorted(worst.items())))
    print()
    if fails:
        print("XCHECK: FAIL")
        return 1
    print("XCHECK: PASS — the independent numpy projection and update reproduce "
          "the C results to within %g relative, and m matches 2k - rank(Hf) in "
          "every track." % tol_num)
    return 0


if __name__ == "__main__":
    sys.exit(main())
