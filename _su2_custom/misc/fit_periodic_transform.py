# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 15:22:11 2025

@author: BriceM
"""

import numpy as np, re, math
from pathlib import Path

SU2_FILE = "C:\\Users\\BriceM\\Documents\\SU2 CFD Data\\3D_Tests\\3d_test.su2"
M1, M2 = "Symmetry1", "Symmetry2"  # receiver, donor order you *intend* to use

# ---------- SU2 reader (nodes + marker node sets) ----------
lines = Path(SU2_FILE).read_text(encoding="utf-8", errors="ignore").splitlines()
i_n = next(i for i,l in enumerate(lines) if l.strip().startswith("NPOIN="))
n = int(re.search(r"NPOIN=\s*(\d+)", lines[i_n]).group(1))
xyz = np.loadtxt(lines[i_n+1:i_n+1+n])[:, :3]  # X Y Z

marker_i = [i for i,l in enumerate(lines) if l.strip().startswith("MARKER_TAG=")]
def nodes_of(tag):
    i = next(i for i in marker_i if lines[i].split("=")[1].strip() == tag)
    ne = int(lines[i+1].split("=")[1])
    ids = set()
    for cl in lines[i+2:i+2+ne]:
        parts = list(map(int, cl.split()))
        ids.update(parts[2:])  # node ids (1-based)
    return np.array(sorted(ids)) - 1

idx1, idx2 = nodes_of(M1), nodes_of(M2)
P = xyz[idx1]   # receiver
Q = xyz[idx2]   # donor

# ---------- Fit planes to each marker to get normals ----------
def fit_plane(points):
    c = points.mean(0)
    U, S, Vt = np.linalg.svd(points - c, full_matrices=False)
    nrm = Vt[-1]  # normal = last right-singular vector
    if nrm[0] < 0:  # orient roughly +x for consistency (optional)
        nrm = -nrm
    return c, nrm/np.linalg.norm(nrm)

cP, nP = fit_plane(P)
cQ, nQ = fit_plane(Q)

# The wedge axis should be parallel to +x. Compute the line that is the
# intersection of the two planes and project it to be parallel to +x.
axis_dir = np.array([1.0, 0.0, 0.0])  # enforce +x axis
# Find best (Yc,Zc): minimize distance between rotated donor and receiver
# when rotating around the axis line x-parallel through (Yc,Zc).

def rotate_about_x(v, theta_deg):
    th = math.radians(theta_deg)
    R = np.array([[1, 0, 0],
                  [0,  math.cos(th), -math.sin(th)],
                  [0,  math.sin(th),  math.cos(th)]], dtype=float)
    return v @ R.T

def residual(Yc, Zc, theta_deg, Tx=0.0, Ty=0.0, Tz=0.0):
    # rotate donor Q around x-axis through (0,Yc,Zc), then translate
    Qshift = Q.copy()
    Qshift[:,1] -= Yc
    Qshift[:,2] -= Zc
    Qrot = rotate_about_x(Qshift, theta_deg)
    Qrot[:,1] += Yc + Ty
    Qrot[:,2] += Zc + Tz
    Qrot[:,0] += Tx  # allow tiny axial shift if needed
    d2 = ((Qrot - P)**2).sum(1)
    return math.sqrt(d2.mean()), Qrot

# 1) coarse scan for theta (deg) around 60 and 120 to see which bucket is right
cands = np.linspace(0, 180, 361)
rmse_theta = []
for th in cands:
    r, _ = residual(cP[1], cP[2], th)  # start around receiver centroid
    rmse_theta.append(r)
th0 = float(cands[int(np.argmin(rmse_theta))])

# 2) refine (Yc,Zc,theta) by local search
Yc, Zc, th = cP[1], cP[2], th0
step = 0.1 * max(np.ptp(P[:,1]), 1.0)
for _ in range(50):
    improved = False
    for dY in (0, +step, -step):
        for dZ in (0, +step, -step):
            for dT in (0, +2.0, -2.0):
                r, _ = residual(Yc+dY, Zc+dZ, th+dT)
                if r + 1e-9 < residual(Yc, Zc, th)[0]:
                    Yc, Zc, th = Yc+dY, Zc+dZ, th+dT
                    improved = True
    step *= 0.5
    if not improved: break

# 3) with Yc,Zc,theta fixed, fit tiny translation (Tx,Ty,Tz) by least squares
_, Qrot = residual(Yc, Zc, th)
Delta = P - Qrot
Tx = float(np.median(Delta[:,0]))
Ty = float(np.median(Delta[:,1]))
Tz = float(np.median(Delta[:,2]))
rmse_final, _ = residual(Yc, Zc, th, Tx, Ty, Tz)

print(f"Estimated center (Yc,Zc) = ({Yc:.6f}, {Zc:.6f})")
print(f"Estimated rotation about +x θ = {th:.3f} deg")
print(f"Estimated translation (Tx,Ty,Tz) = ({Tx:.6e}, {Ty:.6e}, {Tz:.6e})")
print(f"Final RMSE = {rmse_final:.6e} (mesh units)")

# Emit SU2 line (2 names + 9 numbers; center first, then angles [deg], then translation)
print("\nPaste this (verify marker order/sign):")
print(f"MARKER_PERIODIC = ( {M1}, {M2}, "
      f"0.0, {Yc:.9f}, {Zc:.9f}, "
      f"{th:.9f}, 0.0, 0.0, "
      f"{Tx:.9f}, {Ty:.9f}, {Tz:.9f} )")
