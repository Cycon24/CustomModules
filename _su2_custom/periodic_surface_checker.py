# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 15:06:11 2025

@author: BriceM
"""

import numpy as np, re

su2_path = "C:\\Users\\BriceM\\Documents\\SU2 CFD Data\\3D_Tests\\3d_test.su2"
m1, m2 = "Symmetry1", "Symmetry2"   # adjust if swapped

# --- minimal SU2 parser for nodes and marker faces ---
with open(su2_path, "r", encoding="utf-8", errors="ignore") as f:
    lines = f.readlines()

# nodes
n_start = next(i for i,l in enumerate(lines) if l.strip().startswith("NPOIN="))
n = int(re.search(r"NPOIN=\s*(\d+)", lines[n_start]).group(1))
xyz = np.loadtxt(lines[n_start+1:n_start+1+n])[:, :3]  # assume X Y Z

# gather elements to get which nodes belong to markers
# find start of markers
marker_idxs = [i for i,l in enumerate(lines) if l.strip().startswith("MARKER_TAG=")]
def read_marker_nodes(tag):
    # find marker
    i = next(i for i in marker_idxs if lines[i].strip().endswith(tag))
    # read number of elems on this marker
    n_e = int(re.search(r"MARKER_ELEMS=\s*(\d+)", lines[i+1]).group(1))
    conn_lines = lines[i+2:i+2+n_e]
    nodes = set()
    for cl in conn_lines:
        parts = list(map(int, cl.split()))
        etype, nn, ids = parts[0], parts[1], parts[2:]
        # collect node ids (SU2 is 1-based)
        nodes.update(ids)
    return np.array(sorted(nodes)) - 1  # to 0-based

idx1 = read_marker_nodes(m1)
idx2 = read_marker_nodes(m2)
P = xyz[idx1]
Q = xyz[idx2]

# center both around origin (for pure rotation test)
Pc = np.zeros_like(P) + .125 # P - P.mean(0)
Qc = np.zeros_like(P) + 0.125 # Q - Q.mean(0)

# Orthogonal Procrustes: find R minimizing ||R*Pc - Qc||
U, S, Vt = np.linalg.svd(Qc.T @ Pc)
R = U @ Vt
if np.linalg.det(R) < 0:  # ensure a rotation, not reflection
    U[:,-1] *= -1
    R = U @ Vt

# rotation angle about x from R
theta_x = np.degrees(np.arctan2(R[2,1], R[2,2]))  # rotation around x that best explains mapping
rmse = np.sqrt(((Qc - Pc @ R.T)**2).sum(1).mean())

print(f"Best-fit rotation about x ≈ {theta_x:.3f} deg, RMSE ≈ {rmse:.4e} (mesh units)")
