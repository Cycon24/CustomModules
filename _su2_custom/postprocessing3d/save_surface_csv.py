# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 15:30:14 2026

@author: BriceM
"""
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from scipy.spatial.distance import cdist


def import_csv_to_df(filename:str, filepath:str):
    df = pd.read_csv(filepath+"\\" +  filename)
    return df

def extract_points_in_plane(plane_x:float, dataframe, tol:float=1e-6):
    new_df = dataframe[(abs(dataframe['x']-plane_x)) < tol]
    return new_df

def generate_plane_csv(plane_x:float, filename:str, filepath:str, plane_filename:str="plane") -> None:
    extract_points_in_plane(plane_x,
                            import_csv_to_df(filename, filepath)
                            ).save_csv(filepath + "\\" + plane_filename + ".csv", index=False)
    return None
    


def mass_flow_rate_yz(
    df: pd.DataFrame,
    edge_factor: float = 3.0,  # larger = less aggressive masking
    return_diagnostics: bool = False,
    gc:float = 32.174,
):
    """
    Compute mass flow through a surface in the YZ plane:
        mdot = ∬ rho(y,z) * Vx(y,z) dA

    Uses 2D Delaunay triangulation on (y,z) points and sums flux over triangles.

    Parameters
    ----------
    df : DataFrame with columns y, z, Vx, rho (names configurable).
    edge_factor : masks triangles that have any edge longer than edge_factor * h,
                  where h is an estimated local point spacing.
    return_diagnostics : if True, also returns triangulation + per-triangle arrays.

    Returns
    -------
    mdot : float
        Mass flow rate in lbm/s, is converted internally from slug/s.
    (optional) diagnostics dict
    """

    yz = df[["y", "z"]].to_numpy(dtype=float)
    rho = df["Density"].to_numpy(dtype=float)
    Vx  = df["Velocity_x"].to_numpy(dtype=float)

    if yz.shape[0] < 3:
        raise ValueError("Need at least 3 points to define an area.")

    # Integrand at points
    q = rho * Vx  # rho*Vn

    # 2D triangulation in YZ plane
    tri = Delaunay(yz)
    simplices = tri.simplices  # (ntri, 3) indices into points

    # --- Triangle areas in 2D (y,z)
    p0 = yz[simplices[:, 0], :]
    p1 = yz[simplices[:, 1], :]
    p2 = yz[simplices[:, 2], :]

    # area = 0.5 * |det([p1-p0, p2-p0])|
    v1 = p1 - p0
    v2 = p2 - p0
    areas = 0.5 * np.abs(v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0])

    # --- Optional: mask "bad" triangles that span too far (helps with concave regions / gaps)
    # Estimate a characteristic point spacing h using nearest-neighbor distances (fast-ish for ~few k points)
    # If you have huge point sets, swap to KDTree.
    d = cdist(yz, yz)
    np.fill_diagonal(d, np.inf)
    h = np.median(np.min(d, axis=1))  # typical nearest-neighbor spacing

    # Compute triangle edge lengths
    e01 = np.linalg.norm(p1 - p0, axis=1)
    e12 = np.linalg.norm(p2 - p1, axis=1)
    e20 = np.linalg.norm(p0 - p2, axis=1)
    max_edge = np.maximum.reduce([e01, e12, e20])

    good = (areas > 0.0) & (max_edge <= edge_factor * h)

    # --- Triangle-averaged integrand
    q_tri = q[simplices].mean(axis=1)

    mdot = np.sum(q_tri[good] * areas[good]) * gc

    if not return_diagnostics:
        return mdot

    diagnostics = {
        "triangulation": tri,
        "simplices": simplices,
        "areas": areas,
        "q_tri": q_tri,
        "good_mask": good,
        "h_spacing": h,
        "edge_factor": edge_factor,
        "used_triangle_fraction": float(np.mean(good)),
        "total_area_used": float(np.sum(areas[good])),
    }
    return mdot, diagnostics



if __name__=="__main__":
    import matplotlib.pyplot as plt
    
    filename = "entire_surface_restart.csv"
    filelocation = r"C:\Users\BriceM\Documents\SU2 CFD Data\3D_Tests\Pb_sweep_coarse\Pb_1693"
    
    imported = import_csv_to_df(filename, filelocation)
    inlet = extract_points_in_plane(0/12, imported, tol=0.01/12)
    mdot, diag = mass_flow_rate_yz(inlet, return_diagnostics=True)
    
    print(f"mdot_s = {mdot:.4f} lbm/s")
    print(f"mdot_t = {6*mdot:.4f} lbm/s")
    print(f"mdot_t = {6*mdot/2.2:.4f} kg/s")
    plt.figure()
    plt.scatter(inlet["y"], inlet["z"])
    plt.ylabel("z [ft]")
    plt.xlabel("y [ft]")
    
    
    