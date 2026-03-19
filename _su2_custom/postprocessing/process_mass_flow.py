# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 15:30:14 2026

@author: BriceM
"""
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay, cKDTree
from scipy.spatial.distance import cdist
import _aerodynamics.GasDynamics as GD



def import_csv_to_df(filename:str, filepath:str):
    df = pd.read_csv(filepath+"\\" +  filename)
    return df

def extract_points_in_plane(plane_x:float, dataframe, tol:float=1e-6):
    '''
    Extracts a slice (all points) from dataframe at a specified axial position (x) in the units
    of the dataframe (m as of 20260216) within a tolerance of tol. 
    
    1 m = 100/2.54 in 

    Parameters
    ----------
    plane_x : float
        Axial location.
    dataframe : DataFrame
        dataframe of entire mesh and the values at each point.
    tol : float, optional
        Maximum axial distance from which points are pulled from plane_x. The default is 1e-6.

    Returns
    -------
    new_df : DataFrame
        The sliced dataframe containing only points within tol of the plane at plane_x.

    '''
    new_df = dataframe[(abs(dataframe['x']-plane_x)) < tol]
    return new_df

def generate_plane_csv(plane_x:float, filename:str, filepath:str, plane_filename:str="plane") -> None:
    '''
    Save the plane pulled from plane_x to a csv file.

    Parameters
    ----------
    plane_x : float
        DESCRIPTION.
    filename : str
        DESCRIPTION.
    filepath : str
        DESCRIPTION.
    plane_filename : str, optional
        DESCRIPTION. The default is "plane".

    Returns
    -------
    None
        DESCRIPTION.

    '''
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


def total_pressure_slice(
    df: pd.DataFrame,
    Gamma: float = 1.4,
    ):
    """
    Compute mass total pressure in a slice of the mesh. Calculated
        with isentropic relations from mach and pressure at each point

    Parameters
    ----------
    df : DataFrame with columns y, z, Pressyre, Mach 
    Gamma : Float equal to the specific heat ratio

    Returns
    -------
    Ptot : array of floats
        Total pressures in same units as original df
    new_df : DataFrame with a Total Pressure column added
    """
    p_col = "Pressure"
    M_col = "Mach"
    
    
    # for r in range(0, len(Ptots)):
    Ptots = np.array(df[p_col] * GD.Po_P_ratio(df[M_col], Gamma))
    new_df = df.copy() 
    new_df["Total Pressure"] = Ptots
    
    return Ptots, new_df


def swirl_number_ang(angle #df: pd.DataFrame,
                #edge_factor: float = 3.0,  # larger = less aggressive masking
                #return_diagnostics: bool = False,
                #gc:float = 1
                ):
    SN = (2/3) * np.tan(np.deg2rad(angle))
    return SN

def swirl_angle(SN #df: pd.DataFrame,
                #edge_factor: float = 3.0,  # larger = less aggressive masking
                #return_diagnostics: bool = False,
                #gc:float = 1
                ):
    angle = np.rad2deg(np.arctan((3/2) * SN))
    return angle


def swirl_number_yz(
    df: pd.DataFrame,
    edge_factor: float = 3.0,
    return_diagnostics: bool = False,
    center: tuple[float, float] | None = (0,0),   # (y0, z0); if None, uses mean(y), mean(z)
    R_ref: float | None = None,                  # if None, uses max(r) from points
    R_ref_mode: str = "max",                     # "max" or "eq_area" (sqrt(A/pi))
    y_col: str = "y",
    z_col: str = "z",
    rho_col: str = "Density",
    vx_col: str = "Velocity_x",
    vy_col: str = "Velocity_y",
    vz_col: str = "Velocity_z",
):
    """
    Compute swirl number S about the x-axis for a cross-section in the YZ plane.

        S = ( ∬ rho * Vx * Vtheta * r dA ) / ( R_ref * ∬ rho * Vx^2 dA )

    - Vtheta is the tangential component in the YZ plane about (y0,z0):
        e_theta = (-sinθ, cosθ) in (y,z)
        Vtheta  = Vy*(-sinθ) + Vz*(cosθ) = (-z_rel/r)*Vy + (y_rel/r)*Vz

    Integration is done by 2D Delaunay triangulation on (y,z) and summing
    triangle-averaged integrands times triangle areas. Optionally masks "bad"
    triangles with overly long edges.

    Parameters
    ----------
    df : DataFrame
        Must contain columns y, z, Density, Velocity_x, Velocity_y, Velocity_z (names configurable).
    edge_factor : float
        Masks triangles with any edge > edge_factor * h, where h is median NN spacing.
    return_diagnostics : bool
        If True, returns (S, diagnostics dict).
    center : (y0, z0) or None
        Swirl axis location in the YZ plane. If None, uses (mean(y), mean(z)).
    R_ref : float or None
        Reference radius in the swirl definition. If None:
          - if R_ref_mode == "max": uses max(r) from point set
          - if R_ref_mode == "eq_area": uses sqrt(A/pi), where A is integrated area used
    R_ref_mode : str
        "max" or "eq_area" (only used if R_ref is None)

    Returns
    -------
    S : float
        Swirl number (dimensionless).
    (optional) diagnostics : dict
        Contains triangulation, masks, integrated numerator/denominator, etc.
    """

    yz = df[[y_col, z_col]].to_numpy(dtype=float)
    rho = df[rho_col].to_numpy(dtype=float)
    vx  = df[vx_col].to_numpy(dtype=float)
    vy  = df[vy_col].to_numpy(dtype=float)
    vz  = df[vz_col].to_numpy(dtype=float)

    if yz.shape[0] < 3:
        raise ValueError("Need at least 3 points to define an area.")

    # --- choose center (y0,z0)
    if center is None:
        y0 = float(np.mean(yz[:, 0]))
        z0 = float(np.mean(yz[:, 1]))
    else:
        y0, z0 = float(center[0]), float(center[1])

    y_rel = yz[:, 0] - y0
    z_rel = yz[:, 1] - z0
    r = np.sqrt(y_rel**2 + z_rel**2)

    # --- compute Vtheta about x-axis
    # e_theta = (-sinθ, cosθ) with sinθ=z/r, cosθ=y/r
    # Vtheta = Vy*(-z/r) + Vz*(y/r)
    eps = 1e-30
    inv_r = 1.0 / np.maximum(r, eps)
    vtheta = (-z_rel * inv_r) * vy + (y_rel * inv_r) * vz
    # At r ~ 0, the direction is undefined; set Vtheta=0 there.
    vtheta = np.where(r > 0.0, vtheta, 0.0)

    # --- Triangulate in YZ plane
    tri = Delaunay(yz)
    simplices = tri.simplices  # (ntri, 3)

    p0 = yz[simplices[:, 0], :]
    p1 = yz[simplices[:, 1], :]
    p2 = yz[simplices[:, 2], :]

    v1 = p1 - p0
    v2 = p2 - p0
    areas = 0.5 * np.abs(v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0])

    # --- point spacing estimate via KDTree (fast; avoids O(N^2) cdist)
    tree = cKDTree(yz)
    dists, _ = tree.query(yz, k=2)  # k=1 is self (0), k=2 is nearest neighbor
    nn = dists[:, 1]
    h = float(np.median(nn[np.isfinite(nn)]))

    # --- edge-length mask
    e01 = np.linalg.norm(p1 - p0, axis=1)
    e12 = np.linalg.norm(p2 - p1, axis=1)
    e20 = np.linalg.norm(p0 - p2, axis=1)
    max_edge = np.maximum.reduce([e01, e12, e20])

    good = (areas > 0.0) & (max_edge <= edge_factor * h)

    # --- integrands at points
    num_pt = rho * vx * vtheta * r      # rho*Vx*Vtheta*r
    den_pt = rho * vx * vx             # rho*Vx^2

    # triangle-averaged integrands
    num_tri = num_pt[simplices].mean(axis=1)
    den_tri = den_pt[simplices].mean(axis=1)

    # integrate over cross-section
    N = float(np.sum(num_tri[good] * areas[good]))
    D = float(np.sum(den_tri[good] * areas[good]))

    if D == 0.0:
        raise ZeroDivisionError("Axial momentum flux integral is zero; cannot compute swirl number.")

    # --- choose R_ref if not provided
    if R_ref is None:
        if R_ref_mode.lower() == "max":
            R = float(np.max(r))
        elif R_ref_mode.lower() in ("eq_area", "equivalent_area", "equiv_area"):
            A_used = float(np.sum(areas[good]))
            if A_used <= 0.0:
                raise ValueError("Integrated area is non-positive; cannot form equivalent-area radius.")
            R = float(np.sqrt(A_used / np.pi))
        else:
            raise ValueError("R_ref_mode must be 'max' or 'eq_area'.")
    else:
        R = float(R_ref)

    if R <= 0.0:
        raise ValueError("R_ref must be > 0.")

    S = N / (R * D)

    if not return_diagnostics:
        return S

    diagnostics = {
        "triangulation": tri,
        "simplices": simplices,
        "areas": areas,
        "good_mask": good,
        "h_spacing": h,
        "edge_factor": edge_factor,
        "used_triangle_fraction": float(np.mean(good)),
        "total_area_used": float(np.sum(areas[good])),
        "center_y0": y0,
        "center_z0": z0,
        "R_ref_used": R,
        "numerator_integral": N,
        "denominator_integral": D,
    }
    return S, diagnostics


if __name__=="__main__":
    import matplotlib.pyplot as plt
    
    filename = "entire_surface_restart.csv"
    filelocation = r"C:\Users\BriceM\Documents\SU2 CFD Data\3D_Tests\MedSwirl_01\Test03"
    
    imported = import_csv_to_df(filename, filelocation)
    inlet = extract_points_in_plane(0, imported, tol=0.01e-3)
    outlet = extract_points_in_plane(17*2.54/100, imported, tol=0.01e-3)
    TE_slice = extract_points_in_plane(8*2.54/100, imported, tol=0.10e-3)
    
    mdot, diag = mass_flow_rate_yz(inlet, return_diagnostics=True, gc=1)
    mdot_o, diag_o = mass_flow_rate_yz(outlet, return_diagnostics=True, gc=1)
    R_refs = np.linspace(0.2*2.54/100, 2.0*2.54/100, 100)
    # SNs = []
    # for R_ref in R_refs:
    SN, diag_SN = swirl_number_yz(outlet, return_diagnostics=True, R_ref=None) #2.0*2.54/100)
    SN_TE, diag_SN_TE = swirl_number_yz(TE_slice, return_diagnostics=True, R_ref=None) #2.0*2.54/100)
   
        # SNs.append(SN)
    Ptots_i, inlet =  total_pressure_slice(inlet, 1.4)
    Ptots_o, outlet =  total_pressure_slice(outlet, 1.4)
    dPtots = Ptots_o - Ptots_i 
    pi_s = Ptots_o / Ptots_i
    pi_avg = np.average(pi_s)
    
    dmdot = 6*(mdot_o - mdot)
    # print(f"mdot_s = {mdot:.4f} lbm/s")
    # print(f"mdot_t = {6*mdot:.4f} lbm/s")
    print(r"- $\dot{m}$" + f" = {6*mdot:.4f} kg/s")
    print(r"- $\Delta\dot{m}$" + f" = {dmdot:.4f} kg/s")
    print(r"- $\pi_{avg}$" + f" = {pi_avg:.4f} ")
    print(r"- $\pi_{min}$" + f" = {np.min(pi_s):.4f} ")
    print(r"- $\pi_{max}$" + f" = {np.max(pi_s):.4f} ")
    print(f"- SN     = {SN:.4f}")
    print(r"- $\alpha_{SN}$" + f"     = {swirl_angle(SN):.4f}" + r"$^o$")
    print(r"- SN_{TE}   = "+f"{SN_TE:.4f}")
    print(r"- $\alpha_{SN, TE}$" + f"     = {swirl_angle(SN_TE):.4f}" + r"$^o$")
    # plt.figure()
    # plt.scatter(inlet["y"], inlet["z"])
    # plt.ylabel("z [ft]")
    # plt.xlabel("y [ft]")
    
    
    # fig, ax = plt.subplots(figsize=(8, 6))
    # YY, ZZ = np.meshgrid(inlet["y"], inlet["z"], indexing="xy")
    # cf = ax.contourf(YY, ZZ, pi_s)
    
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(inlet["y"], inlet["z"], pi_s)
    ax.set_xlabel('y')
    ax.set_ylabel('z')
    ax.set_zlabel(r'$\pi$')
    
    # plt.figure()
    # plt.plot(R_refs*100/2.54, SNs)
    
    