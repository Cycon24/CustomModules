# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 15:20:14 2025

Author: BriceM (+ helper functions by ChatGPT)
Description:
  - Builds a swirler wedge + blade, extrudes along x, boolean CUT (duct - blade)
  - Automatically sets periodicity between wedge side faces (handles split patches)
  - Automatically detects inlet, outlet, inner/outer walls, and airfoil surfaces
    for physical group tagging.

Notes:
  - Periodic rotation is computed from actual wedge edge angles (no hard-coding).
  - Works even when CUT splits the periodic faces into many patches.
"""

import sys
import math
import numpy as np
import gmsh
from pathlib import Path

# Optional: your module import (unchanged)
import _aerodynamics.AirfoilGenerator as AG


# =============================================================================
# Helpers
# =============================================================================

def _wrap_pm_pi(a: float) -> float:
    """Wrap angle to (-pi, pi]."""
    a = (a + math.pi) % (2 * math.pi) - math.pi
    return a


def _get_external_faces() -> list[int]:
    """
    Return tags of all outward boundary faces (dim=2) of all current 3D volumes.
    Faces that appear exactly once in the union of volume boundaries are external.
    """
    from collections import Counter
    cnt = Counter()
    for (dim, vtag) in gmsh.model.occ.getEntities(3):
        for (d, f) in gmsh.model.getBoundary([(3, vtag)], oriented=False, recursive=False):
            if d == 2:
                cnt[f] += 1
    ext = [f for f, c in cnt.items() if c == 1]
    return ext


def _com(face_tag: int) -> tuple[float, float, float]:
    return gmsh.model.occ.getCenterOfMass(2, face_tag)


def _bbox(face_tag: int) -> tuple[float, float, float, float, float, float]:
    return gmsh.model.occ.getBoundingBox(2, face_tag)


def _x_span(face_tag: int) -> float:
    xmin, ymin, zmin, xmax, ymax, zmax = _bbox(face_tag)
    return abs(xmax - xmin)


def _theta_yz_from_point(p: tuple[float, float, float]) -> float:
    _, y, z = p
    return math.atan2(z, y)


def set_wedge_periodic(arcangle: float, r_hub: float, L_tot: float) -> dict:
    """
    Identify periodic lateral faces of the wedge, compute the correct rotation Δ,
    pair split patches in axial order, and call gmsh.model.mesh.setPeriodic.

    Returns a dict with:
      {
        "slaves": [...],
        "masters": [...],
        "R": [16-element row-major 4x4],
        "delta": float,         # radians
        "theta1": float,        # radians
        "theta2": float,        # radians
        "sideA": [...],         # classified faces near theta1
        "sideB": [...],         # classified faces near theta2
      }
    """
    gmsh.model.occ.synchronize()

    # External faces only
    ext_faces = _get_external_faces()

    # Derive target azimuths from your actual radial edges (consistent with how wedge edges are created)
    # inner-radius points were defined as:
    #   r1p1: (0, -r_hub*cos(a/2), +r_hub*sin(a/2))  -> theta1 ≈  π - a/2
    #   r1p2: (0, +r_hub*cos(a/2), +r_hub*sin(a/2))  -> theta2 ≈ +a/2
    y1 = -r_hub * math.cos(arcangle)
    z1 =  r_hub * math.sin(arcangle)
    y2 =  r_hub * math.cos(arcangle)
    z2 =  r_hub * math.sin(arcangle)
    theta1 = math.atan2(z1, y1)
    theta2 = math.atan2(z2, y2)

    # Precompute COMs and filter to axial faces (exclude inlet/outlet by x-span fraction of L_tot)
    faces_com = {f: _com(f) for f in ext_faces}
    axial_faces = [f for f in ext_faces if _x_span(f) > 0.25 * L_tot]

    # Classify by azimuth proximity to theta1/theta2
    def ang_dist(a, b):
        d = (a - b) % (2 * math.pi)
        if d > math.pi:
            d -= 2 * math.pi
        return abs(d)

    sideA, sideB = [], []
    ang_tol = 0.30  # ~17 deg; widen/narrow if your faces are chopped differently
    for f in axial_faces:
        th = _theta_yz_from_point(faces_com[f])
        dA = ang_dist(th, theta1)
        dB = ang_dist(th, theta2)
        if dA < dB and dA < ang_tol:
            sideA.append(f)
        elif dB < dA and dB < ang_tol:
            sideB.append(f)

    if not sideA or not sideB:
        raise RuntimeError("Periodic classification failed: check wedge angles / geometry.")

    # Sort each side by COM.x to establish stable 1-1 patch order along the axis
    sideA.sort(key=lambda f: faces_com[f][0])
    sideB.sort(key=lambda f: faces_com[f][0])

    # Equalize length (if asymmetric split occurred)
    n = min(len(sideA), len(sideB))
    sideA = sideA[:n]
    sideB = sideB[:n]
    if n == 0:
        raise RuntimeError("No overlapping lateral patches to pair for periodicity.")

    # Compute actual angular separation Δ between the side planes
    delta = _wrap_pm_pi(theta1 - theta2)  # for your build this equals (pi - arcangle)
    cd = float(math.cos(delta))
    sd = float(math.sin(delta))

    # 4x4 transforms about x by ±Δ
    R_d  = [1,0,0,0,  0,cd,-sd,0,  0,sd,cd,0,  0,0,0,1]
    R_md = [1,0,0,0,  0,cd, sd,0,  0,-sd,cd,0, 0,0,0,1]

    def rot_d(p):
        x, y, z = p
        return (x, cd*y - sd*z, sd*y + cd*z)

    def rot_md(p):
        x, y, z = p
        return (x, cd*y + sd*z, -sd*y + cd*z)

    def total_err(slaves, masters, rot_fn):
        err2 = 0.0
        for fs, fm in zip(slaves, masters):
            xs, ys, zs = faces_com[fs]
            xm, ym, zm = faces_com[fm]
            xr, yr, zr = rot_fn((xm, ym, zm))  # rotate master COM into slave frame
            err2 += (xr - xs) ** 2 + (yr - ys) ** 2 + (zr - zs) ** 2
        return err2

    # Try both orientations/signs: (A<-B), (B<-A) with ±Δ; pick minimal error
    configs = [
        ("A<-B +Δ", sideA, sideB, R_d,  rot_d),
        ("B<-A +Δ", sideB, sideA, R_d,  rot_d),
        ("A<-B -Δ", sideA, sideB, R_md, rot_md),
        ("B<-A -Δ", sideB, sideA, R_md, rot_md),
    ]
    best = min(configs, key=lambda cfg: total_err(cfg[1], cfg[2], cfg[4]))
    label, slaves, masters, R_use, rot_fn = best

    # Optional diagnostics
    print(f"[Periodic] theta1={theta1:.6f}, theta2={theta2:.6f}, Δ={delta:.6f} rad "
          f"({math.degrees(delta):.3f} deg), choice={label}")

    # Apply periodicity
    gmsh.model.mesh.setPeriodic(2, slaves, masters, R_use)
    gmsh.model.occ.synchronize()

    return {
        "slaves": slaves,
        "masters": masters,
        "R": R_use,
        "delta": delta,
        "theta1": theta1,
        "theta2": theta2,
        "sideA": sideA,
        "sideB": sideB,
    }


def auto_physical_groups(r_hub: float, r_duct: float, L_tot: float,
                         periodic_sides: tuple[list[int], list[int]]) -> dict:
    """
    Auto-detect inlet, outlet, inner/outer walls, and airfoil surfaces.

    Args:
      r_hub, r_duct, L_tot: geometry scales for thresholds
      periodic_sides: (sideA, sideB) lists of periodic faces to exclude from walls

    Returns dict with lists of face tags:
      {
        "inlet": [...],
        "outlet": [...],
        "inner_walls": [...],
        "outer_walls": [...],
        "airfoils": [...],
      }
    """
    gmsh.model.occ.synchronize()

    # External faces only
    ext = _get_external_faces()

    # Precompute COM, bbox, auxiliary
    faces_com = {f: _com(f) for f in ext}
    faces_xc  = {f: faces_com[f][0] for f in ext}
    faces_xspan = {f: _x_span(f) for f in ext}

    # 1) Inlet & Outlet: smallest & largest x-center among "short" x-span faces
    #    (planar-ish cross-sections). We’ll select faces with x-span < 5% of L_tot.
    flat_thresh = 0.05 * L_tot
    flat_faces  = [f for f in ext if faces_xspan[f] <= flat_thresh]
    # If CUT created multiple patches at inlet/outlet, we collect all near min/max x-center.
    # Use a proximity window (± 1% of L_tot) around global min/max xc.
    if flat_faces:
        xmin_center = min(faces_xc[f] for f in flat_faces)
        xmax_center = max(faces_xc[f] for f in flat_faces)
        win = 0.01 * L_tot
        inlet  = [f for f in flat_faces if abs(faces_xc[f] - xmin_center) <= win]
        outlet = [f for f in flat_faces if abs(faces_xc[f] - xmax_center) <= win]
    else:
        # Fallback: pick absolute min/max among all faces
        xmin_center = min(faces_xc[f] for f in ext)
        xmax_center = max(faces_xc[f] for f in ext)
        win = 0.01 * L_tot
        inlet  = [f for f in ext if abs(faces_xc[f] - xmin_center) <= win]
        outlet = [f for f in ext if abs(faces_xc[f] - xmax_center) <= win]

    # 2) Periodic sides (exclude from walls). Provided as (sideA, sideB)
    sideA, sideB = periodic_sides
    periodic_all = set(sideA) | set(sideB)

    # 3) Inner/Outer walls: axial faces whose radius ~ r_hub or r_duct
    #    Use COM radius; require axial faces (long x-span) and not in periodic/inlet/outlet.
    radial_faces = [f for f in ext
                    if (f not in periodic_all) and (f not in inlet) and (f not in outlet)
                    and faces_xspan[f] > 0.25 * L_tot]

    def rad(p): 
        _, y, z = p
        return math.hypot(y, z)

    tol_inner = max(1e-3, 0.03 * r_hub)   # 3% of inner radius or 0.001"
    tol_outer = max(1e-3, 0.02 * r_duct)  # 2% of outer radius or 0.001"

    inner_walls = [f for f in radial_faces if abs(rad(faces_com[f]) - r_hub)  <= tol_inner]
    outer_walls = [f for f in radial_faces if abs(rad(faces_com[f]) - r_duct) <= tol_outer]

    # 4) Airfoils: whatever external faces remain after removing inlet/outlet/walls/periodic
    used = set(inlet) | set(outlet) | set(inner_walls) | set(outer_walls) | periodic_all
    airfoils = [f for f in ext if f not in used]

    return {
        "inlet": inlet,
        "outlet": outlet,
        "inner_walls": inner_walls,
        "outer_walls": outer_walls,
        "airfoils": airfoils,
    }


# =============================================================================
#  Input Params (YOUR SETTINGS)
# =============================================================================
def GenerateMesh3D(**kwargs):
    mshName = kwargs.get("MeshName", "3d_swirler_mesh")
    filepath = kwargs.get("FileLocation", "")
    r_hub = kwargs.get("r_hub", 0.2)  # in
    r_pipe = kwargs.get("r_pipe", 2)    # in
    r_duct = kwargs.get("r_duct", 1.825) # in
    h_blade = r_duct - r_hub
    
    nBlades = kwargs.get("nBlades", 6)
    
    intersect_tol = 0.05  # in
    
    chord_r = kwargs.get("Chord_root", 1)  # in
    chord_t = kwargs.get("Chord_tip", 2)  # in
    
    AoA_r = kwargs.get("AoA_root", 5.0)  # deg 2.5
    AoA_t = kwargs.get("AoA_tip", 10.0)  # deg 4
    
    NACA_r = kwargs.get("NACA_root", "4406")
    NACA_t = kwargs.get("NACA_tip",  "6406")
    
    L_Upstream = kwargs.get("L_Upstream", 5)   # in, distance from inlet to LE
    L_Downstream = kwargs.get("L_Downstream", 12)  # in, distance from LE to outlet
    L_tot = L_Downstream + L_Upstream
    
    n_af_pts = kwargs.get("n_pts_airfoils", 500)
    af_mesh_size = kwargs.get("mesh_size_airfoils", 0.01)
    pts_msh_size = kwargs.get("mesh_size_general", 0.1)
    
    
    # =============================================================================
    # Geometry & Meshing
    # =============================================================================
    gmsh.initialize()
    gmsh.model.add(mshName)
    geo = gmsh.model.occ
    
    # 1) Airfoil loft (root & tip)
    rotation_angles = [0, 0, -AoA_r]
    dxdydz = [L_Upstream, 0, r_hub - intersect_tol]
    afcurve_r, af_line_tags_r, af_point_tags_r = AG.generateGMSH_NACA4(
        geo, NACA_r, dxdydz[0], dxdydz[1], dxdydz[2],
        c=chord_r, numPoints=n_af_pts, rot_ang=rotation_angles,
        mesh_size=af_mesh_size
    )
    gmsh.model.occ.synchronize()
    
    rotation_angles = [0, 0, -AoA_t]
    dxdydz = [L_Upstream, 0, r_duct + intersect_tol]
    afcurve_t, af_line_tags_t, af_point_tags_t = AG.generateGMSH_NACA4(
        geo, NACA_t, dxdydz[0], dxdydz[1], dxdydz[2],
        c=chord_t, numPoints=n_af_pts, rot_ang=rotation_angles,
        mesh_size=af_mesh_size
    )
    blade = gmsh.model.occ.addThruSections([afcurve_r, afcurve_t], makeSolid=True, makeRuled=True)[0][1]
    gmsh.model.occ.synchronize()
    
    # 2) Build wedge cross-section (in y-z) and extrude along x
    arcangle = 2.0 * math.pi / nBlades
    origin = geo.addPoint(0, 0, 0, meshSize=pts_msh_size)
    
    # Convenience
    half = arcangle / 2.0
    yc_hub  =  r_hub  * math.sin(half)   # +y for both inner endpoints
    zc_hub  =  r_hub  * math.cos(half)   # magnitude only
    
    yc_pipe =  r_duct * math.sin(half)   # +y for both outer endpoints
    zc_pipe =  r_duct * math.cos(half)   # magnitude only
    
    # Inner radius arc endpoints (z = ±, y = +)
    r1p1 = geo.addPoint(0.0, -yc_hub,  zc_hub, pts_msh_size)   # -arcangle/2  (quadrant IV)
    r1p2 = geo.addPoint(0.0, yc_hub,  zc_hub, pts_msh_size)   # +arcangle/2  (quadrant I)
    rad1 = geo.addCircleArc(startTag=r1p1, endTag=r1p2, middleTag=origin, center=True)
    
    # Outer radius arc endpoints (z = ±, y = +)
    r2p1 = geo.addPoint(0.0, -yc_pipe, zc_pipe, pts_msh_size) # -arcangle/2
    r2p2 = geo.addPoint(0.0, yc_pipe,  zc_pipe, pts_msh_size) # +arcangle/2
    rad2 = geo.addCircleArc(startTag=r2p1, endTag=r2p2, middleTag=origin, center=True)
    
    # Radial lines
    l1 = geo.addLine(r1p1, r2p1)
    l2 = geo.addLine(r1p2, r2p2)
    gmsh.model.occ.synchronize()
    
    inletLoop = geo.addCurveLoop([l1, rad2, -l2, -rad1])
    inletSurf = geo.addSurfaceFilling(inletLoop)
    ductVol_tags = geo.extrude([(2, inletSurf)], L_tot, 0, 0, recombine=False)
    gmsh.model.occ.synchronize()
    
    vol_after_extrude = ductVol_tags[1][1]
    
    # 3) Boolean: CUT (duct minus blade)
    geo.cut([(3, vol_after_extrude)], [(3, blade)])
    gmsh.model.occ.synchronize()
    
    # 4) Periodic pairing (auto)
    periodic_info = set_wedge_periodic(arcangle=arcangle, r_hub=r_hub, L_tot=L_tot)
    
    # 5) Auto-detect physical groups (faces)
    pg = auto_physical_groups(r_hub=r_hub, r_duct=r_duct, L_tot=L_tot,
                              periodic_sides=(periodic_info["sideA"], periodic_info["sideB"]))
    
    # 6) Create physical groups (you can rename as you like)
    #    - Periodic faces typically tagged as Symmetry1 / Symmetry2
    sym1 = gmsh.model.addPhysicalGroup(2, periodic_info["masters"], name="Symmetry1")
    sym2 = gmsh.model.addPhysicalGroup(2, periodic_info["slaves"],  name="Symmetry2")
    inlet_pg  = gmsh.model.addPhysicalGroup(2, pg["inlet"],        name="Inlet")
    outlet_pg = gmsh.model.addPhysicalGroup(2, pg["outlet"],       name="Outlet")
    inner_pg  = gmsh.model.addPhysicalGroup(2, pg["inner_walls"],  name="Walls_Inner")
    outer_pg  = gmsh.model.addPhysicalGroup(2, pg["outer_walls"],  name="Walls_Outer")
    foil_pg   = gmsh.model.addPhysicalGroup(2, pg["airfoils"],     name="Airfoils")
    
    # Volume(s) physical group
    vol_tags = [t for (d, t) in gmsh.model.occ.getEntities(3)]
    vol_pg = gmsh.model.addPhysicalGroup(3, vol_tags, name="FlowVolume")
    
    gmsh.model.occ.synchronize()
    
    # 7) Mesh options (optional)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.000001)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.05)
    
    # 8) Generate mesh & export
    gmsh.model.mesh.generate(3)
    
    
# =============================================================================
#     Save mesh
# =============================================================================
    msh_filepath = filepath + "\\" if filepath != "" else ""
    out_path = Path(msh_filepath + mshName + ".su2")
    # Ensure parent directory exists (no-op if '.' or already present).
    out_path.parent.mkdir(parents=True, exist_ok=True)
        
        
    fullname = msh_filepath + mshName
    gmsh.write(fullname + ".msh")
    gmsh.write(fullname + '.su2')
    # elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(dim=2, tag=periodic_info["masters"][0])
    # num_pts_sym1 = len(list(set(nodeTags[0].tolist())))
    # elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(dim=2, tag=periodic_info["slaves"][0])
    # num_pts_sym2 = len(list(set(nodeTags[0].tolist())))
    # print(f"Sym1: {num_pts_sym1}, Sym2: {num_pts_sym2}")
    
    
    # Launch GUI
    OPEN_GMSH_VISUALIZATION = kwargs.get("OPEN_GMSH_VISUALIZATION", False)
    if OPEN_GMSH_VISUALIZATION:
        gmsh.fltk.run()
    
    gmsh.finalize()
    return True
    
if __name__=="__main__":
    filepath = "C:\\Users\\BriceM\\Documents\\SU2 CFD Data\\3D_Tests"
    mshName = "refined_3D_test"
    GenerateMesh3D(MeshName=mshName, FileLocation=filepath, OPEN_GMSH_VISUALIZATION=True)
