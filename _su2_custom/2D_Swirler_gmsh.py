# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 21:09:31 2025

@author: BriceM
"""

import sys
import os
# Add parent import capabilities
parentdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if parentdir not in sys.path:
    sys.path.insert(0, parentdir)
import _aerodynamics.AirfoilGenerator as AG
import gmsh 
import numpy as np
import SU2_Mesh_Orientation_Fix as MOF
import traceback

# =============================================================================
#  Define Mesh Options
# =============================================================================
# Swirler Options
airfoil = "8412" # Airfoil of the blades
num_Blades = 6   # Number of blades on the swirler
chord = 2        # in, chord length of airfoil 
blade_AoA = 10   # deg, angle of attack of airfoils relative to the axial direction
cross_sec_rad = 1 # in

# Tunnel Options
L_swirler = 2 # in
L_upstream = 8 # in 
L_downstream = 12 # in 
L_tot = L_swirler + L_downstream + L_upstream

# Mesh Options
meshName = f"swirler_n{num_Blades}_c{chord}_AoA{blade_AoA}_NACA{airfoil}"
point_mesh_size = 0.1


# =============================================================================
# Make Cross section calculations
# =============================================================================
#  - we are unraveling a circular cut-out at radius cross_sec_rad into a 2d plane
#  - y component of things will be in the tangential direction
#  - Will not be super accurate since radial movement will not be captured
cross_circ = 2*np.pi*cross_sec_rad # Circumference of cross section
LE_dist = cross_circ / num_Blades  # 


# =============================================================================
# Triangle Fizer (brokeN, from Chat GPT)
# =============================================================================
def fix_triangle_orientation():
    elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(2)
    elemTypes = list(elemTypes)  # convert to Python list

    if 2 not in elemTypes:
        return 0

    idx = elemTypes.index(2)
    tris = np.array(nodeTags[idx]).reshape(-1, 3)
    flipped = 0

    node_ids, node_coords, _ = gmsh.model.mesh.getNodes()
    coords = {tag: node_coords[i*3:(i+1)*3] for i, tag in enumerate(node_ids)}

    for i, tri in enumerate(tris):
        p1 = np.array(coords[tri[0]][:2])
        p2 = np.array(coords[tri[1]][:2])
        p3 = np.array(coords[tri[2]][:2])

        area = 0.5 * ((p2[0] - p1[0]) * (p3[1] - p1[1]) -
                      (p2[1] - p1[1]) * (p3[0] - p1[0]))

        if area < 0:
            # flip orientation
            tris[i][1], tris[i][2] = tris[i][2], tris[i][1]
            flipped += 1

    gmsh.model.mesh.setElements(
        2, [2], [elemTags[idx]], [tris.flatten().tolist()]
    )
    return flipped


def periodic_curve_check(sym1, sym2, period):
    """
    sym1: tag of bottom curve
    sym2: tag of top curve
    period: translation distance in y (your 2*pi*radius)
    """
    # 1) Get node ids and coords on each curve (dim=1)
    n1, c1, _ = gmsh.model.mesh.getNodes(dim=1, tag=sym1, includeBoundary=True)
    n2, c2, _ = gmsh.model.mesh.getNodes(dim=1, tag=sym2, includeBoundary=True)

    # 2) Basic count check
    if len(n1) != len(n2):
        print(f"[FAIL] Node count mismatch: {len(n1)} vs {len(n2)}")
        return False

    # 3) Build coordinate arrays
    p1 = np.array(c1).reshape(-1, 3)  # (N,3)
    p2 = np.array(c2).reshape(-1, 3)

    # 4) Sort nodes along x to ensure consistent ordering (robust against tag scrambling)
    i1 = np.lexsort((p1[:,1], p1[:,0]))   # sort by x then y
    i2 = np.lexsort((p2[:,1], p2[:,0]))

    p1s = p1[i1]
    p2s = p2[i2]

    # 5) Apply periodic translation to sym1 nodes and compare to sym2
    p1s_shift = p1s.copy()
    p1s_shift[:,1] += period  # y-translation

    d = np.linalg.norm(p1s_shift - p2s, axis=1)
    print(f"[INFO] periodic max |Δ| = {d.max():.3e}, mean |Δ| = {d.mean():.3e}")

    # 6) Orientation hint: check monotonic direction
    dx1 = np.sign((p1s[-1,0] - p1s[0,0]) + 1e-15)
    dx2 = np.sign((p2s[-1,0] - p2s[0,0]) + 1e-15)
    if dx1 != dx2:
        print("[WARN] Curves sorted opposite along x; consider reversing one.")
        # If you need to reverse in your API version: gmsh.model.mesh.reverse([(1, sym2)])

    ok = d.max() < 1e-10
    print("[PASS]" if ok else "[FAIL]", "Periodic point match")
    return ok

# =============================================================================
# Start GMSH Operations
# =============================================================================
gmsh.initialize()



try:
    gmsh.model.add(meshName)
    
    
    geo = gmsh.model.geo 
    
    # Define the wind tunnel rectangle, the width direction will be z-plane
    # Add tunnel points
    ms = point_mesh_size

    # Add airfoil curves
    af_curves = []
    all_af_lines = []
    for i in range(0, num_Blades):
        # generate airfoil curve
        rotation_angles = [0, 0, -blade_AoA]
        afcurve, af_line_tags, af_point_tags = AG.generateGMSH_NACA4(geo, airfoil, dx=L_upstream, dy=LE_dist*(1/2 + i), dz=0,c=chord,numPoints=40, rot_ang=rotation_angles)
        
        af_curves.append(afcurve)
        all_af_lines.extend(af_line_tags)
        
    # Generate space around the airfoils
    # Add points
    p1 = geo.addPoint(0,0,0, ms)
    p2 = geo.addPoint(L_tot,0,0, ms)
    p3 = geo.addPoint(L_tot, cross_circ,0, ms)
    p4 = geo.addPoint(0, cross_circ, 0, ms)
    
    # Add Inlet/Outlet lines
    outlet = geo.addLine(p2, p3)
    inlet = geo.addLine(p4, p1)
    
    # Make bottom as a true line
    sym1 = geo.addLine(p1, p2)
    sym2 = geo.addLine(p3, p4)
    geo.synchronize()
    
    # Make top as an exact copy, then translate it by the period
    # copied = gmsh.model.geo.copy([(1, sym1)])     # returns [(1, newTag)]
    # sym2 = copied[0][1]
    # gmsh.model.geo.translate([(1, sym2)], 0.0, cross_circ, 0.0)
    # geo.synchronize()
    
  
    
# # =============================================================================
# #     From GPT
# # =============================================================================
#     # Force both symmetry curves to have the same discretization
#     nDiv = 40
#     gmsh.model.mesh.setTransfiniteCurve(sym1, nDiv)
#     gmsh.model.mesh.setTransfiniteCurve(sym2, nDiv)
    
#     # # Get endpoints of both curves
#     # bnd_sym1 = gmsh.model.getBoundary([(1, sym1)], oriented=False)
#     # bnd_sym2 = gmsh.model.getBoundary([(1, sym2)], oriented=False)
    
#     # p1_tag = bnd_sym1[0][1]
#     # p2_tag = bnd_sym2[0][1]
    
#     # x1, y1, z1 = gmsh.model.getValue(0, p1_tag, [])
#     # x2, y2, z2 = gmsh.model.getValue(0, p2_tag, [])
    
#     # # Reverse sym2 if its orientation doesn't match sym1
#     # if abs(x1 - x2) > 1e-12 or abs(y1 - y2) > 1e-12:
#     #     gmsh.model.mesh.reverse([(1, sym2)])
    
#     # gmsh.model.geo.synchronize()
    
#     # Apply periodicity with translation
#     transform = [
#         1, 0, 0, 0,
#         0, 1, 0, cross_circ,
#         0, 0, 1, 0,
#         0, 0, 0, 1
#     ]
#     gmsh.model.mesh.setPeriodic(1, [sym2], [sym1], transform)

# # =============================================================================
# #     
# # =============================================================================
#     geo.synchronize()
    # Make lines into a curve loop
    curve1 = geo.addCurveLoop([sym1, outlet, sym2, inlet])

    
    # Make it into a surface with the airfoil hole
    flowsurf = geo.addPlaneSurface([curve1, *af_curves])#, afcurve1])
    
    #Sync
    geo.synchronize()
    
    # Make physical groups
    gmsh.model.addPhysicalGroup(dim=1, tags=all_af_lines, name="Airfoils")
    gmsh.model.addPhysicalGroup(dim=1, tags=[inlet], name="Inlet")
    gmsh.model.addPhysicalGroup(dim=1, tags=[sym1], name="Symmetry1")
    gmsh.model.addPhysicalGroup(dim=1, tags=[sym2], name="Symmetry2")
    gmsh.model.addPhysicalGroup(dim=1, tags=[outlet], name="Outlet")
    gmsh.model.addPhysicalGroup(dim=2, tags=[flowsurf], name="FlowSurface")
    
    # Synchronize the model
    geo.synchronize()
    
    # Finalizations
    print("Generating Mesh")
    gmsh.model.mesh.generate(2) # Geneerate a 2D mesh
    
    # Check sym
    _ = periodic_curve_check(sym1, sym2, period=2*np.pi*cross_sec_rad)
    
    # Fix mesh (broken)
    # flipped_count = fix_triangle_orientation()
    # print(f"Corrected {flipped_count} flipped triangles before writing SU2.")
    
    print("Writing Mesh")
    gmsh.write(meshName + '.su2') # Save it to disk (Note, .su2 can be specifed here)
    gmsh.write(meshName + '.msh')
    # print("Checking Mesh")
    # infile = meshName + '.su2'
    # outfile = meshName + '.su2'
    # lines = MOF.read_su2_mesh(infile)
    # fixed_lines = MOF.fix_tri_orientation(lines)
    # MOF.write_su2_mesh(outfile, fixed_lines)
    # print(f"Fixed mesh written to {outfile}")
    
    
   # Run the GUI to visualize the created mesh
    # if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
        
    
except Exception as e:
    print("\nError Occured")
    print(e)
    traceback.print_exc()
    # Add any cleanup code here
    
finally:
    print("Finalizing Gmsh...")
    gmsh.finalize()
    sys.exit(0)


