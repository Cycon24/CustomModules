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

# =============================================================================
#  Define Mesh Options
# =============================================================================
# Swirler Options
airfoil = "0012" # Airfoil of the blades
num_Blades = 8   # Number of blades on the swirler
chord = 2        # in, chord length of airfoil 
blade_AoA = 20   # deg, angle of attack of airfoils relative to the axial direction
cross_sec_rad = 2.5 # in

# Tunnel Options
L_swirler = 2 # in
L_upstream = 5 # in 
L_downstream = 5 # in 
L_tot = L_swirler + L_downstream + L_upstream

# Mesh Options
meshName = f"swirler_n{num_Blades}_c{chord}_AoA{blade_AoA}_NACA{airfoil}"
point_mesh_size = 0.5 


# =============================================================================
# Make Cross section calculations
# =============================================================================
#  - we are unraveling a circular cut-out at radius cross_sec_rad into a 2d plane
#  - y component of things will be in the tangential direction
#  - Will not be super accurate since radial movement will not be captured
cross_circ = 2*np.pi*cross_sec_rad # Circumference of cross section
LE_dist = cross_circ / num_Blades  # 

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

    
    # afcurve1, _, _ = AG.generateGMSH_NACA4(geo, airfoil, dx=af_LE_x, dy=af_LE_y, dz=0,c=c,numPoints=25)
    # Add airfoil curves
    af_curves = []
    all_af_lines = []
    for i in range(0, num_Blades):
        # generate airfoil curve
        rotation_angles = [0, 0, blade_AoA]
        afcurve, af_line_tags, af_point_tags = AG.generateGMSH_NACA4(geo, airfoil, dx=L_upstream, dy=LE_dist*(1/2 + i), dz=0,c=chord,numPoints=25, rot_ang=rotation_angles)
        
        af_curves.append(afcurve)
        all_af_lines.extend(af_line_tags)
        
    # Generate space around the airfoils
    # Add points
    p1 = geo.addPoint(0,0,0, ms)
    p2 = geo.addPoint(L_tot,0,0, ms)
    p3 = geo.addPoint(L_tot, cross_circ,0, ms)
    p4 = geo.addPoint(0, cross_circ, 0, ms)
    
    # Add tunnel lines
    sym1 = geo.addLine(p1, p2)
    outlet = geo.addLine(p2, p3)
    sym2 = geo.addLine(p3, p4)
    inlet = geo.addLine(p4, p1)
    
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
    print("Saving Mesh")
    gmsh.model.mesh.generate(2) # Geneerate a 2D mesh
    gmsh.write(meshName + '.su2') # Save it to disk (Note, .su2 can be specifed here)
    
    # Run the GUI to visualize the created mesh
    # if '-nopopup' not in sys.argv:
    #     gmsh.fltk.run()
        
    
except Exception as e:
    print("\nError Occured")
    print(e)
    # Add any cleanup code here
    
finally:
    print("Finalizing Gmsh...")
    gmsh.finalize()
    sys.exit(0)


