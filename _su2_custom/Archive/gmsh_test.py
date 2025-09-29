# -*- coding: utf-8 -*-
"""
Created on Sat Sep  6 15:37:13 2025

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



# Goal of this it to make a wind tunnel with an airfoil across the seciton
# Tunnel dimensions
L = 10 # m x-dir
W = 3 # m z-dir
H = 3 # M y-dir

# Set chord length
c = 1 

# Airfoil locations
af_LE_x = L/2 - c/2
af_LE_y = H/2

# Set airfoil
NACA = "8412"
meshName = "tunnel_NACA" + NACA

# NACApoints = AG.generateNACA4(NACA, c, 25)

# =============================================================================
# Mesh making
# =============================================================================
gmsh.initialize()

try:
    gmsh.model.add(meshName)
    
    
    geo = gmsh.model.geo 
    
    # Define the wind tunnel rectangle, the width direction will be z-plane
    # Add tunnel points
    wt_s = 0.5
    
    wtp1 = geo.addPoint(0,0,0, wt_s)
    wtp2 = geo.addPoint(L,0,0, wt_s)
    wtp3 = geo.addPoint(L,W,0, wt_s)
    wtp4 = geo.addPoint(0,W,0, wt_s)
    
    # Add tunnel lines
    wtl1 = geo.addLine(wtp1, wtp2)
    wtl2 = geo.addLine(wtp2, wtp3)
    wtl3 = geo.addLine(wtp3, wtp4)
    wtl4 = geo.addLine(wtp4, wtp1)
    
    # Make lines into a curve loop
    wtcurve1 = geo.addCurveLoop([wtl1, wtl2, wtl3, wtl4])
    
    afcurve1, _, _ = AG.generateGMSH_NACA4(geo, NACA, dx=af_LE_x, dy=af_LE_y, dz=0,c=c,numPoints=25)
    
    # Make it into a surface with the airfoil hole
    wtwall1 = geo.addPlaneSurface([wtcurve1, afcurve1])
    
    # Extrude it into the rest of the section
    vol1 = geo.extrude([(2,wtwall1)], dx=0,dy=0,dz=W)
    
    # Synchronize the model
    geo.synchronize()
    
    # Finalizations
    gmsh.model.mesh.generate(3) # Geneerate a 2D mesh
    gmsh.write(meshName + '.msh') # Save it to disk (Note, .su2 can be specifed here)
    
    # Run the GUI to visualize the created mesh
    # if '-nopopup' not in sys.argv:
    #     gmsh.fltk.run()
        
    
except KeyboardInterrupt:
    print("\nKeyboard interrupt received. Stopping gracefully...")
    # Add any cleanup code here
    
finally:
    print("Finalizing Gmsh...")
    gmsh.finalize()
    sys.exit(0)
