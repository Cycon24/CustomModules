# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 14:35:23 2025

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
import traceback
from pathlib import Path


class SwirlerMeshGenerator():
    def __init__(self, **kwargs):
        print('Initiating Swirl Mesh Generator')
        
        # Mesh Settings
        self.MeshName = kwargs.get("MeshName", "swirlerMesh")
        self.FileLocation = kwargs.get("FileLocation", "")
        self.AirfoilNACA = kwargs.get("NACA", "0012") # NACA 4 digit airfoil number
        self.FileType = kwargs.get("FileType", "BOTH") # SU2, GMSH, or BOTH currently
        self.point_mesh_size = kwargs.get("GenPointSize", 0) # Passed into point creation, increases mesh density around points
        self.af_mesh_size = kwargs.get("AF_PointSize", 0) # Will be overwritten if BL is enabled
        self.AF_numPoints = kwargs.get("AF_numPoints", 25)
        self.FullMesh = kwargs.get("FullMesh", True) # Defines wether we look at a single blade or all
        
        # Tunnel Settings
        self.numBlades = kwargs.get("numBlades", 1) # number of blades to model
        self.BladeAoA  = kwargs.get("BladeAoA", 0)
        self.Radius    = kwargs.get("Radius", 1) # Radius at which cross-section is made
        self.chord     = kwargs.get("chord", 1) # Airfoil chord length
        self.L_Upstream = kwargs.get("L_Upstream", 5*self.chord) # length of tunnel upstream of LE
        self.L_Downstream = kwargs.get("L_Downstream", 10*self.chord) # length of tunel downstream of TE
        
        # BL Settings
        self.GenBL = kwargs.get("GenBL", False)
        if self.GenBL:
            self.BL_nLayers     = kwargs.get("BL_nLayers", 20)     # number of prism layers
            self.BL_Growth      = kwargs.get("BL_Growth", 1.20)   # layer growth ratio (1.15–1.25)
            self.BL_Quad        = kwargs.get("BL_Quad", 1)       # 1 = request quads in BL
            self.BL_Factor      = kwargs.get("BL_Factor", 1.2)    # total BL thickness ~ BL_factor * laminar BL estimate
            self.BL_h1          = kwargs.get("BL_h1", 1e-5) # First BL height
            self.BL_Thickness    = kwargs.get("BL_Thickness", 1e-1)
        # TE Fan Settings
        self.GenFan_TE = kwargs.get("GenFan_TE", False)
        if self.GenFan_TE:
            self.fanLength = kwargs.get("TEF_Length")
            self.fanMinSize = kwargs.get("TEF_minSize")
             
        
        
    @staticmethod
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
        print(f"[INFO] curve1 = {len(n1):.0f}, curve2 = {len(n2):.0f}")

        # 6) Orientation hint: check monotonic direction
        dx1 = np.sign((p1s[-1,0] - p1s[0,0]) + 1e-15)
        dx2 = np.sign((p2s[-1,0] - p2s[0,0]) + 1e-15)
        if dx1 != dx2:
            print("[WARN] Curves sorted opposite along x; consider reversing one.")
            # If you need to reverse in your API version: gmsh.model.mesh.reverse([(1, sym2)])

        ok = d.max() < 1e-10
        print("[PASS]" if ok else "[FAIL]", "Periodic point match")
        return ok
    
    def GenerateMesh(self, OpenGMSHVisual=False)-> bool:
        
        # =============================================================================
        # Start GMSH Operations
        # =============================================================================
        gmsh.initialize()
        try:
            # Start logger
            log = gmsh.logger()
            log.start()
            
            gmsh.model.add(self.MeshName)
            # Make an alias to speed up typing
            geo = gmsh.model.geo 
            
            # Define the wind tunnel rectangle, the width direction will be z-plane
            # Add tunnel points
            ms = self.point_mesh_size
            
        # =============================================================================
        #     Add airfoil curves
        # =============================================================================
            af_curves = []
            all_af_lines = []
            all_af_pts = []
            
            LE_dist = 2*np.pi*self.Radius / self.numBlades 
            
            # Check if we are rendering all blades or just a single one.
            if self.FullMesh:
                for i in range(0, self.numBlades):
                    # generate airfoil curve
                    rotation_angles = [0, 0, -self.BladeAoA]
                    dxdydz = [self.L_Upstream, (LE_dist*(1/2 + i)), 0] # Add AoA adjustment on dy:  - (chord/2)*np.sin(np.deg2rad(-blade_AoA))
                    afcurve, af_line_tags, af_point_tags = AG.generateGMSH_NACA4(geo,
                                                                                 self.AirfoilNACA,
                                                                                 dxdydz[0], dxdydz[1], dxdydz[2],
                                                                                 c=self.chord, 
                                                                                 numPoints=self.AF_numPoints, 
                                                                                 rot_ang=rotation_angles, 
                                                                                 mesh_size=self.af_mesh_size)
                    # Append curve, lines, and point tags to one list each
                    af_curves.append(afcurve)
                    all_af_lines.extend(af_line_tags)
                    all_af_pts.extend(af_point_tags)
            else:
                rotation_angles = [0, 0, -self.BladeAoA]
                dxdydz = [self.L_Upstream, LE_dist/2, 0] # Add AoA adjustment on dy:  - (chord/2)*np.sin(np.deg2rad(-blade_AoA))
                afcurve, af_line_tags, af_point_tags = AG.generateGMSH_NACA4(geo,
                                                                             self.AirfoilNACA,
                                                                             dxdydz[0], dxdydz[1], dxdydz[2],
                                                                             c=self.chord, 
                                                                             numPoints=self.AF_numPoints, 
                                                                             rot_ang=rotation_angles, 
                                                                             mesh_size=self.af_mesh_size)
                # Append curve, lines, and point tags to one list each
                af_curves.append(afcurve)
                all_af_lines.extend(af_line_tags)
                all_af_pts.extend(af_point_tags)
                
        # =============================================================================
        #         BL and Fan generatiom
        # =============================================================================
            geo.synchronize()
            
            
            if self.GenBL:
                self.GenerateBoundaryLayer(geo, all_af_lines) # af_curves) # 
                
            if self.GenFan_TE:
                self.GenerateFan_TE(all_af_pts)
                
                
            
        # =============================================================================
        #     Tunnel Definition
        # =============================================================================
            # Generate space around the airfoils (Tunnel)
            y_t = 2*np.pi*self.Radius  if self.FullMesh else LE_dist # max y of tunnel
            self.PeriodicOffset = y_t
            L_tot = self.L_Upstream + self.chord + self.L_Downstream # If AoA != 0, isnt quite accurate but oh well dont care
        
            # Add points
            p1 = geo.addPoint(0,0,0, ms)
            p2 = geo.addPoint(L_tot,0,0, ms)
            p3 = geo.addPoint(L_tot, y_t,0, ms)
            p4 = geo.addPoint(0, y_t, 0, ms)
            
            # Add Inlet/Outlet lines
            outlet = geo.addLine(p2, p3)
            inlet = geo.addLine(p4, p1)
            
            # Make bottom as a true line
            sym1 = geo.addLine(p1, p2)
            sym2 = geo.addLine(p4, p3)
            # sym1 = geo.addPolyline([p1,p2])
            # sym2 = geo.addPolyline([p4,p3])
            geo.synchronize()
            
            # Apply periodicity with translation
            transform = [
                1, 0, 0, 0,
                0, 1, 0, 2*np.pi*self.Radius,
                0, 0, 1, 0,
                0, 0, 0, 1
            ]
            gmsh.model.mesh.setPeriodic(1, [sym2], [sym1], transform)
           

        # =============================================================================
        #     Add curve loop and surface
        # =============================================================================
            # Make lines into a curve loop    
            geo.synchronize()
            curve1 = geo.addCurveLoop([sym1, outlet, -sym2, inlet])
            # curve2 = geo.addCurveLoop([pline1, pline2, pline3, pline4])

            # Make it into a surface with the airfoil hole
            flowsurf = geo.addPlaneSurface([curve1, *af_curves])#, afcurve1])
            # probesurf = geo.addPlaneSurface([curve2])
            # Make Rectangular
            geo.synchronize()
            

        # =============================================================================
        #     Final defs and saves
        # =============================================================================
            gmsh.model.geo.mesh.setRecombine(2, flowsurf)
            gmsh.model.mesh.setTransfiniteSurface(flowsurf)
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
            print(" -- Periodic Check -- ")
            _ = self.periodic_curve_check(sym1, sym2, period=y_t)
            
            
            msh_filepath = self.FileLocation + "\\" if self.FileLocation != "" else ""
            out_path = Path(msh_filepath + self.MeshName + ".su2")
            # Ensure parent directory exists (no-op if '.' or already present).
            out_path.parent.mkdir(parents=True, exist_ok=True)
            
            
            
            print("Writing Mesh")
            if self.FileType == "BOTH" or self.FileType == "SU2":
                gmsh.write(msh_filepath + self.MeshName + '.su2') # Save it to disk (Note, .su2 can be specifed here)
            if self.FileType == "BOTH" or self.FileType == "GMSH":
                gmsh.write(msh_filepath + self.MeshName + '.msh')
               
            
            if OpenGMSHVisual:
                gmsh.fltk.run()
            
            success = True
                
            
        except Exception as e:
            print("\nError Occured")
            print(e)
            traceback.print_exc()
            # Add any cleanup code here
            success = False
            
        finally:
            log.stop()
            logs = log.get()
            for l in logs:
                log.write(logs)
            self.Logs = logs
            
            print("Finalizing Gmsh...")
            gmsh.finalize()
            
        return success    
        
            
    def GenerateBoundaryLayer(self, geo, BL_curves) -> None:
    # =============================================================================
    #     Set up Boundary Layer (Didnt seem to actually affect mesh)
    # STILL doesnt seem to actually affect the mesh
    # =============================================================================
        # Get the LE and TE points in all the airfoil points 
        geo.synchronize()
        print('Adding Boundary Layer...')
        
        # Add boundary layer and set options
        bl = gmsh.model.mesh.field.add("BoundaryLayer") # BoundaryLayer isnt a definitive type, just is for us to keep track of
        gmsh.model.mesh.field.setNumbers(bl,"CurvesList", BL_curves) 
        # gmsh.model.mesh.field.setNumber(bl, 'hwall_n', self.BL_h1)
        gmsh.model.mesh.field.setNumber(bl, 'BetaLaw', 1)
        gmsh.model.mesh.field.setNumber(bl, 'Beta', self.BL_Growth)
        # gmsh.model.mesh.field.setNumber(bl, "Thickness",  self.BL_Thickness)  
        gmsh.model.mesh.field.setNumber(bl, "NbLayers", self.BL_nLayers)
        # gmsh.model.mesh.field.setNumber(bl, "Ratio", self.BL_Growth)   # growth ratio between layers
        gmsh.model.mesh.field.setNumber(bl, "Quads", self.BL_Quad)   # if 1, generate recombined elements
        gmsh.model.mesh.field.setNumber(bl, "IntersectMetrics", 1)    # robust BL intersections
        gmsh.model.mesh.field.setNumber(bl, "Size", self.BL_h1)
            
        
        
        gmsh.model.mesh.field.setAsBoundaryLayer(bl)
        geo.synchronize()
        
        self.bl = bl
        
        # Recombination (may not need?)
        # gmsh.option.setNumber("Mesh.RecombineAll", 1)
        
        # Recombine flow surf?
        # gmsh.model.mesh.setRecombine(2, flowsurf)
        return None
    
    
    def GenerateFan_TE(self, all_af_pts)-> None:
        te_core_size    = max(0.25 * self.fanMinSize, 5e-4)        # very small size right at TE
        te_size_near    = max(self.fanMinSize, 1e-3)      # within DistMin
        te_size_far     = self.point_mesh_size      # beyond DistMax
        te_DistMin      = 0.1*self.fanLength        # inner radius of strong refinement
        te_DistMax      = self.fanLength #8.0 * te_size_far        # fade-out radius
        
        
        pts_xyz = {p: gmsh.model.getValue(0, p, []) for p in all_af_pts}
        x_max = max(pts_xyz[p][0] for p in all_af_pts)
        x_min= min(pts_xyz[p][0] for p in all_af_pts)
        eps_x = 1e-12
        te_points = [p for p in all_af_pts if abs(pts_xyz[p][0] - x_max) < eps_x]
        le_points = [p for p in all_af_pts if abs(pts_xyz[p][0] - x_min) < eps_x]
        
        
        # Add TE and LE pts to BL
        le_te_pts = []
        le_te_pts.extend(le_points)
        le_te_pts.extend(te_points)
        
        # Check if bl exists            
        bl = self.bl  if self.GenBL else gmsh.model.mesh.field.add("BoundaryLayer")
        
        gmsh.model.mesh.field.setNumbers(bl, "FanPointsList", le_te_pts)
        
        # Add refinement at TE
        for p in te_points:
            gmsh.model.mesh.setSize([(0, p)], te_core_size)
            
        # Add a distance field for radial refinement around TE
        dist = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(dist, "PointsList", te_points)
        gmsh.model.mesh.field.setNumber(dist,  "Sampling", 100)   # denser sampling along distance rays
    
         # Add threshold field for transition from distance field
        thr = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(thr,  "InField", dist)
        gmsh.model.mesh.field.setNumber(thr,  "SizeMin", te_size_near)  # inside DistMin
        gmsh.model.mesh.field.setNumber(thr,  "SizeMax", te_size_far)   # beyond DistMax
        gmsh.model.mesh.field.setNumber(thr,  "DistMin", te_DistMin)
        gmsh.model.mesh.field.setNumber(thr,  "DistMax", te_DistMax)
      
        
        gmsh.option.setNumber('Mesh.BoundaryLayerFanElements', self.BL_nLayers*6)
        
        # Combine BL and thr/ist fields
        minField = gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.setNumbers(minField, "FieldsList", [bl, thr])

        gmsh.model.mesh.field.setAsBackgroundMesh(minField)
        
        # Mesh.BoundaryLayerFanElements, default = 5

        
        return None
    
    @staticmethod 
    def Exit():
       sys.exit(0)


if __name__ == "__main__":
    mesh_params = {
    # Core
    "MeshName": "swirlerMeshTest",
    "FileLocation": 'TestMeshes',
    "NACA": "8412",
    "FileType": "BOTH",
    "GenPointSize": 0.5,
    "AF_PointSize": 0,
    "AF_numPoints": 100,
    "FullMesh": False,

    # Tunnel Settings
    "numBlades": 6,
    "BladeAoA": 10,
    "Radius": 1,
    "chord": 2,
    "L_Upstream": 5,      # why: default derives from 5 * chord; chord default is 1
    "L_Downstream": 10,   # why: default derives from 10 * chord; chord default is 1

    # Boundary-Layer Settings (used only if GenBL is True)
    "GenBL": True,
    "BL_nLayers": 20,
    "BL_Growth": 1.005,
    "BL_Quad": 1,
    "BL_Factor": 1.1,
    "BL_h1": 1e-4,
    "BL_Thickness": 8,

    # Trailing-Edge Fan Settings (used only if GenFan_TE is True)
    "GenFan_TE": True,
    "TEF_Length": 0.7,   # why: no default in source; explicit sentinel
    "TEF_minSize": 5e-3,
    
    # Optional run settings
    "SystemExit": True
    
    }
    
    meshGen = SwirlerMeshGenerator(**mesh_params)
    meshGen.GenerateMesh(True)
 
    
