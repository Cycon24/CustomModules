# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 14:23:40 2025

@author: BriceM

This code will be interacting with the previously refined mesh generation code,
SU2 runner, and the configuration files. In order to do this, we need interfaces
for each that includes the parameters we will use as inputs. We want to develop a 
consistent naming system for results and file management, as we will be generating
a large number of files for each sweep. Additionally, we want to maintain a log
file that tracks the overall sweep progression. Below we will define the parameters
necessary for each code.

General Flow (used by mesh and config):
- Mi
- Pi
- Ti
- Gas:
    - gamma
    - R
    - mu 
- Other calculated flow vars (total properties, vel, Re, etc)

Mesh
- Pipe L Upstream
- Pipe L Downstream
- Crossection Radius
- Number of Blades
- Chord Length
- Blade shape
- Blade AoA
- Boundary Layer Parameters
- General Flow vars

Config File
- Iterations
- General Flow vars
- File names
- File output names

SU2 Runner
- File names
- File locations/folders

"""

import base_params as params
import base_params_RANS as paramsRANS
import base_params_RANS_3D as paramsRANS_3D

import GMSH_2D_Swirler_Generator as genSwirler 
import GMSH_3D_Swirler_Generator as gen3DSwirler

import SU2_Runner as su2run 
import SU2_cfg_generator as su2cfg 
import Parameter_Adjustor as PA 

import Data_Processor as DP 
import Data_Plotter_Params as DPP 
import Data_Plotter_Sweep as DPS

from meshCopier import copy_mesh

# =============================================================================
# # Options
# =============================================================================
# Toggles
RUN_MESH_ONLY = False
RUN_SWEEP = True 
RUN_DATA_PROCESSING = True 
RANS = True
RANS_3D = True

SU2_RUNS_PATH = "C:\\Users\\BriceM\\Documents\\SU2 CFD Data\\"

# if RUN_SWEEP, then we need to set what parameter to sweep
P_atm = 2116.216 # lbf/ft^2
sweep_param = {"Pb": [0.99*P_atm, 0.95*P_atm]}
''' 
Currently Supported params:
    - Angle of Attack ("AoA")
    - Cross-section Radius ("Radius")
    - Number of Blades ("nBlades")
    - Airfoil Shape ("NACA") (valid inputs are currently only NACA 4-digit airfoils)
    - Chord length ("Chord")
    3D
    - Back Pressure ("Pb") in lbf/ft^2 (1 atm is 2116.216)
    
Whenever a new one is made need to update:
    - Parameter_Adjustor.py case structure
    - label_subs.yaml file
'''
# Set General Params (Single)
single_Folder = SU2_RUNS_PATH + "RANS_Tests_Single"
single_FilePath = "RANS_Test_02"
mshName = "RansMSH" 
cfgName = "Rans_test.cfg"

# Set General Params (Sweep) 
sweep_FilePath = SU2_RUNS_PATH + "3D_Tests"
USE_DEFAULT_MESH = True
default_mesh_location = SU2_RUNS_PATH + "3D_tests"
default_mesh_name = "refined_3D_test"

# =============================================================================
# Parameter Manual Adjustments (will be applied to ALL runs)
# =============================================================================
# Import base parameters
if RANS_3D:
    cfg_params = paramsRANS_3D.cfg_params 
    msh_params = paramsRANS_3D.mesh_params 
elif RANS:
    cfg_params = paramsRANS.cfg_params 
    msh_params = paramsRANS.mesh_params 
else:
    cfg_params = params.cfg_params 
    msh_params = params.mesh_params 

# Change params (if needed)
# cfg_params["ITER"] = 7500


# =============================================================================
# Running code
# =============================================================================
if not RUN_SWEEP:
    filepath = single_Folder + "\\" +  single_FilePath
    msh_params["FileLocation"] = filepath
    msh_params["FileType"] = "BOTH"
    msh_params["MeshName"] = mshName 
    cfg_params["MESH_FILENAME"] = mshName + ".su2"
    
    # Generate CFG Mesh file (need to use values calculated in mesh for cfg inputs
    
    swirler = genSwirler.SwirlerMeshGenerator(**msh_params)
    SWIRLER_GEN_COMPLETE = swirler.GenerateMesh(OpenGMSHVisual= RUN_MESH_ONLY)
    
    su2cfg.dict_to_cfg(msh_params, mshName+"_params.txt", filepath)
    
    # Save cfg file
    cfg_params["MARKER_PERIODIC"] = cfg_params["MARKER_PERIODIC"].format(swirler.PeriodicOffset)
    CFG_GEN_COMPLETE = su2cfg.dict_to_cfg(cfg_params, cfgName, filepath)
    
    
    # Run CFD
    if SWIRLER_GEN_COMPLETE and CFG_GEN_COMPLETE and not RUN_MESH_ONLY:
        print("# =============================================================================")
        print("#   BEGINNING SU2 RUN")
        print("# =============================================================================")
        
        su2run.basic_CFD_run(cfgName, filepath, save_output=True)
        
    if RUN_DATA_PROCESSING and not RUN_MESH_ONLY: 
        DP.main(single_Folder, tol_in=0.1)
        DPP.main(single_Folder, "png", dpi=200)
      
        

else:
    for key in sweep_param: 
        # sweep through each parameter, each parameter has its own sweep
        param_FilePath = sweep_FilePath + "\\" + key +"_sweep" 
        
        
        # a few things to keep in mind. We want to initial params to be kept in
        # between sweeps, aka be able to revert to the previous value. We can do this
        # by either saving it in a variable or by making the new dicts as copies of the old.
        # Storing a single variable seems the most straightforward, however, for more impactful 
        # parameters where multiple values are changing this method becomes detrimental. 
        # Therefore, we should make copies of the dicts rather than alter the base.
        # We only need to copy them prior to each sweep parameter to "reset" from 
        # a previous run
        sweep_cfg_params = cfg_params.copy()
        sweep_msh_params = msh_params.copy()
        
        
        for val in sweep_param[key]:
            # Sweep through all values in the sweep. Each will get its own file as well
            match val:
                case int():
                    endStr = f"_{val}"
                case float():
                    endStr = f"_{val:.0f}"
                case str():
                    endStr = "_" + val 
                case _:
                    endStr = f"_{val}"
                    
            run_FilePath = param_FilePath + "\\" + key + endStr
            mshName = "Mesh_" + key + endStr 
            cfgName = "CFG_"+ key + endStr + ".cfg" 
            
            # Set naming params
            sweep_msh_params["FileLocation"] = run_FilePath
            sweep_msh_params["FileType"] = "BOTH"
            sweep_msh_params["MeshName"] = mshName 
            sweep_cfg_params["MESH_FILENAME"] = mshName + ".su2"
            
            # Make the adjustment
            sweep_cfg_params, sweep_msh_params = PA.AdjustParams(key, val, sweep_cfg_params, sweep_msh_params)
            
            if not RANS_3D:
                # Generate CFG Mesh file (need to use values calculated in mesh for cfg inputs)
                swirler = genSwirler.SwirlerMeshGenerator(**sweep_msh_params)
                SWIRLER_GEN_COMPLETE = swirler.GenerateMesh(OpenGMSHVisual= RUN_MESH_ONLY)
                su2cfg.dict_to_cfg(sweep_msh_params, mshName+"_params.txt", run_FilePath)
                
                # Save cfg file
                sweep_cfg_params["MARKER_PERIODIC"] = "(Symmetry2, Symmetry1) , (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, -{}, 0.0)".format(swirler.PeriodicOffset)
                CFG_GEN_COMPLETE = su2cfg.dict_to_cfg(sweep_cfg_params, cfgName, run_FilePath)
            else:
                # gen 3d mesh
                if USE_DEFAULT_MESH:
                    copy_mesh(default_mesh_location, default_mesh_name, run_FilePath, mshName)
                    SWIRLER_GEN_COMPLETE = True    
                else:
                    SWIRLER_GEN_COMPLETE = gen3DSwirler.GenerateMesh3D(**sweep_msh_params)
                su2cfg.dict_to_cfg(sweep_msh_params, mshName+"_params.txt", run_FilePath) 
                
                # save cfg file
                CFG_GEN_COMPLETE = su2cfg.dict_to_cfg(sweep_cfg_params, cfgName, run_FilePath)
                
            # Run CFD
            if SWIRLER_GEN_COMPLETE and CFG_GEN_COMPLETE and not RUN_MESH_ONLY:
                print("# =============================================================================")
                print("#   BEGINNING SU2 RUN")
                print("# =============================================================================")
                
                su2run.basic_CFD_run(cfgName, run_FilePath, save_output=True)
    
    
        if RUN_DATA_PROCESSING and not RUN_MESH_ONLY: 
            DP.main(param_FilePath, tol_in=0.1)
            DPP.main(param_FilePath, "png", dpi=200)
            DPS.main(param_FilePath, "png", dpi=200)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
 