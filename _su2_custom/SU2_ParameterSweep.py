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
import GMSH_2D_Swirler_Generator as genSwirler 
import SU2_Runner as su2run 
import SU2_cfg_generator as su2cfg 
import Parameter_Adjustor as PA 
import Data_Processor as DP 
import Data_Plotter_Params as DPP 
import Data_Plotter_Sweep as DPS
import myLogger as MyLog

# =============================================================================
# # Options
# =============================================================================
# Toggles
RUN_MESH_ONLY = False
RUN_SWEEP = True 
RUN_DATA_PROCESSING = True 
RANS = True

SU2_RUNS_PATH = "C:\\Users\\BriceM\\Documents\\SU2 CFD Data\\"

# if RUN_SWEEP, then we need to set what parameter to sweep
sweep_param = {"NACA": ["2412", "4412", "6412", "8412"]}
''' 
Currently Supported params:
    - Angle of Attack ("AoA")
    - Cross-section Radius ("Radius")
    - Number of Blades ("nBlades")
    - Airfoil Shape ("NACA") (valid inputs are currently only NACA 4-digit airfoils)
    - Chord length ("Chord")
    
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
sweep_FilePath = SU2_RUNS_PATH + "PARAM_RANS_Sweep_01"

# =============================================================================
# Parameter Manual Adjustments (will be applied to ALL runs)
# =============================================================================
# Import base parameters
if RANS:
    cfg_params = paramsRANS.cfg_params 
    msh_params = paramsRANS.mesh_params 
else:
    cfg_params = params.cfg_params 
    msh_params = params.mesh_params 

# Change params (if needed)
cfg_params["ITER"] = 7500


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
            run_FilePath = param_FilePath + "\\" + key + f"_{val}"
            mshName = "Mesh_" + key + f"_{val}" 
            cfgName = "CFG_"+ key + f"_{val}.cfg" 
            
            # Set naming params
            sweep_msh_params["FileLocation"] = run_FilePath
            sweep_msh_params["FileType"] = "BOTH"
            sweep_msh_params["MeshName"] = mshName 
            sweep_cfg_params["MESH_FILENAME"] = mshName + ".su2"
            
            # Make the adjustment
            sweep_cfg_params, sweep_msh_params = PA.AdjustParams(key, val, sweep_cfg_params, sweep_msh_params)
            
            # Generate CFG Mesh file (need to use values calculated in mesh for cfg inputs)
            swirler = genSwirler.SwirlerMeshGenerator(**sweep_msh_params)
            SWIRLER_GEN_COMPLETE = swirler.GenerateMesh(OpenGMSHVisual= RUN_MESH_ONLY)
            su2cfg.dict_to_cfg(sweep_msh_params, mshName+"_params.txt", run_FilePath)
            
            # Save cfg file
            sweep_cfg_params["MARKER_PERIODIC"] = "(Symmetry2, Symmetry1) , (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, -{}, 0.0)".format(swirler.PeriodicOffset)
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
 