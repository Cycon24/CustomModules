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
import GMSH_2D_Swirler_Generator as genSwirler 
import SU2_Runner as su2run 
import SU2_cfg_generator as su2cfg 


# Set General Params
sweep_FilePath = "sweepTest003"
mshName = "testMesh" 
cfgName = "testCFG.cfg"

# Import base parameters
cfg_params = params.cfg_params 
msh_params = params.mesh_params 

# Change params (if needed)
msh_params["FileLocation"] = sweep_FilePath
msh_params["FileType"] = "BOTH"
msh_params["MeshName"] = mshName 
cfg_params["MESH_FILENAME"] = mshName + ".su2"
cfg_params["ITER"] = 5000

# Generate CFG Mesh file (need to use values calculated in mesh for cfg inputs)
swirler = genSwirler.SwirlerMeshGenerator(**msh_params)
SWIRLER_GEN_COMPLETE = swirler.GenerateMesh(True)
su2cfg.dict_to_cfg(msh_params, mshName+".txt", sweep_FilePath)

# Save cfg file
cfg_params["MARKER_PERIODIC"] = f"(Symmetry2, Symmetry1) , (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, -{swirler.PeriodicOffset}, 0.0)"
CFG_GEN_COMPLETE = su2cfg.dict_to_cfg(cfg_params, cfgName, sweep_FilePath)


# Run CFD
if SWIRLER_GEN_COMPLETE and CFG_GEN_COMPLETE:
    print("# =============================================================================")
    print("#   BEGINNING SU2 RUN")
    print("# =============================================================================")
    
    su2run.basic_CFD_run(cfgName, sweep_FilePath, save_output=True)

