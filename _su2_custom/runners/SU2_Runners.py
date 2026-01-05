# -*- coding: utf-8 -*-
"""
Created on Mon Jan  5 18:44:16 2026

@author: BriceM

This file is to serve as an interaction point for all previously generated files
so that it is easier to run CFD by impporting this file as a module.
"""


# =============================================================================
#  CFD Runners
# =============================================================================
def runSinglePoint_CFD():
    
    
    # From Previous File
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
    return None


def runParameterSweep_CFD():
    
    return None 

def runParameterSweep_CFD():
    
    return None 



# =============================================================================
#  Mesh Generators
# =============================================================================
def generateMesh():
    
    return None


# =============================================================================
#  Data Processors
# =============================================================================
