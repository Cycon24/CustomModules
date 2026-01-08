# -*- coding: utf-8 -*-
"""
Created on Mon Jan  5 18:44:16 2026

@author: BriceM

This file is to serve as an interaction point for all previously generated files
so that it is easier to run CFD by impporting this file as a module.
"""
# Imports (temporary until file structure remedied)
import SU2_cfg_generator as su2cfg 
import SU2_base_runner as su2run

import tools.Parameter_Adjustor as PA 
from   tools.Mesh_Copier import copy_mesh

import GMSH_2D_Swirler_Generator as genSwirler 
import GMSH_3D_Swirler_Generator as gen3DSwirler

import Data_Processor as DP 
import Data_Plotter_Params as DPP 
import Data_Plotter_Sweep as DPS

import postprocessing3d.Data3D_Processor as DP3
import postprocessing3d.Data3D_Plotter_Params as DPP3
import postprocessing3d.Data3D_Plotter_Sweep as DPS3

# =============================================================================
#  CFD Runners
# =============================================================================
def runSinglePoint_CFD(cfg_params:dict, msh_params:dict, **kwargs):
    # STILL NEED A 2D/3D Toggle? Or it can just be supplied in config and mesh param files?
    # Toggles
    USE_DEFAULT_MESH = kwargs.get("USE_DEFAULT_MESH", False)
    
    
    # Values
    filepath = kwargs.get("FilePath", None) # Folder location and folder where config file and all data are stored, does not end in "\\"
    mshName = kwargs.get("MeshName", "Mesh")
    cfgName = kwargs.get("ConfigName","Config")  # Dont need these to be inputted controls?
    default_mesh_location = kwargs.get("default_mesh_location", None)
    default_mesh_name = kwargs.get("default_mesh_name", None)
    
    # update mesh parameters dictionary with inputted information (need here?)
    msh_params["FileLocation"] = filepath
    msh_params["MeshName"] = mshName 
    
    SWIRLER_GEN_COMPLETE = False
    if USE_DEFAULT_MESH:
        # Check if location and name provided:
        if (default_mesh_location != None) and (default_mesh_name != None):
            copy_mesh(default_mesh_location, default_mesh_name, filepath, mshName)
            SWIRLER_GEN_COMPLETE = True  
            periodicOffset = su2cfg.txt_to_dict(default_mesh_name+"_params.txt", default_mesh_location)["PERIODIC_OFFSET"]
        else:
            print('[Warning]\t Default Mesh not found, generating mesh from parameters')
    
    if not SWIRLER_GEN_COMPLETE:
        # Generate Mesh (need to use values calculated in mesh for cfg inputs)
        # Separate this into its own function
        Generated, periodicOffset = generateMesh(msh_params, OpenGMSHVisual=False)        # Mesh generation
        # periodicOffset = swirler.PeriodicOffset                         # Pull periodic offset out to update config file
        # SWIRLER_GEN_COMPLETE = swirler.GENERATED
    
    # Update cfg and mesh parameters
    if periodicOffset != None:
        cfg_params["MARKER_PERIODIC"] = cfg_params["MARKER_PERIODIC"].format(periodicOffset)
    cfg_params["MESH_FILENAME"] = mshName + ".su2"
    msh_params["PERIODIC_OFFSET"] = periodicOffset
    
    # Save cfg and mesh params to file
    CFGtxt_GEN_COMPLETE = su2cfg.dict_to_cfg(cfg_params, cfgName, filepath)
    MSHtxt_GEN_COMPLETE = su2cfg.dict_to_cfg(msh_params, mshName+"_params.txt", filepath)
    
    # Run CFD
    if (SWIRLER_GEN_COMPLETE and MSHtxt_GEN_COMPLETE and CFGtxt_GEN_COMPLETE):
        print("# =============================================================================")
        print("#   BEGINNING SU2 RUN")
        print("# =============================================================================")
        
        su2run.base_CFD_run(cfgName, filepath, save_output=True)
    else:
        print("[Error]\t CFD not ran, previous task failed:")
        print(f"\t SWIRLER_GEN_COMPLETE = {SWIRLER_GEN_COMPLETE}")
        print(f"\t MSHtxt_GEN_COMPLETE  = {MSHtxt_GEN_COMPLETE}")
        print(f"\t CFGtxt_GEN_COMPLETE  = {CFGtxt_GEN_COMPLETE}")
       
    return None


def runParameterSweep_CFD(sweep_params:dict, cfg_params:dict, msh_params:dict, **kwargs):
    sweep_FilePath = kwargs.get("FilePath", "")
    
    # Iterate through each sweep prameter
    for key in sweep_params: 
        # sweep through each parameter, each parameter has its own sweep
        param_FilePath = sweep_FilePath + "\\" + key +"_sweep" 
        
        # Need to copy them prior to each sweep parameter to "reset" from 
        # a previous run and keep the initial parameter values
        base_cfg_params = cfg_params.copy()
        base_msh_params = msh_params.copy()
        
        # Iterate through each value within a particular sweep parameter
        for val in sweep_params[key]:
            # Go through all values in the sweep. Each will get its own file as well
            match val:
                case int():
                    endStr = f"_{val}"
                case float():
                    endStr = f"_{val:.0f}"
                case str():
                    endStr = "_" + val 
                case list():
                    endStr = ""
                    for v in val:
                        endStr += f"_{v:.0f}"
                case _:
                    endStr = f"_{val}"
                    
            run_FilePath = param_FilePath + "\\" + key + endStr
            mshName = "Mesh_" + key + endStr 
            cfgName = "CFG_"+ key + endStr + ".cfg" 
            
            # Set naming params
            base_msh_params["FileLocation"] = run_FilePath
            base_msh_params["MeshName"] = mshName 
            base_cfg_params["MESH_FILENAME"] = mshName + ".su2"
            
            # Make the parameter adjustment
            sweep_cfg_params, sweep_msh_params = PA.AdjustParams(key, val, base_cfg_params, base_msh_params)
            
            # Run the CFD by passing kwargs into the single run, mesh generation handled within single run
            runKwargs = kwargs.copy()
            runKwargs["MeshName"] = mshName
            runKwargs["ConfigName"] = cfgName
            runKwargs["FilePath"] = run_FilePath
            runSinglePoint_CFD(sweep_cfg_params, sweep_msh_params, **runKwargs)
            
            
            
            # if not RANS_3D:
            #     # Generate CFG Mesh file (need to use values calculated in mesh for cfg inputs)
            #     swirler = genSwirler.SwirlerMeshGenerator(**sweep_msh_params)
            #     SWIRLER_GEN_COMPLETE = swirler.GenerateMesh(OpenGMSHVisual= RUN_MESH_ONLY)
            #     su2cfg.dict_to_cfg(sweep_msh_params, mshName+"_params.txt", run_FilePath)
                
            #     # Save cfg file
            #     sweep_cfg_params["MARKER_PERIODIC"] = "(Symmetry2, Symmetry1) , (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, -{}, 0.0)".format(swirler.PeriodicOffset)
            #     CFG_GEN_COMPLETE = su2cfg.dict_to_cfg(sweep_cfg_params, cfgName, run_FilePath)
            # else:
            #     # gen 3d mesh
            #     if USE_DEFAULT_MESH:
            #         copy_mesh(default_mesh_location, default_mesh_name, run_FilePath, mshName)
            #         SWIRLER_GEN_COMPLETE = True    
            #     else:
            #         SWIRLER_GEN_COMPLETE = gen3DSwirler.GenerateMesh3D(**sweep_msh_params)
            #     su2cfg.dict_to_cfg(sweep_msh_params, mshName+"_params.txt", run_FilePath) 
                
            #     # save cfg file
            #     CFG_GEN_COMPLETE = su2cfg.dict_to_cfg(sweep_cfg_params, cfgName, run_FilePath)
                
            
    
    return None 





# =============================================================================
#  Mesh Generators
# =============================================================================
def generateMesh(mesh_parameters:dict, OpenGMSHVisual=True):
    IS_3D = mesh_parameters["IS_3D"]
    
    if IS_3D:
        GENERATED = gen3DSwirler.GenerateMesh3D(OpenGMSHVisual=OpenGMSHVisual,**mesh_parameters)
        periodicOffset = None
    else:
        swirlerObj = genSwirler.SwirlerMeshGenerator(**mesh_parameters)
        GENERATED = swirlerObj.GenerateMesh(OpenGMSHVisual)
        periodicOffset = swirlerObj.PeriodicOffset
        
    return GENERATED, periodicOffset


# =============================================================================
#  Data Processors
# =============================================================================
def processSweepData(filepath:str, IS_3D:bool=False, 
                     cfg_path:str   = r'C:\Users\BriceM\Documents\Modules\_su2_custom\3dpostprocessing\post3d_config.yaml',
                     units_path:str = r'C:\Users\BriceM\Documents\Modules\_su2_custom\3dpostprocessing\plot_units.yaml',
                     dpi:int=250,
                     tol_in = 0.1):
    
    if not IS_3D:
        DP.main(filepath, tol_in=tol_in)
        DPP.main(filepath, "png", dpi=dpi)
        DPS.main(filepath, "png", dpi=dpi)
    else:
        # 3-D pipeline
        try:
            DP3.main(filepath, cfg_path, units_path=units_path)
            DPP3.main(filepath, cfg_path, units_path=units_path, img_format="png", dpi=dpi)
            DPS3.main(filepath, cfg_path, units_path=units_path, img_format="png", dpi=dpi)
        except:
            print("[Error]\t Error in 3d processing")