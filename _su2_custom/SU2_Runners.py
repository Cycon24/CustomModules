# -*- coding: utf-8 -*-
"""
Created on Mon Jan  5 18:44:16 2026

@author: BriceM

This file is to serve as an interaction point for all previously generated files
so that it is easier to run CFD by impporting this file as a module.
"""
# Imports (temporary until file structure remedied)


# import Data_Processor as DP 
# import Data_Plotter_Params as DPP 
# import Data_Plotter_Sweep as DPS

# import postprocessing3d.Data3D_Processor as DP3
# import postprocessing3d.Data3D_Plotter_Params as DPP3
# import postprocessing3d.Data3D_Plotter_Sweep as DPS3

# Solidified
import yaml 
from pathlib import Path
import numpy as np 

import runners.SU2_base_runner as su2run
# from postprocessing.Base_Plotter import process_single_run

import tools.Parameter_Adjustor as PA 
import tools.File_Manipulator as su2cfg 
import tools.Flow_Calculator as flowcalc

from   meshgeneration.Mesh_Copier import copy_mesh
import meshgeneration.GMSH_2D_Swirler_Generator as genSwirler 
import meshgeneration.GMSH_3D_Swirler_Generator as gen3DSwirler

# =============================================================================
#  CFD Runners
# =============================================================================
def runSinglePoint_CFD(cfg_params:dict, msh_params:dict, filepath:Path, flow_params:dict|None = None, **kwargs):
    # STILL NEED A 2D/3D Toggle? Or it can just be supplied in config and mesh param files?
    # Toggles
    USE_DEFAULT_MESH = kwargs.get("USE_DEFAULT_MESH", False)
    
    # Values
    # Folder location and folder data should be stored, does not end in "\\"
    mshName = kwargs.get("MeshName", "Mesh")
    cfgName = kwargs.get("ConfigName","Config.cfg")  # Dont need these to be inputted controls?
    default_mesh_location = kwargs.get("default_mesh_location", None)
    default_mesh_name = kwargs.get("default_mesh_name", None)
    
    # update mesh parameters dictionary with inputted information (need here?)
    msh_params["FileLocation"] = str(filepath)
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
        # If mesh wasnt copied (default mesh), generate mesh with the provided parameters
        SWIRLER_GEN_COMPLETE, periodicOffset = generateMesh(msh_params, OpenGMSHVisual=False)
    
    # Update cfg and mesh parameters
    if periodicOffset != None:
        cfg_params["MARKER_PERIODIC"] = cfg_params["MARKER_PERIODIC"].format(periodicOffset)
    cfg_params["MESH_FILENAME"] = mshName + ".su2"
    msh_params["PERIODIC_OFFSET"] = periodicOffset
    
    # Save cfg and mesh params to file
    CFGtxt_GEN_COMPLETE = su2cfg.dict_to_cfg(cfg_params, cfgName, filepath)
    MSHtxt_GEN_COMPLETE = su2cfg.dict_to_cfg(msh_params, mshName+"_params.txt", filepath)
    if type(flow_params) == dict:
        su2cfg.dict_to_cfg(flow_params, f"flow_params.txt", filepath)
    
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


def runParameterSweep_CFD(sweep_params:dict, cfg_params:dict, msh_params:dict, sweep_dir:Path, **kwargs):
    # Ensure the directory exists
    sweep_dir.mkdir(exist_ok=True)
    
    # Iterate through each sweep prameter
    for key in sweep_params: 
        # sweep through each parameter, each parameter has its own sweep, generate directory
        param_dir = sweep_dir / f"{key}_sweep" 
        param_dir.mkdir(exist_ok=True)
        
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
                    endStr = f"_{val:.2f}".replace(".","-")
                case str():
                    endStr = "_" + val 
                case list():
                    endStr = ""
                    for v in val:
                        endStr += f"_{v:.2f}".replace(".","-")
                case _:
                    endStr = f"_{val}"
                    
            run_dir = param_dir / f"{key}{endStr}"
            mshName = "Mesh_" + key + endStr 
            cfgName = "CFG_"+ key + endStr + ".cfg" 
            
            # Set naming params
            base_msh_params["FileLocation"] = str(run_dir)
            base_msh_params["MeshName"] = mshName 
            base_cfg_params["MESH_FILENAME"] = mshName + ".su2"
            
            # Make the parameter adjustment
            sweep_cfg_params, sweep_msh_params = PA.AdjustParams(key, val, base_cfg_params, base_msh_params)
            
            # Run the CFD by passing kwargs into the single run, mesh generation handled within single run
            runKwargs = kwargs.copy()
            runKwargs["MeshName"] = mshName
            runKwargs["ConfigName"] = cfgName
            # runKwargs["FilePath"] = str(run_dir)
            runSinglePoint_CFD(sweep_cfg_params, sweep_msh_params, filepath=run_dir, **runKwargs)
                
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
        if GENERATED:
            periodicOffset = swirlerObj.get_PeriodicOffset()
        
    return GENERATED, periodicOffset


# =============================================================================
#  Data Processors
# =============================================================================
def processpParamData(filepath:str, IS_3D:bool=True, 
                     cfg_path:str|None   = None,
                     units_path:str|None = None,
                     ):
    process_single_run(filepath, cfg_path, units_path)
    
    
    # if not IS_3D:
    #     DP.main(filepath, tol_in=tol_in)
    #     DPP.main(filepath, "png", dpi=dpi)
    #     DPS.main(filepath, "png", dpi=dpi)
    # else:
    #     # 3-D pipeline
    #     try:
    #         DP3.main(filepath, cfg_path, units_path=units_path)
    #         DPP3.main(filepath, cfg_path, units_path=units_path, img_format="png", dpi=dpi)
    #         DPS3.main(filepath, cfg_path, units_path=units_path, img_format="png", dpi=dpi)
    #     except:
    #         print("[Error]\t Error in 3d processing")
    return None
            


# =============================================================================
#  yaml imports
# =============================================================================
def importBaseParams_raw(filepath:str = None,
                     msh_filename:str = "mesh_base_params.yaml",
                     cfg_filename:str = "cfg_base_params.yaml") -> (dict, dict):
    '''
    Imports the .yaml configuration/parameter files in their raw form, a nested
    dictionary containing parameters and their values for each solver and mesh type. 

    Parameters
    ----------
    filepath : str, optional
        Filepath to the folder containing the .yaml files. The default is None which utilizes base parameters within module.
    msh_filename : str, optional
        filename (including extention) for the mesh parameters. The default is "mesh_base_params.yaml".
    cfg_filename : str, optional
        filename (including extention) for the su2 config parameters. The default is "cfg_base_params.yaml".

    Returns
    -------
    (msh_params: dict, cfg_params: dict)
        Returns the two dictionaries of nested parameters and their values..

    '''
    
    if filepath == None:
        # Get the absolute path of the current file
        current_file_path = Path(__file__).resolve()
        # If you only need the directory containing the file
        filepath = str(current_file_path.parent) + "\\config"
    
    # Combine the filepaths
    msh_filename = filepath + "\\" + msh_filename 
    cfg_filename = filepath + "\\" + cfg_filename
    
    try:
        with open(msh_filename, 'r') as file:
            msh_params = yaml.safe_load(file)
        
        print("[Info]\t Mesh parameters loaded.")
    except FileNotFoundError:
        print(f"[Error]\t The file {msh_filename} was not found.")
    except yaml.YAMLError as e:
        print(f"[Error]\t parsing YAML file: {e}")
        
        
    try:
        with open(cfg_filename, 'r') as file:
            cfg_params = yaml.safe_load(file)
        
        print("[Info]\t Config parameters loaded.")
    except FileNotFoundError:
        print(f"[Error]\t The file {cfg_filename} was not found.")
    except yaml.YAMLError as e:
        print(f"[Error]\t parsing YAML file: {e}")
        
    # Makes flow updates:
    cfg_params = updateFlow(cfg_params)
           
        
    return msh_params, cfg_params

def importBaseParams_filtered(solverType:str="RANS",
                              meshType:str  ="3D",
                              filepath:str = None,
                              msh_filename:str = "mesh_base_params.yaml",
                              cfg_filename:str = "cfg_base_params.yaml") -> (dict, dict):
    '''
    Imports the yaml configuration files from filepath and filteres the mesh and configuration
    dictionaries according to solverType and meshType

    Parameters
    ----------
    solverType : str, optional
        Either "RANS" or "Navier". The default is "RANS".
    meshType : str, optional
        Either "2D" or "3D". The default is "3D".
    filepath : str, optional
        Filepath to the folder containing the .yaml files. The default is None which utilizes base parameters within module.
    msh_filename : str, optional
        filename (including extention) for the mesh parameters. The default is "mesh_base_params.yaml".
    cfg_filename : str, optional
        filename (including extention) for the su2 config parameters. The default is "cfg_base_params.yaml".

    Raises
    ------
    ValueError
        When solverType or mshType does not match available options.

    Returns
    -------
    (msh_params: dict, cgf_params: dict, flow_params: dict)
        Returns the two dictionaries of unnested parameters and their values.

    '''
    
    
    
    # Import dictionaries
    msh_params_raw, cfg_params_raw = importBaseParams_raw(filepath, msh_filename, cfg_filename)
    
    # Begin Building base dictionaries
    msh_params = {**msh_params_raw["General"]}
    cfg_params = {**cfg_params_raw["General"],
                  **cfg_params_raw["GasProperties"],
                  **cfg_params_raw["BoundaryConditions"],
                  **cfg_params_raw["Convergence"],
                  **cfg_params_raw["Outputs"]}
    flow_params = {**cfg_params_raw["FlowConditions"]}
    
    # Add specific based on solver and mesh:
    match solverType:
        case "RANS":
            cfg_params.update(cfg_params_raw["Solver_RANS"])
            cfg_params.update(cfg_params_raw["Convergence_RANS"])
            cfg_params.update(cfg_params_raw["Outputs_RANS"])
            print("[Info]\t RANS solver type selected.")
        case "Navier":
            cfg_params.update(cfg_params_raw["Solver_Navier"])
            cfg_params.update(cfg_params_raw["Outputs_Navier"])
            print("[Info]\t Navier Stokes solver type selected.")
        case _: 
            # Default Case
            raise ValueError(f"[Error]\t Incorrect solverType imported: {solverType}\n\t Must be 'RANS' or 'NAVIER'")
    
    match meshType:
        case "3D":
            cfg_params.update(cfg_params_raw["BC_3D"])
            msh_params.update(msh_params_raw["mesh3D"])
            msh_params["IS_3D"] = True
            print("[Info]\t 3D mesh type selected.")
        case "2D":
            cfg_params.update(cfg_params_raw["BC_2D"])
            msh_params.update(msh_params_raw["mesh2D"])
            msh_params["IS_3D"] = False
            print("[Info]\t 2D mesh type selected.")
        case _: 
            # Default Case
            raise ValueError(f"[Error]\t Incorrect meshType imported: {meshType}\n\t Must be '2D' or '3D'")
    
    
    return msh_params, cfg_params, flow_params

def updateFlow(cfg_params:dict):
    T = cfg_params["FlowConditions"]["Inlet_Static_Temperature"]
    P = cfg_params["FlowConditions"]["Inlet_Static_Pressure"]
    # M = cfg_params["FlowConditions"]["Inlet_Target_Mach"]
    gam = cfg_params["GasProperties"]["GAMMA_VALUE"]
    Pb = cfg_params["FlowConditions"]["Back_Pressure"]
    R = cfg_params["GasProperties"]["GAS_CONSTANT"]
    rh = cfg_params["FlowConditions"]["r_hub_in"] * 2.54 / 100
    rt = cfg_params["FlowConditions"]["r_tip_in"] * 2.54 / 100
    mdot = cfg_params["FlowConditions"]["Inlet_Target_mdot"]         
         
    # Calculate density and velocity from static conditions and target
    A = np.pi*(rt**2 - rh**2)
    # mass flow rate               
    rho = P / (R*T)
    V = mdot / (rho * A)
    M = V / np.sqrt(gam * R * T)
    
    Ttot = flowcalc.Stag_Temperature(M, T, gam)
    Ptot = flowcalc.Stag_Pressure(M, P, gam)
    
    # print(f"V1 = {V:.3f} m/s\t V2 = {V2:.3f}m/s")
    # print(f"rho1 = {rho:.3f} kg/m^3\t rho2 = {rho2:.3f}kg/m^3")
    # print(f"P1 = {P:.3f} Pa\t P2 = {P2:.3f} Pa")
    
    # Make updates:
    if cfg_params["BoundaryConditions"]["INLET_TYPE"] == "MASS_FLOW":
        cfg_params["BoundaryConditions"]["MARKER_INLET"] = cfg_params["BoundaryConditions"]["MARKER_INLET"].format(rho, V)
    
    else:
        cfg_params["BoundaryConditions"]["MARKER_INLET"] = cfg_params["BoundaryConditions"]["MARKER_INLET"].format(Ttot, Ptot)
    
    cfg_params["BoundaryConditions"]["MARKER_OUTLET"] = cfg_params["BoundaryConditions"]["MARKER_OUTLET"].format(Pb)
    
    
    return cfg_params
   

# =============================================================================
#  Testing
# =============================================================================
if __name__=="__main__":
    msh_params, cfg_params, flow_params = importBaseParams_filtered("RANS", "3D")
    cfg_params["ITER"] = 10000 
  
    
  
    filepath=Path(r"C:\Users\BriceM\Documents\SU2 CFD Data\3D_Tests\MinThickTests\Test06")
    
    USE_DEFAULT_MESH = False
    default_mesh_name = "MedSwirl_CutTE_Refined04"
    default_mesh_location = Path(r"C:\Users\BriceM\Documents\SU2 CFD Data\3D_Tests\Meshes")
    
    # For saving a default mesh
    msh_params["FileLocation"] = str(default_mesh_location)
    msh_params["MeshName"] = str(default_mesh_name)
    su2cfg.dict_to_cfg(msh_params, default_mesh_name+".txt", default_mesh_location)
    generateMesh(msh_params, OpenGMSHVisual=True)
    
    
    
    
    # runSinglePoint_CFD(cfg_params, msh_params, filepath, flow_params=flow_params,
    #                    USE_DEFAULT_MESH=USE_DEFAULT_MESH,
    #                    default_mesh_location=default_mesh_location,
    #                    default_mesh_name=default_mesh_name)
    
    
    # runParameterSweep_CFD({"AoA_rt": [[5.0, 12.0], [5.0, 16.0], [5.0, 20.0]]}, cfg_params, msh_params, filepath)
    
    # processpParamData(filepath)
    # processSweepData(filepath, IS_3D=True)