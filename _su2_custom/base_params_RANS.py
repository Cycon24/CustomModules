# -*- coding: utf-8 -*-
"""
Created on Thu Sep 25 14:33:53 2025

@author: BriceM
"""

mesh_params = {
# Core
"MeshName": "swirlerMeshTest",
"FileLocation": 'meshs',
"NACA": "6412",
"FileType": "BOTH",
"GenPointSize": 0.05,
"AF_PointSize": 0,
"AF_numPoints": 800,
"FullMesh": False,

# Tunnel Settings
"numBlades": 6,
"BladeAoA": 5,
"Radius": 2,
"chord": 2,
"L_Upstream": 5,      # why: default derives from 5 * chord; chord default is 1
"L_Downstream": 10,   # why: default derives from 10 * chord; chord default is 1

# Boundary-Layer Settings (used only if GenBL is True)
"GenBL": True,
"BL_nLayers": 30,
"BL_BetaCoef": 1.005,
"BL_Quad": 1,
# "BL_Factor": 1.1,
"BL_h1": 1e-4,
"BL_Thickness": 5,

# Trailing-Edge Fan Settings (used only if GenFan_TE is True)
"GenFan_TE": True,
"TEF_Length": 4,   # why: no default in source; explicit sentinel
"TEF_FanMultFactor": 2,
"TEF_minSize": 5e-3

}

cfg_params: dict[str, str] = {
    # Mesh definition
    "MESH_FORMAT": "SU2",
    "MESH_FILENAME": "swirler_n6_c2_AoA10_NACA8412.su2",
    # Solver Type
    "MATH_PROBLEM": "DIRECT",
    "SOLVER": "RANS",
    "KIND_TURB_MODEL": "SA", # Or SST, SA or SA_NEG
    "PRANDTL_TURB": "0.9", # degault is 0.9
    # # May need to add:
    # "FREESTREAM_TURBULENCEINTENSITY": 0.05,
    # "FREESTREAM_TURB2LAMVISCRATIO": 10.0,
    # # For SA turb model
    # "FREESTREAM_NU_FACTOR": 3.0,
    # Restart
    "RESTART_SOL": "NO",
    "WRT_RESTART_COMPACT": "NO",
    
    # Gas and Reference
    "SYSTEM_MEASUREMENTS": "US",
    "FLUID_MODEL": "STANDARD_AIR",
    "GAMMA_VALUE": "1.4",
    "GAS_CONSTANT": "1716.49",
    "MACH_NUMBER": "0.30",
    "AOA": "0.0",
    "FREESTREAM_TEMPERATURE": "1660.0",
    "FREESTREAM_PRESSURE": "2116.216",
    "REYNOLDS_NUMBER": "8.981e+04",
    "REYNOLDS_LENGTH": "2",

    "REF_LENGTH": "2.0",
    "REF_AREA": "1.0",
    
    # Viscocity and thermal
    "VISCOSITY_MODEL": "SUTHERLAND",
    "MU_REF": "3.737e-7",
    "MU_T_REF": "518.67",
    "SUTHERLAND_CONSTANT": "198.72",
    "PRANDTL_LAM": "0.72",
    
    # Boundary Conditions
    "MARKER_PERIODIC": "(Symmetry2, Symmetry1) , (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, -{}, 0.0)",
    "MARKER_INLET": "( Inlet, 1689.88, 2252.564, 1.0, 0.0, 0.0 )",
    "MARKER_OUTLET": "( Outlet, 2116.216 )",
    "MARKER_HEATFLUX": "( Airfoils, 0.0)",
    # Plotting Markers
    "MARKER_PLOTTING": "( Airfoils, Outlet )",
    "MARKER_MONITORING": "( Airfoils )",

    # Residuals
    "NUM_METHOD_GRAD": "WEIGHTED_LEAST_SQUARES", # GREEN_GAUSS or WEIGHTED_LEAST_SQUARES
    "MUSCL_FLOW": "YES",
    "SLOPE_LIMITER_FLOW": "VENKATAKRISHNAN",
    "VENKAT_LIMITER_COEFF": "0.05", # default is 0.05, may need to increase to reduce interpolation slope limiting
    
    # Convergence Params
    "CONV_NUM_METHOD_FLOW": "ROE", # May need to change to JST
    "ENTROPY_FIX_COEFF": "0.05",
    "LOW_MACH_PREC": "NO",
    "ROE_KAPPA": "0.5",
    "MIN_ROE_TURKEL_PREC": "0.01",
    "MAX_ROE_TURKEL_PREC": "0.2",
    # For JST
    "JST_SENSOR_COEFF": "(0.5, 0.02)", # default (0.5, 0.02)
    
    # Turbulence Convergence Params
    "CONV_NUM_METHOD_TURB": "SCALAR_UPWIND",
    "MUSCL_TURB": "NO", # other rans used NO
    "SLOPE_LIMITER_TURB": "VENKATAKRISHNAN",
    "TIME_DISCRE_TURB": "EULER_IMPLICIT",
    
    # Time Convergence Params
    "TIME_DISCRE_FLOW": "EULER_IMPLICIT", # Should be okay
    "LINEAR_SOLVER": "FGMRES",
    "LINEAR_SOLVER_PREC": "ILU",
    "LINEAR_SOLVER_ITER": "20", # Reduced from 120 to 20, rans sim used 10
    "LINEAR_SOLVER_ERROR": "1e-10", # Increased from e-6 to e-10 from other rans sim
    
    
    # CFL Params
    "CFL_NUMBER": "0.5", # May need ot increase CFL Vals, other rans sim used 100 no adapt and 100000 max CFL, 15 min
    "CFL_ADAPT": "YES",
    "CFL_ADAPT_PARAM": "( 0.5, 1.2, 0.01, 3.0,1e-3, 0 )", # Try keeping CFL low (3 rather than 30)
    "RK_ALPHA_COEFF": "( 0.66667, 0.66667, 1.000000 )",
    
    # Overall Convergence Params
    "ITER": "10000",
    "CONV_STARTITER": "20",
    "CONV_RESIDUAL_MINVAL": "-9",
    "CONV_CAUCHY_ELEMS": "100",
    "CONV_CAUCHY_EPS": "1e-6",
    # Time convergence monitoring
    # "WINDOW_CAUCHY_CRIT"   : "YES",
    "CONV_FIELD"           : "DRAG",
    # "CONV_WINDOW_FIELD"    : "(TAVG_DRAG, TAVG_LIFT)",# % List of time convergence fields
    # "CONV_WINDOW_STARTITER"  : "0",
    # "CONV_WINDOW_CAUCHY_EPS" : "1E-3", #% Epsilon to control the series convergence
    # "CONV_WINDOW_CAUCHY_ELEMS" : "10", #  % Number of elements to apply the criteria

    # Outputs
    "OUTPUT_FILES": "( RESTART, PARAVIEW, SURFACE_CSV, CSV, SURFACE_PARAVIEW )",
    "VOLUME_OUTPUT": "(COORDINATES, SOLUTION, PRIMITIVE, CUSTOM)",
    "OUTPUT_WRT_FREQ": "250, 250, 250, 250, 250",
    "TABULAR_FORMAT": "CSV",

    # Output Filenames
    "SOLUTION_FILENAME": "solution_flow",
    "SOLUTION_ADJ_FILENAME": "solution_adj",
    "CONV_FILENAME": "history",
    "RESTART_FILENAME": "entire_surface_restart",
    "RESTART_ADJ_FILENAME": "restart_adj",

    "VOLUME_FILENAME": "flow",
    "VOLUME_ADJ_FILENAME": "adjoint",
    "GRAD_OBJFUNC_FILENAME": "of_grad",

    "SURFACE_FILENAME": "surface_flow_markers",
    "SURFACE_ADJ_FILENAME": "surface_adjoint",
    # Screen and mesh output
    "SCREEN_OUTPUT": "(INNER_ITER, WALL_TIME, RMS_DENSITY, RMS_MOMENTUM-X, RMS_MOMENTUM-Y, RMS_ENERGY, RMS_NU_TILDE, LIFT, DRAG)",
    "MESH_OUT_FILENAME": "mesh_out",
}