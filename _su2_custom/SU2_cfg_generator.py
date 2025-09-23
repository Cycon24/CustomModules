# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 22:17:18 2025

@author: BriceM
"""
config_params = {
    "MESH_FILENAME": "mesh.su2",
    "SOLVER": "UNSTEADY_NAVIER_STOKES",
    "FREESTREAM_MACH": 0.8,
    "FREESTREAM_TEMPERATURE": 288.15,
    "FREESTREAM_PRESSURE": 101325.0,
    "MARKER_FARFIELD": "( Farfield, 1.0, 0.0, 0.0, 100000.0, 288.15 )",
    "MARKER_EULER": "( Airfoil, 0.0, 0.0, 0.0, 0.0, 0.0 )",
    "TIME_ITER": 1000,
    "CFL_NUMBER": 1.0,
    "CONV_RESIDUAL_MINVAL": -12,
    # ... add more parameters as needed
}

'''
Need to convert all of this into a dictionary, and then alter the variables needed.
% --- Mesh / I-O
MESH_FORMAT= SU2 
MESH_FILENAME = swirler_n6_c2_AoA10_NACA8412.su2

MATH_PROBLEM  = DIRECT
SOLVER        = NAVIER_STOKES
KIND_TURB_MODEL = NONE            % (set SA later if you want RANS)

RESTART_SOL= NO
WRT_RESTART_COMPACT= NO

% --- Gas & reference (Imperial)
SYSTEM_MEASUREMENTS= US 
FLUID_MODEL = STANDARD_AIR        %, FLUID_MIXTURE or COOLPROP
GAMMA_VALUE    = 1.4             % This gets overwritten to 1.4 anyways
GAS_CONSTANT   = 1716.49          % ft*lbf/(slug*R)
MACH_NUMBER    = 0.30
AOA            = 0.0 
FREESTREAM_TEMPERATURE = 1660.0   % R (use 518.67 for a cool smoke test)
FREESTREAM_PRESSURE    = 2116.216 % lbf/ft^2 (1 atm)
REYNOLDS_NUMBER = 8.981e+04
REYNOLDS_LENGTH = 2 % 

% Reference length for pitching, rolling, and yawing non-dimensional moment
% Reference area for force coefficients (0 implies automatic calculation)
REF_LENGTH= 2.0
REF_AREA= 1.0


% --- Viscosity & thermal models
VISCOSITY_MODEL     = SUTHERLAND
MU_REF              = 3.737e-7     % slug/(ft*s) at MU_T_REF
MU_T_REF            = 518.67       % R
SUTHERLAND_CONSTANT = 198.72       % R
PRANDTL_LAM         = 0.72

% --- Boundaries
% Periodic (slave, master, rot_cx, rot_cy, rot_cz, rot_ax, rot_ay, rot_az, tx,ty,tz)
MARKER_PERIODIC = (Symmetry2, Symmetry1) , (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, -6.283185307179586, 0.0)

% Inlet - total T [R], total P [lbf/ft^2], dir vec (unit)
MARKER_INLET   = ( Inlet, 1689.88, 2252.564, 1.0, 0.0, 0.0 )

% Outlet - static pressure (subsonic)
MARKER_OUTLET  = ( Outlet, 2116.216 )

% Walls - viscous no-slip, adiabatic BC
MARKER_HEATFLUX = ( Airfoils, 0.0) % , Symmetry1, 0.0, Symmetry2, 0.0)

MARKER_PLOTTING = ( Airfoils, Outlet )
MARKER_MONITORING = ( Airfoils, Outlet )

% --- Multigrid
%MGLEVEL         = 3
%MG_CYCLE        = V_CYCLE
%MG_PRE_SMOOTH   = ( 1, 1, 1, 1 )
%MG_POST_SMOOTH  = ( 0, 0, 0, 0 )

% --- Space discretization & numerics
NUM_METHOD_GRAD         = WEIGHTED_LEAST_SQUARES
MUSCL_FLOW              = YES
SLOPE_LIMITER_FLOW      = VENKATAKRISHNAN
VENKAT_LIMITER_COEFF    = 0.02

% Robust flux for internal subsonic flow
CONV_NUM_METHOD_FLOW    = ROE 
ENTROPY_FIX_COEFF       = 0.05
LOW_MACH_PREC           = NO
ROE_KAPPA               = 0.5
MIN_ROE_TURKEL_PREC     = 0.01
MAX_ROE_TURKEL_PREC     = 0.2

% --- Time (steady implicit), linear solver
TIME_DISCRE_FLOW       = EULER_IMPLICIT
LINEAR_SOLVER          = FGMRES
LINEAR_SOLVER_PREC     = ILU
LINEAR_SOLVER_ITER     = 120
LINEAR_SOLVER_ERROR    = 1e-6

% --- CFL (use the 6-parameter adaptation?)
CFL_NUMBER       = 0.5
CFL_ADAPT        = YES
CFL_ADAPT_PARAM  = ( 0.5, 1.2, 0.01, 30.0,1e-3, 0 )
%                = (down, up, CFL_min, CFL_max, acceptable_lin_res, start_iter)


% --- Convergence settings 
ITER                 = 1250
CONV_STARTITER       = 20
CONV_RESIDUAL_MINVAL = -8
CONV_CAUCHY_ELEMS    = 100
CONV_CAUCHY_EPS      = 1e-6



% --- Output 
% CUSTOM_OUTPUTS      = 'probe_VY : Probe{VELOCITY_Y}[0.005, 0.005]'
OUTPUT_FILES        = ( RESTART, PARAVIEW, SURFACE_CSV, CSV, SURFACE_PARAVIEW )
VOLUME_OUTPUT       = (COORDINATES, SOLUTION, PRIMITIVE, CUSTOM)
OUTPUT_WRT_FREQ     = 250, 250, 250, 250, 250
TABULAR_FORMAT      = CSV


% Need these? 
SOLUTION_FILENAME= solution_flow
SOLUTION_ADJ_FILENAME= solution_adj
CONV_FILENAME= history
RESTART_FILENAME= entire_surface_data % restart
RESTART_ADJ_FILENAME= restart_adj

% vtu output
VOLUME_FILENAME= flow
VOLUME_ADJ_FILENAME= adjoint

% Objective func
GRAD_OBJFUNC_FILENAME= of_grad

% Surface outputs
SURFACE_FILENAME= surface_flow_markers
SURFACE_ADJ_FILENAME= surface_adjoint

SCREEN_OUTPUT  = (INNER_ITER, WALL_TIME, RMS_DENSITY, RMS_MOMENTUM-X, RMS_MOMENTUM-Y, RMS_ENERGY)
MESH_OUT_FILENAME= mesh_out

'''

config_filename = "my_su2_config.cfg"

with open(config_filename, "w") as f:
    for key, value in config_params.items():
        f.write(f"{key}= {value}\n")

print(f"SU2 configuration file '{config_filename}' created successfully.")