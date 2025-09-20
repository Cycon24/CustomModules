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


config_filename = "my_su2_config.cfg"

with open(config_filename, "w") as f:
    for key, value in config_params.items():
        f.write(f"{key}= {value}\n")

print(f"SU2 configuration file '{config_filename}' created successfully.")