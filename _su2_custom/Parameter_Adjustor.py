# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 13:26:17 2025

@author: BriceM
"""

def AdjustParams(adj_param, adj_value, cfg_dict, msh_dict):
    # Check if adj_param and adj_value are lists
    # if type(adj_param) != list: 
    #     # Both are not lists, make lists
    #     adj_param = [adj_param]
    # if  type(adj_value) != list:
    #     adj_value = [adj_value]
    # # If either are not a list now, then one must be in list form and the other not,
    # # return an error
    # if type(adj_param) != list or type(adj_value) != list:
    #     raise ValueError("Parameter and Value inputs must both be lists or both be singular values.")
    
    # # # Check list lengths
    # if len(adj_param) != len(adj_value):
    #     raise ValueError("Parameter and Value lists must be the same length")
    
    # Now we can proceed;
    # for i, param in enumerate(adj_param): 
    param = adj_param
    val = adj_value
    
    # Now enter switch statement
    match param:
        case "AoA":
            # Angle of Attack Adjustments
            #  - mesh AoA, in degrees
            msh_dict["BladeAoA"] = val 
        case "nBlades":
            # Number of Blades in swirler
            #  - mesh numBlades, int
            msh_dict["numBlades"] = val 
            
        case "Radius":
            # Radius of cross-section in swirler
            #  - mesh Radius, float
            msh_dict["Radius"] = val
            
        case "NACA":
            # Airfoil shape as NACA 4-digit airfoil
            #  - mesh NACA, string
            if len(val) != 4:
                raise ValueError("NACA airfoil must be a 4-digit string")
            msh_dict["NACA"] = val
            
        case "Chord":
            msh_dict['chord'] = val
# =============================================================================
#                 3D Meshes
# =============================================================================
        case "Pb":
            # Backpressure
            cfg_dict["MARKER_OUTLET"] = f"( Outlet, {val} )"
            
        case "AoA_rt":
            # Angle of attack of root a nd tip profiles
            msh_dict["BladeAoA_root"] = val[0]
            msh_dict["BladeAoA_tip"] = val[1]
            
        case _: 
            # Default Case
            raise ValueError("Adjust Parameter not Recognized")
    
    return cfg_dict, msh_dict 