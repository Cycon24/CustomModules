# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 11:57:36 2026

@author: BriceM
"""
import _aerodynamics.GasDynamics as GD 
import cantera as ct 

def Stag_Temperature(M:float, T:float, gamma:float=1.4):
    return T * GD.To_T_ratio(M, gamma)

def Stag_Pressure(M:float, P:float, gamma:float=1.4):
    return P * GD.Po_P_ratio(M, gamma)

def combustion_products(phi:float, Ti:float=None)->dict:
    
    gas1 = ct.Solution('gri30.yaml')  
    # Define fuel (propane) and oxidizer (air)
    fuel = "C3H8"
    oxidizer = "O2:1, N2:3.76"
    
    # species = ["CO2", "CO", "H2O", "H2", "H", "OH", "O2", "O", "NO", "N2", "N"]
    
    # set equivs
    # phi1 = 0.25
    # phi2 = 1.2 
    
    # Set TP
    Ti = 300 if Ti == None else Ti
    
    # Set the mixture
    gas1.set_equivalence_ratio(phi, fuel=fuel, oxidizer=oxidizer)
    # Set propterties, assuming T~300K and P as 1 atm
    gas1.TP = Ti, ct.one_atm # (Pressure not affecting combustion)
    
    # Run the reaciton, keep enthalpy and pressure constant
    gas1.equilibrate('HP')
    # print(gas1())
    
    cp = gas1.cp_mass 
    cv = gas1.cv_mass
    density = gas1.density 
    T = gas1.T
    P = gas1.P
    
    gamma = cp / cv 
    R = cp - cv 
    
    return {"Gamma": gamma, "R": R, "T":T, "P":P}

if __name__=="__main__":
    print(combustion_products(0.235, 290))