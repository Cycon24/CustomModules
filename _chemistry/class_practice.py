# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 14:30:10 2025

@author: cycon
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Practice 1
# =============================================================================
# gas = ct.Solution('gri30.yaml')
# gas.TP = 500, ct.one_atm 
# gas.X =  "H2:1, O2:1, N2:3.76" # Set molar concentrations of each
# # Ensure properties have been set
# gas()


# rates = gas.net_rates_of_progress

# for i in range(gas.n_reactions):
#     if rates[i] != 0.0:
#         if rates[i] > 1e-60:
#             print(rates[i])
    

# =============================================================================
# Practice 2        
# =============================================================================
gas1 = ct.Solution('gri30.yaml')  
# She added
fuel = "C2H4"
oxidizer = "O2:1, N2:3.76"
#
# set up range of equivs
phis = np.arange(0.5, 2.0 + 0.05, 0.05)
Temps = np.zeros_like(phis)

for i, phi in enumerate(phis):
    # Dont need to compute manually, there is a funciton
    # X = (1-phi)/phi 
    # nO2 = X*(2+4/4)
    # nN2 = X*(2+4/4)*3.76
    # gas1.X = f"C2H4:1, O2:{nO2}, N2:{nN2}"
    gas1.set_equivalence_ratio(phi, fuel=fuel, oxidizer=oxidizer)
    # Now set propterties, need to set AFTER mixture set
    gas1.TP = 300, ct.one_atm 
    
    # Run the reaciton, keep enthalpy and pressure constant
    gas1.equilibrate('HP')
    Temps[i] = gas1.T
    
plt.figure()
plt.plot(phis, Temps)
plt.xlabel("Equiv Ratio")
plt.ylabel("Flame Temp")
plt.grid()
plt.show()
    
    
            
            