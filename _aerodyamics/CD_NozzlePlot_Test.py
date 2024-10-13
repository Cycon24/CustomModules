# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 14:18:50 2024

@author: cycon
"""
import numpy as np
import matplotlib.pyplot as plt
import GasDynamics as GD

def x_pos(A_ft2, BeforeThroat):
    #𝐴(𝑥) = (𝜋/4)[(0.02246 ⋅ 𝑋) + 0.2504]^2  − 0.01327 where where X = distance from throat
    A_in2 = A_ft2 * (12/1)**2
    x_in = (-0.2504+np.sqrt((A_in2 + 0.01327)*(4/np.pi)))/0.02246
    
    return -x_in/8 if BeforeThroat else x_in


class CD_Nozzle():
    def __init__(self, units='SI'):
        print('CD Nozzle Class')
        # We want a general form CD nozzle which can fully explain the theoretical flow through the nozzle
        # Want to have the option for self-made x vs area functions, 
        self.Po = None
        self.Patm = None
        
        # Nozzle Geometry
        self.Ai = None 
        self.At = None
        self.Ae = None 
        self.A_star = None
        
        # Gas Constants
        self.gamma = 1.4 # 1.4 for diatomic, 5/3 for monotomic, 1.3 for other
        self.R_u = 10.73158     # [psi⋅ft3/lbmol⋅°R]  - Universal Gas Constant (Imperial)
        self.R_u = 8314         # [J/kmol*K]
        self.MW  = 28.966       # [kg/kmol] [lbm/lbmol]
        self.R_air = 287.1   # [J/kg*K)]
        
        # Unit coversions
        self.gc = [1, 32.174] # [NonDim] [𝑙𝑏𝑚−𝑓𝑡/𝑙𝑏𝑓−𝑠^2] - Gravitational Constant
        
    def set_gas_constants( MW=None, R=287.1, gamma=1.4):
        '''
        Set the gas constants of the gas flowing throug the nozzle. Can input 
        R if known or can utilize the gas Molecular Weight.If no inputs given 
        the gas is assumed to be air and in SI units.

        Parameters
        ----------
        R : Float, optional
            Gas Constant. The default is 287.1 J/kg*K.
        gamma : Float, optional
            Specific Heat Ratio. The default is 1.4.
        MW : Float, optional
            Molecular weight of gas [kg/kmol] or [lbm/lbmol]. The default is None.

        Returns
        -------
        None.

        '''
        # do things
        
    def set_area_function(self,function):
        print('NEED CREATED')
        self.f_A_x = function
        
    def calculate_critical_pressure_ratios(self):
        print("NEEDS CREATED")
        
    def find_shock_locations(self, PR_bs):
        print("NEEDS CREATED")


# Inputs
Po = 99.36135 # Psi
P_atm = 14.361353041951999  # psi

Pbs_g = np.array([0.000000000000000000e+00,
                3.000000000000000000e+01,
                4.000000000000000000e+01,
                5.000000000000000000e+01,
                6.000000000000000000e+01,
                6.500000000000000000e+01,
                7.000000000000000000e+01,
                7.500000000000000000e+01,
                8.000000000000000000e+01]) #86.03-P_atm
Pbs_g.sort()
Pbs = Pbs_g + P_atm
PR_bs = Pbs/Po

Ae = 0.00035676 # ft^2
At = 0.00024980 # ft^2
Ai_At = 10 # APproximately stag condition
Ae_At = Ae/At 
A_rats = np.append(np.linspace(Ai_At, 1,num=50), np.linspace(1,Ae_At,num=50,endpoint=True))  

R_univ = 10.731577089016    # [psi⋅ft3/lbmol⋅°R]  - Universal Gas Constant
MW_air = 28.966             # [lbm/lbmol]
R      = R_univ / MW_air    # [psi⋅ft3/lbm⋅°R] Gas Constant if  air
g_c   = 32.174              # [𝑙𝑏𝑚−𝑓𝑡/𝑙𝑏𝑓−𝑠^2] - Gravitational Constant
gamma = 1.4

# Get isentropic ratio
M_isen_sub, M_isen_sup = GD.Mach_at_A(Ae_At, Gamma=1.4)
PR_isen_sup = 1/GD.Po_P_ratio(M_isen_sup, Gamma=1.4) # Pe/Po
PR_isen_sub = 1/GD.Po_P_ratio(M_isen_sub, Gamma=1.4) # Pe/Po
PR_exit_shock = GD.P2_P1_n(M_isen_sup, gamma)*(1/GD.Po_P_ratio(M_isen_sup, gamma))
PR_crit = 0.5283 # DEPENDENT ON GAMMA

# Find Shock Area Locations from back pressure ratios
As_Ats = np.empty(PR_bs.shape)
for iPb, Pb_Po in enumerate(PR_bs):
    if (Pb_Po < PR_isen_sub) and (Pb_Po > PR_exit_shock):
        # Flow is choked and in between sub isen solutions and exit shock -> shock in noz
        # Find Location of shock
        As_Ats[iPb] = GD.As_At_n(Ae_At, Pb_Po, gamma,)
    else:
        As_Ats[iPb] = None
        

# Setup for solving for pressure ratios
# Will be len(A_rats) rows by len(PR_bs) cols
Px_Pos = np.empty([len(A_rats),len(PR_bs)],dtype=float) 
BeforeThroat = True     # Define toggle for before and after throat
tol = 1e-5

for iA, Ax_At in enumerate(A_rats):
    if abs(Ax_At - 1) < tol:
        BeforeThroat=False # Toggle swapped so now everything is after throat
    
    # if not AfterThroat:
    #     # Pressure ratio will  be subsonic solution of A ratio 
    #     M_sub, M_sup = GD.Mach_at_A(Ax_At, gamma) # only for choked
    #     Px_Pos[iA,:] = 1/GD.Po_P_ratio(M_sub, gamma)
    # else:
        # Cycle through the back pressure ratios for after throat cases to find pressure
    for iPb, Pb_Po in enumerate(PR_bs):
        if Pb_Po > PR_isen_sub:
            # Not Choked Flow, find Mach number at exit plane with back pressure, (pe=pb)
            # Subsonic Cases (Assuming Pe=Pb)
            # Works for all Aratios when not choked
            
            # Find theoretical mach at exit
            Me_the = GD.Mach_at_PR(1/Pb_Po, gamma)
            # Find theoretical throat area ratio that would need to be choked
            Ae_Astar = GD.A_ratio(Me_the, gamma)
            
            # Get ratio of Ax/A* to find pressure at x loction
            Ax_Astar = Ax_At * (1/Ae_At) * Ae_Astar
            
            # Find Mach at x location 
            '''  THIS NEEDS A CHECK FOR WHEN APPROACHING CHOKED FLOW 
            IT WILL NOT SOLVE '''
            M_sub, M_sup = GD.Mach_at_A(Ax_Astar, gamma)
            
            # Pressure Ratio at x location
            Px_Pos[iA, iPb] = 1/GD.Po_P_ratio(M_sub, gamma)
            # The above can be adapted to  find PR throught nozzle
            
        elif Pb_Po > PR_exit_shock:
            # choked flow, But shock occurs in nozzle (not ful expand, pb /= pe)
            # Will follow supersonic isentropic solution until As_At is reached
            if BeforeThroat:
                # Need to use subsonic solution
                M_sub, M_sup = GD.Mach_at_A(Ax_At, gamma)
                # Pressure Ratio at x location
                Px_Pos[iA, iPb] = 1/GD.Po_P_ratio(M_sub, gamma)
                
            else: 
                # After throat
                if Ax_At < As_Ats[iPb]: # Checks if before shock
                    # BEFORE SHOCK LOCATION 
                    # Need to use sonic solution
                    M_sub, M_sup = GD.Mach_at_A(Ax_At, gamma)
                    # Pressure Ratio at x location
                    Px_Pos[iA, iPb] = 1/GD.Po_P_ratio(M_sup, gamma)
                
                else:
                    # AFTER SHOCK LOCATION
                    # Need to get Ax_Astar after shock and use subsonic solution
                    # Get shock conditions
                    M_empty, Ms1 = GD.Mach_at_A(As_Ats[iPb],gamma)
                    Ms2 = GD.Mach2_n(Ms1, gamma)
                    Po2_Po1 = GD.Po2_Po1_n(Ms1,gamma)
                    
                    As_Astar = GD.A_ratio(Ms2, gamma) 
                    Ax_Astar = Ax_At * (1/As_Ats[iPb]) * As_Astar
                  
                    
                    M_sub, M_empty = GD.Mach_at_A(Ax_Astar, gamma)
                    
                    # Pressure Ratio at x location
                    Px_Po2 = 1/GD.Po_P_ratio(M_sub, gamma)
                    Px_Pos[iA, iPb] = Px_Po2 * Po2_Po1
                    if iPb == 3 and (iA  == 95 or iA == 96):
                        print('M1 = {:.4f} M2 = {:.4f} P1/Po = {:.4f} P2/Po = {:.4f}'.format(Ms1, Ms2,Px_Pos[iA-1, iPb],Px_Pos[iA, iPb]))
                    
                    
        elif Pb_Po > PR_isen_sup:
            # Flow is overexpanded so shock outside of nozzle
            # Will follow the supersonic  isentropic solution up to exit plane
            description = "Flow Overexpanded: Shock outside of nozzle"
            M_sub, M_sup = GD.Mach_at_A(Ax_At, gamma)
            if BeforeThroat:
                # Need to use subsonic solution
                # Pressure Ratio at x location
                Px_Pos[iA, iPb] = 1/GD.Po_P_ratio(M_sub, gamma)
            else:
                # Need to use sonic solution
                # Pressure Ratio at x location
                # if Pb_Po==PR_bs[-1]:
                #     Px_Pos[iA, iPb] = Pb_Po
                # else:
                Px_Pos[iA, iPb] = 1/GD.Po_P_ratio(M_sup, gamma) 
            Px_Pos[-1,iPb] = Pb_Po
        else:
            # Choked and fully expanded flow (pe = pb)
            # Will follow supersonic isentropic solution
            M_sub, M_sup = GD.Mach_at_A(Ax_At, gamma)
            if BeforeThroat:
                # Need to use subsonic solution
                # Pressure Ratio at x location
                Px_Pos[iA, iPb] = 1/GD.Po_P_ratio(M_sub, gamma)
            else:
                # Need to use sonic solution
                # Pressure Ratio at x location
                Px_Pos[iA, iPb] = 1/GD.Po_P_ratio(M_sup, gamma)
        
BeforeThroat=True
x_locs = np.empty(A_rats.shape)
for iA, Ax_At in enumerate(A_rats):
    if abs(Ax_At - 1) < tol:
        BeforeThroat=False

    x_locs[iA] = x_pos(Ax_At*At, BeforeThroat)
# =============================================================================
# Plotting
# =============================================================================
plt.close('all')
plt.figure(5)

stretch = 1.1/1.625
x_locs = x_locs*stretch +  0.9
# 0->1.62479: dx = 1.625
#into 
# 0.9 -> 2: dx = 1.1

for i, Pb in enumerate(PR_bs):
    plt.plot(x_locs,Px_Pos[:,i],label='Pb={:.4f}'.format(Pbs_g[i]))
plt.plot([min(x_locs), max(x_locs)], [0.5283, 0.5283], '--', label="Pb=P_crit")
plt.plot([min(x_locs), max(x_locs)], [PR_isen_sup, PR_isen_sup],'--', label="Pb=P_fe")
plt.plot([min(x_locs), max(x_locs)], [PR_isen_sub, PR_isen_sub],'--', label="Pb=P_ue")
plt.plot([min(x_locs), max(x_locs)], [PR_exit_shock, PR_exit_shock],'--', label="Pb=P_e_shock")
plt.xlabel('x loc [in]')
plt.ylabel('P/Po')
plt.xlim(-0.5, 2.5); plt.ylim(-0.25, 1.1)
plt.legend()
plt.grid('on')