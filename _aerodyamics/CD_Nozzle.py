# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 14:38:18 2024

@author: cycon
"""
import numpy as np
import matplotlib.pyplot as plt
import GasDynamics as GD
# CHANGE TO INHERIT FROM STAGE
class CD_Nozzle():
    def __init__(self, units='SI'):
        self.SI_units = True if units=='SI' else False #?? Better way of handling?
        
        # We want a general form CD nozzle which can fully explain the theoretical flow through the nozzle
        # Want to have the option for self-made x vs area functions, 
        self.Po = None
        self.Poe =None
        self.Patm = None
        
        # Nozzle Geometry
        self.Ai = None 
        self.At = None
        self.Ae = None 
        self.A_star = None
        
        # Gas Constants
        self.gamma = None # 1.4 for diatomic, 5/3 for monotomic, 1.3 for other
        self.MW  = None         # [kg/kmol] [lbm/lbmol]
        self.R = None           # [J/kg*K)] [psi⋅ft3/lbm⋅°R]
        
        # Unit coversions
        self.gc = [1, 32.174] # [NonDim] [𝑙𝑏𝑚−𝑓𝑡/𝑙𝑏𝑓−𝑠^2] - Gravitational Constant
        
    def set_gas_constants(self, MW=None, R=287.1, gamma=1.4):
        '''
        Set the gas constants of the perfect gas flowing throug the nozzle. Can 
        input R if known or can utilize the gas Molecular Weight.
        
        If MW is not input, then air is assumed. The units will be SI unless the units
        have been changed prior to the function call.
            - Units can  be changed manually with self.SI_units = True/False or when
            creating a class object through inputs.
            
        Parameters
        ----------
        MW : Float, optional
            Molecular weight of gas [kg/kmol] or [lbm/lbmol]. The default is 28.966 for air.
        R : Float, optional
            Gas Constant. The default is 287.1 J/kg*K.
        gamma : Float, optional
            Specific Heat Ratio. The default is 1.4.
        
        Returns
        -------
        None.

        '''
        R_u = 8314 if self.SI_units else 10.73158  # [J/kmol*K] [psi⋅ft3/lbmol⋅°R]  - Universal Gas Constant (Imperial)
        
        self.gamma = gamma # 1.4 for diatomic, 5/3 for monotomic, 1.3 for other
        if MW == None:
            self.R = R
        else:
            self.R = R_u/MW     # [J/kg*K)] [psi⋅ft3/lbm⋅°R]
            self.MW = MW        # [kg/kmol] [lbm/lbmol]
                 
        
        
    def set_area_function(self,function):
        print('NEED CREATED')
        self.f_A_x = function
        
    def calculate_critical_pressure_ratios(self):
        '''
        Calculates all of the critical pressure rations associated with a CD Nozzle:
            - 1st Critical PR: Subsonic Isentropic Pressure Ratio
                self.PR_isen_sub
            - 2nd Critical PR: Normal Shock at Nozzle Exit Pressure Ratio
                self.PR_exit_shock
            - 3rd Critical PR: Throat Critical PR for Choked Flow
                self.PR_crit
            - 4th Critical PR: Supersonic Isentropic Pressure Ratio.
                self.PR_isen_sup

        Raises
        ------
        TypeError
            If Exit to Throat Area is not yet set and if gamma is not defined.

        Returns
        -------
        None.

        '''
        if self.Ae_At == None:
            raise TypeError("Exit to Throat Area Ratio not defined.")
        if self.gamma == None:
            raise TypeError("Specific Heat Ratio (γ) not defined.")
        
        M_isen_sub, M_isen_sup = GD.Mach_at_A(self.Ae_At, Gamma=self.gamma)
        self.PR_isen_sup = 1/GD.Po_P_ratio(M_isen_sup, Gamma=self.gamma) # Pe/Po
        self.PR_isen_sub = 1/GD.Po_P_ratio(M_isen_sub, Gamma=self.gamma) # Pe/Po
        self.PR_exit_shock = GD.P2_P1_n(M_isen_sup, self.gamma)*(1/GD.Po_P_ratio(M_isen_sup, self.gamma))
        self.PR_crit = GD.Pcrit_Po(self.gamma)
        
    def find_shock_locations(self, BackPressureRatios):
        '''
        Finds the area ratios of the shock to physical throat based
        on the given backpressures. If no normal shock occurs within 
        the nozzle at a given backpressure, the area ratio will be None.

        Parameters
        ----------
        BackPressureRatios : 1D Array-like or Float
            A list/np array of back pressure ratios (Pb/Po) (Both in absolute Pressure).

        Returns
        -------
        A numpy array the same size of input array that contains the corresponding area ratios (As/At).

        '''
        PR_bs =  np.array(BackPressureRatios)
        As_Ats = np.empty(PR_bs.shape)
        for iPb, Pb_Po in enumerate(PR_bs):
            if (Pb_Po < self.PR_isen_sub) and (Pb_Po > self.PR_exit_shock):
                # Flow is choked and in between sub isen solutions and exit shock -> shock in noz
                # Find Location of shock
                As_Ats[iPb] = GD.As_At_n(self.Ae_At, Pb_Po, self.gamma)
            else:
                As_Ats[iPb] = None
        self.BackPressureRatios = PR_bs
        self.ShockAreaRatios = As_Ats
        return As_Ats
    
    def calculate_pressure_through_nozzle(self, AreaRatios=None, BackPressureRatios=None,tol=1e-5): 
        if BackPressureRatios == None and self.BackPressureRatios == None:
            raise TypeError("No back pressures provided or defined")
        if BackPressureRatios == None:
            PR_bs = self.BackPressureRatios
        else:
            # Inputted Pb/Pos take priority over stored values
            PR_bs = BackPressureRatios
        
        A_rats = None #''' NEEDS RESOLVED''' Either a given array or auto generated from Ai/At/Ae
        As_Ats = [None] #''' NEEDS RESOLVED''' Check  if shock locs ran yet
        
        
        # Will be len(A_rats) rows by len(PR_bs) cols
        Px_Pos = np.empty([len(A_rats),len(PR_bs)],dtype=float) 
        BeforeThroat = True     # Define toggle for before and after throat
        descriptions = ['' for i in PR_bs]
        
        for iA, Ax_At in enumerate(A_rats):
            if abs(Ax_At - 1) < tol:
                BeforeThroat=False # Toggle swapped so now everything is after throat
            
            # Cycle through the back pressure ratios for after throat cases to find pressure
            for iPb, Pb_Po in enumerate(PR_bs):
                if Pb_Po > self.PR_isen_sub:
                    descriptions[iPb] = 'Subsonic Unchoked Flow'
                    # Not Choked Flow, find Mach number at exit plane with back pressure, (pe=pb)
                    # Subsonic Cases (Assuming Pe=Pb)
                    # Works for all Aratios when not choked
                    
                    # Find theoretical mach at exit
                    Me_the = GD.Mach_at_PR(1/Pb_Po,self.gamma)
                    # Find theoretical throat area ratio that would need to be choked
                    Ae_Astar = GD.A_ratio(Me_the,self.gamma)
                    
                    # Get ratio of Ax/A* 
                    Ax_Astar = Ax_At * (1/self.Ae_At) * Ae_Astar
                    
                    # Find Mach at x location 
                    M_sub, M_sup = GD.Mach_at_A(Ax_Astar,self.gamma)
                    
                    # Pressure Ratio at x location
                    Px_Pos[iA, iPb] = 1/GD.Po_P_ratio(M_sub,self.gamma)
                    
                elif Pb_Po > self.PR_exit_shock:
                    if (Pb_Po - self.PR_isen_sub) < tol:
                        descriptions[iPb] = 'Choked subsonic isentropic flow'
                    else: 
                        descriptions[iPb] = 'Choked flow, normal shock in nozzle'
                    
                    # choked flow, But shock occurs in nozzle (not ful expand, pb /= pe)
                    # Will follow supersonic isentropic solution until As_At is reached
                    if BeforeThroat:
                        # Need to use subsonic solution
                        M_sub, M_sup = GD.Mach_at_A(Ax_At,self.gamma)
                        # Pressure Ratio at x location
                        Px_Pos[iA, iPb] = 1/GD.Po_P_ratio(M_sub,self.gamma)
                        
                    else: 
                        # After throat
                        if Ax_At < As_Ats[iPb]: # Checks if before shock
                            # BEFORE SHOCK LOCATION 
                            # Need to use sonic solution from  area ratio
                            M_sub, M_sup = GD.Mach_at_A(Ax_At,self.gamma)
                            # Pressure Ratio at x location
                            Px_Pos[iA, iPb] = 1/GD.Po_P_ratio(M_sup,self.gamma)
                        
                        else:
                            # AFTER SHOCK LOCATION
                            # Need to get Ax_Astar after shock and use subsonic solution
                            # Get shock conditions
                            M_empty, Ms1 = GD.Mach_at_A(As_Ats[iPb],self.gamma)
                            Ms2 = GD.Mach2_n(Ms1,self.gamma)
                            Po2_Po1 = GD.Po2_Po1_n(Ms1,self.gamma)
                            
                            As_Astar = GD.A_ratio(Ms2,self.gamma) 
                            Ax_Astar = Ax_At * (1/As_Ats[iPb]) * As_Astar
                            
                            M_sub, M_empty = GD.Mach_at_A(Ax_Astar,self.gamma)
                            
                            # Pressure Ratio at x location
                            Px_Po2 = 1/GD.Po_P_ratio(M_sub,self.gamma)
                            Px_Pos[iA, iPb] = Px_Po2 * Po2_Po1
                           
                            
                elif Pb_Po > self.PR_isen_sup:
                    if (Pb_Po - self.PR_exit_shock) < tol:
                        descriptions[iPb] = 'Choked flow, normal shock at exit of nozzle'
                    else: 
                        descriptions[iPb] = 'Choked flow, overexpanded'
                  
                    # Flow is overexpanded so shock outside of nozzle
                    # Will follow the supersonic  isentropic solution up to exit plane
                    M_sub, M_sup = GD.Mach_at_A(Ax_At,self.gamma)
                    if BeforeThroat:
                        # Need to use subsonic solution
                        # Pressure Ratio at x location
                        Px_Pos[iA, iPb] = 1/GD.Po_P_ratio(M_sub,self.gamma)
                    else:
                        # Need to use sonic solution
                        # Pressure Ratio at x location
                        # if Pb_Po==PR_bs[-1]:
                        #     Px_Pos[iA, iPb] = Pb_Po
                        # else:
                        Px_Pos[iA, iPb] = 1/GD.Po_P_ratio(M_sup,self.gamma) 
                    Px_Pos[-1,iPb] = Pb_Po
                else:
                    if (Pb_Po - self.PR_isen_sup) < tol:
                        descriptions[iPb] = 'Choked supersonic isentropic flow'
                    else: 
                        descriptions[iPb] = 'Choked flow, underexpanded'
                    
                    # Choked and fully expanded flow (pe = pb)
                    # Will follow supersonic isentropic solution
                    M_sub, M_sup = GD.Mach_at_A(Ax_At,self.gamma)
                    if BeforeThroat:
                        # Need to use subsonic solution
                        # Pressure Ratio at x location
                        Px_Pos[iA, iPb] = 1/GD.Po_P_ratio(M_sub,self.gamma)
                    else:
                        # Need to use sonic solution
                        # Pressure Ratio at x location
                        Px_Pos[iA, iPb] = 1/GD.Po_P_ratio(M_sup,self.gamma)
                        
    def plot_pressure_ratios(self, PR_bs):
        plt.figure()
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