# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 15:06:55 2024

@author: cycon
"""
from GasDynamics import *
import numpy as np
import matplotlib.pyplot as plt

class FannoDuct():
    def __init__(self, **kwargs):
        self.Gamma = kwargs.get('Gamma',1.4)
        
        '''
        Process for solving duct conditions
        
        
        
        
        NOTE:How to solve for initial conditions from only the back pressure (=Pe)
        
        '''
        
        
        
        
    def stationProps(self, fL_D, fLmax_D_ref):
        fLmax_D_stat = fLmax_D_ref - fL_D 
        fanno = Fanno_Flow(fLmax_D=fLmax_D_stat)
        return [fanno['Mach'], fanno['P_Pstar'], fanno['T_Tstar']]
    
    def CalculateDuctFlow(self, Mi, Pi_Po, Ti_To, L_D_duct, f, dL_D):
        # For now, assuming Pb low enoughand L/D of duct is low enough
        # Need to add shock locations, conditional checks, etc
        
        # Setup L/D 
        L_Ds = np.arange(0,L_D_duct+dL_D, dL_D)
        P_Pos= np.empty(L_Ds.shape)
        T_Tos= np.empty(L_Ds.shape)
        Ms   = np.empty(L_Ds.shape)
        
        # Get reference (Lmax/D for duct at these inlet condition)
        MachKey = 'Mach_sup' if Mi > 1 else 'Mach'
        fanno_i = Fanno_Flow(Mach=Mi, Gamma=self.Gamma)
        fLmax_D_ref = fanno_i['fLmax_D']
        print('Ref/Max: ', fLmax_D_ref)
        Pstar_Poi = Pi_Po/fanno_i['P_Pstar']
        Tstar_Toi = Ti_To/fanno_i['T_Tstar']
        
        for j, L_D in enumerate(L_Ds):
            fLmax_D_stat = fLmax_D_ref - f*L_D 
            print(fLmax_D_stat)
            fanno_j = Fanno_Flow(fLmax_D=fLmax_D_stat, Gamma=self.Gamma)
            P_Pos[j] = fanno_j['P_Pstar']*Pstar_Poi
            T_Tos[j] = fanno_j['T_Tstar']*Tstar_Toi
            Ms[j]    = fanno_j[MachKey]
            
        return L_Ds, P_Pos, T_Tos, Ms
    
    def PlotDuctFlow(self, Mi, Pi_Po, Ti_To, L_D_duct, f, dL_D):
        L_Ds, P_Pos, T_Tos, Ms = self.CalculateDuctFlow(Mi, Pi_Po, Ti_To, L_D_duct, f, dL_D)
        plt.figure()
        plt.plot(L_Ds, Ms,label='Mach')
        plt.plot(L_Ds, P_Pos,label='P/Po')
        plt.plot(L_Ds, T_Tos,label='T/To')
        plt.legend()
        plt.xlabel('L/D')
        plt.ylabel('M, PR, TR')
        plt.grid()
        plt.show()
        
    def CalculateDuctFlow_Shock(self, Mi, Pi_Po, Ti_To, L_D_duct, f, dL_D, dS):
        # For now, assuming Pb low enoughand L/D of duct is low enough
        # Need to add shock locations, conditional checks, etc
        
        # Setup L/D 
        L_Ds = np.arange(0,L_D_duct+dL_D,dL_D,)
        L_Dsh = L_D_duct*dS
        tol = dL_D
        
        P_Pos= np.empty(L_Ds.shape)
        T_Tos= np.empty(L_Ds.shape)
        Ms   = np.empty(L_Ds.shape)
        
        # Get reference (Lmax/D for duct at these inlet condition)
        MachKey = 'Mach_sup' if Mi > 1 else 'Mach'
        fanno_i = Fanno_Flow(Mach=Mi, Gamma=self.Gamma)
        fLmax_D_ref = fanno_i['fLmax_D']
        print('Ref/Max: ', fLmax_D_ref)
        Pstar_Poi = Pi_Po/fanno_i['P_Pstar']
        Tstar_Toi = Ti_To/fanno_i['T_Tstar']
        AfterShock = False 
        
        for j, L_D in enumerate(L_Ds):
            if abs(L_D - L_Dsh) < tol and not AfterShock:
                # Shock at this location
                # Get props right before shock
                fLmax_D_stat = fLmax_D_ref - f*L_D 
                fanno_j = Fanno_Flow(fLmax_D=fLmax_D_stat, Gamma=self.Gamma,ForceSupersonic=MachKey=='Mach_sup')
                P1_Pstar_j = fanno_j['P_Pstar']
                T1_Tstar_j = fanno_j['T_Tstar']
                
                # Get Mach right after the shock
                Mj_s = Mach2_n(fanno_j['Mach'],Gamma=self.Gamma)
                P2_P1 = P2_P1_n(fanno_j['Mach'],Gamma=self.Gamma)
                T2_T1 = T2_T1_n(fanno_j['Mach'],Gamma=self.Gamma)
                MachKey = 'Mach'# change Mach  key to retrieve subsonic version
               
                # Find fanno flow cond. after shck
                fanno_js = Fanno_Flow(Mach=Mj_s,Gamma=self.Gamma)
                # Reassign reference fLmax/D for duct after shock
                fLmax_D_ref = fanno_js['fLmax_D']
                P2_Pstar_js = fanno_js['P_Pstar']
                T2_Tstar_js = fanno_js['T_Tstar']
                
                # Reassign Pstar_Poi since Pstar changes
                Pstar_Poi = Pstar_Poi * P1_Pstar_j * (1/P2_Pstar_js) * P2_P1
                # P*2 / Poi = (P*1/Poi) * (P1/P*1) * (P*2/P2) * (P2/P1)  
                
                # Reassign Tstar_Toi since Tstar changes
                Tstar_Toi = Tstar_Toi * T1_Tstar_j * (1/T2_Tstar_js) * T2_T1
                AfterShock = True 
            
            # Need to re-adjust L/D after shock since before shock doesnt matter
            L_D = L_D - L_Dsh if AfterShock else L_D 
            fLmax_D_stat = fLmax_D_ref - f*L_D 
            print(AfterShock, fLmax_D_stat)
            
            if fLmax_D_stat < 0.001:
                P_Pos[j:] = np.zeros(P_Pos[j:].shape)
                T_Tos[j:] = np.zeros(P_Pos[j:].shape)
                Ms[j:]    = np.zeros(P_Pos[j:].shape)
                return L_Ds, P_Pos, T_Tos, Ms 
            
            fanno_j = Fanno_Flow(fLmax_D=fLmax_D_stat, Gamma=self.Gamma, ForceSupersonic=MachKey=='Mach_sup')
            P_Pos[j] = fanno_j['P_Pstar']*Pstar_Poi
            T_Tos[j] = fanno_j['T_Tstar']*Tstar_Toi
            Ms[j]    = fanno_j[MachKey]
        self.Lmax_D_final = fLmax_D_ref/f    
        self.L_D_after_shock = L_D_duct - L_Dsh
        return L_Ds, P_Pos, T_Tos, Ms
    
    def PlotDuctFlow_Shock(self, Mi, Pi_Po, Ti_To, L_D_duct, f, dL_D, dS):
        L_Ds, P_Pos, T_Tos, Ms = self.CalculateDuctFlow_Shock(Mi, Pi_Po, Ti_To, L_D_duct, f, dL_D, dS)
        plt.figure()
        plt.plot(L_Ds, Ms,label='Mach')
        plt.plot(L_Ds, P_Pos,label='P/Po')
        plt.plot(L_Ds, T_Tos,label='T/To')
        plt.legend()
        plt.xlabel('L/D')
        plt.ylabel('M, PR, TR')
        plt.grid()
        plt.show()
    
        
    
    
    
    
    
    
    
    
    
    
if __name__=='__main__':
    duc = FannoDuct()
    A_rat = 2
    f = 0.02
    # L_D = 15/1
    frac_max = 1.5
    
    Mi = Mach_at_A(A_rat)[1]
    Pi_Po = 1/Po_P_ratio(Mi)
    Ti_To = 1/To_T_ratio(Mi)
    
    L_D = frac_max*Fanno_Flow(Mach=Mi)['fLmax_D']/f
    # duc.PlotDuctFlow_Shock(Mi,Pi_Po, Ti_To,L_D, f,0.05,0.5)
    print('Duct: ', f*L_D)
    
    # Calculate exit conditions for sweep of shock locations
    shoc_locs = np.linspace(0,1, 200, endpoint=True)
    Pe_Pos= np.empty(shoc_locs.shape)
    Te_Tos= np.empty(shoc_locs.shape)
    Mes   = np.empty(shoc_locs.shape)
    Lmax_D_finals= np.empty(shoc_locs.shape)
    L_D_after_shock= np.empty(shoc_locs.shape)
    
    for i, dS in enumerate(shoc_locs):
        print(i)
        L_Ds, P_Pos, T_Tos, Ms = duc.CalculateDuctFlow_Shock(Mi, Pi_Po, Ti_To, L_D, f, 0.1, dS)
        Pe_Pos[i] = P_Pos[-1]
        Te_Tos[i] = T_Tos[-1]
        Mes[i]    = Ms[-1]
        Lmax_D_finals[i] = duc.Lmax_D_final
        L_D_after_shock[i] = duc.L_D_after_shock
        
        
    plt.figure()
    plt.plot(shoc_locs, Pe_Pos, label='Pe/Po')
    plt.plot(shoc_locs, Te_Tos, label='Te/To')
    plt.plot(shoc_locs, Mes, label='Me')
    plt.legend()
    plt.xlabel('Shock Location % L/D')
    plt.ylabel('Exit Conditions')
    plt.grid()
    plt.show()
    
    plt.figure()
    plt.plot(shoc_locs, Lmax_D_finals, label='Lmax/D')
    plt.plot(shoc_locs, L_D_after_shock, label='L/D after shock')
    plt.plot([0,1],[L_D,L_D],label='L/D duct')
    plt.legend()
    plt.xlabel('Shock Location % L/D')
    plt.ylabel('Exit Conditions')
    plt.grid()
    plt.show()
    
    
'''
Discovered Observations:
For supersonic flow in a duct where L/D < Lmax/D
    - The flow will remain supersonic if Pe > Pb assuming no shock occurs
        - A constant regime, resistant to ΔPb as long as Pe > Pb
    - As Pe increases, as shock will form at the exit plane first
        - A shock at the exit plane will produce the lowest Pe for this duct
    - As Pe increases further, it will reach max Pe for shock in duct
        - A shock at the inlet plane will produce the highest Pe for this duct 
        (without affecting conditions before duct)
    - As Pe increases even more, the conditions prior to the inlet will be affected
        - This could result in:
            1: A shock forming in geometry upstream (a CD nozzle for example)
            2: A mass flow reduction 
        - Determining which  happens will require  knowing upstream information
            1: If CD nozzle is upstream of duct, the shock will move from CD 
            exit plane towards the throat in a manner that will ensure Pe=Pb for
            the duct.
            2: If the Pb increases further, the nozzle will unchoke and mass flow
            rate will decrease.
            
    Therefore, we find a few regimes by looking at changing back pressure:
        1. 0 <= Pb < Pe1 found assuming no shock and supersonic flow through entire duct
            - changes in Pb  will not affect the conditions in the duct at all
            - flow is under-expanded and expansion fans will exist
            on duct outlet
            
        2. Pb = Pe1
            - flow is fully expanded, no shock or expansion waves
            
        3. Pe1 < Pb < Pe2 from a normal shock at the exit of the duct 
            - the flow is over-expanded and obliue shock waves will form at the
            end of the duct, as Pb -> Pe2, they become stronger and closer to a 
            normalshock
            
        4. Pb = Pe2
            - a normal shock exists at the exit of the duct
            - No fanno flow occurs after shock
            
        5. Pe2 < Pb < Pe3 found from assuming a normal shock at duct inlet
            - a shock will form somewhere in the duct
            - the shock will form so that Pe_duct = Pb
        
        6. Pb = Pe3
            - a normal shock exists at the inlet of the duct
            
        7. Pe3 < Pb
            - The shock moves upstream of the duct and either:
                1. Shocks in upstream geometry
                2. Unchokes the flow and decreases mass flow rate
'''
'''
Discovered Observations:
For supersonic flow in a duct where L/D < Lmax/D
 - A shock must form somewhere in the duct so that the flow 
 at the exit is sonic
 - If Pb > Pe for sonic condition, the shock will move towards
 the inlet of the duct
 - A shock at the inlet of the duct generates the highest exit
 pressure and the lowest exit mach number
 
So the flow regimes are:
    1. 0 < Pb < Pe1 found by iteration for Me = 1
        - the pressure change will not change the location of 
        the shock
    2. Pb = Pe1 
        - the shock is still at the location producing the
    shock condition at exit
    3. Pe1 < Pb < Pe2 found by assuming shock at the duct inlet
        - The shock will form somewhere between the inlet and
        the shock location to satisfy Me = 1
        - The exit Mach number will be subsonic
        - Increases in the back pressure will move the duct further
        backwards
        - Location of shock found through iteration to satisfy Pe=Pb
    4. Pb = Pe2 
        - The shock forms at the inlet of the duct
    5. Pe2 < Pb
        - The shock moves upstream of the duct and either:
            1. Shocks in upstream geometry
            2. Unchokes the flow and decreases mass flow rate

'''