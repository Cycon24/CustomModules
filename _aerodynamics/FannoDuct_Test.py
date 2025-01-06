# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 11:03:05 2024

@author: cycon
"""

from GasDynamics import *


# Get Inlet mach number and actual fL/D 
# - Determine outlet conditions
# - Determine if a shock occurs and where


fL_D = 20*0.02/1
A = np.pi*(0.01/2)**2 # m^2
Gamma = 1.4
Mi = Mach_at_A(2, Gamma)[1]
Pb = 0
Po = 650e3
To = 1000
Pb_Po  = Pb/Po
tol = 2e-5
R=287

#def Fanno_Duct(FL_D, Mi, Po, Pb?, tol=1e-6):
fanno_i = Fanno_Flow(Mach=Mi, Gamma=Gamma)
Pi_Pstar = fanno_i['P_Pstar']
Ti_Tstar = fanno_i['T_Tstar']
fLmax_D_i = fanno_i['fLmax_D']
Pi_Po = Isentropic_Flow(Mach=Mi, Gamma=Gamma)['P_Po']
Ti_To = Isentropic_Flow(Mach=Mi, Gamma=Gamma)['T_To']

'''
Convention
Inlet      b4 shock { after shock                   exit    back(atm)
_________________________________________________________
                    {
i                1  {  2                                e b
                    {
____________________{____________________________________

'''
for i in range(0,1):
    if fL_D > fLmax_D_i:
        if Mi >= 1:
            # A shock forms in the duct such that Pe = Pb,or if duct is choked
            # We need to guess where fL_D_s is (where shock is located at L_s)
            # then find the conditions before the shock and after the shock
            # Initial Guess
            Len_frac = 0.01
            k = -0.1
            k_min = 1e-6
            error = -1
            while abs(error) > tol:
                fL_D_i1 = Len_frac*fLmax_D_i
                fLmax_D_1 = fLmax_D_i - fL_D_i1
                
                # Station 1: Conditions before shock
                fanno_1 = Fanno_Flow(fLmax_D=fLmax_D_1, Gamma=Gamma)
                M1 = fanno_1['Mach_sup']
                P1_Pstar = fanno_1['P_Pstar']
                T1_Tstar = fanno_1['T_Tstar']
                
                # Station 2: Conditions after shock
                M2 = Mach2_n(M1,Gamma=Gamma)
                P2_P1 = P2_P1_n(M1, Gamma=Gamma)
                T2_T1 = T2_T1_n(M1, Gamma=Gamma)
                
                fanno_2 = Fanno_Flow(Mach=M2, Gamma=Gamma)
                fLmax_D_2 = fanno_2['fLmax_D']
                P2_Pstar2 = fanno_2['P_Pstar']
                fL_D_2e = fL_D - fL_D_i1
                
                fLmax_D_e = fLmax_D_2 - fL_D_2e
                if fLmax_D_e < 0:
                    print('Need a new guess, duct length after shock too long')
                    Len_frac += 2*k*error 
                    k /= 5 if abs(k) > k_min else 1
                # elif fLmax_D_e < tol:
                #     print('Close to choked at exit')
                #     Len_frac += k*error 
                #     k /= 10 if abs(k) > k_min else 1
                else:
                    try:
                        fanno_e = Fanno_Flow(fLmax_D=fLmax_D_e, Gamma=Gamma)
                        Pe_Pstar2 = fanno_e['P_Pstar']
                        Te_Tstar2 = fanno_e['T_Tstar']
                        Me = fanno_e['Mach']
                        Pe_Po = Pe_Pstar2 * (1/P2_Pstar2) * P2_P1 * P1_Pstar * (1/Pi_Pstar) * Pi_Po
                        Pe = Pe_Po*Po
                        error = 1 - Me
                        print(error)
                        
                    except:
                        print('Flow choked at exit')
                        Me = 1
                        Pe_Pstar2 = 1 
                        Te_Tstar2 = 1 
                        Pe_Po = Pe_Pstar2 * (1/P2_Pstar2) * P2_P1 * P1_Pstar * (1/Pi_Pstar) * Pi_Po
                        Pe = Pe_Po*Po
                        error = 1 - Me
                        print(error)
                  
                Len_frac += - k*error
                if  Len_frac > 1:
                    print('WARNING: Len Frac > 1, changing sign of k')
                    k *= -1
                
        else:
            # Subsonic flow, mass flow rate decreases
            print('Mass flow decreases')
    else:
        # Standard fanno flow
        print('standard fanno flow')
        fanno_e = Fanno_Flow(fLmax_D=fanno_i['fLmax_D']-fL_D, Gamma=Gamma)
        Me = fanno_e['Mach_sup'] if Mi > 1 else fanno_e['Mach']
        Pe_Po = fanno_e['P_Pstar']*(1/Pi_Pstar)*(Pi_Po)
        Te_To = fanno_e['T_Tstar']*(1/Ti_Tstar)*(Ti_To)
        Pe = Pe_Po * Po
        Te = Te_To * To 
        mdot = mdot_s(Pe, Te, A, Mach=Me, Gamma=Gamma, R=R)
        if Pe < Pb:
            print('Back pressure higher than the exit pressure: Shock forms in duct and must be re-evaluted')
        


'''
For L/D > Lmax/D|i, we can check if at a given shock location
if the L/D after shock (length of duct after shock) = Lmax/D|s

'''
