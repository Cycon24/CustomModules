# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 11:34:05 2023

@author: cycon
"""
# Imports
import numpy as np
import matplotlib.pyplot as plt
import math
import RootAlgorithms as rt
import Interpolator as intrp

def shockDeflection(Mach=None,Beta_deg=None,Theta_deg=None, Gamma = 1.4): # Beta and Theta should be entered into fct as deg.
    # Mach - The Mach number of the incoming flow (freestream)
    # Theta - Wedge Angle
    # Beta - Shock Angle
    Beta = None;Theta=None;
    if Beta_deg != None: Beta = np.radians(Beta_deg)
    if Theta_deg != None: Theta = np.radians(Theta_deg)
    
    # -------------- Find Mach ------------------------------------------------------------------    
    if Mach == None and (Beta != None and Theta!= None):
        # Define functions internally to avoid the functions using globals. 
        # Allows for modularization of the overall functions
        def F_M(x):
            M = x
            theta = Theta
            beta = Beta
            gam = Gamma
            F = (2/np.tan(beta))*((M**2 * np.sin(beta)**2 - 1)/(M**2 * (gam + np.cos(2*beta)) + 2)) - np.tan(theta)
            return F
        
        def dF_M(x):
            M = x
            beta = Beta
            gam = Gamma
            dF_M = (4*(gam + 1)*M*(1/np.tan(beta)))/(M**2 * np.cos(2*beta) + gam * M**2 + 2)**2
            return dF_M
        
        # Find Mach Number
        # Root search - Returns the lower and upper bounds of a root between start and end
        low_r, high_r = rt.rootsearch(F_M, 0, 8, 0.01)

        # Newton Raphson
        root = rt.newtonRaphson(F_M,dF_M,low_r,high_r)
        
        
        return  {'Mach': root, 'Beta': np.degrees(Beta), 'Theta': np.degrees(Theta)}
    
    # -------------- Find Beta ------------------------------------------------------------------
    elif Beta == None and (Mach != None and Theta!= None):
        # Find Beta 
        def F_Beta(x):
            M = Mach
            theta = Theta
            beta = x
            gam = Gamma
            F = (2/np.tan(beta))*((M**2 * np.sin(beta)**2 - 1)/(M**2 * (gam + np.cos(2*beta)) + 2)) - np.tan(theta)
            return F

        def dF_Beta(x):
            M = Mach
            beta = x
            gam = Gamma
            dF_B = (2*(gam*(M**4)*np.cos(2*beta) + 1/np.sin(beta)**2 *\
                       ((gam + 1)*M**2 + 2) + M**4 - 4*M**2))/((M**2)*np.cos(2*beta) + gam*M**2 + 2)**2
            return dF_B
        
        # Find Beta
        # Root search
        low_r, high_r = rt.rootsearch(F_Beta, 0,2*np.pi, .01)
        
        # Newton Raphson
        root1 = rt.newtonRaphson(F_Beta,dF_Beta,low_r,high_r)
        
        # Find second shock angle
        low_r, high_r = rt.rootsearch(F_Beta, high_r, 2*np.pi, 0.01)
        
        root2 = rt.newtonRaphson(F_Beta,dF_Beta,low_r,high_r)
    
        return {'Mach': Mach, 'Weak Beta': np.degrees(root1), 'Strong Beta': np.degrees(root2),'Theta': np.degrees(Theta)}
    
    # -------------- Find Theta ------------------------------------------------------------------
    elif Theta == None and (Beta != None and Mach != None):
        # Find Theta 
        
        def F_Theta(x):
            M = Mach
            theta = x
            beta = Beta
            gam = Gamma
            F = (2/np.tan(beta))*((M**2 * np.sin(beta)**2 - 1)/(M**2 * (gam + np.cos(2*beta)) + 2)) - np.tan(theta)
            return F

        def dF_Theta(x):
            theta = x
            dF_T = -(1/np.cos(theta))**2
            return dF_T
        
        # Find Theta
        # Root search
        low_r, high_r = rt.rootsearch(F_Theta,  0,2*np.pi, .01)
        
        # Newton Raphson
        root = rt.newtonRaphson(F_Theta,dF_Theta,low_r,high_r)
        
        return {'Mach': Mach, 'Beta': np.degrees(Beta), 'Theta': np.degrees(root)}
    
    # -------------- Unsolvable  ------------------------------------------------------------------
    else:
        # More than one entry is None, shockDeflection unsolvable
        print('\nError: Two values needed to solve the equation, less than two were supplied.\n')
        return None


def plotShockDeflection(M=None):
    if M == None:
        M = [0]
        k = 0
        incr = 0.01
        step = incr
        while M[k] <=15:
            M.append(M[k]+step)
            step += incr + incr*k/3
            k += 1
            
    if type(M) != list and type(M) != np.ndarray:
        M= np.array([M])

    M = np.array(M)
    B = np.arange(1,90,0.5)
    
    Matrix = np.zeros((np.size(B),2,np.size(M)))
    
    # Solve for theta for all of the Mach numbers from Beta range 0-90
    for i, m in enumerate(M):
        for j, b in enumerate(B):
            res = shockDeflection(m,b,None)
            if res['Theta'] < 80:
                Matrix[j,0,i] = res['Theta']
            elif Matrix[j-1,0,i] != 0:
                Matrix[j,0,i] = 0
            else:
                Matrix[j,0,i] = None
            Matrix[j,1,i] = b
    
    # Plotting
    plt.plot(Matrix[:,0,:],Matrix[:,1,:])
    plt.title('Shock Deflection Angle')
    plt.xlabel('Theta (°)')
    plt.ylabel('Beta (°)')
    
    return None


# _______ Isentropic Exapnsion/Compression _______________
# To/T
def To_T_ratio(Mach,Gamma=1.4):
    return 1 + (Gamma-1)*(Mach**2)/2

# Po/P
def Po_P_ratio(Mach, Gamma =1.4):
    return (1 + (Gamma-1)*(Mach**2)/2)**(Gamma/(Gamma-1))

# Po/P
def Po_P_ratio2(To, T, Gamma=1.4):
    return (To/T)**(Gamma/(Gamma-1))

def RHO_o_Rho(Mach, Gamma=1.4):
    return (1 + (Gamma-1)*(Mach**2)/2)**(1/(Gamma-1))


def mdot(Po, To, A, Mach=1, Gamma=1.4, R=287):
    first = (Po*A/(np.sqrt(R*To))) *Mach*np.sqrt(Gamma)
    second = 1 + ((Gamma-1)/2)*(Mach**2)
    power = -(Gamma+1)/(2*(Gamma-1))
    return first*np.power(second, power)

# _______ Prendlt Meyer Expansion Fan ?? _________________
# Define Prendtl-Meyer Function to pass into approximations
def nu_PM(Mach, Gamma=1.4):
    return np.sqrt((Gamma+1)/(Gamma-1)) * \
            np.arctan(np.sqrt((Gamma-1)*(np.power(Mach,2) - 1)/(Gamma+1))) - \
            np.arctan(np.sqrt(np.power(Mach,2) - 1))

# Exact Solutions
def dv_PM(Mach, Gamma=1.4):
    Num = 2*np.sqrt(np.power(Mach,2) -1)
    Den = Mach*(Gamma*np.power(Mach,2) - np.power(Mach,2) + 2)
    return np.divide(Num,Den)

# Full Prendtly-Meyer Solver to get Mach from input
def ExpansionFan(Mach1, Theta, Gamma=1.4):
    def PM(Mach2):
        return nu_PM(Mach2) - nu_PM(Mach1) - np.radians(Theta)
    
    # Root Search
    low_r, high_r = rt.rootsearch(PM, Mach1, 20, 0.01)
    # Newton Raphson
    M2 = rt.newtonRaphson(PM, dv_PM, low_r,high_r)
    return M2



# ________________________________________________________

# _______ Normal Shock Relations Properties ______________
# velocity ratio across normal shock: v2/v1
def vR_n(Mach1, Mach2, Temp1, Temp2):
    # Temperatures are static, in absolute units
    return (Mach2/Mach1)*np.sqrt(Temp2/Temp1)

# density ratio across normal shock: rho2/rho1
def rhoR_n(Mach1, Mach2, Temp1, Temp2):
    # Temperatures are static, in absolute units
    return (Mach1/Mach2)*np.sqrt(Temp1/Temp2)

# static temperature ratio across normal shock: T2/T1
def TR_n(Mach1, Mach2, Gamma=1.4):
    return (1+(Mach1**2)*(Gamma-1)/2)/(1+(Mach2**2)*(Gamma-1)/2)

# static pressure ratio across normal shock: p2/p1
def PR_n(Mach1, Mach2, Gamma=1.4):
    return (1+Gamma*Mach1**2)/(1+Gamma*Mach2**2)

# mach number after a shock wave
def Mach2_n(Mach1, Gamma=1.4):
    num = Mach1**2 + 2/(Gamma-1)
    den = -1 + (2*Gamma/(Gamma-1))*Mach1**2   
    return np.sqrt(num/den)

# ---- Stagnation Properties ----
# stagnation pressure ratio across normal shock: po2/po1
def stag_PR_n(Pres1, Pres2, Temp1, Temp2, Gamma=1.4):
    return (Pres2/Pres1)*(Temp1/Temp2)**(Gamma/(Gamma-1))

# stagnation density ratio across normal shock: ρo2/ρo1
def stag_rhoR_n(stagPres1, stagPres2):
    return stagPres1/stagPres2

# throat area across normal shock? A2*/A1*
def throatAreaR_n(stagPres1, stagPres2):
    return stagPres1/stagPres2

#_________________________________________________________

# ____________ Turning Flow - Oblique Shock ______________

# Notes: Oblique shock at angle Beta relative to Vinf.
#   Resultant flow will be at angle Theta.
#   M1n = M1*sin(Beta)        M1t = M1*cos(Beta)
#   M2n = M2*sin(Beta-Theta)  M2t = M2*cos(Beta-Theta)
#   Tangential velocity does not change across shock (Mt does)
#   Normal velocity changes same as machΔ across Normal shock
#   To find property changes across oblique shock, use M1n and M2n
#   - Stag Temperature stays constant

def Mach2_ob(Mach1, Beta, Theta, Gamma=1.4):
    B_rad = np.radians(Beta)
    th_rad = np.radians(Theta)
    
    num = (np.sin(B_rad)**2)*(Mach1**2) + 2/(Gamma-1)
    den = (2*Gamma/(Gamma-1))*(np.sin(B_rad)**2)*(Mach1**2) -1
    
    m2 = np.sqrt(num/(den*np.sin(B_rad - th_rad)))/np.sin(B_rad-th_rad)
    return m2
    
def Mach_ob_to_n(Mach1, Beta, Theta):
    '''

    Parameters
    ----------
    Mach1 : Float
        Mach number of flow upstream of oblique shock.
    Beta : Float
        In Degrees, the oblique shock angle.
    Theta : Float
        In Degrees, the turn angle of the flow (due to geometry).

    Returns
    -------
    m1n : Float
        Normal component of upstream flow relative to shock.
    m1t : Float
        Tangential component of upstream flow relative to shock.
    m2n : Float
        Normal component of downstream flow relative to shock.
    m2t : Float
        Tangential component of downstream flow relative to shock.

    '''
    m1n = Mach1*np.sin(np.radians(Beta))
    m2n = Mach2_n(m1n)
    Mach2 = m2n/np.sin(np.radians(Beta-Theta))
    m1t = Mach1*np.cos(np.radians(Beta))
    m2t = Mach2*np.cos(np.radians(Beta-Theta))
    return m1n, m1t, m2n, m2t

def Mach_Tang(Mach,Beta):
    return Mach*np.sin(np.radians(Beta))

# After swapping from oblique to normal Machs, can use them
#   to find the properties after the oblique shock


# __________ Nozzle Functions _____________

def Mach_at_A(Ai_At, Gamma=1.4):
    def FA_ratio(M):
        return -Ai_At + (1/M)*((1 + M**2*(Gamma-1)/2)/((Gamma+1)/2))**((Gamma+1)/(2*Gamma-2))
    
    def dAR_M(M):
        first_num = - (2 + (M**2)*(Gamma-1))**((Gamma+1)/(2*Gamma-2))
        first_den = (M**2)*(Gamma+1)**((Gamma+1)/(2*Gamma-2))
        second = ((2 + (M**2)*(Gamma-1))/(Gamma+1))**((-Gamma+3)/(2*Gamma-2))
        return first_num/first_den + second
    
    # Solve for first root (subsonic solution)
    low_r, high_r = rt.rootsearch(FA_ratio, 0.001, 8, 0.01)
    Msub = rt.newtonRaphson(FA_ratio,dAR_M,low_r,high_r)
    
    # Solve for second root (sonic solution)
    low_r, high_r = rt.rootsearch(FA_ratio, Msub+0.01, 15, 0.01)
    Msup = rt.newtonRaphson(FA_ratio,dAR_M, low_r,high_r)
    
    return Msub, Msup

# Area ratio: A/A*
def A_ratio(M, Gamma=1.4):
    return (1/M)*((1 + M**2*(Gamma-1)/2)/((Gamma+1)/2))**((Gamma+1)/(2*Gamma-2))


def A_throat(m_dot, Po, To, Gamma=1.4, R=287.1):
    return m_dot/((Po/np.sqrt(R*To))*np.sqrt(Gamma)*np.power(1+(Gamma-1)/2, (Gamma+1)/(2-2*Gamma)))


def Mach_at_PR(Po_P, Gamma=1.4):
    return np.sqrt(((Po_P)**((Gamma-1)/Gamma) - 1)*2/(Gamma-1))
# __________________________________________




if __name__ == "__main__":
    # tests = np.array([[2,None,20]])
    
    # for i in range(0,1):
    #     try:
    #         res = shockDeflection(tests[i,0],tests[i,1],tests[i,2])
    #         print('Test {}:'.format(i))
    #         for key, value in res.items():
    #             print("{:>15s} = {:8.4f}".format(key,value))
    #         print('\n')
    #     except:
    #         print('No Solution')
    print(shockDeflection(3.2,None,17))
    mach2n = 0.630
    m2 = mach2n/np.sin(np.radians(32.98-17))
    print(shockDeflection(m2,None,17)['Weak Beta'] - 17)
    m3n = m2*np.sin(np.radians(42.2)) # Normal component: ??? look at cass notes
    m3 = m3n*np.sin(np.radians(25.2)) # normal component time sin(Beta_reflected)
    print(m3)
    
    plotShockDeflection()

   