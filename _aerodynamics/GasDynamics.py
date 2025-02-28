# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 11:34:05 2023

@author: cycon
"""
# Imports
import sys
sys.path.append("C:\\Users\\cycon\\Documents\\Modules")
import _tools.Interpolator as intrp
import _tools.RootAlgorithms as rt
import numpy as np
import matplotlib.pyplot as plt
import math
'''
Planned to-do log:
    2-20-2025: (Unresolved) Need to incorperate root search params in IsentropicFlow function
    2-20-2025: (Unresolved) Need to add ability to find Mach from MFP function and incorp. in Isentropic FlowFunction
    2-20-2025: (Unresolved) Need to add docstrings for total Flow Functions (Fanno/Rayleigh)

Done:
    2-27-2025: (Resolved) Renamed all Normal Shock Equations to start with 'NS' to find all related functions easier
    2-27-2025: (Resolved) Made headway in adding docstrings for all functions, still need:
                                - All Fanno Flow functions
                                - Full Fanno and Rayleigh Flow docstring finish

'''


def calorically_imperfect_gas(Temp, Cv_perf, Cp_perf, gam_p):
    # Standin for now  just to save these equations/details
    Theta = 3055.556  # 5500 [Deg R] or 3055.556 [K]
    T = Temp
    e_frac = (np.exp(Theta/T)) / (np.exp(Theta/T)-1)**2
    T_frac = (Theta/T)**2
    
    Cv = Cv_perf*(1 + (gam_p - 1)*(T_frac*e_frac))
    Cp = Cp_perf*(1 + (gam_p - 1)*(T_frac*e_frac)/gam_p)
    Gamma = 1  + (gam_p-1)/(1 + (gam_p - 1)*(T_frac*e_frac))
    return Gamma, Cp, Cv

# =============================================================================
#               Shock Deflection
# =============================================================================
def shockDeflection(Mach=None, Beta_deg=None, Theta_deg=None, Gamma=1.4):
    '''
    Outputs a dictionary containing the Mach number, Beta (shock angle), and Theta (deflection angle) for a shock deflection.
    Any two iputs are needed to solve. Angles should be input in degrees

    Parameters
    ----------
    Mach : Float, optional
        Mach Number of Incoming Freestream Flow. The default is None.
    Beta_deg : Float, optional
        Shock angle [deg]. The default is None.
    Theta_deg : TYPE, optional
        Wedge angle (aka exit flow angle) [deg]. The default is None.
    Gamma : Float, optional
        Specific Heat Ratio. The default is 1.4.

    Returns
    -------
    Dictionary
        returns {'Mach'=M, 'Beta'=beta, 'Theta'=theta} in degrees.

    '''
    Beta = None
    Theta = None
    if Beta_deg != None:
        Beta = np.radians(Beta_deg)
    if Theta_deg != None:
        Theta = np.radians(Theta_deg)

    # -------------- Find Mach ------------------------------------------------------------------
    if Mach == None and (Beta != None and Theta != None):
        # Define functions internally to avoid the functions using globals.
        # Allows for modularization of the overall functions
        def F_M(x):
            M = x
            theta = Theta
            beta = Beta
            gam = Gamma
            F = (2/np.tan(beta))*((M**2 * np.sin(beta)**2 - 1) /
                                  (M**2 * (gam + np.cos(2*beta)) + 2)) - np.tan(theta)
            return F

        def dF_M(x):
            M = x
            beta = Beta
            gam = Gamma
            dF_M = (4*(gam + 1)*M*(1/np.tan(beta))) / \
                (M**2 * np.cos(2*beta) + gam * M**2 + 2)**2
            return dF_M

        # Find Mach Number
        # Root search - Returns the lower and upper bounds of a root between start and end
        low_r, high_r = rt.rootsearch(F_M, 0, 8, 0.01)

        # Newton Raphson
        root = rt.newtonRaphson(F_M, dF_M, low_r, high_r)

        return {'Mach': root, 'Beta': np.degrees(Beta), 'Theta': np.degrees(Theta)}

    # -------------- Find Beta ------------------------------------------------------------------
    elif Beta == None and (Mach != None and Theta != None):
        # Find Beta
        def F_Beta(x):
            M = Mach
            theta = Theta
            beta = x
            gam = Gamma
            F = (2/np.tan(beta))*((M**2 * np.sin(beta)**2 - 1) /
                                  (M**2 * (gam + np.cos(2*beta)) + 2)) - np.tan(theta)
            return F

        def dF_Beta(x):
            M = Mach
            beta = x
            gam = Gamma
            dF_B = (2*(gam*(M**4)*np.cos(2*beta) + 1/np.sin(beta)**2 *
                       ((gam + 1)*M**2 + 2) + M**4 - 4*M**2))/((M**2)*np.cos(2*beta) + gam*M**2 + 2)**2
            return dF_B

        # Find Beta
        # Root search
        low_r, high_r = rt.rootsearch(F_Beta, 0, 2*np.pi, .01)

        # Newton Raphson
        root1 = rt.newtonRaphson(F_Beta, dF_Beta, low_r, high_r)

        # Find second shock angle
        low_r, high_r = rt.rootsearch(F_Beta, high_r, 2*np.pi, 0.01)

        root2 = rt.newtonRaphson(F_Beta, dF_Beta, low_r, high_r)

        return {'Mach': Mach, 'Weak Beta': np.degrees(root1), 'Strong Beta': np.degrees(root2), 'Theta': np.degrees(Theta)}

    # -------------- Find Theta ------------------------------------------------------------------
    elif Theta == None and (Beta != None and Mach != None):
        # Find Theta

        def F_Theta(x):
            M = Mach
            theta = x
            beta = Beta
            gam = Gamma
            F = (2/np.tan(beta))*((M**2 * np.sin(beta)**2 - 1) /
                                  (M**2 * (gam + np.cos(2*beta)) + 2)) - np.tan(theta)
            return F

        def dF_Theta(x):
            theta = x
            dF_T = -(1/np.cos(theta))**2
            return dF_T

        # Find Theta
        # Root search
        low_r, high_r = rt.rootsearch(F_Theta,  0, 2*np.pi, .01)

        # Newton Raphson
        root = rt.newtonRaphson(F_Theta, dF_Theta, low_r, high_r)

        return {'Mach': Mach, 'Beta': np.degrees(Beta), 'Theta': np.degrees(root)}

    # -------------- Unsolvable  ------------------------------------------------------------------
    else:
        # More than one entry is None, shockDeflection unsolvable
        print('\nError: Two values needed to solve the equation, less than two were supplied.\n')
        return None


def plotShockDeflection(M=None):
    '''
    Plots the shock deflection map at a given Mach number M
    or for a range of Mach numbers from 0 to 15 (if M=None)

    Parameters
    ----------
    M : Float, optional
        Mach Number of Freestream Flow. The default is None.

    Returns
    -------
    Matplotlib Figure.

    '''
    if M == None:
        M = [0]
        k = 0
        incr = 0.01
        step = incr
        while M[k] <= 15:
            M.append(M[k]+step)
            step += incr + incr*k/3
            k += 1
    print(M)
    if type(M) != list and type(M) != np.ndarray:
        M = np.array([M])

    M = np.array(M)
    B = np.arange(1, 90, 0.5)

    Matrix = np.zeros((np.size(B), 2, np.size(M)))

    # Solve for theta for all of the Mach numbers from Beta range 0-90
    for i, m in enumerate(M):
        for j, b in enumerate(B):
            res = shockDeflection(m, b, None)
            if res['Theta'] < 80:
                Matrix[j, 0, i] = res['Theta']
            elif Matrix[j-1, 0, i] != 0:
                Matrix[j, 0, i] = 0
            else:
                Matrix[j, 0, i] = None
            Matrix[j, 1, i] = b

    # Plotting
    fig = plt.figure()
    plt.plot(Matrix[:, 0, :], Matrix[:, 1, :])
    plt.title('Shock Deflection Angle')
    plt.xlabel('Theta (°)')
    plt.ylabel('Beta (°)')
    plt.grid()
    return fig


# =============================================================================
#               Isentropic Exapnsion/Compression
# =============================================================================
# To/T
def To_T_ratio(Mach, Gamma=1.4):
    '''
    Calculates the Stagnation to Static Temperature ratio at a Mach number
    To/T

    Parameters
    ----------
    Mach : Float
        Mach of gas flow.
    Gamma : Float, optional
        Specific Heat Ratio. The default is 1.4.

    Returns
    -------
    Float
        To/T [NonDim].

    '''
    return 1 + (Gamma-1)*(Mach**2)/2

# Po/P


def Po_P_ratio(Mach, Gamma=1.4):
    '''
    Calculates the Stagnation to Static Pressure ratio at a Mach number
    Po/P

    Parameters
    ----------
    Mach : Float
        Mach of gas flow.
    Gamma : Float, optional
        Specific Heat Ratio. The default is 1.4.

    Returns
    -------
    Float
        Po/P [NonDim].

    '''
    return (1 + (Gamma-1)*(Mach**2)/2)**(Gamma/(Gamma-1))

# Po/P


def Po_P_ratio2(To, T, Gamma=1.4):
    '''
    Calculates the Stagnation to Static Pressure ratio at a given stagnation and static Temperature 
    Po/P = f(To, T)

    Parameters
    ----------
    Mach : Float
        Mach of gas flow.
    Gamma : Float, optional
        Specific Heat Ratio. The default is 1.4.

    Returns
    -------
    Float
        Po/P [NonDim].

    '''
    return (To/T)**(Gamma/(Gamma-1))


def RHO_o_Rho(Mach, Gamma=1.4):
    '''
    Calculates the Stagnation to Static Density ratio at a Mach number
    Rho_o/Rho

    Parameters
    ----------
    Mach : Float
        Mach of gas flow.
    Gamma : Float, optional
        Specific Heat Ratio. The default is 1.4.

    Returns
    -------
    Float
        Rho_o/Rho [NonDim] (Stag/Static).

    '''
    return (1 + (Gamma-1)*(Mach**2)/2)**(1/(Gamma-1))


def mdot(Po, To, A, Mach=1, Gamma=1.4, R=287, gc=1):
    '''
    Calculates the mass flow rate from stagnation conditions

    Parameters
    ----------
    Po : Float
        Stagnation Pressure [Pa] or [lbf/ft^2].
    To : Float
        Stagnation Temperature [K] or [R].
    A : Float
        Area of channel [m^2] or [ft^2].
    Mach : Float, optional
        Mach of gas. The default is 1.
    Gamma : Float, optional
        DESCRIPTION. The default is 1.4.
    R : Float, optional
        Gas Constant [J/kg*K]. The default is 287.
        Or [ft*lbf/R*lbm], 53.535 for air.
    gc : Float, optional
        Unit conversion parameter. The default is 1 [kg·m/N·s^2]
        For English Units: gc = 32.174 [lbm·ft/lbf·s^2]
    Returns
    -------
    Float
        mdot [kg/s].

    '''
    first = (Po*A/(np.sqrt(R*To/gc))) * Mach*np.sqrt(Gamma)
    second = 1 + ((Gamma-1)/2)*(Mach**2)
    power = -(Gamma+1)/(2*(Gamma-1))
    return first*np.power(second, power)


def Pcrit_Po(Gamma=1.4):
    '''
    Calculates the critical pressure ratio (static/total) for the sonic condition.

    Parameters
    ----------
    Gamma : TYPE, optional
        Specific Heat Ratio of gas. The default is 1.4 for air.

    Returns
    -------
    Float
        Ratio P_crit/Po of the gas at the sonic condition.

    '''
    return (1 + 0.5*(Gamma-1))**(-Gamma/(Gamma-1))

def Mach_at_PR(Po_P, Gamma=1.4, tol=1e-9):
    '''
    Calculates the Mach number of isentropic flow
    At a pressure ratio of Po/P

    Parameters
    ----------
    Po_P : Float
        Pressure ratio.
    Gamma : Float, optional
        Gas dependent. The default is 1.4 (air).

    Returns
    -------
    Float
        Mach number at the point where P = Po / (Po/P).

    '''
    if abs(Po_P - 1) < tol:
        return 0
    else:   
        return np.sqrt(((Po_P)**((Gamma-1)/Gamma) - 1)*2/(Gamma-1))

def Mach_at_TR(To_T, Gamma=1.4, tol=1e-9):
    '''
    Calculates the Mach number of isentropic flow
    At a pressure ratio of To/T

    Parameters
    ----------
    Po_P : Float
        Pressure ratio.
    Gamma : Float, optional
        Gas dependent. The default is 1.4 (air).

    Returns
    -------
    Float
        Mach number at the point where P = Po / (Po/P).

    '''
    if abs(To_T - 1) < tol:
        return 0
    else:   
        return np.sqrt(((To_T) - 1)*2/(Gamma-1))
# =============================================================================
#                Nozzle Relations (Isentropic Flow)
# =============================================================================

def Mach_at_A(Ai_At, Gamma=1.4):
    '''
    Calcultes the MachNumberat the area ratio A/A*
    Returns subsonic and supersonic results

    Parameters
    ----------
    Ai_At : Float
        Area ratio: A/A*.
    Gamma : Float, optional
        Gas Dependent. The default is 1.4 (air).

    Returns
    -------
    M_subsonic, M_supersonic : Float
        Subsonic and Supersonic Mach number solution.

    '''
    if (Ai_At - 1.0) < 1e-5:
        return 1.0, 1.0

    def FA_ratio(M):
        return -Ai_At + (1/M)*((1 + M**2*(Gamma-1)/2)/((Gamma+1)/2))**((Gamma+1)/(2*Gamma-2))

    def dAR_M(M):
        first_num = - (2 + (M**2)*(Gamma-1))**((Gamma+1)/(2*Gamma-2))
        first_den = (M**2)*(Gamma+1)**((Gamma+1)/(2*Gamma-2))
        second = ((2 + (M**2)*(Gamma-1))/(Gamma+1))**((-Gamma+3)/(2*Gamma-2))
        return first_num/first_den + second

    # Solve for first root (subsonic solution)
    low_r, high_r = rt.rootsearch(FA_ratio, 0.001, 8, 0.01)
    Msub = rt.newtonRaphson(FA_ratio, dAR_M, low_r, high_r)

    # Solve for second root (sonic solution)
    low_r, high_r = rt.rootsearch(FA_ratio, Msub+0.01, 15, 0.01)
    Msup = rt.newtonRaphson(FA_ratio, dAR_M, low_r, high_r)

    return Msub, Msup

# Area ratio: A/A*


def A_ratio(M, Gamma=1.4):
    '''
    Calculates the ratio of Area / Throat Area from Mach number and Gamma

    Parameters
    ----------
    M : Float
        Mach Number.
    Gamma : Float, optional
        The default is 1.4.

    Returns
    -------
    A_At : Float
        A/A*: Area ratio.

    '''
    A_At = (1/M)*((1 + M**2*(Gamma-1)/2) /
                  ((Gamma+1)/2))**((Gamma+1)/(2*Gamma-2))
    return A_At


def A_throat(m_dot, Po, To, Gamma=1.4, R=287.1):
    '''
    Calculates the Throat Area [m^2] for choked flow.

    Parameters
    ----------
    m_dot : Float
        Mass flow rate in kg/s.
    Po : Float
        Total Pressure in Pa.
    To : Float
        Total Temperature in K.
    Gamma : Float, optional
        Gas dependent. The default is 1.4 (air).
    R : Float, optional
        Gas Constant J/kg*K. The default is 287.1 (air).

    Returns
    -------
    Float
        Area of the throat of choked flow.

    '''
    return m_dot/((Po/np.sqrt(R*To))*np.sqrt(Gamma)*np.power(1+(Gamma-1)/2, (Gamma+1)/(2-2*Gamma)))



def MassFlowParam_norm(Mach, Gamma=1.4):
    '''
    Calculates the mass flow parameter from the Mach number and spec. heat ratio.
    Returns the value:
        MFP*sqrt(R/gc) [NonDim] = mdot * sqrt(To) / (Po * A)

    Parameters
    ----------
    Mach : Float
        Mach number of flow.
    Gamma : Float, optional
        Specific heat ratio of the gas. The default is 1.4.

    Returns
    -------
    MFPsqrtRgc : Float
        Normalized Mass Flow Parameter.

    '''
    MFPsqrtRgc = Mach*np.sqrt(Gamma)*(1 + (Mach**2)*(Gamma-1)/2)**(-(Gamma+1)/(2*Gamma-2))
    
    return MFPsqrtRgc


# _________________________________________________________


def Isentropic_Flow(**kwargs):
    '''
    Calculates all related ratios/values for isentropic compressible flow
    from a single input (Currently does not support input of MFP*sqrt(R/gc)).

    Parameters
    ----------
    **kwargs : Dictionary of inputs
        Gas Properties:
            'Gamma' - Float, default of 1.4 for air
        Flow Properties (enter one):
            'Mach' - Mach number of flow
            'T_To' - Static to Stagnation Temperature ratio
            'P_Po' - Static to Stagnation Pressure ratio
            'A_At' - Current Area to Throat Area ratio
        Function Properties:
            'IntegralPoints' - Defines the number of points used for root solves (not used currently)
            'root_dx' - Step size of dx (Mach) for root solves (not used currently)
            'root_max_search' - Maximum x (Mach) for root solves (not used currently)
            

    Raises
    ------
    ValueError
        Raised when not enough inputs to solve flow equations.

    Returns
    -------
    dict
        Contains:
            'Mach': Mach number 
            'T_To': Static to Stagnation Temperature Ratio
            'P_Po': Static to Stagnation Pressure Ratio
            'A_At': Current Cross-section Area to Throat Area Ratio
            'MFPsqrtR_gc': Normalized Mass Flow Parameter (MFP*sqrt(R/g_c))

    '''
    Mach = kwargs.get('Mach',None)
    T_To = kwargs.get('T_To',None)
    P_Po = kwargs.get('P_Po',None)
    A_At = kwargs.get('A_At',None)
    Gamma = kwargs.get('Gamma',1.4)
    # IntegralPoints = kwargs.get('IntegralPoints',1000)
    root_dx = kwargs.get('root_dx',0.1)
    root_max_search = kwargs.get('root_max_search',8)
    
    if Mach == None: 
        Props = [T_To, P_Po, A_At]
        PropsF = [lambda x: Mach_at_TR(x,Gamma), lambda x: Mach_at_PR(x,Gamma), lambda x: Mach_at_A(x,Gamma)]
        MachFound = False
        for i,  prop in enumerate(Props):
            if prop != None:
                # Perfom root search to get mach number
                Mach = PropsF[i](prop)
                MachFound=True
        if not MachFound:
            raise ValueError('Isentropic Flow: Insufficient inputs')
            return None
        
    T_To = 1/To_T_ratio(Mach,Gamma)
    P_Po = 1/Po_P_ratio(Mach,Gamma)
    A_At = A_ratio(Mach,Gamma)
    MFPsqrtR_gc = MassFlowParam_norm(Mach, Gamma)

    return {'Mach':Mach, 'T_To':T_To,'P_Po':P_Po,'A_At':A_At, 'MFPsqrtR_gc':MFPsqrtR_gc}
    

def mdot_s(P, T, A, Mach=1, Gamma=1.4, R=287):
    '''
    Calculates the mass flow from the static properties

    Parameters
    ----------
    P : Float
        Static Pressure [Pa].
    T : Flaot
        Static Temperature [Pa].
    A : Float
        Area of channel [m^2].
    Mach : Float, optional
        Mach of airfow. The default is 1.
    Gamma : Float, optional
        DESCRIPTION. The default is 1.4.
    R : Float, optional
        Gas constant [J/kg*K]. The default is 287.

    Returns
    -------
    Float
        mass flow rate [kg/s].

    '''
    return (P/(R*T))*A*Mach*np.sqrt(Gamma*R*T)
# =============================================================================
#               Prendlt Meyer Expansion Fan
# =============================================================================
# Define Prendtl-Meyer Function to pass into approximations

def mu_MachWave(Mach):
    '''
    Calculates μ which is the mach wave angle for supersonic flow.

    Parameters
    ----------
    Mach : Float
        Mach number >= 1.

    Returns
    -------
    mu : Float
        Angle μ in radians.

    '''
    mu = np.arcsin(1/Mach)
    return mu 

def nu_PM(Mach, Gamma=1.4):
    '''
    Calculates the angle ν (nu) for supersonic flow that is required for sonic flow to turn
    in order to reach the inputted Mach number

    Parameters
    ----------
    Mach : Float
        Mach number of flow.
    Gamma : Float, optional
        Specific Heat Ratio of the gas. The default is 1.4 for air.

    Returns
    -------
    nu : Float
        The angle sonic flow has turned to reach the inputted Mach number in Radians.

    '''
    nu = np.sqrt((Gamma+1)/(Gamma-1)) * \
        np.arctan(np.sqrt((Gamma-1)*(np.power(Mach, 2) - 1)/(Gamma+1))) - \
        np.arctan(np.sqrt(np.power(Mach, 2) - 1))
    return nu


# Full Prendtly-Meyer Solver to get Mach from input


def Expansion_Fan(Mach1, Theta, Gamma=1.4):
    '''
    Calculates the Mach after expansion through a turn of Theta degrees.

    Parameters
    ----------
    Mach1 : Float
        Mach number of flow before expansion wave.
    Theta : Float
        The turn angle in degrees of the flow (Angle between M1 and M2).
    Gamma : Float, optional
        Specific heat ratio. The default is 1.4.

    Returns
    -------
    Mach2: Float
        Mach number of flow after expansion fan.

    '''
    def PM(Mach2):
        return nu_PM(Mach2) - nu_PM(Mach1) - np.radians(Theta)

    def dv_PM(Mach, Gamma=1.4):
        Num = 2*np.sqrt(np.power(Mach, 2) - 1)
        Den = Mach*(Gamma*np.power(Mach, 2) - np.power(Mach, 2) + 2)
        return np.divide(Num, Den)
    
    # Root Search
    low_r, high_r = rt.rootsearch(PM, Mach1, 20, 0.01)
    # Newton Raphson
    M2 = rt.newtonRaphson(PM, dv_PM, low_r, high_r)
    return M2

# ________________________________________________________


# =============================================================================
#            Normal Shock Relations Properties
# =============================================================================
# velocity ratio across normal shock: v2/v1
def NS_Vel_Ratio(Mach1, Mach2, Temp1, Temp2):
    '''
    Calculates the velocity ratio across a normal shock.
    v2/v1

    Parameters
    ----------
    Mach1 : Float
        Mach of freestream before shock.
    Mach2 : Float
        Mach of freestream after shock.
    Temp1 : Float
        Static Temperature (abs) of gas before shock.
    Temp2 : Float
        Static Temperature (abs) of gas after shock.

    Returns
    -------
    Float
        Velocity Ratio v2/v1.

    '''
    # Temperatures are static, in absolute units
    return (Mach2/Mach1)*np.sqrt(Temp2/Temp1)


def NS_Mach2(Mach1, Gamma=1.4):
    '''
    Calculates the Mach number of flow after a normal shock

    Parameters
    ----------
    Mach1 : Float
        Mach number of flow upstream of shock.
    Gamma : Float, optional
        Specific Heat Ratio of gas. The default is 1.4 for air.

    Returns
    -------
    Mach2: Float
        Mach number of flow downstream of the shock.

    '''
    num = Mach1**2 + 2/(Gamma-1)
    den = -1 + (2*Gamma/(Gamma-1))*Mach1**2
    Mach2 = np.sqrt(num/den)
    return Mach2

# throat area across normal shock? A2*/A1*
def NS_Throat_A_Ratio(Mach, Gamma=1.4):
    '''
    Calculates the area ratios of the throat condition before and after a shock:
        A2*/A1*
    MAY NEED VERIFIED 
    
    Parameters
    ----------
    Mach: Float
        Mach number of flow upstream of shock.

    Returns
    -------
    At2_At1: Float
        Area ratio of throats for flow downstream (2) and upstream (1) of shock.

    '''
    At2_At1 = 1/NS_Po2_Po1(Mach, Gamma)
    return At2_At1

# Properties after shock from Mach number
# static temperature ratio across normal shock: T2/T1
def NS_T2_T1(Mach1, Gamma=1.4):
    '''
    Calculates the static temperature ratio T2/T1 across a normal shock from the Mach number.

    Parameters
    ----------
    Mach1 : Float
        Mach number of freestream prior to shock.
    Gamma : Float, optional
        Specific heat ratio. The default is 1.4.

    Returns
    -------
    Float
        Ratio of T2/T1.

    '''
    Mach2 = NS_Mach2(Mach1, Gamma)
    return (1+(Mach1**2)*(Gamma-1)/2)/(1+(Mach2**2)*(Gamma-1)/2)


def NS_P2_P1(Mach1, Gamma=1.4):
    '''
    Calculates the static pressure ratio P2/P1 across a normal shock from the Mach number.

    Parameters
    ----------
    Mach1 : Float
        Mach number of freestream prior to shock.
    Gamma : Float, optional
        Specific heat ratio. The default is 1.4.

    Returns
    -------
    Float
        Ratio of P2/P1.

    '''
    Mach2 = NS_Mach2(Mach1, Gamma)
    return (1+Gamma*Mach1**2)/(1+Gamma*Mach2**2)

def NS_Po2_Po1(Mach1, Gamma=1.4):
    '''
    Calculates the stagnation pressure ratio P2/P1 across a normal shock from the Mach number.

    Parameters
    ----------
    Mach1 : Float
        Mach number of freestream prior to shock.
    Gamma : Float, optional
        Specific heat ratio. The default is 1.4.

    Returns
    -------
    Float
        Ratio of Po2/Po1.

    '''
    p2_p1 = NS_P2_P1(Mach1, Gamma)
    T1_T2 = 1/NS_T2_T1(Mach1, Gamma)
    return p2_p1*T1_T2**(Gamma/(Gamma-1))


def NS_Rho2_Rho1(Mach1, Gamma=1.4):
    '''
    Calculates the Static Density Ratio across a normal shock:
        ρ2/ρ1

    Parameters
    ----------
    Mach1 : Float
        Mach number of flow upstream of shock.
    Gamma : Float, optional
        Specific heat ratio. The default is 1.4.

    Returns
    -------
    ρ2/ρ1: Float
        Static Density Ratio (1 - Upstrea, 2 - Downstream).

    '''
    Mach2 = NS_Mach2(Mach1, Gamma)
    T1_T2 = 1/NS_T2_T1(Mach1, Gamma)
    return (Mach1/Mach2)*np.sqrt(T1_T2)

def NS_Rho_o2_Rho_o1(Mach1, Gamma=1.4):
    '''
    Calculates the Stagnation Density Ratio across a normal shock:
        ρo2/ρo1

    Parameters
    ----------
    Mach1 : Float
        Mach number of flow upstream of shock.
    Gamma : Float, optional
        Specific heat ratio. The default is 1.4.

    Returns
    -------
    ρo2/ρo1: Float
        Stagnation Density Ratio (1 - Upstrea, 2 - Downstream).

    '''
    return 1/NS_Po2_Po1(Mach1, Gamma)




def As_At_n(Ae_At, Pb_Po, Gamma=1.4, dx=1e-3, max_iter=30):
    '''
    Calculates the Area of Shock / Throat area for a CD nozzle in which the back pressure
    induces a normal shock within the nozzle.
    Procedure:
        We will take an iterative approach to finding the area ratio. This
        means that we will start by guessing the As_At ratio, solving for Pb/Po,
        and then comparing it to the given Pb_Po. We then adjust the As/At until 
        we reach the Pb/Po ratio within a tolerance.

    Definitions:
        At = The area of the real throat
        As = The area of the shock
        Ae = The area of the exit
        A* = The area of the hypotetical shock (aka when choked)

        Mt  = 1 (choked)
        Ms1 = Mach number right before the shock
        Ms2 = Mach numebr right after the shock
        Me  = Mach number at the exit

        Po1, To1 = Stagnation properties before the shock
        Po2, To2 = Stagnation properties after the shock

    Assume:
        Pe = Pb
        To1 = To2
        As1 = As2 (Shock thickness is negligible)
        As/At not= As/A2*

    '''

    Pb_Po1_target = Pb_Po
    # As_At = 1 + (1-Ae_At)/2  # Initial guess, at area halfway between throat and shock

    def F_shock_area(As_At):
        # Find mach at shock from guessed area ratio As/At:
        M_sub, M1 = Mach_at_A(As_At, Gamma)

        # Find properties after shock
        Ms2 = NS_Mach2(M1, Gamma)
        Po2_Po1 = NS_Po2_Po1(M1, Gamma)

        # Check if the pressure right after shock is = Pb
        P2_Po2 = 1/Po_P_ratio(Ms2, Gamma)
        P2_Po1 = P2_Po2 * Po2_Po1
        # if abs(P2_Po1 - Pb_Po1_target) < 1e-4:
        #     print('Exit condition found {:.4f}'.format(Pb_Po1_target))
        #     return 0

        # Find new area ratio after shock: As/A2*
        As_Astar = A_ratio(Ms2, Gamma)
        Ae_Astar = Ae_At*(1/As_At)*As_Astar

        # Find Mach at exit from new area ratio Ae/A2*
        Me, M_sup = Mach_at_A(Ae_Astar, Gamma)

        # Assuming Pe=Pb
        Pe_Po2 = 1/Po_P_ratio(Me, Gamma)
        Pe_Po1 = Pe_Po2 * Po2_Po1
        return Pb_Po1_target - Pe_Po1

    low_r, high_r = rt.rootsearch(F_shock_area, 1.0001, Ae_At+dx, dx)
    try:
        # print('PR_b = {:6.5f} r1 = {:7.4f}  r2 = {:7.4f}'.format(Pb_Po1_target,low_r,high_r))
        As_At = (low_r+high_r)/2
        return As_At
    except:
        print('FAILED: PR_b = {:6.5f} r1 = '.format(
            Pb_Po1_target), low_r, '  r2 = ', high_r)



# =============================================================================
#                Turning Flow - Oblique Shock
# =============================================================================

# Notes: Oblique shock at angle Beta relative to Vinf.
#   Resultant flow will be at angle Theta.
#   M1n = M1*sin(Beta)        M1t = M1*cos(Beta)
#   M2n = M2*sin(Beta-Theta)  M2t = M2*cos(Beta-Theta)
#   Tangential velocity does not change across shock (Mt)
#   Normal velocity changes same as machΔ across Normal shock
#   To find property changes across oblique shock, use M1n and M2n
#   - Stag Temperature stays constant

def Mach2_ob(Mach1, Beta, Theta, Gamma=1.4):
    '''
    Calculates the magnitude of M2 from Mach,
    Beta (shock angle), and Theta (deflecction angle)
      
    Parameters
    ----------
    Mach1 : Float
        Mach number of flow before shock.
    Beta : Float, degrees
        Oblique shock angle relative to line parallel to Vinf.
    Theta : Float, degrees
        Deflection angle relative to line parallel to Vinf.
    Gamma : Float, optional
        Gas constant ratio. The default is 1.4.

    Returns
    -------
    Mach2 : Float
        Mach numberof flow after shock.

    '''
    Beta = np.radians(Beta)
    Theta = np.radians(Theta)

    M1n = Mach1*np.sin(Beta)
    num = M1n**2 + 2/(Gamma-1)
    den = (2*Gamma*M1n**2)/(Gamma-1) - 1
    Mach2 = (1/np.sin(Beta-Theta))*np.sqrt(num/den)
    return Mach2


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
    m2n = NS_Mach2(m1n)
    Mach2 = m2n/np.sin(np.radians(Beta-Theta))
    m1t = Mach1*np.cos(np.radians(Beta))
    m2t = Mach2*np.cos(np.radians(Beta-Theta))
    return m1n, m1t, m2n, m2t


def Mach_Tang(Mach, Beta):
    return Mach*np.sin(np.radians(Beta))

# After swapping from oblique to normal Machs, can use them
#   to find the properties after the oblique shock

def Shock_Angle_ob(Mach, P2_P1, gamma=1.4):
    '''
    Calculates the shock angle Beta from the pressure ratio
    before (P1) and after (P2) the shock. 

    Parameters
    ----------
    Mach : Float
        Mach number of flow before shock.
    P2_P1 : Float
        Static Pressure ratio of before (P1) and after (P2) shock.
    gamma : Float, optional
        Specifc Heat Ratio. The default is 1.4.

    Returns
    -------
    Float
        Shock angle Beta.

    '''
    return np.rad2deg(np.arcsin(np.sqrt((P2_P1*(gamma+1) + (gamma-1))/(2*gamma*Mach**2))))


# __________________________________________

# =============================================================================
#           Fanno Flow
# =============================================================================

def Fanno_T_Tstar(Mach, Gamma=1.4, IntegralPoints=1000):
    '''
    Sonic Reference Temperature Ratio for Fanno Flow:
        T/T*

    Parameters
    ----------
    Mach : Float
        Mach number of flow.
    gamma : Float, optional
        Specifc Heat Ratio. The default is 1.4.
    IntegralPoints : TYPE, optional
        DESCRIPTION. The default is 1000.

    Returns
    -------
    Float
        T/T*.

    '''
    def M_internal(M):
        return -(Gamma-1)*(M) / (1 + np.power(M,2)*(Gamma-1)/2)
    
    def MachIntegral(M):
        Ms = np.linspace(M,1,IntegralPoints)
        Ms_int = M_internal(Ms)
        Integral = np.trapz(Ms_int,Ms)
        return Integral
    
    T_Tstar  = 1/np.exp(MachIntegral(Mach))
    return T_Tstar

def Fanno_P_Pstar(Mach, Gamma=1.4, IntegralPoints=1000):
    def M_internal(M):
        return -(1/M)*(1 + (Gamma-1)*np.power(M,2) ) / (1 + np.power(M,2)*(Gamma-1)/2)
    
    def MachIntegral(M):
        Ms = np.linspace(M,1,IntegralPoints)
        Ms_int = M_internal(Ms)
        Integral = np.trapz(Ms_int,Ms)
        return Integral
    
    P_Pstar  = 1/np.exp(MachIntegral(Mach))
    return P_Pstar

def Fanno_Po_Postar(Mach, Gamma=1.4, IntegralPoints=1000):
    Po_P = Po_P_ratio(Mach,Gamma)
    Pstar_Postar = 1/Po_P_ratio(1, Gamma)
    Po_Postar = Po_P * Fanno_P_Pstar(Mach,Gamma,IntegralPoints) * Pstar_Postar 
    return Po_Postar 

def Fanno_fLmax_D(Mach,Gamma=1.4,IntegralPoints=1000):
    def M_internal(M):
        return (2/M)*(1 - np.power(M,2) ) / ((1 + np.power(M,2)*(Gamma-1)/2)*Gamma*np.power(M,2))
    
    def MachIntegral(M):
        Ms = np.linspace(M,1,IntegralPoints)
        Ms_int = M_internal(Ms)
        Integral = np.trapz(Ms_int,Ms)
        return Integral 
    
    fLmax_D  = MachIntegral(Mach)
    return fLmax_D
    
def Fanno_I_Istar(Mach, Gamma=1.4):
    I_Istar = (1 + Gamma*Mach**2) / (Mach*np.sqrt( (2*(Gamma+1)) * (1 + (Mach**2)*(Gamma-1)/2)  ))    
    return I_Istar




def Fanno_Flow(**kwargs):
    # Define functions for use later
    def T_Tstar_DF(M):
        return -(Gamma-1)*(M) / (1 + np.power(M,2)*(Gamma-1)/2)
    
    def T_Tstar_F(Mach):
        def MachIntegral(M):
            Ms = np.linspace(M,1,IntegralPoints)
            Ms_int = T_Tstar_DF(Ms)
            Integral = np.trapz(Ms_int,Ms)
            return Integral
        
        T_Tstar  = 1/np.exp(MachIntegral(Mach))
        return T_Tstar
    
    def P_Pstar_DF(M):
        return -(1/M)*(1 + (Gamma-1)*np.power(M,2) ) / (1 + np.power(M,2)*(Gamma-1)/2)

    def P_Pstar_F(Mach):        
        def MachIntegral(M):
            Ms = np.linspace(M,1,IntegralPoints)
            Ms_int = P_Pstar_DF(Ms)
            Integral = np.trapz(Ms_int,Ms)
            return Integral
        
        P_Pstar  = 1/np.exp(MachIntegral(Mach))
        return P_Pstar

    def fLmax_D_DF(M):
        return (2/M)*(1 - np.power(M,2) ) / ((1 + np.power(M,2)*(Gamma-1)/2)*Gamma*np.power(M,2))

    def fLmax_D_F(Mach):
        def MachIntegral(M):
            Ms = np.linspace(M,1,IntegralPoints)
            Ms_int = fLmax_D_DF(Ms)
            Integral = np.trapz(Ms_int,Ms)
            return Integral 
        
        fLmax_D  = MachIntegral(Mach)
        return fLmax_D
    
    # Pull out the kwargs
    Mach = kwargs.get('Mach',None)
    T_Tstar = kwargs.get('T_Tstar',None)
    P_Pstar = kwargs.get('P_Pstar',None)
    fLmax_D = kwargs.get('fLmax_D',None)
    Gamma = kwargs.get('Gamma',1.4)
    
    IntegralPoints = kwargs.get('IntegralPoints',1000)
    root_dx = kwargs.get('root_dx',0.1)
    root_max_search = kwargs.get('root_max_search',8)
    root_min_search = kwargs.get('root_min_search',1e-6)
    forceSupersonic = kwargs.get('forceSupersonic',False)
    Mach_sup = None        
        
    if Mach == None: 
        Props = [T_Tstar, P_Pstar, fLmax_D]
        PropsF = [T_Tstar_F, P_Pstar_F, fLmax_D_F]
        PropsDF = [T_Tstar_DF, P_Pstar_DF, fLmax_D_DF]
        MachFoud = False
        for i,  prop in enumerate(Props):
            if prop != None:
                # Perfom root search to get mach number
                try:
                    low_r, high_r = rt.rootsearch(lambda M: PropsF[i](M) - prop, root_min_search, root_max_search, root_dx)
                    Mach = rt.newtonRaphson(lambda M: PropsF[i](M) - prop, PropsDF[i], low_r, high_r)
                    
                    if fLmax_D != None and fLmax_D > 0.7681818527269877:
                        root_min_search = 7 
                        root_max_search = 100
                        root_dx = 1
                     
                    if fLmax_D != None and fLmax_D > 0.8192125375987388:
                        Mach_sup = None
                    else:
                        low_r, high_r = rt.rootsearch(lambda M: PropsF[i](M) - prop, 1 + root_min_search, root_max_search, root_dx)
                        Mach_sup = rt.newtonRaphson(lambda M: PropsF[i](M) - prop, PropsDF[i], low_r, high_r)
                    
                    MachFound =True 
                except:
                    MachFound = False
                    raise ValueError('Could not sovle for Mach number')
        if not MachFound:
            raise ValueError('Not enough inputs')
            return None
    
    if forceSupersonic:
        if Mach < 1:
            Mach= Mach_sup if Mach_sup != None else Mach
            # print('Warning: Could not force supersonic condition')
            
    T_Tstar = Fanno_T_Tstar(Mach,Gamma)
    P_Pstar = Fanno_P_Pstar(Mach,Gamma)
    Po_Postar = Fanno_Po_Postar(Mach,Gamma)
    fLmax_D = Fanno_fLmax_D(Mach,Gamma)
    I_Istar = Fanno_I_Istar(Mach, Gamma)
    
    return {'Mach':Mach,'Mach_sup':Mach_sup, 'T_Tstar':T_Tstar,'P_Pstar':P_Pstar,'Po_Postar':Po_Postar,'fLmax_D':fLmax_D, 'I_Istar':I_Istar}
    

# =============================================================================
#  Rayleigh Flow - Heat Transfer in a Duct 
# =============================================================================

def Rayleigh_P_Pstar(Mach, Gamma=1.4):
    '''
    Calculates the Static Pressure ratio for Rayleigh flow relative to sonic condition:
        P/P*

    Parameters
    ----------
    Mach : Float
        Current Mach number of flow.
    Gamma : Float, optional
        Specific heat ratio. The default is 1.4.

    Returns
    -------
    P_Pstar : Float
        Static Reference Pressure Ratio.

    '''
    P_Pstar = (1 + Gamma) / (1 + Gamma*Mach**2)
    return P_Pstar

def Rayleigh_T_Tstar(Mach, Gamma=1.4):
    '''
    Calculates the Static Temperature ratio for Rayleigh flow relative to sonic condition:
        T/T*

    Parameters
    ----------
    Mach : Float
        Current Mach number of flow.
    Gamma : Float, optional
        Specific heat ratio. The default is 1.4.

    Returns
    -------
    T_Tstar : Float
        Static Reference Temperature Ratio.

    '''
    T_Tstar = (Mach**2)*(1 + Gamma)**2 / (1 + Gamma*Mach**2)**2
    return T_Tstar

def Rayleigh_V_Vstar(Mach, Gamma=1.4):
    '''
    Calculates the Velocity ratio for Rayleigh flow relative to sonic condition:
        V/V*

    Parameters
    ----------
    Mach : Float
        Current Mach number of flow.
    Gamma : Float, optional
        Specific heat ratio. The default is 1.4.

    Returns
    -------
    V_Vstar : Float
        Reference Velocity Ratio.

    '''
    P_Pstar = Rayleigh_P_Pstar(Mach,Gamma)
    T_Tstar = Rayleigh_T_Tstar(Mach,Gamma)
    V_Vstar = T_Tstar/P_Pstar 
    return V_Vstar

def Rayleigh_To_Tostar(Mach,Gamma=1.4):
    '''
    Calculates the Stagnation Temperature ratio for Rayleigh flow relative to sonic condition:
        To/To*

    Parameters
    ----------
    Mach : Float
        Current Mach number of flow.
    Gamma : Float, optional
        Specific heat ratio. The default is 1.4.

    Returns
    -------
    To_Tostar : Float
        Stagnation Reference Temperature Ratio.

    '''
    To_Tostar = ((Mach**2)*(1 + Gamma)**2 / (1 + Gamma*Mach**2)**2) * (1 + (Mach**2)*(Gamma-1)/2)/(1 + (Gamma-1)/2)
    return To_Tostar 

def Rayleigh_Po_Postar(Mach,Gamma=1.4):
    '''
    Calculates the Stagnation Pressure ratio for Rayleigh flow relative to sonic condition:
        To/To*

    Parameters
    ----------
    Mach : Float
        Current Mach number of flow.
    Gamma : Float, optional
        Specific heat ratio. The default is 1.4.

    Returns
    -------
    Po_Postar : Float
        Stagnation Reference Pressure Ratio.

    '''
    Po_Postar = Rayleigh_P_Pstar(Mach,Gamma)*((1+(Mach**2)*(Gamma-1)/2)**(Gamma/(Gamma-1))) / (1+(Gamma - 1)/2)**(Gamma/(Gamma-1))
    return Po_Postar

def Rayleigh_phi_MS(Mach,Gamma=1.4):
    '''
    Calculates the Phi(M^2) for Rayleigh flow from the Mach number and Gamma.

    Parameters
    ----------
    Mach : Float
        Mach number of flow.
    Gamma : Float, optional
        Specific heat ratio of gas. The default is 1.4.

    Returns
    -------
    Phi_MS : Float
        Phi(M^2) of flow.

    '''
    Phi_MS = (Mach**2)*(1+(Mach**2)*(Gamma-1)/2) / (1+Gamma*Mach**2)**2
    return Phi_MS



def Rayleigh_Flow(**kwargs):
    '''
    Calculates all related flow parameters for Rayleigh flow (Inviscid, constant area duct with heat exchange).
    Can calcualte all parameters from a single input.

    Parameters
    ----------
    **kwargs : Dict,
        Flow Params (One Necessary)
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    
    def Mach_at_Rayleigh(Prop_ratio, Prop_F,Gamma=1.4):
        low, high = rt.rootsearch(lambda P: Prop_ratio - Prop_F(P,Gamma), 1e-6, 5, 1e-6)
        Mach = (low + high)/2 
        return Mach
    
    Mach = kwargs.get('Mach',None)
    T_Tstar = kwargs.get('T_Tstar',None)
    P_Pstar = kwargs.get('P_Pstar',None)
    To_Tostar = kwargs.get('To_Tostar',None)
    Po_Postar = kwargs.get('Po_Postar',None)
    
    Gamma = kwargs.get('Gamma',1.4)
    IntegralPoints = kwargs.get('IntegralPoints',1000)
    root_dx = kwargs.get('root_dx',0.1)
    root_max_search = kwargs.get('root_max_search',8)
    
    Props = [T_Tstar, P_Pstar, To_Tostar, Po_Postar]
    PropsF = [Rayleigh_T_Tstar, Rayleigh_P_Pstar, Rayleigh_To_Tostar, Rayleigh_Po_Postar]
    
    if Mach == None: 
        MachFound = False
        for i,  prop in enumerate(Props):
            if prop != None:
                # Perfom root search to get mach number
                Mach = Mach_at_Rayleigh(Props[i],PropsF[i],Gamma)
                MachFound=True
        if not MachFound:
            raise ValueError('Not enough inputs')
            return None
        
    for i, F in enumerate(PropsF):
        Props[i] = F(Mach,Gamma)
        
    Phi_MS = Rayleigh_phi_MS(Mach, Gamma)

    return {'Mach':Mach, 'T_Tstar':Props[0],'P_Pstar':Props[1],'To_Tostar':Props[2], 'Po_Postar':Props[3], 'Phi_MS':Phi_MS}
    





# def IsoThermal(fL_D, Mach, Gamma=1.4, dM = 1e-6, tol=1e-9):
#     M1 = Mach
#     M2 = M1
#     error = -1
    
#     while prevsign==currentsign:
#         M2 += dM
#         error = -fL_D + (1 - Gamma*M1**2)/(Gamma*M1**2) - (1 - Gamma*M2**2)/(Gamma*M2**2) + np.log(M1**2/M2**2)
#         print(error, M2)
#     return M2
# =============================================================================
#           Main
# =============================================================================

if __name__ == "__main__":
    # def ob_shock_p():
    #     M1 = 2.0
    #     shock_angle = 40 # deg
    #     P1 = 20      # kPa
    #     T1 = 263.15  # K
        
    #     Theta = shockDeflection(M1, shock_angle, None)['Theta']
    #     m1n, m1t, m2n, m2t = Mach_ob_to_n(M1, shock_angle, Theta)
        
    #     P2 = P1*PR_n(m1n, m2n)
    #     T2 = T1*TR_n(m1n, m2n)
    #     M2 = np.sqrt(m2n**2 + m2t**2)
    #     # Find static pressure and temperature after shock
    #     # find shock angle 
    #     # find mach 2
    
    def ex_windtunnel():
        T2 = 216.7
        P2 = 5.53e3
        M2 = 2.0
        A = np.pi*(0.25/2)**2
        
        Po = P2*Po_P_ratio(M2)
        To = T2*To_T_ratio(M2)
        
        m_dot = mdot(Po, To, A, M2)
        V2 = M2*np.sqrt(1.4*287*T2)
        cp = 1.005 # J/kg*K
        
        h2 = cp*T2
        h1 = cp*To
        
  
    # # Problem 1 exam
    # Po = 200e3
    # To = 500
    # M = 1 
    # A = 15*(1/100)**2
    # gam = 1.4
    # R = 287
    # print(mdot(Po, To, A,M,gam,R))
    
    
    # PR_exit_shock = P2_P1_n(2.6374, 1.4)*(1/Po_P_ratio(2.6374, 1.4))

    # flow_ang = Shock_Angle_ob(2.0,0.5)
    # M2 = ExpansionFan(2.0, flow_ang)

    f = 0.02
    A_At = 2 
    Po = 500
    L_D = 25
    Pb = 0
    Me = 1 
    
    Mi = Mach_at_A(A_At)[1]
    i = Isentropic_Flow(Mach=Mi)
    Pi_Po  = i['P_Po']
    
    i_f = Fanno_Flow(Mach=Mi)
    Pi_Pstar = i_f['P_Pstar'] 
    fLmax_D_i = i_f['fLmax_D']
    
   
    # P1 = s1['P_Po']*500
    
    # s1_fanno = Fanno_Flow(Mach=s1['Mach'])
    
    # fLmax2_D = s1_fanno['fLmax_D'] - 0.2 
    
    # s2 = Fanno_Flow(fLmax_D=fLmax2_D)
    # P2_P1 = s2['P_Pstar']*s1_fanno['P_Pstar']
    
    
    # M3 = NS_Mach2(s2['Mach_sup'])
    # P3_P2 = P2_P1_n(s2['Mach_sup'])
    # P3 = P3_P2*P2_P1*P1
    # IsoThermal(62.5,0.03375,1.32,1e-9)

    
    
    
    
    
    
    
    