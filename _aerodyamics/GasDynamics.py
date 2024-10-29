# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 11:34:05 2023

@author: cycon
"""
# Imports
import sys
sys.path.append("..")
import _tools.Interpolator as intrp
import _tools.RootAlgorithms as rt
import numpy as np
import matplotlib.pyplot as plt
import math




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
    Outputs a dictionary containing the Mach number, Beta angle, and Theta angle for a shock deflection.
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


def mdot(Po, To, A, Mach=1, Gamma=1.4, R=287):
    '''
    Calculates the mass flow rate from stagnation conditions

    Parameters
    ----------
    Po : Float
        Stagnation Pressure [Pa].
    To : Float
        Stagnation Temperature [K].
    A : Float
        Area of channel [m^2].
    Mach : Float, optional
        Mach of gas. The default is 1.
    Gamma : Float, optional
        DESCRIPTION. The default is 1.4.
    R : TYPE, optional
        Gas Constant [J/kg*K]. The default is 287.

    Returns
    -------
    Float
        mdot [kg/s].

    '''
    first = (Po*A/(np.sqrt(R*To))) * Mach*np.sqrt(Gamma)
    second = 1 + ((Gamma-1)/2)*(Mach**2)
    power = -(Gamma+1)/(2*(Gamma-1))
    return first*np.power(second, power)


def Pcrit_Po(Gamma=1.4):
    return (1 + 0.5*(Gamma-1))**(-Gamma/(Gamma-1))


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


def nu_PM(Mach, Gamma=1.4):
    return np.sqrt((Gamma+1)/(Gamma-1)) * \
        np.arctan(np.sqrt((Gamma-1)*(np.power(Mach, 2) - 1)/(Gamma+1))) - \
        np.arctan(np.sqrt(np.power(Mach, 2) - 1))

# Exact Solutions


def dv_PM(Mach, Gamma=1.4):
    Num = 2*np.sqrt(np.power(Mach, 2) - 1)
    Den = Mach*(Gamma*np.power(Mach, 2) - np.power(Mach, 2) + 2)
    return np.divide(Num, Den)

# Full Prendtly-Meyer Solver to get Mach from input


def ExpansionFan(Mach1, Theta, Gamma=1.4):
    def PM(Mach2):
        return nu_PM(Mach2) - nu_PM(Mach1) - np.radians(Theta)

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
def vR_n(Mach1, Mach2, Temp1, Temp2):
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

# density ratio across normal shock: rho2/rho1


def rhoR_n(Mach1, Mach2, Temp1, Temp2):
    # Temperatures are static, in absolute units
    return (Mach1/Mach2)*np.sqrt(Temp1/Temp2)

# static temperature ratio across normal shock: T2/T1
def TR_n(Mach1, Mach2, Gamma=1.4):
    return (1+(Mach1**2)*(Gamma-1)/2)/(1+(Mach2**2)*(Gamma-1)/2)

# static pressure ratio across normal shock: p2/p1


def PR_n(Mach1, Mach2, Gamma=1.4):
    '''
    Returns P2/P1

    Parameters
    ----------
    Mach1 : TYPE
        DESCRIPTION.
    Mach2 : TYPE
        DESCRIPTION.
    Gamma : TYPE, optional
        DESCRIPTION. The default is 1.4.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
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

# Properties after shock from Mach number
# static temperature ratio across normal shock: T2/T1


def T2_T1_n(Mach1, Gamma=1.4):
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
    Mach2 = Mach2_n(Mach1, Gamma)
    return (1+(Mach1**2)*(Gamma-1)/2)/(1+(Mach2**2)*(Gamma-1)/2)


def P2_P1_n(Mach1, Gamma=1.4):
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
    Mach2 = Mach2_n(Mach1, Gamma)
    return (1+Gamma*Mach1**2)/(1+Gamma*Mach2**2)


def Rho2_Rho1_n(Mach1, Gamma=1.4):
    Mach2 = Mach2_n(Mach1, Gamma)
    T1_T2 = 1/T2_T1_n(Mach1, Gamma)
    return (Mach1/Mach2)*np.sqrt(T1_T2)


def Po2_Po1_n(Mach1, Gamma=1.4):
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
    p2_p1 = P2_P1_n(Mach1, Gamma)
    T1_T2 = 1/T2_T1_n(Mach1, Gamma)
    return p2_p1*T1_T2**(Gamma/(Gamma-1))


def As_At_n(Ae_At, Pb_Po, Gamma=1.4, dx=1e-3, max_iter=30):
    '''
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
        Ms2 = Mach2_n(M1, Gamma)
        Po2_Po1 = Po2_Po1_n(M1, Gamma)

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


# _________________________________________________________

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
    B_rad = np.radians(Beta)
    th_rad = np.radians(Theta)

    num = (np.sin(B_rad)**2)*(Mach1**2) + 2/(Gamma-1)
    den = (2*Gamma/(Gamma-1))*(np.sin(B_rad)**2)*(Mach1**2) - 1

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
    P1_P0 : Float
        Static Pressure ratio of before (P1) and after (P2) shock.
    gamma : Float, optional
        Specifc Heat Ratio. The default is 1.4.

    Returns
    -------
    Float
        Shock angle Beta.

    '''
    return np.rad2deg(np.arcsin(np.sqrt((P2_P1*(gamma+1) + (gamma-1))/(2*gamma*Mach**2))))

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


def Mach_at_PR(Po_P, Gamma=1.4):
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
    return np.sqrt(((Po_P)**((Gamma-1)/Gamma) - 1)*2/(Gamma-1))
# __________________________________________


# =============================================================================
#           Main
# =============================================================================

if __name__ == "__main__":
    # T1 = 300 # K
    # Vs1 = 1035.853 # m/s
    # Vg1 = 25 # m/s
    # Vg2 = 50 # m/s
    # gam = 5/3
    # R = 2076.9423
    # V1 = Vs1 
    # V2 = Vs1 - Vg1 
    
    # Ms1 = V1/np.sqrt(gam*R*T1)
    # T2 = T1*T2_T1_n(Ms1, gam)
    
    # a2 = np.sqrt(gam*R*T2)
    # M2 = Mach2_n(Ms1, gam)
    
    # def f(Vs2):
    #     return ((gam+1)*((Vs2-Vg1)/a2)**2)/((gam-1)*((Vs2-Vg1)/a2)**2 + 2) - (Vs2-25)/(Vs2-50)
    
    # def f2(Vs2):
    #     return 
   
    # # find roots: 
    # print(rt.rootsearch(f, 0, 2000, 0.01))
    # Po = 5000 # kpa
    # A_At = A_ratio(1.75)
    # M_sub, M_sup = Mach_at_A(A_At)
    # Pb_max = Po/Po_P_ratio(M_sub)
    # Pb_sup = Po/Po_P_ratio(M_sup)
    
    # Ms2 = Mach2_n(M_sup)
    # P2 = Pb_sup*PR_n(M_sup, Ms2)
    
    Po = 3000e3#kpa
    Pb1 = 101e3 # kpa
    Pb2 = 0# kpa
    gam = 1.4
    R = 400
    To = 1600
    Me = Mach_at_PR(Po/Pb1,gam)
    Ae_At = A_ratio(Me,gam)
    mdotA = mdot(Po, To, A=1,Mach=Me,Gamma=gam,R=R)
    
    T_rat = 1/To_T_ratio(Me,gam)
    Te = To*T_rat
    Ve = Me*np.sqrt(gam*R*Te)
    Rho_e = Pb1/(R*Te)
    
    mdotA_2 = mdot_s(Pb1,Te,1,Me,gam,R)
    
    Ft1_Ft2 = mdotA_2*Ve/(mdotA_2*Ve+Pb1)
    
    Ft1_Ft2_2 = (Rho_e*Ve**2)/(Rho_e*Ve**2 + Pb1)
    
    
    # M = 2.75 + ((0.03367 - 0.037425)/(0.024663-0.037425))*(3-2.75)
    # 2.75  0.037425 
    # 3.00  0.024663
    
    
    
    
    
    
   
