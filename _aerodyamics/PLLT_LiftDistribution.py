# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 01:23:34 2023

@author: cycon
"""
import numpy as np
import scipy.integrate as spint 
import matplotlib.pyplot as plt
import MatrixManipulator as MM
import RootAlgorithms as RT
import time

# z(x) = 4*h*x(1-x) - the camber line
def TAFT5(NACA, c=1):
    # 5 Series: NACA LPQXX
    #   L - amount of camber: Design CL = 3*l/2, l = L*10%
    #   P - Designator for pos of max camber, xf = p/2, p = P*10%
    #   Q - 0 = Standard 5 digit airfoil camber
    #       1 = "Reflexed" camber
    #   XX - Max thickess of airfoil as % chord: h = XX*100%
    
    L = int(NACA[0])*0.10
    p = int(NACA[1])*0.10
    Q = NACA[2]
    h = int(NACA[3:])*0.01
    
    # Locaiton of max camber:
    xf = p/2
    
    def F_m(m):
      return m*(1-np.sqrt(m/3))
    
    def dF_m(m):
       return 1-np.sqrt(m/3)-np.sqrt(m)/(2*np.sqrt(3))
    

    # Root search - Returns the lower and upper bounds of a root between start and end
    low_r, high_r = RT.rootsearch(F_m, 0, 2, 0.01)
    # Newton Raphson
    root = RT.newtonRaphson(F_m,dF_m,low_r,high_r)
    
    m = root
    # Cacluate K1
    CLi = L*3/2
    Q = (3*m - 7*m**2 + 8*m**3 - 4*m**4)/(np.sqrt(m*(1-m))) - 1.5*(1-2*m)*(np.pi/2 - np.arcsin(1-2*m))
    K1 = 6*CLi/Q
    
    # Mean Camber line functions
    lim = m
    if Q == '0':
        def dz_dx(x_c):
            if x_c <= lim:
                return (K1/6)*( 3*(x_c**2) - 6*m*x_c + (m**2)*(3-m) )
            else:
                return -K1*(m**3)/6
    elif Q == '1':
        def dz_dx(x_c):
            return None
    else:
        print('Error: Invalid airfoil number. Q != 0 or 1')
        return None
    
    
def TAFT(NACA, c=1):
    # 4 Series: NACA XYZZ
    #   X - Max camber as percent of chord: m = X*100% <br>
    #   ZZ - Mach Thickness of airfoil as % chord: h = ZZ*100%
    
    # z is the mean camber line 
    
    # Cant check NACA for correct input yet
    # Set up parameters 
    m = int(NACA[0])*0.01*c
    p = int(NACA[1])*0.1*c
    h = int(NACA[2:])*0.01*c
    
    # Define Piece wise functions
    lim_th = np.arccos(1-2*p/c)
    def dz_dx(theta):
        if theta <= lim_th:
            return ((2*m/p)-(m*c/p**2)*(1-np.cos(theta)))
        else:
            return ((2*p*m - m*c +m*c*np.cos(theta))/((c-p)**2))

    def dz_dx_Cl(theta):
        return dz_dx(theta)*(np.cos(theta)-1)

    def dz_dx_A1(theta):
        return dz_dx(theta)*np.cos(theta)

    def dz_dx_A2(theta):
        return dz_dx(theta)*np.cos(2*theta)

    # Compute Integrals
    a_L0 = -(1/np.pi)*spint.quad(dz_dx_Cl,0,np.pi)[0]
    A_1  = (2/np.pi)*spint.quad(dz_dx_A1,0,np.pi)[0]
    A_2  = (2/np.pi)*spint.quad(dz_dx_A2,0,np.pi)[0]

    # Compute Coefficients
    def CL(alpha):
        return 2*np.pi*(alpha - a_L0)

    def CM_LE(alpha): 
        return - (CL(alpha)/4 + np.pi*(A_1-A_2)/4)
    
    CM_c4 = (np.pi/4)*(A_2-A_1)
    
    def x_CP(alpha):
        return (c/4)*(1+ np.pi*(A_1-A_2)/CL(alpha))
    
    return np.degrees(a_L0), CM_c4, CL, CM_LE, x_CP


# Function assuming a constant airfoil acros cross section, sym about midpoint
def LiftDistribution(a, a0, b, c, Vinf, S=None, N=50, N_p=100): 
    ## !!!!
    # Want to add a variable AoA to account for geometric twist (angle of incidence)
    
    # If c is passed as a scalar, covnert to a function
    if not callable(c):
        chord_length = c;
        def def_c(y):
            return chord_length
        c = def_c
        
    if not callable(a0):
        AoA_L0 = a0;
        def def_a0(y):
            return AoA_L0
        a0 = def_a0
    
    if S==None:
        AR = b*c(1) # Rectangular Wing
    else: 
        AR = b**2/S 
        
    # Convert to radians for calculations
    a_rad = np.radians(a)
    
    def a0_rad(y):
        return np.radians(a0(y))
    
    # Set up our thetas
    theta = np.linspace(0,np.pi,N)
    theta_p = np.linspace(0,np.pi,N_p)
    
    # NEED TO USE THETA TO CALCULATE y CUZ SPACING IS GOING TO DIFFER
    y = -(b/2)*np.cos(theta)
    y_p = -(b/2)*np.cos(theta_p)
    
    # Initialize arrays/Matrices
    Gamma = np.zeros(N_p)
    A = np.zeros((N,N),dtype=float)
    sol = np.empty(N)
    for i, yi in enumerate(y):
        sol[i] = a_rad - a0_rad(yi)
    
    # Calculate the matrix to find the coefficients An
    for i, th in enumerate(theta):
        for n in range(1,N+1):
            A[i, n-1] = 2*b /(np.pi*c(y[i])) * np.sin(n*th)
            if np.sin(th) != 0:
                A[i, n-1] += n*np.sin(n*th)/np.sin(th)
            else: 
                A[i, n-1] += n*n*np.cos(n*th)/np.cos(th)

    # Get the coefficient vector
    An = MM.gaussPivot(A.copy(), sol.copy())
   
    # Get error in case it is needed
    res = MM.accuracy(A, An, sol)
    
    # Calculate Gamma
    for i, th in enumerate(theta_p):
        sig=0
        for n in range(1,N+1):
            sig += An[n-1]*np.sin(n*th)
        Gamma[i] = 2*b*Vinf*sig
    
    if S != None:
        # Calculate Lift Coefficient
        CL1 = 2*np.trapz(Gamma, x=y_p)/(Vinf*S)
        CL2 = An[0]*np.pi*AR
        
        CL = (CL1 + CL2)/2
        
        # Calculate induced angle of attack?
        a_i = np.empty(N_p)
        
        for i, yo in enumerate(y_p):
            a_i[i] = a_rad - a0_rad(yo) - Gamma[i]/(np.pi*Vinf*c(yo))
        
        Cdi1 = 2*np.trapz(Gamma*a_i, x=y_p)/(Vinf*S)
        
        sig = 0
        for n in range(2, N+1):
            sig += n*(An[n-1]/An[0])**2
            
        Cdi2 = np.pi*AR*An[0]**2 * (1 + sig)
        
        Cdi = (Cdi1 + Cdi2)/2
    
        # print('Cdi1 = {} \nCdi2 = {}'.format(Cdi,Cdi2))
        # print('CL1 = {} \nCL2 = {}'.format(CL,CL2))
        
        # Seems to be more accurate if we average out the CL values between each method
    return [theta_p, y_p, Gamma] if S == None else [theta_p, y_p, Gamma, CL, Cdi]

def chord_Lin_func_maker(max_c, min_c, b):
    # will work according to 
    
    def chord_funct_y(y):
        c = abs(y)*(-min_c/(2*b)) + max_c
        return c
    
    def chord_funct_theta(th):
        c = (max_c-min_c)*np.sin(th) + min_c
        return c
    
    return chord_funct_theta, chord_funct_y


def plotPlaniform(b,c_func):
    # Assuming c is callable for now
    # Works on theta:
    th = np.linspace(0,np.pi,100,endpoint=True)
    # y = -b*np.cos(th)/2
    y = np.linspace(-b/2, b/2, 100, endpoint=True)
    c = c_func(y) # Chord lengths
    
    upperC = c/2
    lowerC = -c/2
    leftTip = np.array([[y[0], upperC[0]],[y[0], lowerC[0]]])
    rightTip = np.array([[y[-1], upperC[-1]],[y[-1], lowerC[-1]]])
    
    plt.figure()
    plt.plot(y,upperC, y, lowerC, leftTip[:,0], leftTip[:,1], rightTip[:,0], rightTip[:,1])
    c_max = np.max(c)
    plt.ylim([-c_max,c_max])
    # plt.xlim([-3*b/4,3*b])
    plt.grid()
    
    
if __name__ == '__main__':
    Airfoil = '6412'
    
    # Set up some standard constants
    a = 0           # deg - angle of attack of interest
    Vinf = 10       # m/s - freestream velocity
    b = 30          # m - span length
    c = 3           # m - chord length 
    a0 = TAFT(Airfoil , c)[0]  # deg - zero lift angle of attack
    S = b*c         # m^2 - Ref Area (different if chord varies, this is for rect. wing)
    AR = b**2 / S   # Nodim - Aspect Ratio
    
    
    t, y, g, Cl, Cdi = LiftDistribution(a, a0, b, c, Vinf,S)
    
    # ---- Plotting Gamma -----
    if True:
        # Plotting
        plt.figure()
        plt.plot(y, g)
        # plt.yticks(np.arange(0,15,2.5))
        plt.title('NACA{}\n b = {}, c = {}\n $C_L$ = {:.3f} $C_{{{}}}$ = {:.5f}'.format(Airfoil, b, c, Cl,'Di',Cdi))
        plt.xlabel('Span (m)')
        plt.ylabel('Gamma')
        plt.grid()
    
    if False:
        for n in range(5,50):
            t, g, Cl = LiftDistribution(a, a0, b, 1, Vinf,S,n, 100)
        
            # ---- Plotting Gamma -----
            if True:
                # Plotting
                plt.plot(t, g)
                # plt.yticks(np.arange(0,15,2.5))
                plt.xlabel('Theta (rad)')
                plt.ylabel('Gamma')
                plt.title('$C_L$ = {}\nN={}'.format(Cl, n))
                plt.grid()
                plt.pause(0.5)
                plt.clf()
        
        
# # --- Chord Testing ----
# if False:
#     print(chord(np.arange(-5,6,1)))
    
#     chord = chord_Lin_func_maker(1,0.5,5)
#     theta = np.linspace(0,np.pi,10,endpoint=True)
#     print(chord(theta))
    
#     plt.plot(theta, chord(theta))
#     plt.ylim(0,2)
