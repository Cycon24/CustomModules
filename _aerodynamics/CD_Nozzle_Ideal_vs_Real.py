# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 12:58:30 2024

@author: cycon
"""

import numpy as np
import matplotlib.pyplot as plt
import GasDynamics as GD
import pandas as pd


def x_pos(A_ft2, BeforeThroat):
    #𝐴(𝑥) = (𝜋/4)[(0.02246 ⋅ 𝑋) + 0.2504]^2  − 0.01327 where where X = distance from throat
    A_in2 = A_ft2 * (12/1)**2
    x_in = (-0.2504+np.sqrt((A_in2 + 0.01327)*(4/np.pi)))/0.02246
    
    return -x_in/8 if BeforeThroat else x_in

def A_pos(dx_t):
    # 𝐴(𝑥) = (𝜋/4)[(0.02246 ⋅ 𝑋) + 0.2504]^2  − 0.01327
    x = min(abs(dx_t - 0.3), 2) 
    A = (np.pi/4)*(0.02246*x + 0.2504)**2 - 0.01327 # in^2
    A_ft2 = A*(1/12)**2
    return A_ft2


def fL_D(M, A, dA, dM, dx, L, Gamma):
    fLD= L*(2/(Gamma*M**2))*((1/A)*(dA/dx) + (dM/dx)*(((1-M**2)*(1/M))/(1+(M**2)*(Gamma-1)/2)))
    D = np.sqrt(4*A/np.pi) # A = pi/4 d^2 -> D = sqrt(4A/pi)
    f = fLD * D/L
    return fLD, f

# =============================================================================
# Retrieve and organize data
# ============================================================================+
# r'Lab_02/
df_1 = pd.read_csv('AeroLab2_DataSheet.csv').to_numpy()
df_2 = pd.read_csv('AeroLab2_DataSheet_02.csv').to_numpy()

ProbePosition = np.array(df_1[1:,0], dtype=int)     # Integer
ProbePressure = np.array(df_1[1:,1:], dtype=float)  # kPa

BackPressure    = np.array(df_2[0,1:], dtype=float) # psig
InletTemp       = np.array(df_2[1,1:], dtype=float) # C
OutletTemp      = np.array(df_2[2,1:], dtype=float) # C
deltaP          = np.array(df_2[3,1:], dtype=float) # inH2O - Orifice dP
BarometricP     = float(df_2[4,1])                # inHg

# =============================================================================
# Convert Values
# =============================================================================
# Convert into Imperial Units and absolute units (if needed)
InletTemp =  InletTemp *9/5 + 491.67    # Covnert to R
OutletTemp = OutletTemp*9/5 + 491.67    # Covnert to R
deltaP *= 0.0360912                     # Convert to psi [lbf/in^2]
deltaP_ft2 = deltaP*(12/1)**2           # Convert to [lbf/ft^2]
BarometricP *= 0.4911543448             # Convert t0 psia
BackPressure_abs  = BackPressure + BarometricP # Convert to psia
ProbePressure_abs = ProbePressure+ BarometricP # Convert to psia


# =============================================================================
#   Inputs
# =============================================================================
locs = np.arange(-0.5, 2.31, 0.1) # Experimental
# Nozzle Diameters
d_th  = 0.2504              # [in] - Diameter of Nozzle Throat
d_probe = 0.130             # [in] - Diameter of Probe
d_e   = 0.2869              # [in] - Diameter of Nozzle Exit

# Areas
A_probe = (np.pi*(d_probe/2)**2)*(1/12)**2
At  = np.pi*((d_th/2)**2)*(1/12)**2 -A_probe  # [ft^2] - Throat Area 
Ae  =  np.pi*((d_e/2)**2  - (d_probe/2)**2)*(1/12)**2  # [ft^2] - Exit Area 
Ai_At = 10 # APproximately stag condition
Ae_At = Ae/At 
A_rats = np.array([A_pos(x) for x in locs])/At#np.append(np.linspace(Ai_At, 1,num=20), np.linspace(1,Ae_At,num=80,endpoint=True))  
A_rats[0:7] = [Ai_At, Ai_At,Ai_At,Ai_At, Ai_At/2, Ai_At/3, Ai_At/4]
# Constants
F_a   = 1.0                 # [NonDim] - Thermal Expansion Factor
C_d   = 0.623               # [NonDim] - Discharge Coefficient
Beta  = 1.6255/3.1875       # [NonDim] - d/D - Diameter ratio
g_c   = 32.174              # [𝑙𝑏𝑚−𝑓𝑡/𝑙𝑏𝑓−𝑠^2] - Gravitational Constant

R_univ = 10.731577089016    # [psi⋅ft3/lbmol⋅°R]  - Universal Gas Constant
MW_air = 28.966             # [lbm/lbmol]
R      = R_univ / MW_air    # [psi⋅ft3/lbm⋅°R] Gas Constant if  air
gamma = 1.4



# Get isentropic ratios
M_isen_sub, M_isen_sup = GD.Mach_at_A(Ae_At, Gamma=1.4)
PR_isen_sup = 1/GD.Po_P_ratio(M_isen_sup, Gamma=1.4) # Pe/Po
PR_isen_sub = 1/GD.Po_P_ratio(M_isen_sub, Gamma=1.4) # Pe/Po
PR_exit_shock = GD.P2_P1_n(M_isen_sup, gamma)*(1/GD.Po_P_ratio(M_isen_sup, gamma))
PR_crit = GD.Pcrit_Po(gamma) # DEPENDENT ON GAMMA


# =============================================================================
# Get Pressure Ratios Experimentsal
# =============================================================================

# Calculate density (NOT SURE IF RIGHT)
P_f = BarometricP #ProbePressure_abs[28,:] # lbf/in^2
rho_f = P_f / (R*OutletTemp)# [lbm/ft^3] - Fluid Density (Needs verified)

# Calculate Pressure Ratios
Po = 85+BarometricP # psia, Total pressure was 85 psig 
PRs =   ProbePressure_abs/Po # TODO: is ProbePressure total? is this eqtn right idfk
PRs_e = ProbePressure_abs[25,:]/Po
PRs_b = ProbePressure_abs[28,:]/Po


# =============================================================================
# Theoretical Pressure Calculation 
# =============================================================================
# Find Shock Area Locations from back pressure ratios
As_Ats = np.empty(PRs_b.shape)
for iPb, Pb_Po in enumerate(PRs_b):
    if (Pb_Po < PR_isen_sub) and (Pb_Po > PR_exit_shock):
        # Flow is choked and in between sub isen solutions and exit shock -> shock in noz
        # Find Location of shock
        As_Ats[iPb] = GD.As_At_n(Ae_At, Pb_Po, gamma,)
    else:
        As_Ats[iPb] = None
        

# Setup for solving for pressure ratios
# Will be len(A_rats) rows by len(PRs_b) cols
Px_Pos = np.empty([len(A_rats),len(PRs_b)],dtype=float) 
BeforeThroat = True     # Define toggle for before and after throat
tol = 1e-5

for iA, Ax_At in enumerate(A_rats):
    if abs(Ax_At - min(A_rats)) < tol:
        BeforeThroat=False # Toggle swapped so now everything is after throat
    
    
        # Cycle through the back pressure ratios for after throat cases to find pressure
    for iPb, Pb_Po in enumerate(PRs_b):
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
                # if Pb_Po==PRs_b[-1]:
                #     Px_Pos[iA, iPb] = Pb_Po
                # else:
                Px_Pos[iA, iPb] = 1/GD.Po_P_ratio(M_sup, gamma) 
            
            # Px_Pos[-1,iPb] = Pb_Po
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
        
# =============================================================================
#  Calculate Mach Number for Theoretical and Experimental 
# =============================================================================
# Seperate pressures and xlocations for just between throat and exit
ie = 25
it = 0

Px_Po_th = Px_Pos[it:ie+1, :]
Px_Po_ex = PRs[it:ie+1, :]
x_short = locs[it:ie+1]
As = (A_rats[it:ie+1]*At)*(12/1)**2  # in^2
M_th = np.zeros(Px_Po_th.shape)
M_ex = np.zeros(Px_Po_th.shape)
fL_D_th = np.zeros(Px_Po_th.shape)
fL_D_ex = np.zeros(Px_Po_th.shape)
f_th = np.zeros(Px_Po_th.shape)
f_ex = np.zeros(Px_Po_th.shape)

for j in range(0, len(Px_Po_th[0,:])):
    for i in range(0, len(Px_Po_th[:,0])):
        M_th[i,j] = GD.Mach_at_PR(1/Px_Po_th[i,j], Gamma=gamma)
        M_ex[i,j] = GD.Mach_at_PR(1/Px_Po_ex[i,j], Gamma=gamma)
        
        if i > 8:
            fL_D_th[i,j], f_th[i,j] = fL_D(M_th[i,j], As[i], As[i]-As[i-1], M_th[i,j]-M_th[i-1,j],x_short[i]-x_short[i-1],x_short[-1]-x_short[0],gamma)
            fL_D_ex[i,j], f_ex[i,j] = fL_D(M_ex[i,j], As[i], As[i]-As[i-1], M_ex[i,j]-M_ex[i-1,j],x_short[i]-x_short[i-1],x_short[-1]-x_short[0],gamma)
            
# =============================================================================
# Plot Theoretical Results
# =============================================================================
plt.close('all')
plt.figure()

x_locs = locs
for i, Pb in enumerate(PRs_b):
    plt.plot(x_locs,Px_Pos[:,i],label='Pb={:.0f} psig'.format(BackPressure[i]))
plt.plot([min(x_locs), max(x_locs)], [PR_crit, PR_crit], '--', label="Pb=P_crit")
plt.plot([min(x_locs), max(x_locs)], [PR_isen_sup, PR_isen_sup],'--', label="Pb=P_fe")
plt.plot([min(x_locs), max(x_locs)], [PR_isen_sub, PR_isen_sub],'--', label="Pb=P_ue")
plt.plot([min(x_locs), max(x_locs)], [PR_exit_shock, PR_exit_shock],'--', label="Pb=P_e_shock")
plt.plot([locs[8], locs[8]], [0,1], '-k',label='Throat')
plt.plot([locs[25], locs[25]], [0,1], '-k',label='Exit')
plt.xlabel('x loc [in]')
plt.ylabel('P/Po')
plt.legend(bbox_to_anchor=(1.1,0.8))
plt.grid()
plt.tight_layout()


# # =============================================================================
# # # Plot Experimental Results
# plt.clf()
# # =============================================================================\

#     # TODO: can use locations or probe positions, unclear which he wants

def normPressProfile_single(Pb, Poss, PRs, inParent=False):
    if inParent:
        plt.plot(Poss, PRs, label="Pb={:.0f} psig".format(Pb))
        plt.xlabel("Position (in.)"); plt.ylabel(r"$\dfrac{P}{P_0}$")
    else:
        plt.figure(f"Normalized Pressure Profile (1a): Pb={Pb:.0f} psig")
        plt.title(f"Normalized Pressure Profile, Pb={Pb:.0f} psig")
        plt.plot(Poss, PRs)
        plt.xlim(-0.5, 2.5); plt.ylim(0, 6)
        plt.xlabel("Position (in.)"); plt.ylabel(r"$\dfrac{P}{P_0}$")
        plt.grid()
    
def normPressProfile_all(BackPressureArray, locs, PRsArray):
    plt.figure(f"Normalized Pressure Profile (1a)")
    plt.title(f"Normalized Pressure Profile")
    for i in range(9):
        normPressProfile_single(BackPressureArray[i], locs, PRsArray[:,i],True)
    # Plot Critical Pressure
    plt.plot([min(x_locs), max(x_locs)], [PR_crit, PR_crit], '--', label="Pb=P_crit")
    plt.plot([min(x_locs), max(x_locs)], [PR_isen_sup, PR_isen_sup],'--', label="Pb=P_fe")
    plt.plot([min(x_locs), max(x_locs)], [PR_isen_sub, PR_isen_sub],'--', label="Pb=P_ue")
    plt.plot([min(x_locs), max(x_locs)], [PR_exit_shock, PR_exit_shock],'--', label="Pb=P_e_shock")
    plt.plot([locs[8], locs[8]], [0,1], '-k',label='Throat')
    plt.plot([locs[25], locs[25]], [0,1], '-k',label='Exit')
    plt.legend(bbox_to_anchor=(1.1,0.8))
    plt.grid()
    plt.tight_layout()

# # can generate multiple plots here (comment out if working on dif. section)
# # for i in range(8):
# #     normPressProfile(BackPressure[i], locs, PRs[:,i])
normPressProfile_all(BackPressure, locs, PRs)
plt.show()

# =============================================================================
# Plot Mach numbers Exp vs Theoretical
# =============================================================================
plt.figure()
colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive']
for i in range(0,len(M_ex[0,:])):
    plt.plot(x_short, M_ex[:,i], color=colors[i],label='Pb={:.0f} psig'.format(BackPressure[i]))
    plt.plot(x_short, M_th[:,i], '--', color=colors[i],label='Pb={:.0f} psig'.format(BackPressure[i]))
plt.legend(bbox_to_anchor=(1.1,0.9))
plt.xlabel('x loc [in]')
plt.ylabel('Mach Number')
plt.grid()
plt.tight_layout()

# =============================================================================
# Plot fL/D Exp vs Theoretical
# =============================================================================
plt.figure()
colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive']
for i in range(0,len(M_ex[0,:])):
    plt.plot(M_ex[:,i], fL_D_ex[:,i], color=colors[i],label='Pb={:.0f} psig'.format(BackPressure[i]))
    plt.plot(M_th[:,i], fL_D_th[:,i], '--', color=colors[i],label='Pb={:.0f} psig'.format(BackPressure[i]))
plt.legend(bbox_to_anchor=(1.1,0.9))
plt.xlabel('Mach Number')
plt.ylabel('fL/D [NonDim]')
plt.grid()
plt.tight_layout()

# =============================================================================
# Plot f Exp vs Theoretical
# =============================================================================
plt.figure()
colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive']
for i in range(0,len(M_ex[0,:])):
    plt.plot(x_short, f_ex[:,i], color=colors[i],label='Pb={:.0f} psig'.format(BackPressure[i]))
    plt.plot(x_short, f_th[:,i], '--', color=colors[i],label='Pb={:.0f} psig'.format(BackPressure[i]))
plt.legend(bbox_to_anchor=(1.1,0.9))
plt.xlabel('Mach Number')
plt.ylabel('f [NonDim]')
plt.grid()
plt.tight_layout()

# =============================================================================
# Plot fL/D * M: Exp vs Theoretical
# =============================================================================
plt.figure()
colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive']
for i in range(0,len(M_ex[0,:])):
    plt.plot(x_short, np.multiply(fL_D_ex[:,i],M_ex[:,i]), color=colors[i],label='Pb={:.0f} psig'.format(BackPressure[i]))
    plt.plot(x_short, np.multiply(fL_D_th[:,i],M_th[:,i]), '--', color=colors[i],label='Pb={:.0f} psig'.format(BackPressure[i]))
plt.legend(bbox_to_anchor=(1.1,0.9))
plt.xlabel('x loc[in]')
plt.ylabel('fL/D * M [NonDim]')
plt.grid()
plt.tight_layout()




# =============================================================================
# Surf Mach Number
# =============================================================================
fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
Y = x_short
X = PRs_b
X,Y = np.meshgrid(X,Y)

surf = ax.plot_surface(X, Y, M_th, cmap='viridis')

# Add labels and title
ax.set_xlabel('PR')
ax.set_ylabel('x loc')
ax.set_zlabel('Mach')

ax = fig.add_subplot(122, projection='3d')
Y = x_short
X = PRs_b
X,Y = np.meshgrid(X,Y)

surf = ax.plot_surface(X, Y, M_ex, cmap='viridis')

# Add labels and title
ax.set_xlabel('PR')
ax.set_ylabel('x loc')
ax.set_zlabel('Mach')


plt.title('Mach Number Surface Plots')
fig.colorbar(surf)
plt.show()


# =============================================================================
# Surf fL/D * M
# =============================================================================
fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
Y = x_short
X = PRs_b
X,Y = np.meshgrid(X,Y)

surf = ax.plot_surface(X, Y, np.multiply(M_th,fL_D_th), cmap='viridis')

# Add labels and title
ax.set_xlabel('PR')
ax.set_ylabel('x loc')
ax.set_zlabel('fL/D * M')
ax.title.set_text('Theoretical')

ax = fig.add_subplot(122, projection='3d')
Y = x_short
X = PRs_b
X,Y = np.meshgrid(X,Y)

surf = ax.plot_surface(X, Y, np.multiply(M_th,fL_D_ex), cmap='viridis')

# Add labels and title
ax.set_xlabel('PR')
ax.set_ylabel('x loc')
ax.set_zlabel('fL/D * M')
ax.title.set_text('Experimental')
# ax.set_ylim([x_short[8],x_short[-1]])


fig.suptitle('Mach Number Surface Plots')
fig.colorbar(surf)
plt.show()


# =============================================================================
# Surf f
# =============================================================================
fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
Y = x_short
X = PRs_b
X,Y = np.meshgrid(X,Y)

surf = ax.plot_surface(X, Y, f_th, cmap='viridis')
fig.colorbar(surf)
# Add labels and title
ax.set_xlabel('PR')
ax.set_ylabel('x loc')
ax.set_zlabel('f')
ax.title.set_text('Theoretical')

ax = fig.add_subplot(122, projection='3d')
Y = x_short
X = PRs_b
X,Y = np.meshgrid(X,Y)

surf = ax.plot_surface(X, Y, f_ex, cmap='viridis')
fig.colorbar(surf)

# Add labels and title
ax.set_xlabel('PR')
ax.set_ylabel('x loc')
ax.set_zlabel('f')
ax.title.set_text('Experimental')
# ax.set_ylim([x_short[8],x_short[-1]])


fig.suptitle('Frictional Factor Plots')
plt.show()


# =============================================================================
# Surf f and Mach experimental
# =============================================================================
fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
Y = x_short
X = PRs_b
X,Y = np.meshgrid(X,Y)

surf = ax.plot_surface(X, Y, M_ex, cmap='viridis')
# fig.colorbar(surf)
# Add labels and title
ax.set_xlabel('PR')
ax.set_ylabel('x loc')
ax.set_zlabel('Mach')
ax.title.set_text('Experimental Mach Number')

ax = fig.add_subplot(122, projection='3d')
Y = x_short
X = PRs_b
X,Y = np.meshgrid(X,Y)

surf = ax.plot_surface(X, Y, f_ex, cmap='viridis')
# fig.colorbar(surf)

# Add labels and title
ax.set_xlabel('P/Po')
ax.set_ylabel('x loc')
ax.set_zlabel('f')
ax.title.set_text('Experimental Frictional Factor')
# ax.set_ylim([x_short[8],x_short[-1]])


fig.suptitle('Experimental Frictional vs Mach')
plt.show()


# =============================================================================
# Surf f and Mach Theoretical
# =============================================================================
fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
Y = x_short
X = PRs_b
X,Y = np.meshgrid(X,Y)

surf1 = ax.plot_surface(X, Y, M_th, cmap='viridis')
# fig.colorbar(surf1)
# Add labels and title
ax.set_xlabel('PR')
ax.set_ylabel('x loc')
ax.set_zlabel('Mach')
ax.title.set_text('Theoretical Mach Number')

ax = fig.add_subplot(122, projection='3d')
Y = x_short
X = PRs_b
X,Y = np.meshgrid(X,Y)

surf2 = ax.plot_surface(X, Y, f_th, cmap='viridis')
# fig.colorbar(surf2)

# Add labels and title
ax.set_xlabel('P/Po')
ax.set_ylabel('x loc')
ax.set_zlabel('f')
ax.title.set_text('Theoretical Frictional Factor')
# ax.set_ylim([x_short[8],x_short[-1]])


fig.suptitle('Theoretical Frictional vs Mach')
plt.show()