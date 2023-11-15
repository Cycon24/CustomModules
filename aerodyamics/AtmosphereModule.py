'''
Atmosphere class that allows for the access of atmoshpheric properties at any height within the atmosophere
'''
import numpy as np
import math
import matplotlib.pyplot as plt

class Atmosphere():
    
    # Selft initialization
    def __init__(self):
        # Create all atmosphereic properties
        Altitude = np.arange(0,65000,1000)
        
        g = 32.1740484 # ft/s2
        NM = 6076.4   # 1 Nautical Mile = ft
        
        # Define Sea Level values
        T_SL_R = 518.67 # °R
        P_SL = 2116.22  # lbs/ft^2
        rho_SL = 0.00237688  # slugs/ft^3
        a_SL = 1116.45  # ft/s
    
        self.Name = "Plane"
    
        # Atmosheric Array definitions
        self.Alt = Altitude
        self.Temp_R = np.zeros(np.size(Altitude))  # Temperature in Rankin
        self.Temp_F = np.zeros(np.size(Altitude))  # Temperature in Fahrenheit
        
        self.rho = np.zeros(np.size(Altitude))     # Density - slugs/ft^3
        self.P = np.zeros(np.size(Altitude))       # Pressure - lbs/ft^2
        
        self.TR = np.zeros(np.size(Altitude))      # θ Temperature Ratio
        self.PR = np.zeros(np.size(Altitude))      # δ Pressure Ratio
        self.DR = np.zeros(np.size(Altitude))      # σ Density Ratio
        
        self.sqrtDR = np.zeros(np.size(Altitude))  # Square root of Density Ratio
        self.QMS = np.zeros(np.size(Altitude))     # Dynamic Pressure / Mach^2
        self.specW = np.zeros(np.size(Altitude))   # Specific Weight of Air
        self.a = np.zeros(np.size(Altitude))       # Speed of sound ft/s
        self.VELA = np.zeros(np.size(Altitude))    # Speed of sound knots
        self.VRKIN = np.zeros(np.size(Altitude))   #
        
        # Array Calculations
        for i, h in enumerate(Altitude):
            if h <= 36000:
                self.Temp_R[i] = T_SL_R - 3.566*h/1000
                
                self.TR[i] = self.Temp_R[i] / T_SL_R
                self.PR[i] = self.TR[i]**5.2562

            else:
                self.Temp_R[i] = 389.99
                
                self.TR[i] = 0.7519
                self.PR[i] = 0.223361 * math.exp((-0.0481/1000)*(h-36089))
                
            self.DR[i] = self.PR[i] / self.TR[i]
            self.a[i] = a_SL*np.sqrt(self.TR[i])
            
            self.rho[i] = self.DR[i] * rho_SL
            self.P[i] = self.PR[i] * P_SL
            self.sqrtDR[i] = math.sqrt(self.DR[i])
            self.QMS[i] = 1481.4354 * self.PR[i]
            self.specW[i] = self.rho[i] * g
            self.VELA[i] = 3600 * self.a[i] / NM
            
        self.Temp_F = self.Temp_R - 459.67
        
    
    
    def linterp_h(self, h_range, target_property):
        '''

        Parameters
        ----------
        h_range : Float/Float Array
            Values of altitude in ft to calculate target property at.
        target_property : Atmospheric property array
            The property who's values are needed at the given altitudes.

        Returns
        -------
        Float/Float Array
            A value/array of the target property associated with the altitude inputs.

        '''
        # Linearly interpolates the value of the target property(y) at h(x)
        # Linterp: y = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        # Need to get the heights above and below inputted height\
        # can handle vectors
        if type(h_range) == float or type(h_range) == int:
            temp_h = np.empty((1,1), dtype='float')
            temp_h[0] = h_range
            h_range = temp_h
            
            
        propArray = np.zeros(np.size(h_range))
        
        for i_h, h in enumerate(h_range):
            found = False
            for idx, h_test in enumerate(self.Alt):
                if h_test == h:
                    # Values are equal so just return value at h
                    propArray[i_h] = target_property[idx]
                elif h_test >= h:
                    # We have passed the target altitude, so need to break from
                    # loop after saving values
                    h1 = self.Alt[idx-1]
                    h2 = h_test
                    y1 = target_property[idx-1]
                    y2 = target_property[idx]
                    
                    propArray[i_h] = y1 + (y2 - y1) * (h - h1) / (h2 - h1)
                    found = True
                    break
            if not found:
                print('\nError: Target Alitutde not found within Range. At {} ft and {} ft\n'.format(h, h_test))
        
        return propArray if np.size(h_range) > 1 else propArray[0]
    
    def linterp(x, x_array, y_array):
        # Linearly interpolates the value of the target property(y) at x
        # Linterp: y = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        
        for idx, x_test in enumerate(x_array):
            if x_test == x:
                # Values are equal so just return value at h
                return y_array[idx]
            elif x_test >= x:
                x1 = x_array[idx-1]
                x2 = x_test
                y1 = y_array[idx-1]
                y2 = y_array[idx]
                
                return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        print('\nError: Target Value not found within Range.\n')
                
if __name__=="__main__":
    atm = Atmosphere()
    print(atm.linterp_h(10000, atm.rho))
    # Po / Ps = (1 + gamma - 1)/2 * M^2