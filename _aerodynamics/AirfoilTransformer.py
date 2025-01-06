# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 15:07:42 2023

@author: cycon
"""
import sys
import os
# Add parent import capabilities
parentdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if parentdir not in sys.path:
    sys.path.insert(0, parentdir)

import _tools.CoordinateTransformer as CT
import numpy as np
import matplotlib.pyplot as plt


class AirfoilTransformer():
    def __init__(self, airfoil_file_location, airfoil_file_name):
        self.airfoil_name = airfoil_file_name
        self.file_loc = airfoil_file_location
        self.airfoil = np.genfromtxt(airfoil_file_location+airfoil_file_name+'.txt')
        self.airfoil_trans = self.airfoil.copy()
        
    def Transform(self, scale=1, translations=[0,0,0], rotations=[0,0,0]):
        '''
        Transforms the airfoil from initialization by scaling by scale, translating by [dx, dy, dz]
        and rotating about [ax, ay, az] (in degrees). Saves transformed airfoil into self.airfoil_trans
        
        Parameters
        ----------
        scale : Float, optional
            Scales the Airfoil. The default is 1.
        translations : Float List, optional
            Translations for the airfoil in order [dx, dy, dz]. The default is [0,0,0].
        rotations : Float List, optional
            Rotations for the airfoil in order [ax, ay, az] in degrees. The default is [0,0,0].
        
        Returns
        -------
        None.
        
        '''
        # Translation
        [dx, dy, dz] = translations
        # Rotation
        [ax, ay, az] = rotations
        # Shape
        [R,C] = np.shape(self.airfoil)
        for r in range(0,R):
            v = self.airfoil[r,:]
            # Scale to needed size
            v = CT.Scale(v.copy(), scale,scale,scale)
            # Rotate
            v = CT.Rotate(v.copy(), ax, ay, az)
            # Translate
            v = CT.Translate(v.copy(), dx,dy,dz)
            
            self.airfoil_trans[r,:] = v.copy()

    def PlotAirfoils(self, margin=0.0):
        '''
        Plots the original and transformed airfoils on a 3D Plot to compair/verify correct rotaions

        Parameters
        ----------
        margin : Float, optional
            The axis limit margin as a percent of the max/min values of the airfoil. The default is 0.0.

        Returns
        -------
        None.

        '''
        fig = plt.figure("Airfoil Transformer")
        ax = fig.add_subplot(121, projection='3d')
        
        plt.title(self.airfoil_name)
        ax.plot(self.airfoil[:,0],self.airfoil[:,1],self.airfoil[:,2], label='Original')
        ax.plot(self.airfoil_trans[:,0],self.airfoil_trans[:,1],self.airfoil_trans[:,2], label='Transformed')
        plt.legend()
        
        max1 = max(np.max(self.airfoil), np.max(self.airfoil_trans))
        min1 = min(np.min(self.airfoil), np.min(self.airfoil_trans))
        
        margin += 1
        ax.set_xlim((margin*min1, margin*max1))            
        ax.set_ylim((margin*min1, margin*max1))
        ax.set_zlim((margin*min1, margin*max1))
        
    def SaveAsTXT(self, NewFN_tag="_Transformed"):
       with open(self.file_loc+NewFN_tag+'.txt', "w") as my_output_file:
               [ my_output_file.write("{:.4f} {:.4f} {:.4f}".format(row[0], row[1], row[2])+'\n') for row in self.airfoil_trans]
       my_output_file.close()
       
    def NormalizeAifoil(self, numPoints, detection_tolerance=1e-5):
        if numPoints%2 == 0:
            numPoints += 1 
            print('Warning: Number of points must be odd, automatically added one point')
        num_top = numPoints//2
        num_bot = num_top
        
        tol = detection_tolerance
        tol_v = np.array([tol, tol, tol])
        # Get indices of LE, TE and second LE
        [R,C] = np.shape(self.airfoil)
        for r in range(0,R):
            v = self.airfoil[r,:]
            if v + tol_v <= [1,0,0] and v - tol_v >= [1,0,0]:
                # LE
                LE_i = r
                break
        
       
'''       
# Operation Toggles
PlotAirfoils = True    # Plots the original and rotated point cloud on 3d Plot
SaveNewAirfoil = False # Saves pointcloud to a txt file

# Scale Uniformly:
currentChord = 1
neededChord = 5
scale = neededChord/currentChord # m

# Translate
dx = 0
dy = 0
dz = 0

# rotate
ax = 5
ay = 10
az = 0


FileLoc = 'C:\\Users\\cycon\\Documents\\Airfoils\\'

FileName = 'NACA 4212'
NewFN_tag = FileName + '-Root'

PointCloud = np.genfromtxt(FileLoc+FileName+'.txt')
[R,C] = np.shape(PointCloud)
new_PointCloud = PointCloud.copy()

for r in range(0,R):
    v = PointCloud[r,:]
    # Scale to needed size
    v = CT.Scale(v.copy(), scale,scale,scale)
    # Rotate
    v = CT.Rotate(v.copy(), ax, ay, az)
    # Translate
    v = CT.Translate(v.copy(), dx,dy,dz)
    
    
    new_PointCloud[r,:] = v.copy()

if SaveNewAirfoil:
    with open(FileLoc+NewFN_tag+'.txt', "w") as my_output_file:
            [ my_output_file.write("{:.4f} {:.4f} {:.4f}".format(row[0], row[1], row[2])+'\n') for row in new_PointCloud]
    my_output_file.close()


if PlotAirfoils:
    fig = plt.figure("Airfoil Transformer")
    ax = fig.add_subplot(121, projection='3d')
    
    ax.plot(PointCloud[:,0],PointCloud[:,1],PointCloud[:,2], label='Original')
    ax.plot(new_PointCloud[:,0],new_PointCloud[:,1],new_PointCloud[:,2], label='Transformed')
    plt.legend()
    
    max1 = max(np.max(PointCloud), np.max(new_PointCloud))
    min1 = min(np.min(PointCloud), np.min(new_PointCloud))
    
    margin = 1 + 0.2
    ax.set_xlim((margin*min1, margin*max1))            
    ax.set_ylim((margin*min1, margin*max1))
    ax.set_zlim((margin*min1, margin*max1))
'''







# Elipse Equation:
#   (x/a)^2 + (y/b)^2 = 1
#   a = radius in x dir
#   b = radius in y dir
# Angle of the Line from origin to point on curve:
#   theta = arctan(y/x)

# Leading edge and trailing edge will be at (x,y) along curve
# We need to scale each airfoil by the same amount so chord lengths are equal
# Need to rotate each airfoil by their respective theta based on chords
# Then translate to their final location

# Need total width and height of Ovular Wing:
# Span = 5     # m
# Height = 2.5 # m

# a = Span/2
# b = Height/2 

# def topCurve(x):
#     return np.sqrt(b*(1-(x/a)**2))
         
# def botCurve(x):
#     return -np.sqrt(b*(1-(x/a)**2))
  

# x = np.linspace(-a,a,100)
# y1 = topCurve(x)
# y2 = botCurve(x)

# # plt.plot(x,y1,x,y2)

# x_AfLocs = np.array([-a,-a/2,0,a/2,a])
# y1_AfLocs = topCurve(x_AfLocs)
# y2_AfLocs = botCurve(x_AfLocs)

# TopAirfoils = '4212'
# BotAirfoils = '4212'
