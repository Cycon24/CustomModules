# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 15:07:42 2023

@author: cycon
"""

import CoordinateTransformer as CT
import numpy as np
import matplotlib.pyplot as plt

# Scale Uniformly:
currentChord = 10
neededChord = 5
scale = neededChord/currentChord # m

# Translate
dx = 0
dy = 0
dz = 0

# rotate
ax = 0
ay = 0
az = 0


FileLoc = 'C:\\Users\\cycon\\OneDrive - University of Cincinnati\\Spring \
2023\\Aerodynamics\\Assignments\\Team Project\\AirfoilCurves\\'

FileName = 'NACA 2709'
NewFN_tag = FileName + '-Root'

PointCloud = np.genfromtxt(FileLoc+FileName+'.txt')
[R,C] = np.shape(PointCloud)

for r in range(0,R):
    v = PointCloud[r,:]
    # Scale to needed size
    v = CT.Scale(v.copy(), scale,scale,scale)
    # Rotate
    v = CT.Rotate(v.copy(), ax, ay, az)
    # Translate
    v = CT.Translate(v.copy(), dx,dy,dz)
    
    
    PointCloud[r,:] = v.copy()


with open(FileLoc+NewFN_tag+'.txt', "w") as my_output_file:
        [ my_output_file.write("{:.4f} {:.4f} {:.4f}".format(row[0], row[1], row[2])+'\n') for row in PointCloud]
my_output_file.close()


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
Span = 5     # m
Height = 2.5 # m

a = Span/2
b = Height/2 

def topCurve(x):
    return np.sqrt(b*(1-(x/a)**2))
         
def botCurve(x):
    return -np.sqrt(b*(1-(x/a)**2))
  

x = np.linspace(-a,a,100)
y1 = topCurve(x)
y2 = botCurve(x)

# plt.plot(x,y1,x,y2)

x_AfLocs = np.array([-a,-a/2,0,a/2,a])
y1_AfLocs = topCurve(x_AfLocs)
y2_AfLocs = botCurve(x_AfLocs)

TopAirfoils = '4212'
BotAirfoils = '4212'
