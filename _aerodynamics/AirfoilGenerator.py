# -*- coding: utf-8 -*-
"""
Created on Sat Sep  6 14:30:00 2025

@author: BriceM
"""

'''
This module will be built using GMsh to generate NACA 4 and 5 digit airfoils
'''
import sys
import os
# Add parent import capabilities
parentdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if parentdir not in sys.path:
    sys.path.insert(0, parentdir)

import _tools.CoordinateTransformer as CT
import numpy as np
import matplotlib.pyplot as plt


def generateNACA4(NACA4, c=1, numPoints=100):
    '''
    Generates a list of points [x, y] that defines the full airfoil curve.

    Parameters
    ----------
    NACA4 : String
        The four digit string for a NACA airfoil.
        "0012" would be an airfoil with max camber of 0, at 0, with a thickness of 12% camber
    c : Float
        The cord length of the airfoil.
    numPoints : Int
        The number of points along which the airfoil will be generated.
        The entire curve will be defined by numPoints*2-2  points.
        
    Returns
    -------
    Points : npArray.
        The points in x-y pairs starting at the LE and ending back at the LE, but avoids double points
        

    '''
    M = int(NACA4[0]) /100
    P = int(NACA4[1]) / 10
    XX = int(NACA4[2:]) /100
    T = XX
    
    a0 = 0.2969 
    a1 = -0.126
    a2 = -0.3516 
    a3 = 0.2843
    a4 = -0.1036
    
    # Camber position
    def yc(x):
        # Before max camber position
        yc = np.zeros_like(x)
        
        for i, xi in enumerate(x):
            if 0 <= xi and xi < P:
                yc[i] = (M/P**2)*(2*P*xi - np.power(xi,2))
            elif P <= xi and xi <= 1:
                yc[i] = (M/(1-P)**2)*(1 - 2*P +2*P*xi - np.power(xi,2))
            else:
                raise ValueError("x-value exceeds the bounds [0,1]")
        return yc
    
    # Gradient
    def dyc_dx(x):
        # Before max camber position
        dyc_dx = np.zeros_like(x)
        
        for i, xi in enumerate(x):
            if 0 <= xi and xi < P:
                dyc_dx[i] = (2*M/P**2)*(P - xi)
            elif P <= xi and xi <= 1:
                dyc_dx[i] = (2*M/(1-P)**2)*(P - xi)
            else:
                raise ValueError("x-value exceeds the bounds [0,1]")
        return dyc_dx
    
    # Thickenss distribution (half the total thickness and needs to be applied to both
    # sides of the camber line)
    def yt(x):
        yt = (T/0.2)*(a0*np.power(x,0.5) + a1*x + a2*np.power(x,2) + a3*np.power(x,3) + a4*np.power(x,4))
        return yt
    
    
    beta = np.linspace(0, np.pi, 100, endpoint=True)
    x = (1-np.cos(beta))/2 
    
    theta = np.arctan(dyc_dx(x))
    
    xu = c*(x - np.multiply(yt(x), np.sin(theta)))
    xl = c*(x + np.multiply(yt(x), np.sin(theta)))
    
    yu = c*(yc(x) + np.multiply(yt(x), np.cos(theta)))
    yl = c*(yc(x) - np.multiply(yt(x), np.cos(theta)))
    
    points = np.zeros((len(beta)*2-2, 3))
    
    # for i, p in enumerate(points[:,0]):
    #     points[i, :] = np.array([xu[i]])
    
    # Set up the points on upper surface
    for i in range(0, len(xu)):
        points[i, :] = xu[i], yu[i], 0
        
    for i in range(1, len(xl)-1):
        points[i+len(xu)-1] = xl[len(xu)-1-i], yl[len(xu)-1-i], 0
        
    # reconnect to LE
    #points[-1, :] = xu[0], yu[0]
    
    return points #, xu,xl,yu,yl 

def generateGMSH_NACA4(geo, NACA4, dx=0, dy=0, dz=0, c=1, numPoints=100, rot_ang=None, mesh_size=0):
    '''
    Utilizes a gmsh object to generate the points, lines, and curve for a
    NACA 4-digit airfoil. 

    Parameters
    ----------
    geo : gmsh obj
        Pass in the following "gmsh.model.geo" after gmsh has been initialized.
    NACA4 : String
        The four digit string for a NACA airfoil..
    dx : Float, optional
        x-location of the LE of the airfoil. The default is 0.
    dy : Float, optional
        y-location of the LE of the airfoil. The default is 0.
    dz : Float, optional
        z-location of the LE of the airfoil. The default is 0.
    c : Float, optional
        Chord length. The default is 1.
    numPoints : Int, optional
        Number of points to define the airfoil points. The default is 100.

    Returns
    -------
    afcurve : Int
        Tag ID that relates to the curve of the airfoil.
    af_line_tags : List of Ints
        List of tag IDs that relate to the lines of the airfoil.
    af_point_tags : List of Ints
        List of tag IDs that relate to the points of the airfoil.

    '''
    af_point_tags = []
    af_line_tags = [] 
    
    pts = generateNACA4(NACA4, c, numPoints)
    if  rot_ang != None:
        # Rotate each row point about its own (0,0,0) aka LE
        for r in range(0,len(pts[:,0])):
            pts[r,:] = CT.Rotate(pts[r,:], *rot_ang)
        
            
            
    # make all the point objects
    for i in range(0, len(pts[:,0])):
        af_point_tags.append(geo.addPoint(pts[i, 0] + dx, pts[i,1] + dy, pts[i,2] + dz, mesh_size))
    
    # Make all of the line objects 
    for i in range(0, len(af_point_tags) -1):
        af_line_tags.append(geo.addLine(af_point_tags[i], af_point_tags[i+1]))
    # Need the final line to combine the first and last points
    af_line_tags.append(geo.addLine(af_point_tags[-1], af_point_tags[0]))    
        
    # Make them into a curve loop
    afcurve = geo.addCurveLoop(af_line_tags)
    
    return afcurve, af_line_tags, af_point_tags
    
    

if __name__=="__main__":
    pts = generateNACA4("8412", c=1)
    pts2 = generateNACA4("8412", c=2)
    print("plotting")
    plt.figure()
    # plt.plot(xU, yU)
    # plt.plot(xL, yL)
    plt.plot(pts[:,0], pts[:,1])
    plt.plot(pts2[:,0], pts2[:,1])
    plt.ylim([-0.6,0.6])
    plt.xlim([-0.1,2.1])
    plt.grid()
    plt.show()
    
    
    
    
    
    
    
    
    