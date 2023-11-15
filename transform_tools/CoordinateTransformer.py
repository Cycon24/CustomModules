# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 14:31:47 2023

@author: cycon
"""
import numpy as np 

def Rotate(Vector, Alpha, Beta, Phi): # Angles: X Y Z

    a = np.radians(Alpha)
    b = np.radians(Beta)
    p = np.radians(Phi)    

    Rx = np.array([[1, 0, 0],\
                   [0, np.cos(a), -np.sin(a)],\
                   [0, np.sin(a), np.cos(a)]])
    Ry = np.array([[np.cos(b), 0, np.sin(b)],\
                   [0, 1, 0],\
                   [-np.sin(b), 0, np.cos(b)]])
    Rz = np.array([[np.cos(p), -np.sin(p),0],\
                   [np.sin(p), np.cos(p),0],\
                   [0, 0, 1]])
        
    R = np.matmul(np.matmul(Rx,Ry),Rz)
    
    return np.matmul(R, Vector)

def Translate(Vector, dx, dy, dz):
    NewVec = Vector.copy()
    NewVec[0] += dx
    NewVec[1] += dy
    NewVec[2] += dz
    return NewVec

def Scale(Vector, sx=1, sy=1, sz=1):
    NewVec = Vector.copy()
    NewVec[0] *= sx
    NewVec[1] *= sy
    NewVec[2] *= sz
    return NewVec

if __name__=='__main__':
    v = np.array([1,0,0])
    print(Scale(v, 2,1,6))
    
    
    
    