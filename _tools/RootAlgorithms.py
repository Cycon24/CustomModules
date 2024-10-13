# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 12:58:21 2023

@author: cycon
"""

import numpy as np
import math
import error
from random import random
import cmath
from MatrixManipulator import gaussPivot

# --------- Root Search ------------------
'''
x1,x2 = rootsearch(f,a,b,dx).
Searches the interval (a,b) in increments dx for
the bounds (x1,x2) of the smallest root of f(x).
Returns x1 = x2 = None if no roots were detected.
'''
def rootsearch(f,a,b,dx):
    x1 = a; f1 = f(a)
    x2 = a + dx; f2 = f(x2)
    while np.sign(f1) == np.sign(f2):
        if x1 >= b: 
            return 'x1={}'.format(x1),'b={}'.format(b)# None,None
        x1 = x2; f1 = f2
        x2 = x1 + dx; f2 = f(x2)
    return x1,x2


# --------- Bisection --------------------
'''
root = bisection(f,x1,x2,switch=0,tol=1.0e-9).
Finds a root of f(x) = 0 by bisection.
The root must be bracketed in (x1,x2).
Setting switch = 1 returns root = None if
f(x) increases upon bisection.
'''
def bisection(f,x1,x2,switch=1,tol=1.0e-9):
    f1 = f(x1)
    if f1 == 0.0: return x1
    
    f2 = f(x2)
    if f2 == 0.0: return x2
    
    if np.sign(f1) == np.sign(f2):
        error.err('Root is not bracketed')
    n = int(math.ceil(math.log(abs(x2 - x1)/tol)/math.log(2.0)))
    
    for i in range(n):
        x3 = 0.5*(x1 + x2); f3 = f(x3)
        if (switch == 1) and (abs(f3) > abs(f1)) and (abs(f3) > abs(f2)):
            return None
        if f3 == 0.0: return x3
        if np.sign(f2)!= np.sign(f3): x1 = x3; f1 = f3
        else: x2 = x3; f2 = f3
    return (x1 + x2)/2.0

# --------- Newton Raphson ---------------
'''
root = newtonRaphson(f,df,a,b,tol=1.0e-9).
Finds a root of f(x) = 0 by combining the Newton-Raphson
method with bisection. The root must be bracketed in (a,b).
Calls user-supplied functions f(x) and its derivative df(x).
'''

def newtonRaphson(f,df,a,b,tol=1.0e-9):
    fa = f(a)
    if fa == 0.0: return a
    fb = f(b)
    if fb == 0.0: return b
    if np.sign(fa) == np.sign(fb): error.err('Root is not bracketed')
    x = 0.5*(a + b)
    
    for i in range(30):
        fx = f(x)
        if fx == 0.0: return x
        
        # Tighten the brackets on the root
        if np.sign(fa) != np.sign(fx): b = x
        else: a = x
        
        # Try a Newton-Raphson step
        dfx = df(x)
        
        # If division by zero, push x out of bounds
        try: dx = -fx/dfx
        except ZeroDivisionError: dx = b - a
        x = x + dx
        
        # If the result is outside the brackets, use bisection
        if (b - x)*(x - a) < 0.0:
            dx = 0.5*(b - a)
            x = a + dx
            
        # Check for convergence
        if abs(dx) < tol*max(abs(b),1.0): return x
    print('Too many iterations in Newton-Raphson')

# --------- Newton Raphson 2 --------------
''' soln = newtonRaphson2(f,x,tol=1.0e-9).
Solves the simultaneous equations f(x) = 0 by
the Newton-Raphson method using {x} as the initial
guess. Note that {f} and {x} are vectors.
'''
def newtonRaphson2(f,x,tol=1.0e-9):
    def jacobian(f,x):
        h = 1.0e-4
        n = len(x)
        jac = np.zeros((n,n))
        f0 = f(x)
        
        for i in range(n):
            temp = x[i]
            x[i] = temp + h
            f1 = f(x)
            x[i] = temp
            jac[:,i] = (f1 - f0)/h
        return jac,f0
    
    for i in range(30):
        jac,f0 = jacobian(f,x)
        if math.sqrt(np.dot(f0,f0)/len(x)) < tol: return x
        dx = gaussPivot(jac,-f0)
        x = x + dx
        if math.sqrt(np.dot(dx,dx)) < tol*max(max(abs(x)),1.0):
            return x
    print('Too many iterations in Newton-Raphson 2')


# --------- Poly Roots -----------------------
''' p,dp,ddp = evalPoly(a,x).
Evaluates the polynomial
p = a[0] + a[1]*x + a[2]*xˆ2 +...+ a[n]*xˆn
with its derivatives dp = p’ and ddp = p’’
at x.
'''
def evalPoly_roots(a,x):
    n = len(a) - 1
    p = a[n]
    dp = 0.0 + 0.0j
    ddp = 0.0 + 0.0j
    for i in range(1,n+1):
        ddp = ddp*x + 2.0*dp
        dp = dp*x + p
        p = p*x + a[n-i]
    return p,dp,ddp

def deflPoly(a,root): # Deflates a polynomial
        n = len(a)-1
        b = [(0.0 + 0.0j)]*n
        b[n-1] = a[n]
        for i in range(n-2,-1,-1):
            b[i] = a[i+1] + root*b[i+1]
        return b

def polyRoots(a,tol=1.0e-12):
    def laguerre(a,tol):
        x = random() # Starting value (random number)
        n = len(a) - 1
        for i in range(30):
            p,dp,ddp = evalPoly_roots(a,x)
            if abs(p) < tol: return x
            g = dp/p
            h = g*g - ddp/p
            f = cmath.sqrt((n - 1)*(n*h - g*g))
            if abs(g + f) > abs(g - f): dx = n/(g + f)
            else: dx = n/(g - f)
            x = x - dx
            if abs(dx) < tol: return x
        print('Too many iterations')
    
    
    n = len(a) - 1
    roots = np.zeros((n),dtype=complex)
    for i in range(n):
        x = laguerre(a,tol)
        if abs(x.imag) < tol: x = x.real
        roots[i] = x
        a = deflPoly(a,x)
    return roots


# -------- Testing ----------------------------
if __name__ == '__main__':
    print('')
    e = 1.42
    Mh = 7.3284
    
    def f(F): return e*np.sinh(F)-F-Mh
    def df(F): return e*np.cosh(F) - 1
    
    low_r, high_r = rootsearch(f,0.001, 4,0.001)
    F = newtonRaphson(f, df, low_r, high_r)
    