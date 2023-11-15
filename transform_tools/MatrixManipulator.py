# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 14:56:18 2023

@author: cycon
"""
import numpy as np
import math

## module error 
''' err(string).
Prints 'string' and terminates program.
'''
import sys
def err(string):
    print(string)
    input('Press return to exit')
    sys.exit()

## module swap
''' swapRows(v,i,j).
    Swaps rows i and j of a vector or matrix [v].
    
    swapCols(v,i,j).
    Swaps columns of matrix [v].
    
    swapCramer(a, b, i).
    Swaps i-th column of matrix [a] with array [b].
'''
def swapRows(v,i,j):
    if len(v.shape) == 1:
        v[i],v[j] = v[j],v[i]
    else:
        v[[i,j],:] = v[[j,i],:]
    return

def swapCols(v,i,j):
    v[:,[i,j]] = v[:,[j,i]]
    return

def swapCramer(a, b, i):
    import numpy as np
    ai = a.copy()
    ai[:, i] = np.transpose(b)
    return ai

# --------------- Print Function --------------------------- 
def printMat(mat, form='{:^12.2e}',StrReturn=False, **kwargs):
    mat = np.array(mat)
    if len(np.shape(mat)) == 1:
        newShape = (1,np.shape(mat)[0])
    for key, value in kwargs.items():
        if key == 'vecType':
            if value == 'col':
                newShape = (np.shape(mat)[0],1)
            elif value == 'row':
                newShape = (1,np.shape(mat)[0])

    if len(np.shape(mat)) == 1: # So data is in shape (nrows,) to change to (1,ncols=nrows) since we view [x1, x2, ... xn] as a 1row vec
        mat = np.reshape(mat, newShape)
    
    [R,C] = np.shape(mat)
    printStr = ''
    for r in range(0,R):
        printStr += '|'
        for c in range(0,C):
            printStr += form.format(mat[r,c])
        printStr +='|\n'
        
    if StrReturn:
        return printStr
    else:
        print(printStr)
        return None
        
# ------------ Books Guass Pivot ----------------
def gaussPivot(a_in, b_in,tol=1.0e-12):
    a = np.copy(a_in)
    b = np.copy(b_in)
    
    n = len(b)
    
    # Set up scale factors
    s = np.zeros(n)
    for i in range(n):
        s[i] = max(np.abs(a[i,:]))
    
    
    for k in range(0,n-1):
        # Row interchange, if needed
        p = np.argmax(np.abs(a[k:n,k])/s[k:n]) + k
        
        if abs(a[p,k]) < tol:  err("Matrix is singular")
        
        if p != k:
            swapRows(b,k,p)
            swapRows(s,k,p)
            swapRows(a,k,p)
            
        # Elimination
        for i in range(k+1,n):
            if a[i,k] != 0.0:
                lam = a[i,k]/a[k,k]
                a[i,k+1:n] = a[i,k+1:n] - lam*a[k,k+1:n]
                b[i] = b[i] - lam*b[k]
    
    if abs(a[n-1,n-1]) < tol:  err("Matrix is singular")
    
    # Back substitution
    b[n-1] = b[n-1]/a[n-1,n-1]
    for k in range(n-2,-1,-1):
        b[k] = (b[k] - np.dot(a[k,k+1:n],b[k+1:n]))/a[k,k]
    
    x = b
    
    return x

# ------------- Lu Pivot book ------------------
def LUdecomp(A,tol=1.0e-9):
    a = np.copy(A)
    n = len(a)
    seq = np.array(range(n))
    
    # Set up scale factors
    s = np.zeros((n))
    for i in range(n):
        s[i] = max(abs(a[i,:]))
        
    for k in range(0,n-1):
        # Row interchange, if needed
        p = np.argmax(np.abs(a[k:n,k])/s[k:n]) + k
        if abs(a[p,k]) < tol:  err('Matrix is singular')
        if p != k:
             swapRows(s,k,p)
             swapRows(a,k,p)
             swapRows(seq,k,p)
   
    
        # Elimination
        for i in range(k+1,n):
            if a[i,k] != 0.0:
                lam = a[i,k]/a[k,k]
                a[i,k+1:n] = a[i,k+1:n] - lam*a[k,k+1:n]
                a[i,k] = lam
                
    return a,seq

def LUsolve(A,B,seq):
    a = np.copy(A)
    b = np.copy(B)
    
    n = len(a)
    # Rearrange constant vector; store it in [x]
    x = b.copy()
    
    for i in range(n):
        x[i] = b[seq[i]]
        
    # Solution
    for k in range(1,n):
        x[k] = x[k] - np.dot(a[k,0:k], x[0:k])
    
    x[n-1] = x[n-1]/a[n-1,n-1]
    
    for k in range(n-2,-1,-1):
        x[k] = (x[k] - np.dot(a[k,k+1:n], x[k+1:n]))/a[k,k]
        
    return x

def LUpivot(A,b):
    a, seq = LUdecomp(A)
    x = LUsolve(a,b,seq)
    return x


## module LUdecomp3
''' c,d,e = LUdecomp3(c,d,e).
LU decomposition of tridiagonal matrix [c\d\e]. On output
{c},{d} and {e} are the diagonals of the decomposed matrix.
x = LUsolve(c,d,e,b).
Solves [c\d\e]{x} = {b}, where {c}, {d} and {e} are the
vectors returned from LUdecomp3.
'''
def LUdecomp3(c,d,e):
    n = len(d)
    for k in range(1,n):
        lam = c[k-1]/d[k-1]
        d[k] = d[k] - lam*e[k-1]
        c[k-1] = lam
    return c,d,e

def LUsolve3(c,d,e,b):
    n = len(d)
    for k in range(1,n):
        b[k] = b[k] - c[k-1]*b[k-1]
    
    b[n-1] = b[n-1]/d[n-1]
    for k in range(n-2,-1,-1):
        b[k] = (b[k] - e[k]*b[k+1])/d[k]
    return b


# --------------- My own Cramer ------------------
def Cramer(A, b):
    # Need to get det of A_k and A and x_k = det(A_k)/det(A)
    # Where A_k is A but with colummn k substituted with b
    n = len(b)
    x = np.zeros((n,1))
    
    det_A = np.linalg.det(A)
    
    for k in range(0,n):
        # Insert b into column k of A
        A_k = swapCramer(A, b, k)
        
        #Calculate det of new matrix
        det_Ak = np.linalg.det(A_k)
        
        # Calculate x by det(Ak) / det(A)
        x[k] = det_Ak / det_A
    return x

# ---------------- Book's gaussSeidel -----------------
'''
x,numIter,omega = gaussSeidel(iterEqs,x,tol = 1.0e-9)
Gauss-Seidel method for solving [A]{x} = {b}.
The matrix [A] should be sparse. User must supply the
function iterEqs(x,omega) that returns the improved {x},
given the current {x} (’omega’ is the relaxation factor).
'''

# Needs user to provide the function iterEqs that computes improved x
def gaussSeidel(iterEqs,xi,tol = 1.0e-9):
    # xi is initial solution guess
    x = xi.copy()
    
    omega = 1.0
    k = 10
    p=1
    
    for i in range(1,501):
        # Store last iteration value
        xOld = x.copy()
        
        # Create value and calculate dx
        x = iterEqs(x,omega)
        dx = math.sqrt(np.dot(x-xOld,x-xOld))
        
        # checking to see if we are within tolerance -> converging on sol
        if dx < tol: return x,i,omega
        
        # Compute relaxation factor after k+p iterations
        if i == k: dx1 = dx
        if i == k + p:
            dx2 = dx
            omega = 2.0/(1.0 + math.sqrt(1.0 - (dx2/dx1)**(1.0/p)))
    print('Gauss-Seidel failed to converge')
    return None, None, None



# ---------------- Book's conjGrad ---------------------

'''
 x, numIter = conjGrad(Av,x,b,tol=1.0e-9)
Conjugate gradient method for solving [A]{x} = {b}.
92 Systems of Linear Algebraic Equations
The matrix [A] should be sparse. User must supply
the function Av(v) that returns the vector [A]{v}.
'''
def conjGrad(Av,xi,bi,tol=1.0e-9):
    x = xi.copy()
    b = bi.copy()
    
    n = len(b)
    r = b - Av(x)
    s = r.copy()
    
    for i in range(n):
        u = Av(s)
        alpha = np.dot(s,r)/np.dot(s,u)
        x = x + alpha*s
        r = b - Av(x)
        if(math.sqrt(np.dot(r,r))) < tol:
            break
        else:
            beta = -np.dot(r,u)/np.dot(s,u)
            s = r + beta*s
    return x,i


def accuracy(A, x, b):
    acc = np.linalg.norm(np.dot(A,x) - b)
    return acc



# ------------ Matrix builders ----------------------------------------

def matBuilder(n, mainD, subD, n1s=0, pattern=None):
    PATTERN_IS_CONSTANT = 0
    # Creates a sparese matrix with a values at A_n1 and A_1n, the same value
    # along A_ij when i=j, and the same value above the diagnal and another below 
    # the diagonal
    
    # Setting up for diagonal values
    if type(subD) == np.ndarray or type(subD) == list:
        SUBD_ARRAY = True
        numSDs = len(subD)
    else:
        subD = np.array([subD])
        numSDs = 1
        
    # Setting up for patterns for diag values
    if type(pattern) == np.ndarray or type(pattern) == list:
        NO_PATTERN = False
        numPats = len(pattern)
        
        # Have to fill pattern with values 
        while numPats < numSDs:
            pattern.append(PATTERN_IS_CONSTANT)
            numPats = len(pattern)
    elif pattern == None:
        NO_PATTERN = True
    else:
        NO_PATTERN = False
        pattern = np.array([pattern])
        numPats = 1
        
    # Create empty matrix A nxn
    mat = np.zeros((n,n),dtype=float)
    
    # Populate the matrix
    for i in range(0,n): # Row
        for j in range(0,n): # Col
            #Main Diagonal
            if i==j:
                mat[i,j] = mainD
                
            # Iterate through for sub diagonals
            for k in range(1, numSDs + 1):
                if j == i + k or j == i - k:
                    if NO_PATTERN:
                         # subDiagnals diagonal
                        mat[i,j] = subD[k-1]
                    else:
                        # above diag check, below diag check
                        if pattern[k-1] == PATTERN_IS_CONSTANT:
                            # subDiagnals diagonal
                            mat[i,j] = subD[k-1]
                        elif (j > i and (i+1)%pattern[k-1] == 0) or (i > j and (j+1)%pattern[k-1] == 0): 
                            mat[i,j] = 0
                        else:
                            # subDiagnals diagonal
                            mat[i,j] = subD[k-1]
                   
    mat[n-1, 0] = n1s
    mat[0, n-1] = n1s
    return mat


# -------------- Testing --------------------------
if __name__ == "__main__":
    print('Beebop')
    
    A = np.array([[3, 7, -2],
                  [8, 0, 0],
                  [5, -1, -1]],dtype=float)
    
    b = np.array([[0], [16], [7]],dtype=float)
    # X should be [2, 0, 3]
    
    A = np.array([[1,1],[2,-3]],dtype=float)
    
    b = np.array([[5],[-4]],dtype=float)

    print("A:\n{}".format(A))
    print("b:\n{}".format(b))
    
    x_G = gaussPivot(A,b)
    
    A_L, seq = LUdecomp(A)
    x_L = LUsolve(A_L, b, seq)
    
    x_C = Cramer(A, b)
    
    print("Gauss x: \n{}".format(x_G))
    print("LU x: \n{}".format(x_L))
    print('Cramer:\n{}'.format(x_C))
    