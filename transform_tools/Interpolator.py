import numpy as np
import matplotlib.pyplot as plt
import MatrixManipulator as MM
import math

# ----------- Newton Poly -----------
''' p = evalPoly(a,xData,x).
Evaluates Newton’s polynomial p at x. The coefficient
vector {a} can be computed by the function ’coeffts’.
a = coeffts(xData,yData).
Computes the coefficients of Newton’s polynomial.
p(x) = c[0] + c[1]x + c[2]xˆ2 +...+ c[m]xˆm
'''
def newtonPoly_eval(a,xData,x):
    n = len(xData) - 1 # Degree of polynomial
    p = a[n]
    for k in range(1,n+1):
        p = a[n-k] + (x -xData[n-k])*p
    return p

def newtonPoly_coefs(xD,yD):
    xData = np.array(xD)
    yData = np.array(yD)
    m = len(xData) # Number of data points
    a = yData.copy()
    for k in range(1,m):
        a[k:m] = (a[k:m] - a[k-1])/(xData[k:m] - xData[k-1])
    return a

# Possibly doesnt work?
def newtonPoly_str(xData, yData, cfs=None):
    if cfs.any() == None:
        # returns a string of the polynomial
        cfs = newtonPoly_coefs(xData,yData)
        
    cfs_str = ''
    n = len(cfs)
    cfs_str += '{} '.format(cfs[0]) 
    
    for i in range(1,n):
        if cfs[i] == 1:
            cfs_str += '+ x^{} '.format(i) 
        elif cfs[i] != 0:
            cfs_str += '+ {}x^{} '.format(cfs[i],i) 
            
    return cfs_str

def poly_str(cfs, x='x',y='f(x) =',formatting='({:.3E})', PRINT=False,):
    cfs_str = '{} '.format(y)
    n = len(cfs)
    post_form_str = '{} '.format(formatting)
    cfs_str += post_form_str.format(cfs[0]) 
    
    for i in range(1,n):
        if cfs[i] == 1:
            cfs_str += '+ {}^{} '.format(x,i)
        elif cfs[i] != 0:
            temp = '+ {}'.format(formatting)
            temp += '{}^{} '
            cfs_str += temp.format(cfs[i],x,i) 
    if PRINT:
        print(cfs_str)
    return cfs_str

# ----------- Rational --------------
''' p = rational(xData,yData,x)
Evaluates the diagonal rational function interpolant p(x)
that passes through the data points
'''
def rational(xData,yData,x):
    m = len(xData)
    r = yData.copy()
    rOld = np.zeros(m)
    for k in range(m-1):
        for i in range(m-k-1):
            if abs(x - xData[i+k+1]) < 1.0e-9:
                return yData[i+k+1]
            else:
                c1 = r[i+1] - r[i]
                c2 = r[i+1] - rOld[i+1]
                c3 = (x - xData[i])/(x - xData[i+k+1])
                r[i] = r[i+1] + c1/(c3*(1.0 - c1/c2) - 1.0)
                rOld[i+1] = r[i+1]
    return r[0]


# ----------- Cubic Spline ----------
'''
k = curvatures(xData,yData).
Returns the curvatures of cubic spline at its knots.
y = evalSpline(xData,yData,k,x).
Evaluates cubic spline at x. The curvatures k can be
computed with the function ’curvatures’.
'''

def cubic_curvatures(xData,yData):
    n = len(xData) - 1
    c = np.zeros(n)
    d = np.ones(n+1)
    e = np.zeros(n)
    k = np.zeros(n+1)
    
    c[0:n-1] = xData[0:n-1] - xData[1:n]
    d[1:n] = 2.0*(xData[0:n-1] - xData[2:n+1])
    e[1:n] = xData[1:n] - xData[2:n+1]
    k[1:n] =6.0*(yData[0:n-1] - yData[1:n]) \
            /(xData[0:n-1] - xData[1:n]) \
            -6.0*(yData[1:n] - yData[2:n+1]) \
            /(xData[1:n] - xData[2:n+1])
    c,d,e = MM.LUdecomp3(c,d,e)
    k = MM.LUsolve3(c,d,e,k)
    return k


def cubic_evalSpline(xData,yData,k,x):
    def findSegment(xData,x):
        iLeft = 0
        iRight = len(xData)- 1
        while 1:
            if (iRight-iLeft) <= 1: return iLeft
            i = (iLeft + iRight)//2        # Have to add // in order to cast to int
            if x < xData[i]: iRight = i
            else: iLeft = i
    
    i = findSegment(xData,x)
    h = xData[i] - xData[i+1]
    y = ((x - xData[i+1])**3/h - (x - xData[i+1])*h)*k[i]/6.0 \
            - ((x - xData[i])**3/h - (x - xData[i])*h)*k[i+1]/6.0 \
            + (yData[i]*(x - xData[i+1]) \
            - yData[i+1]*(x - xData[i]))/h
    return y



# -----------Poly Fit ---------------
''' c = polyFit(xData,yData,m).
Returns coefficients of the polynomial
p(x) = c[0] + c[1]x + c[2]xˆ2 +...+ c[m]xˆm
that fits the specified data in the least
squares sense.
sigma = stdDev(c,xData,yData).
Computes the std. deviation between p(x)
and the data.
'''

def polyFit(xData,yData,m):
    a = np.zeros((m+1,m+1))
    b = np.zeros(m+1)
    s = np.zeros(2*m+1)
    for i in range(len(xData)):
        temp = yData[i]
        for j in range(m+1):
            b[j] = b[j] + temp
            temp = temp*xData[i]
        temp = 1.0
        for j in range(2*m+1):
            s[j] = s[j] + temp
            temp = temp*xData[i]
    for i in range(m+1):
        for j in range(m+1):
            a[i,j] = s[i+j]
    return MM.gaussPivot(a,b)

def evalPoly(c,x):
        m = len(c) - 1
        p = c[m]
        for j in range(m):
            p = p*x + c[m-j-1]
        return p
    
def stdDev(c,xData,yData):
    n = len(xData) - 1
    m = len(c) - 1
    sigma = 0.0
    for i in range(n+1):
        p = evalPoly(c,xData[i])
        sigma = sigma + (yData[i] - p)**2
        
    sigma = math.sqrt(sigma/(n - m))
    return sigma

# ----------- Plot Poly -------------

def plotPoly(xData,yData,coeff,xlab='x',ylab='y', title=''):
    m = len(coeff)
    x1 = min(xData)
    x2 = max(xData)
    dx = (x2 - x1)/20.0
    x = np.arange(x1,x2 + dx/10.0,dx)
    y = np.zeros((len(x)))*1.0
    for i in range(m):
        y = y + coeff[i]*x**i
        
    plt.plot(xData,yData,'o',x,y,'-')
    plt.xlabel(xlab); plt.ylabel(ylab)
    plt.title(title)
    plt.grid (True)
    plt.show()


# ------ Main Testing ---------


def main():
    # # Test Cubic Spline
    # xData = np.array([1,2,3,4,5],float)
    # yData = np.array([0,1,0,1,0],float)
    # k = cubic_curvatures(xData,yData)
    # while True:
    #     try: x = eval(input("\nx ==> "))
    #     except SyntaxError: break
    #     print("y =",cubic_evalSpline(xData,yData,k,x))
    # input("Done. Press return to exit")
    xD = [0,1,2,3]
    yD = [1,2,5,10]
    
    cfs = newtonPoly_coefs(xD,yD)
    cfs2 = polyFit(xD,yD,2)
    
    polStr = newtonPoly_str(0,0,cfs)
    polStr2 = newtonPoly_str(0,0,cfs2)
    
    poly_str(cfs, 'T','Cp(T)', True)
    
    plotPoly(xD,yD,cfs2)
    return None

    

if __name__ == "__main__":
    main()

