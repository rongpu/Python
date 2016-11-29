from __future__ import division, print_function
import numpy as np
import statsmodels.api as sm
from scipy.interpolate import interp1d

def hlmean(data,nsamp=-1):
    '''
    PURPOSE:

     Calculate Hodges-Lehmann estimator of mean, using nsamp
      bootstraps from the data

    The Hodges-Lehmann estimator is, formally, the median value of
     (x_i+x_j)/2 over all pairs of indices i,j .

    If there are less than 200 samples, we compute the formal Hodge-
    Lehmann mean. Otherwise, we estimate that quantity using nsamp 
    randomly chosen values of i & j, rather than using all possible 
    values.

    Although it has much of the robustness of an ordinary median, the
    H-L estimator yields much smaller errors (equivalent to the mean of
    >90% as much data, while the median has errors equivalent to the
    standard error of 64% as much data).
    CALLING SEQUENCE:

     result=hlmean(data [,nsamp=nsamp])

    INPUTS:

     data: array of values to calculate the H-L mean of

    OPTIONAL KEYWORD PARAMETERS:

     nsamp : if set, hlmean will use this number of bootstrap
     samples to do the calculation.  If not set, it will use 200*the
     number of elements of the data array random values of i&j.

    OUTPUTS:

     result: bootstrap estimate of the H-L mean estimator

    EXAMPLE:
       test=[1,2,0,1,2,20.]
       print hlmean(test)
    '''

    ndata=len(data)
    if ndata==0:
        print('H-L mean: empty array!!!')
    if ndata < 200:
        pairmean=np.zeros(int(ndata*(ndata+1)/2))
        index=0
        for i in range(ndata):
            for j in range(i,ndata):
                pairmean[index]=(data[i]+data[j])/2
                index+=1
    else:
        if nsamp < 0:
            nsamp=200*ndata
        idx=np.floor(np.random.rand(nsamp,2)*ndata)
        idx=idx.astype(np.int64,copy=False)
        pairmean=np.sum(data[idx],axis=1)/2.
    return(np.median(pairmean))

# ---------------------------------------------------

def mad(a):
    '''
    Calculate the median absolute deviation (MAD).
    
    INPUT:
    1D array of values to calculate MAD of.
    
    OUTPUT:
    MAD of the 1D array.
    '''
    
    medianvalue=np.median(a)
    return np.median(np.abs(a-medianvalue))

# ---------------------------------------------------

# median absolute value (mav)
def mav(a):
    return np.median(np.abs(a))
    
# ---------------------------------------------------

def rlm_fit1d(x, y, t=1.5, order=1):
    '''
    1D robust polynomial fit.
    
    Given x array and y array, calculate the 1D robust 
    polynomial fit of arbitrary order. Huber weight
    function is used. 
    
    See also poly_val1d.py
    
    INPUT:
    1D arrays of x and y values; tunning parameter t; 
    order of the polynomial fit.
    
    OUTPUT:
    Array of parameters of the polynomial [a0, a1, a2 ...] 
    so that y = a0 + a1*x + a2*x**2 + ...
    '''
    
    ncols = order+1
    a = np.zeros((x.size,ncols))
    for i in range(order+1):
        a[:,i] = x**i
    res = sm.RLM(y, a, M=sm.robust.norms.HuberT(t=t)).fit()
    m = res.params
    return(m)

# ---------------------------------------------------

def poly_val1d(x, m):
    '''
    Evaluate the 1D polynomial from x values and polynomial parameters
    
    See also rlm_fit1d.py
    
    INPUT:
    1D array of x values; 
    1D array of polynomial parameters (for example generated by 
    rlm_fit1d.py).
    
    OUTPUT:
    1D array of the evaluated values of the polynomial.
    '''
    
    order = len(m)-1
    z = np.zeros_like(x)
    for i in range(order+1):
        z += m[i] * x**i
    return z

# ---------------------------------------------------

def rlm_fit1d_more(x, y, t=1.5, order=1):
    '''
    Same as rlm_fit1d_more, except that instead of only returning the 
    fitting parameters, the full regression result is return. 
    '''
    
    ncols = order+1
    a = np.zeros((x.size,ncols))
    for i in range(order+1):
        a[:,i] = x**i
    res = sm.RLM(y, a, M=sm.robust.norms.HuberT(t=t)).fit()
    return(res)

# ---------------------------------------------------

def poly_val1d_more(x, res):
    '''
    Same as rlm_fit1d_more, except that instead of only taking the 
    fitting parameters as argument, the full regression result taken. 
    '''
    
    m = res.params
    order = len(m)-1
    z = np.zeros_like(x)
    for i in range(order+1):
        z += m[i] * x**i
    return z

# ---------------------------------------------------

def rlm_fit2d(x, y, z, t=1.5, order=2):
    '''
    2D robust polynomial fit.
    
    Given x, y and z arrays, calculate the 2D robust 
    polynomial fit of arbitrary order. Huber weight
    function is used. 
    
    See also poly_val2d.py
    
    INPUT:
    1D arrays of x, y and z values; tunning parameter t; 
    order of the polynomial fit. 
    Note that x*y is a second order term.
    
    OUTPUT:
    Array of parameters of the polynomial [a0, a1, a2 ...] 
    so that for an n-th order fit
    z = a0 + a1*y + a2*y**2 + ... + a_n*y**n + a_(n+1)*x 
    + a_(n+2)*x*y + a_(n+3)*x*y**2 + ... + a_(2n)*x*y**(n-1) 
    + a_(2n+1)*x**2 + a_(2n+2)*x**2*y + ... ... 
    
    For example 2nd order fit gives:
    z = a0 + a1*y + a2*y**2 + a3*x + a4*x*y + a5*x**2
    
    New values can be evaluated using poly_val2d.py
    '''
    
    ncols = (order+2)*(order+1)/2
    a = np.zeros((x.size,ncols))
    k=0
    for i in range(order+1):
        for j in range(order-i+1):
            a[:,k] = x**i * y**j
            k+=1
    res = sm.RLM(z, a, M=sm.robust.norms.HuberT(t=t)).fit()
    m = res.params
    return(m)

# ---------------------------------------------------

def poly_val2d(x, y, m):
    '''
    Evaluate the 2D polynomial from x and yvalues and polynomial 
    parameters
    
    See also rlm_fit2d.py
    
    INPUT:
    1D array of x and y values; 
    1D array of polynomial parameters (generated by rlm_fit2d.py). 
    
    OUTPUT:
    1D array of the evaluated values of the polynomial.    
    '''
    
    order = int((np.sqrt(8*len(m)+1)-3)/2)
    z = np.zeros_like(x)
    k=0
    for i in range(order+1):
        for j in range(order-i+1):
            z += m[k] * x**i * y**j
            k+=1
    return z

# ---------------------------------------------------


def extrap1d(x, y, kind='linear', extrap='constant'):
    '''
    Interpolation and extrapolation in one function. 
    Constant (default) and linear extrapolation are supported. 

    INPUT:  
    1D array of x and y values;
    kind: interpolation method in scipy.interpolate.interp1d
    extrap: extrapolation method - "linear" or "constant"

    OUTPUT:
    Extrapolated function f(x). 

    '''
    interp = interp1d(x,y, kind=kind)

    if (extrap!='constant') & (extrap!='linear'):
        print('Invalid extrapolation method "'+extrap+'". The default linear extrapolation is used.')
    if extrap=='linear':
        def pointwise(xx):
            if xx < x[0]:
                return y[0]+(xx-x[0])*(y[1]-y[0])/(x[1]-x[0])
            elif xx > x[-1]:
                return y[-1]+(xx-x[-1])*(y[-1]-y[-2])/(x[-1]-x[-2])
            else:
                return interp(xx)
    else:
        def pointwise(xx):
            if xx < x[0]:
                return y[0]
            elif xx > x[-1]:
                return y[-1]
            else:
                return interp(xx)

    def ufunclike(x):
        if isinstance(x, np.ndarray):
            return np.array(map(pointwise, np.array(x)))
        elif isinstance(x, list):
            return np.array(map(pointwise, np.array(x)))
        else:
            return np.array(map(pointwise, [x]))

    return ufunclike
