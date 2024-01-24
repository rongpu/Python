from __future__ import division, print_function
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


def hlmean(data, multiplier=None, verbose=True):
    '''
    Hodges-Lehmann estimator.
    '''

    ndata = len(data)
    if ndata==0 and verbose:
        print('H-L mean: empty array!!!')
    if ndata==1:
        return data[0]
    if ndata < 200:
        pairmean = np.zeros(int(ndata*(ndata+1)/2))
        index = 0
        for i in range(ndata):
            for j in range(i, ndata):
                pairmean[index] = (data[i]+data[j])/2
                index += 1
    else:
        if multiplier is None:
            nsamp = 200 * ndata
        else:
            nsamp = multiplier * ndata
        idx = np.floor(np.random.rand(nsamp, 2)*ndata)
        idx = idx.astype(np.int64, copy=False)
        pairmean = np.sum(data[idx], axis=1)/2.
    return(np.median(pairmean))


def hlmean_v2(data, maxpairs=1e8, random_seed=None, verbose=True):
    '''
    Hodges-Lehmann estimator.
    '''

    import itertools

    maxpairs = int(maxpairs)
    ndata = len(data)

    if ndata==0:
        if verbose:
            print('H-L mean: empty array!!!')
        return None
    if ndata==1:
        return data[0]
    if ndata*(ndata-1)/2 <= maxpairs:
        # only non-identical indices are included
        idx1, idx2 = np.array(list(itertools.combinations(np.arange(ndata), 2))).transpose()
        pairmean1 = (data[idx1]+data[idx2])/2.
        # the pairwise mean of identical indices
        pairmean2 = data
        pairmean = np.concatenate([pairmean1, pairmean2])
        hlmean_value = np.median(pairmean)
    else:
        if verbose:
            print('Too many pairs; only computing {} pairs'.format(maxpairs))
        if random_seed is not None:
            np.random.seed(random_seed)
        idx1, idx2 = np.random.choice(ndata, size=(maxpairs, 2)).transpose()
        pairmean = (data[idx1]+data[idx2])/2.
        hlmean_value = np.median(pairmean)

    return(hlmean_value)


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


# median absolute value (mav)
def mav(a):
    return np.median(np.abs(a))


def poly_fit_1d(x, y, order=1, rlm=False, t=1.5):
    '''
    1D polynomial fit.

    Given x array and y array, calculate the 1D robust polynomial fit of arbitrary order.

    See also poly_val_1d

    INPUT:
    1D arrays of x and y values; tunning parameter t;
    order of the polynomial fit;
    rlm: use robust linear model (HuberT).

    OUTPUT:
    Array of parameters of the polynomial [a0, a1, a2 ...]
    so that y = a0 + a1*x + a2*x**2 + ...
    '''

    import statsmodels.api as sm

    ncols = order+1
    a = np.zeros((x.size, ncols))
    for i in range(order+1):
        a[:, i] = x**i
    if rlm:
        res = sm.RLM(y, a, M=sm.robust.norms.HuberT(t=t)).fit()
    else:
        res = sm.OLS(y, a).fit()
    m = res.params
    return(m)


def poly_val_1d(x, m):
    '''
    Evaluate the 1D polynomial from x values and polynomial parameters

    See also poly_fit_1d

    INPUT:
    1D array of x values;
    1D array of polynomial parameters (for example generated by
    poly_fit_1d).

    OUTPUT:
    1D array of the evaluated values of the polynomial.
    '''

    order = len(m)-1
    z = np.zeros(x.shape)
    for i in range(order+1):
        z += m[i] * x**i
    return z


def poly_fit_2d(x, y, z, order=1, rlm=False, t=1.5):
    '''
    2D polynomial fit.

    Given x, y and z arrays, calculate the 2D polynomial fit of arbitrary order.

    See also poly_val_2d

    INPUT:
    1D arrays of x, y and z values; tunning parameter t;
    order of the polynomial fit;
    rlm: use robust linear model (HuberT);
    Note that x*y is a second order term.

    OUTPUT:
    Array of parameters of the polynomial [a0, a1, a2 ...]
    so that for an n-th order fit
    z = a0 + a1*y + a2*y**2 + ... + a_n*y**n + a_(n+1)*x
    + a_(n+2)*x*y + a_(n+3)*x*y**2 + ... + a_(2n)*x*y**(n-1)
    + a_(2n+1)*x**2 + a_(2n+2)*x**2*y + ... ...

    For example 2nd order fit gives:
    z = a0 + a1*y + a2*y**2 + a3*x + a4*x*y + a5*x**2

    New values can be evaluated using poly_val_2d
    '''

    import statsmodels.api as sm

    ncols = (order+2)*(order+1)//2
    a = np.zeros((x.size, ncols))
    k=0
    for i in range(order+1):
        for j in range(order-i+1):
            a[:, k] = x**i * y**j
            k+=1
    if rlm:
        res = sm.RLM(z, a, M=sm.robust.norms.HuberT(t=t)).fit()
    else:
        res = sm.OLS(z, a).fit()
    m = res.params
    return(m)


def poly_val_2d(x, y, m):
    '''
    Evaluate the 2D polynomial from x and y values and polynomial coefficients

    See also poly_fit_2d

    INPUT:
    1D array of x and y values;
    1D array of polynomial parameters (generated by poly_fit_2d).

    OUTPUT:
    1D array of the evaluated values of the polynomial.
    '''

    order = int((np.sqrt(8*len(m)+1)-3)/2)
    z = np.zeros(x.shape)
    k=0
    for i in range(order+1):
        for j in range(order-i+1):
            z += m[k] * x**i * y**j
            k+=1
    return z


def poly_fit_nd(X, y, rlm=False, order=1, t=1.5):
    '''
    N-dimensional polynomial fit.

    Given X(M, N), where M is number of data points and N is the number
    of dimensions, calculate the robust polynomial fit of arbitrary
    order. Huber weight function is used.

    INPUT:
    X: 2D arrays of the N-dimensional data points
    y: array to be fit
    order: order of the polynomial fit

    OUTPUT:
    m: polynomial coefficients;
    powers_arr: 2D arrays of the power of each variable
    '''

    import statsmodels.api as sm
    from itertools import product

    m, n = X.shape
    powers_all = np.tile(np.arange(order+1), (n, 1))

    a = []
    powers_arr = []

    for powers in product(*powers_all):
        if np.sum(powers)<=order:
            powers_arr.append(powers)
            column = np.ones(m)
            for index, p in enumerate(powers):
                column *= X[:, index]**p
            a.append(column)

    powers_arr = np.array(powers_arr)
    a = np.array(a).T

    if rlm:
        res = sm.RLM(y, a, M=sm.robust.norms.HuberT(t=t)).fit()
    else:
        res = sm.OLS(y, a).fit()
    coeffs = res.params
    return coeffs, powers_arr


def poly_val_nd(X, coeffs, powers_arr):
    '''
    Evaluate the N-dimensional polynomial given input X(M, N), with M being
    the number of data points and N being number of dimensions) and the
    polynomial coefficients from poly_fit_nd
    '''

    m, n = X.shape
    a = np.ones((m, len(powers_arr)))
    for ii, powers in enumerate(powers_arr):
        for jj, p in enumerate(powers):
            a[:, ii] *= X[:, jj]**p
    y = np.dot(a, coeffs)
    return y


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
    interp = interp1d(x, y, kind=kind)

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
            return np.array(list(map(pointwise, np.array(x))))
        elif isinstance(x, list):
            return np.array(list(map(pointwise, np.array(x))))
        else:
            return np.array(list(map(pointwise, [x])))

    return ufunclike


def create_image(data, cmap='gray', dpi=80, vmin=None, vmax=None, origin=None, norm=None):
    '''
    Create an image with exactly the same pixel dimension as the data.
    Example:
        x, y = np.arange(0, 10), np.arange(0, 10)
        xx, yy = np.meshgrid(x, y)
        img = np.array((xx + yy)%2==0, dtype=int)
        ax = create_image(img)
        plt.savefig('img.png')
        plt.close()
    '''
    xpixels, ypixels = data.shape[0], data.shape[1]
    figsize = ypixels / dpi, xpixels / dpi
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(data, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax, origin=origin, norm=norm)
    plt.axis('off')
    return ax

