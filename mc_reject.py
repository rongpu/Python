from __future__ import division, print_function
import numpy as np

def random(f, fmax, xmin, xmax, n, random_seed=None):
    '''
    Generate random numbers from a distribution using the Monte Carlo 
    rejection method. We only use uniform distribution. 

    Inputs
    ------
    f: function of the probability distribution; it should be always 
       positive;
    fmax: maximum value of f, or upper bound of f if its maximum value 
          is unkown;
    xmin, xmax: range of the random variable;
    n: number of random numbers;
    random_seed: random seed, automatic seed by default. 

    Output
    ------
    Numpy array of the random numbers
    '''

    if fmax<0 or xmin>xmax or n<=0:
        raise ValueError('Bad inputs')
        
    count = 0
    randarray = np.zeros(n)
    frac = -1 # fraction of accepted values 

    if random_seed!=None:
        np.random.seed(random_seed)
    nprand1 = np.random.uniform(xmin, xmax, n)
    nprand2 = np.random.uniform(0, 1, n)

    while count<n: 
        mask = f(nprand1)>nprand2*fmax
        count_add = min(n-count, np.sum(mask))
        randarray[count:(count+count_add)] = nprand1[mask][:count_add]
        count = count + count_add
        if frac==-1:
            frac = count/n
        if count!=n:
            nprand1 = np.random.uniform(xmin, xmax, int(np.ceil((n-count)/frac)))
            nprand2 = np.random.uniform(0, 1, int(np.ceil((n-count)/frac)))
            
    return(randarray)