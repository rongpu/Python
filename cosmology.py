from __future__ import division, print_function
import numpy as np
from astropy.cosmology import Planck15 as cosmo
from scipy.integrate import quad

def growth_factor_at_z(z):
    '''
    Linear growth factor normalized to z=0. The mathematics follows from
    Equation 5 of arxiv.org/pdf/1309.5385.pdf
    '''

    gamma = 0.55
    f = lambda z: 1/(1+z)*(cosmo.Odm(z))**0.55

    growth_factor = np.exp(-quad(f, 0, z)[0])

    return growth_factor

