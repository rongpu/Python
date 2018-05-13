from __future__ import division, print_function
import numpy as np
from astropy.cosmology import Planck15 as cosmo
from scipy.integrate import quad

def growth_factor_at_z(z, cosmo):
    '''
    Linear growth factor normalized to z=0. The mathematics follows from
    Equation 5 of arxiv.org/pdf/1309.5385.pdf

    Inputs:
    z: redshift (single value);
    cosmo: astropy comology object, e.g.:
        from astropy.cosmology import Planck15 as cosmo

    Output:
    Linear growth factor
    '''

    gamma = 0.55
    f = lambda z: 1/(1+z)*(cosmo.Om(z))**0.55

    growth_factor = np.exp(-quad(f, 0, z)[0])

    return growth_factor

