from __future__ import division, print_function
import numpy as np

def solid_angle(ra1, dec1, ra2, dec2, return_steradian=False):
    '''
    Calculate the solid angle of a rectangle in the sky. The rectangle 
    is defined by its two corners (ra1, dec1) and (ra2, dec2), 
    where theta is the polar angle, in units of degrees. 

    If return_steradian=False, returns solid angle in square degrees.

    '''
    omega = np.abs(np.pi/180.*(ra2 - ra1)*(np.sin(dec2*np.pi/180) - np.sin(dec1*np.pi/180)))
    if return_steradian:
        return omega
    else:
        return omega*(180/np.pi)**2
