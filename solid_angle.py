from __future__ import division, print_function
import numpy as np

def solid_angle(phi1, theta1, phi2, theta2, return_steradian=True):
    '''
    Calculate the solid angle of a rectangle in the sky. The rectangle 
    is defined by its two corners (phi1, theta1) and (phi2, theta2), 
    where theta is the polar angle, in units of degrees. 

    Returns solid angle in steradian by default. 
    '''
    omega = np.abs(np.pi/180.*(phi2 - phi1)*(np.sin(theta2*np.pi/180) - np.sin(theta1*np.pi/180)))
    if return_steradian:
        return omega
    else:
        return omega*(180/np.pi)**2
