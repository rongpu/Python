from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np

def gc_unit_vector(ra1, dec1, ra2, dec2):
    '''
    Compute the unit vector(s) of a great circle from two points on the circle.
    Inputs are in degrees.
    '''

    ra1 = np.array(ra1)
    dec1 = np.array(dec1)
    ra2 = np.array(ra2)
    dec2 = np.array(dec2)

    # Convert (RA, Dec) to Cartesian coordinates
    x1 = np.cos(ra1/180*np.pi)*np.cos(dec1/180*np.pi)
    y1 = np.sin(ra1/180*np.pi)*np.cos(dec1/180*np.pi)
    z1 = np.sin(dec1/180*np.pi)
    x2 = np.cos(ra2/180*np.pi)*np.cos(dec2/180*np.pi)
    y2 = np.sin(ra2/180*np.pi)*np.cos(dec2/180*np.pi)
    z2 = np.sin(dec2/180*np.pi)

    cc1 = np.array([x1, y1, z1]).transpose()
    cc2 = np.array([x2, y2, z2]).transpose()

    # Find the unit vectors defining the great circle
    v = np.cross(cc1, cc2)
    if len(v.shape)==1:
        v /= (np.sqrt(np.sum(v**2)))
    else:
        v /= (np.sqrt(np.sum(v**2, axis=1))[:, None])

    return v

def distance_to_gc(ra, dec, v):
    '''
    Compute the angular distance between a point (or points) and a great circle.
    
    Inputs
    ------
    ra, dec: spherical coordinates (degrees)
    v: unit vector of the great circle in cartesian coordinates

    Returns
    -------
    d: angular distance (degrees)
    '''

    # Convert (RA, Dec) to unit vector in cartesian coordinates
    x = np.cos(ra/180*np.pi)*np.cos(dec/180*np.pi)
    y = np.sin(ra/180*np.pi)*np.cos(dec/180*np.pi)
    z = np.sin(dec/180*np.pi)
    cc = np.array([x, y, z]).transpose()

    if len(cc.shape)==1:
        d = np.abs(np.arccos(np.sum(cc*v))/np.pi*180. - 90.)
    else:
        d = np.abs(np.arccos(np.sum(cc*v, axis=1))/np.pi*180. - 90.)

    return d

def distance_to_line(ra0, dec0, ra1, dec1, ra2, dec2):
    '''
    Compute the angular distance between a point (or points) and a line (great circle)
    that is defined by two points.
    
    Inputs
    ------
    ra0, dec0: (float or array) the point(s) of interest
    ra1, dec1, ra2, dec2: (float) the two points defining the line

    Returns
    -------
    d: angular distance (degrees)
    '''

    ra0 = np.array(ra0)
    dec0 = np.array(dec0)
    v = gc_unit_vector(ra1, dec1, ra2, dec2)
    return distance_to_gc(ra0, dec0, v)
