# Self-matching and find duplicates within some search radius

from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

def match_self(ra, dec, search_radius=1., return_index=False):
    '''
    Find objects that has a neighbor within search_radius arcsec. 

    Return: 
    Number of such objects. 
    (Optional) Indices of such objects. 
    '''

    ra = np.array(ra)
    dec = np.array(dec)
    skycat = SkyCoord(ra*u.degree,dec*u.degree, frame='icrs')
    idx, d2d, _ = skycat.match_to_catalog_sky(skycat, nthneighbor=2)

    mask = d2d<(search_radius*u.arcsec)

    if return_index:
        idx_dup = np.arange(len(ra))[mask]
        return np.sum(mask), idx_dup
    else:
        return np.sum(mask)