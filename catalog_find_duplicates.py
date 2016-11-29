# Self-matching and find duplicates within some search radius

from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.io import fits
from astropy.coordinates import SkyCoord
from matplotlib.ticker import NullFormatter
import sys
from os.path import expanduser
home = expanduser("~")+'/'
sys.path.append(home+'git/Python/user_modules/')
from catalog_matching_scatter_plot import scatter_plot

# search radius in arcsec
search_radius = 1.

path = 'somewhere'

#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---

cat = fits.getdata(path, ext=1)
ra = np.array(cat['ra'])
dec = np.array(cat['dec'])
skycat = SkyCoord(ra*u.degree,dec*u.degree, frame='icrs')
idx, d2d, _ = skycat.match_to_catalog_sky(skycat)

mask = d2d<(search_radius*u.arcsec)

print('%d duplicates found'%np.sum(mask))
