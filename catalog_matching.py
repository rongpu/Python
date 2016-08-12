# Matching between Terapix+DEEP2/3 and Subaru
# Duplicates are not kept
# Plot scatter-histogram plot

from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from matplotlib.ticker import NullFormatter

# Whether to keep all objects in cat2
# If False, cat1 and cat2 will has the same exact objects
keep_cat2_q = False

cat1=Table.read('data/subaru_28_29_final.fits')
cat2=Table.read('data/DEEP2_Terapix_combined.fits')

# give each object in cat2 a unique number
cat2['foo'] = np.arange(len(cat2))

search_radius = 1.

ra1=np.array(cat1['ra'])
dec1=np.array(cat1['dec'])
ra2=np.array(cat2['RA'])
dec2=np.array(cat2['DEC'])

# Matching catalogs
skycat1=SkyCoord(ra1*u.degree,dec1*u.degree, frame='icrs')
skycat2=SkyCoord(ra2*u.degree,dec2*u.degree, frame='icrs')
idx, d2d, d3d = skycat2.match_to_catalog_sky(skycat1)
# This find a match for each object in cat2 catalog. Not all objects in cat1 catalog is included in the result. 

matchlist=d2d<(search_radius*u.arcsec)
cat2['idx']=idx
cat2['d2d']=d2d

if keep_cat2_q:
    for i in range(len(cat2)):
        mask = ~matchlist
        cat2['idx'][mask] = -99
else:
    cat2=cat2[matchlist]

print('Number of initial matches = %d' %np.sum(matchlist))

#------------------------------removing doubly matched points------------------------------

if keep_cat2_q==False:

    cat2.sort('idx')
    i=0
    while i<=len(cat2)-2:
        if cat2['idx'][i]>=0 and cat2['idx'][i]==cat2['idx'][i+1]:
            end=i+1
            while end+1<=len(cat2)-1 and cat2['idx'][i]==cat2['idx'][end+1]:
                end=end+1
            findmin=np.argmin(cat2['d2d'][i:end+1])
            for j in range(i,end+1):
                if j!=i+findmin:
                    cat2['idx'][j]=-99
            i=end+1
        else:
            i=i+1

    mask = cat2['idx']>=0
    cat2 = cat2[mask]

    cat2.sort('foo')

    print('Final matched objects = %d' %len(cat2))
# -----------------------------------------------------------------------------------------

# This part generates a cat1 catalog that matches the cat2 by index.
if keep_cat2_q==False:
    cat1['matchid']=-99
    for i in np.arange(len(cat2)):
        cat1['matchid'][cat2['idx'][i]]=cat2['foo'][i]
    cat1.sort('matchid')
    mask = cat1['matchid']>=0
    cat1=cat1[mask]

# plotting
if keep_cat2_q:
    mask = cat2['idx']>=0
    d_ra_hist=(cat2['ra'][mask]-cat1[cat2['idx'][mask]]['ra'])*3600    # in arcsec
    d_dec_hist=(cat2['dec'][mask]-cat1[cat2['idx'][mask]]['dec'])*3600 # in arcsec
else:
    d_ra_hist=(cat2['RA']-cat1['ra'])*3600    # in arcsec
    d_dec_hist=(cat2['DEC']-cat1['dec'])*3600 # in arcsec


# ------------ Plot scatter-histogram plot ---------------

nullfmt = NullFormatter()         # no labels

# definitions for the axes
left, width = 0.1, 0.85
bottom, height = 0.1, 0.85

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom, width, 0.3]
rect_histy = [left, bottom, 0.3, height]

# start with a rectangular Figure
plt.figure(figsize=(8,8))

axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)

# # no labels
# axHistx.xaxis.set_major_formatter(nullfmt)
# axHisty.yaxis.set_major_formatter(nullfmt)

# the scatter plot
mask = np.logical_and(np.abs(d_ra_hist)<1.66, np.abs(d_dec_hist)<1.)
axScatter.plot(d_ra_hist[mask], d_dec_hist[mask], 'k.', markersize=1)

axHistx.hist(d_ra_hist, bins=100, histtype='step', color='k', linewidth=2)
axHisty.hist(d_dec_hist, bins=100, histtype='step', color='k', linewidth=2, orientation='horizontal')

axHistx.set_xlim(axScatter.get_xlim())
axHisty.set_ylim(axScatter.get_ylim())

axHistx.axis('off')
axHisty.axis('off')

axScatter.axhline(0, color='k', linestyle='--', linewidth=1.2)
axScatter.axvline(0, color='k', linestyle='--', linewidth=1.2)
axScatter.set_xlabel(('$\\mathbf{RA_{DEEP2} - RA_{Subaru}(arcsec)}$'))
axScatter.set_ylabel(('$\\mathbf{dec_{DEEP2} - dec_{Subaru}(arcsec)}$'))


#--------------- ra dec histogram ---------------------

plt.figure()
plt.hist(d_ra_hist,bins=50)
plt.title('RA difference between DEEP2/3 and Terapix catalog')
plt.xlabel('RA_DEEP2 - RA_Subaru (arcsec)')
plt.grid()

plt.figure()
plt.hist(d_dec_hist,50)
plt.title('Dec difference between DEEP2/3 and Terapix catalog')
plt.xlabel('Dec_DEEP2 - Dec_Subaru (arcsec)')
plt.grid()


plt.show()

cat2.remove_columns(['foo', 'idx', 'd2d'])
