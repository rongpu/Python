# Matching between two catalogs
# Use cat1 data to create a catalog that is line-matched to cat2

# Creating line-matched catalogs (all objects in cat2 are kept):
# 1. Find matched objects
# 2. Flag duplicates (objects in cat1 that have more than one match to cat2)
# 3. Rearrange cat1 so that it is line-matched to cat2 and there's no duplicate
# 4. Save the new cat1
# Note:
# - 1 arcsec search radius
# - Objects with no match are assigned '' for string and 0 for numbers

from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from matplotlib.ticker import NullFormatter
import sys
import fitsio


# search radius in arcsec
search_radius = 1.

fitsio_q = False # use fitsio to read cat2
save_q = False # save catalog
region_q = False # region selection - WARNING: MAY NOT WORK!!!!!!!
correct_offset_q = True
plot_q = True
verbose = True

cat1_path = '/Users/roz18/Documents/Data/LSST_photo-z_testbed/Cross-identification/alldeep.egs.uniq.2012jun13.fits'
cat2_path = '/Users/roz18/Documents/Data/LSST_photo-z_testbed/Cross-identification/Moffat v0.6/3D-HST_Terapix_Subaru_30_31_combined_v0.6b.fits'
output_path = 'whatever.fits'

#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---

if verbose:
    print("Loading data...")
cat1, hdr1 = fits.getdata(cat1_path, ext=1, header=True)
if fitsio_q:
    # Not case sensitive here
    cat2 = fitsio.read(cat2_path, columns=['ra', 'dec'])
else:
    cat2, hdr2 = fits.getdata(cat2_path, ext=1, header=True)

ra1 = np.array(cat1['ra'])
dec1 = np.array(cat1['dec'])
# Case sensitive if fitsio is used
ra2 = np.array(cat2['ra'])
dec2 = np.array(cat2['dec'])

# Remove cat1 objects outside the cat2 square
if region_q:
    ra2max = ra2.max()
    ra2min = ra2.min()
    dec2max = dec2.max()
    dec2min = dec2.min()
    if (ra2min-ra2max)%360<1.:  # If the cat2 region crosses with RA=0
        print('Warning: cat2 crosses with RA=0')
        ra2min = (((ra2+180)%360.).min()-180)%360.
        ra2max = (((ra2+180)%360.).max()-180+360)%360.
        mask = (ra1>ra2min/60.) & (ra1<ra2max+1/60.) & (dec1<dec2max+1/360.) & (dec1>dec2min-1/360.)
    else:  # normal case
        mask = (ra1<ra2max+1/60.) & (ra1>ra2min-1/60.) & (dec1<dec2max+1/360.) & (dec1>dec2min-1/360.)
    cat1 = cat1[mask]
    ra1 = ra1[mask]
    dec1 = dec1[mask]

# Matching catalogs
if verbose:
    print("Matching...")
skycat1 = SkyCoord(ra1*u.degree,dec1*u.degree, frame='icrs')
skycat2 = SkyCoord(ra2*u.degree,dec2*u.degree, frame='icrs')
idx, d2d, _ = skycat2.match_to_catalog_sky(skycat1)
# For each object in cat2, a closest match to cat1 is found. Thus not all cat1 objects are included. 
# idx is the cat1 index for cat2 -- idx[0] is the cat1 index that the first cat2 object matched.
# Similarly d2d is the distance between the matches. 

mask2 = d2d<(search_radius*u.arcsec)
if np.sum(mask2)==0:
    print('0 match!')
    sys.exit()
count1 = np.sum(mask2)
print('Number of initial matches = %d'%count1)

# Correct for systematic offsets
if correct_offset_q:
    ra_offset = np.median(ra2[mask2] - ra1[idx[mask2]])
    dec_offset = np.median(dec2[mask2] - dec1[idx[mask2]])
    ra1 = ra1 + ra_offset
    dec1 = dec1 + dec_offset
    if verbose:
        print('RA  offset = %f arcsec'%(ra_offset*3600))
        print('Dec offset = %f arcsec'%(dec_offset*3600))
    skycat1 = SkyCoord(ra1*u.degree,dec1*u.degree, frame='icrs')
    idx, d2d, _ = skycat2.match_to_catalog_sky(skycat1)
    mask2 = d2d<(search_radius*u.arcsec)
    count1 = np.sum(mask2)
    print('Number of matches after offset correction = %d'%count1)

notmask2 = ~mask2
idx[notmask2] = -99


#------------------------------removing doubly matched points------------------------------

# foo keeps track of cat2 index for idx and d2d
foo = np.arange(len(cat2))
# Sort by idx to find duplicates
sort_index = idx.argsort()
idx.sort()
d2d = d2d[sort_index]
foo = foo[sort_index]

# Find the first non-zero idx
i = np.argmax(idx>=0)
# Find duplicates, keep only the nearest match
while i<=len(idx)-2:
    if idx[i]>=0 and idx[i]==idx[i+1]:
        end = i+1
        while end+1<=len(idx)-1 and idx[i]==idx[end+1]:
            end = end+1
        findmin = np.argmin(d2d[i:end+1])
        for j in range(i,end+1):
            if j!=i+findmin:
                idx[j]=-99
        i = end+1
    else:
        i = i+1
count2 = np.sum(idx>=0)

print('Number of duplicates = %d'%(count1 - count2))
print('Number of final matches = %d'%count2)

# Restore idx and d2d to original order
sort_index = foo.argsort()
foo.sort()
idx = idx[sort_index]
d2d = d2d[sort_index]

# -----------------------------------------------------------------------------------------

# This part remove objects from cat1 that have no successful match, 
# but keep extra space so that there are enough rows for line-matchign
mask2 = idx>=0
notmask2 = ~mask2
# cat1_matchid is the cat2 index for cat1 objects
cat1_matchid = -99.*np.ones(len(cat1))
cat1_matchid[idx[mask2]] = foo[mask2]
# add extra rows for line-matching
extra_ct = len(cat2)-count2
index_list = np.concatenate((np.where(cat1_matchid>=0)[0], np.zeros(extra_ct, dtype=int)))
cat1 = cat1[index_list]
cat1_matchid = cat1_matchid[index_list]
cat1_matchid[-extra_ct:] = foo[notmask2]

#--------------------------- create line-matched cat1 catalog -----------------------------

sort_index = cat1_matchid.argsort()
for index in range(len(cat1.columns)):
    colformat = cat1.columns[index].format
    colname = cat1.columns[index].name        
    cat1[colname] = cat1[colname][sort_index]
    # Assign '' to strings and 0 to numbers
    if colformat.find('A')>=0:
        cat1[colname][notmask2] = ''
    else:
        cat1[colname][notmask2] = 0

if save_q:
    fits.writeto(output_path, cat1, hdr1)

# -----------------------------------------------------------------------------------------

if plot_q:

    mask = idx>=0
    d_ra_hist=(cat2['ra'][mask]-cat1['ra'][mask])*3600    # in arcsec
    d_dec_hist=(cat2['dec'][mask]-cat1['dec'][mask])*3600 # in arcsec

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
    axScatter.set_xlabel(('$\\mathbf{RA_{cat2} - RA_{cat1}(arcsec)}$'))
    axScatter.set_ylabel(('$\\mathbf{dec_{cat2} - dec_{cat1}(arcsec)}$'))

    axScatter.plot(ra_offset*3600, dec_offset*3600, 'r.', markersize=9)

    #--------------- ra dec histogram ---------------------

    plt.figure()
    plt.hist(d_ra_hist,bins=50)
    plt.title('RA difference between cat2/3 and Terapix catalog')
    plt.xlabel('RA_cat2 - RA_cat1 (arcsec)')
    plt.grid()

    plt.figure()
    plt.hist(d_dec_hist,50)
    plt.title('Dec difference between cat2/3 and Terapix catalog')
    plt.xlabel('Dec_cat2 - Dec_cat1 (arcsec)')
    plt.grid()

    plt.show()
