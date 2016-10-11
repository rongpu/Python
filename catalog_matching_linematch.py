# Matching between two catalogs.
# Use cat1 data to create a catalog that is line-matched to cat2.
# Has the option of removing cat1 objects outside cat2 region before matching to make matching faster. 

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
from os.path import expanduser
home = expanduser("~")+'/'
sys.path.append(home+'git/Python/user_modules/')
from catalog_matching_scatter_plot import scatter_plot


# search radius in arcsec
search_radius = 1.

fitsio_q = False # use fitsio to read cat2
save_q = False # save catalog
region_q = True # region selection - WARNING: MAY NOT WORK!!!!!!!
correct_offset_q = True
plot_q = True

# cat1_path = '/Users/roz18/Documents/Data/LSST_photo-z_testbed/Cross-identification/alldeep.egs.uniq.2012jun13.fits'
# cat2_path = '/Users/roz18/Documents/Data/LSST_photo-z_testbed/Cross-identification/Moffat v0.6/3D-HST_Terapix_Subaru_30_31_combined_v0.6b.fits'
cat1_path = '/Users/roz18/Documents/Data/LSST_photo-z_testbed/Cross-identification/Moffat v0.6/3D-HST_Terapix_Subaru_30_31_combined_v0.6b.fits'
cat2_path = '/Users/roz18/Documents/Data/LSST_photo-z_testbed/Cross-identification/alldeep.egs.uniq.2012jun13.fits'
output_path = '/Users/roz18/whatever.fits'

#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---

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


if region_q:

    # Remove cat1 objects outside the cat2 square
    ra2max = ra2.max()
    ra2min = ra2.min()
    dec2max = dec2.max()
    dec2min = dec2.min()
    if (ra2min-ra2max)%360<1.:  # If the cat2 region crosses with RA=0
        print('Warning: cat2 crosses with RA=0')
        ra2min = (((ra2+180)%360.).min()-180)%360.
        ra2max = (((ra2+180)%360.).max()-180+360)%360.
        mask = (ra1>ra2min-1/60.) & (ra1<ra2max+1/60.) & (dec1<dec2max+1/360.) & (dec1>dec2min-1/360.)
    else:  # otherwise
        mask = (ra1<ra2max+1/60.) & (ra1>ra2min-1/60.) & (dec1<dec2max+1/360.) & (dec1>dec2min-1/360.)
    cat1 = cat1[mask]
    ra1 = ra1[mask]
    dec1 = dec1[mask]
    print('%d objects removed from cat1'%(np.sum(~mask)))

    # Remove cat2 objects outside the cat1 square
    ra1max = ra1.max()
    ra1min = ra1.min()
    dec1max = dec1.max()
    dec1min = dec1.min()
    if (ra1min-ra1max)%360<1.:  # If the cat1 region crosses with RA=0
        print('Warning: cat1 crosses with RA=0')
        ra1min = (((ra1+180)%360.).min()-180)%360.
        ra1max = (((ra1+180)%360.).max()-180+360)%360.
        mask = (ra2<ra1max+1/60.) & (ra2>ra1min-1/60.) & (dec2<dec1max+1/360.) & (dec2>dec1min-1/360.)
    else:  # otherwise
        mask = (ra2<ra1max+1/60.) & (ra2>ra1min-1/60.) & (dec2<dec1max+1/360.) & (dec2>dec1min-1/360.)
    # bar keeps track of cat2 original index
    bar = np.arange(len(cat2))
    # effectively reduction of cat2
    ra2 = ra2[mask]
    dec2 = dec2[mask]
    bar = bar[mask]
    print('%d objects removed from cat2'%(np.sum(~mask)))

# Matching catalogs
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

# foo keeps track of cat2 (reduced) index for idx and d2d
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

#--------------------------- create line-matched cat1 catalog -----------------------------
# See Evernote for explanation

mask2 = idx>=0
if region_q:
    mask2full = np.zeros(len(cat2), dtype=bool)
    mask2full[bar[mask2]] = True
    notmask2full = ~mask2full
else:
    mask2full = mask2.copy()
    notmask2full = ~mask2full

# cat1_matchid is the cat2 (original) index for cat1 objects
cat1_matchid = -99.*np.ones(len(cat1))
if region_q:
    cat1_matchid[idx[mask2]] = bar[foo[mask2]]
else:
    cat1_matchid[idx[mask2]] = foo[mask2]
mask1 = cat1_matchid>=0
cat1_index = np.arange(len(cat1))
cat1_matchid = cat1_matchid[mask1]
cat1_index = cat1_index[mask1]

sort_index = cat1_matchid.argsort()
cat1_index = cat1_index[sort_index]
index_list = np.zeros(len(cat2), dtype=int)
index_list[mask2full] = cat1_index

cat1 = cat1[index_list]

# For objects with no match, assign '' to strings and 0 to numbers
for index in range(len(cat1.columns)):
    colformat = cat1.columns[index].format
    colname = cat1.columns[index].name        
    if colformat.find('A')>=0:
        cat1[colname][notmask2full] = ''
    else:
        cat1[colname][notmask2full] = 0

if save_q:
    fits.writeto(output_path, cat1, hdr1)

# -----------------------------------------------------------------------------------------

if plot_q:

    d_ra_hist=(cat2['ra'][mask2full]-cat1['ra'][mask2full])*3600    # in arcsec
    d_dec_hist=(cat2['dec'][mask2full]-cat1['dec'][mask2full])*3600 # in arcsec

    ax = scatter_plot(d_ra_hist, d_dec_hist)
    ax.plot(ra_offset*3600, dec_offset*3600, 'r.', markersize=9)
    
    plt.show()
