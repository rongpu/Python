from __future__ import print_function, division
import numpy as np
from astropy import units as u
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from matplotlib.ticker import NullFormatter


def match_coord(ra1, dec1, ra2, dec2, search_radius=1., nthneighbor=1, plot_q=True, verbose=True):
    '''
    Inputs: RA and Dec of two catalogs;

    Outputs: 
        indices of matched objects in the two catalogs;
        distances;
        the differences (in arcsec) in RA and Dec; 
    '''
    
    t1 = Table()
    t2 = Table()
    
    # protect the global variables from being changed by np.sort
    ra1, dec1, ra2, dec2 = map(np.copy, [ra1, dec1, ra2, dec2])

    t1['ra'] = ra1
    t2['ra'] = ra2
    t1['dec'] = dec1
    t2['dec'] = dec2
    
    t1['foo'] = np.arange(len(t1))
    t2['foo'] = np.arange(len(t2))
    
    # Matching catalogs
    sky1 = SkyCoord(ra1*u.degree,dec1*u.degree, frame='icrs')
    sky2 = SkyCoord(ra2*u.degree,dec2*u.degree, frame='icrs')
    idx, d2d, d3d = sky2.match_to_catalog_sky(sky1, nthneighbor=nthneighbor)
    # This find a match for each object in t2. Not all objects in t1 catalog is included in the result. 
    
    matchlist = d2d<(search_radius*u.arcsec)
    t2['idx'] = idx
    t2['d2d'] = d2d
    t2 = t2[matchlist]
    
    init_count = np.sum(matchlist)

    #------------------------------removing doubly matched objects------------------------------

    t2.sort('idx')
    i = 0
    while i<=len(t2)-2:
        if t2['idx'][i]>=0 and t2['idx'][i]==t2['idx'][i+1]:
            end = i+1
            while end+1<=len(t2)-1 and t2['idx'][i]==t2['idx'][end+1]:
                end = end+1
            findmin = np.argmin(t2['d2d'][i:end+1])
            for j in range(i,end+1):
                if j!=i+findmin:
                    t2['idx'][j]=-99
            i = end+1
        else:
            i = i+1

    mask_match = t2['idx']>=0
    t2 = t2[mask_match]

    t2.sort('foo')

    if verbose:
        print('Doubly matched objects = %d'%(init_count-len(t2)))
        print('Final matched objects = %d'%len(t2))

    # -----------------------------------------------------------------------------------------

    # This rearranges t1 to match t2 by index.
    t1['matchid'] = -99
    for i in np.arange(len(t2)):
        t1['matchid'][t2['idx'][i]] = t2['foo'][i]
    t1.sort('matchid')
    mask_in_t2 = t1['matchid']>=0
    t1 = t1[mask_in_t2]

    d_ra = (t2['ra']-t1['ra'])*3600    # in arcsec
    d_dec = (t2['dec']-t1['dec'])*3600 # in arcsec

    if plot_q:
        scatter_plot(d_ra, d_dec)

    return np.array(t1['foo']), np.array(t2['foo']), np.array(t2['d2d']), np.array(d_ra), np.array(d_dec)


def find_neighbor(ra1, dec1, search_radius=1., nthneighbor=1):
    '''
    Find the n-th nearest neighbor. 
    nthneighbor: the n-th neighbor; the nthneighbor=1 is the first neighbor other than itself. 
    '''
    
    t1 = Table()
    t1['ra'] = ra1
    t1['dec'] = dec1
    t1['foo'] = np.arange(len(t1))

    # Matching catalogs
    sky1 = SkyCoord(ra1*u.degree,dec1*u.degree, frame='icrs')
    idx, d2d, d3d = sky1.match_to_catalog_sky(sky1, nthneighbor=(nthneighbor+1))
    # This find a match for each object in t2. Not all objects in t1 catalog is included in the result. 
    
    matchlist = d2d<(search_radius*u.arcsec)
    t1['idx'] = idx
    t1['d2d'] = d2d
    t1 = t1[matchlist]
    
    return np.array(t1['foo']), np.array(t1['idx']), np.array(t1['d2d'])


def scatter_plot(d_ra, d_dec, title='', x_label='$\\mathbf{RA_{cat2} - RA_{cat1}(arcsec)}$', y_label='$\\mathbf{dec_{cat2} - dec_{cat1}(arcsec)}$'):
    '''
    INPUTS:

     d_ra, d_dec: array of RA and Dec difference in arcsec
    
    OUTPUTS:

     axScatter: scatter-histogram plot
    '''

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
    # mask = np.logical_and(np.abs(d_ra)<1.66, np.abs(d_dec)<1.)
    # axScatter.plot(d_ra[mask], d_dec[mask], 'k.', markersize=1)
    axScatter.plot(d_ra, d_dec, 'k.', markersize=1)

    axHistx.hist(d_ra, bins=100, histtype='step', color='k', linewidth=2)
    axHisty.hist(d_dec, bins=100, histtype='step', color='k', linewidth=2, orientation='horizontal')

    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())

    axHistx.axis('off')
    axHisty.axis('off')

    axScatter.axhline(0, color='r', linestyle='--', linewidth=1.2)
    axScatter.axvline(0, color='r', linestyle='--', linewidth=1.2)
    axScatter.set_xlabel(('$RA_{cat2} - RA_{cat1}(arcsec)$'))
    axScatter.set_ylabel(('$dec_{cat2} - dec_{cat1}(arcsec)}$'))

    return(axScatter)
    # #--------------- ra dec histogram ---------------------

    # plt.figure()
    # plt.hist(d_ra,bins=50)
    # plt.title('RA difference between cat2/3 and Terapix catalog')
    # plt.xlabel('RA_cat2 - RA_cat1 (arcsec)')
    # plt.grid()

    # plt.figure()
    # plt.hist(d_dec,50)
    # plt.title('Dec difference between cat2/3 and Terapix catalog')
    # plt.xlabel('Dec_cat2 - Dec_cat1 (arcsec)')
    # plt.grid()

    plt.show()


def match_self(ra, dec, search_radius=1., return_indices=False, plot_q=False, verbos=True):
    '''
    Find objects that has a neighbor within search_radius arcsec. 

    Return: 
    Number of suspected duplicates. 
    (Optional) idx1, idx2: arrays of indices of suspected duplicates. 
    '''

    ra = np.array(ra)
    dec = np.array(dec)
    skycat = SkyCoord(ra*u.degree,dec*u.degree, frame='icrs')
    idx, d2d, _ = skycat.match_to_catalog_sky(skycat, nthneighbor=2)

    mask = d2d<(search_radius*u.arcsec)
    print(np.sum(mask), "objects with a nearby neighbor")
    n_duplicates = np.sum(mask)
    idx1 = np.arange(len(ra))[mask]
    idx2 = idx[mask]

    if plot_q and (n_duplicates!=0):
        d_ra = ra[idx1] - ra[idx2]
        d_dec = dec[idx1] - dec[idx2]
        scatter_plot(d_ra*3600., d_dec*3600.)

    if return_indices:
        idx_dup = np.arange(len(ra))[mask]
        return np.sum(mask), idx1, idx2
    else:
        return np.sum(mask)