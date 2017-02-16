from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

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

    axScatter.axhline(0, color='k', linestyle='--', linewidth=1.2)
    axScatter.axvline(0, color='k', linestyle='--', linewidth=1.2)
    axScatter.set_xlabel(('$\\mathbf{RA_{cat2} - RA_{cat1}(arcsec)}$'))
    axScatter.set_ylabel(('$\\mathbf{dec_{cat2} - dec_{cat1}(arcsec)}$'))

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

    # plt.show()
