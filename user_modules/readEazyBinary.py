from __future__ import print_function, division
import os
import numpy as np

def readtempfilt(MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT'):
    """
    Read redshift grid, filter, template and flux information from the .templfilt 
    binary file.
        
    """

    root = os.path.join(OUTPUT_DIRECTORY, MAIN_OUTPUT_FILE)
    
    ###### .tempfilt
    file_path = root+'.tempfilt'
    
    if os.path.exists(file_path) is False:
        raise ValueError('File, %s, not found.' %(file_path))

    with open(file_path,'rb') as f:
        # summary data
        s = np.fromfile(file=f,dtype=np.int32, count=4)
        NFILT=s[0]  # number of filters
        NTEMP=s[1]  # number of templates
        NZ=s[2]     # number points on the redshift grid
        NOBJ=s[3]   # number of objects
        # (?) template SED convolved with filter transmission at each redshift
        tempfilt = np.fromfile(file=f,dtype=np.double,count=NFILT*NTEMP*NZ).reshape((NZ,NTEMP,NFILT)).transpose()
        # filter pivot wavelengths
        lc = np.fromfile(file=f,dtype=np.double,count=NFILT)
        # redshift grid
        zgrid = np.fromfile(file=f,dtype=np.double,count=NZ)
        # observed flux
        fnu = np.fromfile(file=f,dtype=np.double,count=NFILT*NOBJ).reshape((NOBJ,NFILT)).transpose()
        # (?) error in observed flux
        efnu = np.fromfile(file=f,dtype=np.double,count=NFILT*NOBJ).reshape((NOBJ,NFILT)).transpose()
        
    tempfilt  = {'NFILT':NFILT,'NTEMP':NTEMP,'NZ':NZ,'NOBJ':NOBJ,\
                 'tempfilt':tempfilt,'lc':lc,'zgrid':zgrid,'fnu':fnu,'efnu':efnu}

    return tempfilt

def readpz(MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', APPLY_PRIOR=True):
    """
    Get Eazy p(z), prior and chi2 from binary data. 

    Inputs
    ------
    APPLY_PRIOR: if True, the prior will be included in the output and it will be applied 
    on p(z);

    Outputs
    -------
    pz_dict: dictionary of
        zgrid: redshift grid
        pz: 2d array (idx, zgrid) of the probability distribution
        prior: 2d array (idx, zgrid) of the prior
        chi2: 2d array (idx, zgrid) of chi2
    
    """

    ###### Read tempfilt data to get zgrid
    tempfilt = readtempfilt(MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, OUTPUT_DIRECTORY=OUTPUT_DIRECTORY)
    zgrid = tempfilt['zgrid']

    ###### Read .pz binary file
    root = os.path.join(OUTPUT_DIRECTORY, MAIN_OUTPUT_FILE)
    file_path = root+'.pz'

    if os.path.exists(file_path) is False:
        raise ValueError('File, %s, not found.' %(file_path))

    with open(file_path,'rb') as f:
        s = np.fromfile(file=f,dtype=np.int32, count=2)
        NZ=s[0]     # number points on the redshift grid
        NOBJ=s[1]   # number of objects
        # (object, redshift) Chi2 of the fit at each redshift
        chi2 = np.fromfile(file=f,dtype=np.double,count=NZ*NOBJ).reshape((NOBJ,NZ))

        if APPLY_PRIOR:
            # prior information
            s = np.fromfile(file=f,dtype=np.int32, count=1)
            # load prior data
            if len(s) > 0:
                NK = s[0] # number of magnitude bins in the prior(redshift, mag) matrix
                # list of magnitude bins
                kbins = np.fromfile(file=f,dtype=np.double,count=NK)
                # prior (mag, redshift) array
                priorkz = np.fromfile(file=f, dtype=np.double, count=NZ*NK).reshape((NK,NZ))
                # list indices of magnitude bins that each object belongs to
                kidx = np.fromfile(file=f,dtype=np.int32,count=NOBJ)
            else:
                raise ValueError('Prior data not found.')

    ###### Convert Chi2 to p(z)
    # [:, None] extends 1D array to 2D
    chi2min = np.ones_like(chi2)*np.min(chi2, axis=1)[:, None]
    if APPLY_PRIOR:
        # prior array for each object
        # kidx==0 for brightest objects and objects with failed photo-z measurements
        prior = priorkz[kidx]
        pz = np.exp(-0.5*(chi2-chi2min))*prior
    else:
        pz = np.exp(-0.5*(chi2-chi2min))

    norm = np.trapz(pz, zgrid)[:, None]
    if np.all(norm > 0):
        pz/=norm
    else:
        raise ValueError('Integrated p(z) must be larger than zero')

    ###### Done
    if APPLY_PRIOR:
        pz_dict = {'zgrid':zgrid, 'pz':pz, 'prior':prior, 'chi2':chi2}
    else:
        pz_dict = {'zgrid':zgrid, 'pz':pz, 'chi2':chi2}
    return pz_dict