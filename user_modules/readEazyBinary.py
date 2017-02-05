import os
import numpy as np

def readtempfilt(MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT'):
    """
    Read Eazy p(z) from BINARY_OUTPUTS files.

    Outputs
    -------
    tempfilt
        
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
    Inputs
    ------
    get_prior: if True, the prior will be applied;
    get_chi2: if True, returns chi2;

    Outputs
    -------
    pz_dict: dictionary of zgrid, p(z), prior(optional) and chi2
        zgrid
        p(z): matrix (idx, zgrid) of the probability distribution with or without a prior

    Get Eazy p(z) from binary data. 
    
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
        chi2fit = np.fromfile(file=f,dtype=np.double,count=NZ*NOBJ).reshape((NOBJ,NZ))

        if APPLY_PRIOR:
            ##################### This will break if APPLY_PRIOR No #####################
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
    if APPLY_PRIOR:
        # [:, None] extends 1D array to 2D
        chi2min = np.ones_like(chi2fit)*np.min(chi2fit, axis=1)[:, None]
        # prior for each object
        prior = priorkz[kidx]
        pz = np.exp(-0.5*(chi2fit-chi2min))*prior
    else:
        pz = np.exp(-0.5*(chi2fit-chi2min))

    norm = np.trapz(pz, zgrid)[:, None]
    if np.all(norm > 0):
        pz/=norm
    else:
        raise ValueError('Integrated p(z) must be larger than zero')

    ###### Done
    if APPLY_PRIOR:
        pz_dict = {'zgrid':zgrid, 'pz':pz, 'prior':prior, 'chi2':chi2fit}
    else:
        pz_dict = {'zgrid':zgrid, 'pz':pz, 'chi2':chi2fit}
    return pz_dict