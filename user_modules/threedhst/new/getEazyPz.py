def getEazyPz(idx, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same', binaries=None, get_prior=False, get_chi2=False):
    """
    zgrid, pz = getEazyPz(idx, \
                      MAIN_OUTPUT_FILE='photz', \
                      OUTPUT_DIRECTORY='./OUTPUT', \
                      CACHE_FILE='Same', binaries=None)
                      
    Get Eazy p(z) for object #idx.
    
    To avoid re-reading the binary files, supply binaries = (tempfilt, pz)
    
    """
    if binaries is None:
        tempfilt, coeffs, temp_seds, pz = readEazyBinary(MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                                                    OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                                                    CACHE_FILE = CACHE_FILE)
    else:
        tempfilt, pz = binaries
        
    if pz is None:
        return None, None
    
    ###### Get p(z|m) from prior grid
    kidx = pz['kidx'][idx]
    #print kidx, pz['priorzk'].shape
    if (kidx > 0) & (kidx < pz['priorzk'].shape[1]):
        prior = pz['priorzk'][:,kidx]
    else:
        prior = np.ones(pz['NZ'])
    
    if get_chi2:
        if get_prior:
            if get_prior:
                return tempfilt['zgrid'], pz['chi2fit'][:,idx], prior
            else:
                return tempfilt['zgrid'], pz['chi2fit'][:,idx]
            
    ###### Convert Chi2 to p(z)
    pzi = np.exp(-0.5*(pz['chi2fit'][:,idx]-min(pz['chi2fit'][:,idx])))*prior#*(1+tempfilt['zgrid'])
    
    if np.sum(pzi) > 0:
        pzi/=np.trapz(pzi, tempfilt['zgrid'])
    
    ###### Done
    if get_prior:
        return tempfilt['zgrid'], pzi, prior
    else:
        return tempfilt['zgrid'], pzi