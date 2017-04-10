def getEazySED(idx, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same', scale_flambda=True, verbose=False, individual_templates=False):
    """
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
     getEazySED(idx, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same')
    
    Get best-fit Eazy template for object number 'idx' from the specified Eazy output files. 

    Output variables are as follows:
        
        lambdaz: full best-fit template (observed) wavelength, interpolated at WAVELENGTH_GRID
        temp_sed:          "        "              flux (F_lambda)
        lci: filter pivot wavelengths
        fobs: observed fluxes, including zeropoint offsets if used, F_lambda
        efobs: observed flux errors,    "            "        "        "
    """
    tempfilt, coeffs, temp_seds, pz = readEazyBinary(MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, CACHE_FILE = CACHE_FILE)
    
    ##### Apply zeropoint factors
    param = EazyParam(PARAM_FILE=OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.param')
    fnumbers = np.zeros(len(param.filters), dtype=np.int)
    for i in range(len(fnumbers)):
        fnumbers[i] = int(param.filters[i].fnumber)
    
    zpfile = OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zeropoint'
    if os.path.exists(zpfile):
        zpfilts, zpf_file = np.loadtxt(zpfile, unpack=True, dtype=np.str)                                    
        zpf = np.ones(tempfilt['NFILT'])
        for i in range(len(zpfilts)):
            match = fnumbers == int(zpfilts[i][1:])
            zpf[match] = np.float(zpf_file[i])
    else:
        zpf = np.ones(tempfilt['NFILT'])

    zpfactors = np.dot(zpf.reshape(tempfilt['NFILT'],1),\
                       np.ones(tempfilt['NOBJ']).reshape(1,tempfilt['NOBJ']))

    if verbose:
        print(zpf)
        
    tempfilt['fnu'] *= zpfactors
    tempfilt['efnu'] *= zpfactors
    
    lci = tempfilt['lc'].copy()
    
    params = EazyParam(PARAM_FILE=OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.param')
    abzp = np.float(params['PRIOR_ABZP'])
        
    # fobs = tempfilt['fnu'][:,idx]/(lci/5500.)**2*flam_factor
    # efobs = tempfilt['efnu'][:,idx]/(lci/5500.)**2*flam_factor
    ### Physical f_lambda fluxes, 10**-17 ergs / s / cm2 / A
    if scale_flambda:
        flam_factor = 10**(-0.4*(params['PRIOR_ABZP']+48.6))*3.e18/1.e-17
    else:
        flam_factor = 5500.**2
    
    missing = (tempfilt['fnu'][:,idx] < -99) | (tempfilt['efnu'][:,idx] < 0)
    fobs = tempfilt['fnu'][:,idx]/lci**2*flam_factor
    efobs = tempfilt['efnu'][:,idx]/lci**2*flam_factor
    fobs[missing] = -99
    efobs[missing] = -99
    #print lci, tempfilt['fnu'][:,idx], tempfilt['efnu'][:,idx]
    
    ##### Broad-band SED
    obs_sed = np.dot(tempfilt['tempfilt'][:,:,coeffs['izbest'][idx]],\
                     coeffs['coeffs'][:,idx])/(lci)**2*flam_factor
    
    zi = tempfilt['zgrid'][coeffs['izbest'][idx]]
    
    ###### Full template SED, observed frame
    lambdaz = temp_seds['templam']*(1+zi)
    temp_sed = np.dot(temp_seds['temp_seds'],coeffs['coeffs'][:,idx])
    if individual_templates:
        temp_sed = temp_seds['temp_seds']*coeffs['coeffs'][:,idx]
    
    temp_sed /= (1+zi)**2
    
    temp_sed *= (1/5500.)**2*flam_factor
    
    ###### IGM absorption
    lim1 = np.where(temp_seds['templam'] < 912)
    lim2 = np.where((temp_seds['templam'] >= 912) & (temp_seds['templam'] < 1026))
    lim3 = np.where((temp_seds['templam'] >= 1026) & (temp_seds['templam'] < 1216))
    
    if lim1[0].size > 0: temp_sed[lim1] *= 0.
    if lim2[0].size > 0: temp_sed[lim2] *= 1.-temp_seds['db'][coeffs['izbest'][idx]]
    if lim3[0].size > 0: temp_sed[lim3] *= 1.-temp_seds['da'][coeffs['izbest'][idx]]
        
    ###### Done
    return lambdaz, temp_sed, lci, obs_sed, fobs, efobs