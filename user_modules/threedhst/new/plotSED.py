def plotExampleSED(idx=20, writePNG=True, MAIN_OUTPUT_FILE = 'photz', OUTPUT_DIRECTORY = 'OUTPUT', CACHE_FILE = 'Same', lrange=[3000,8.e4], axes=None, individual_templates=False, fnu=False):
    """
    PlotSEDExample(idx=20)

    Plot an example Eazy best-fit SED.
    """

    #zout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
    zout = catIO.Readfile(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
    #qz = np.where(zout.z_spec > 0)[0]
    print zout.filename
    qz = np.arange(len(zout.id))
    
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
        getEazySED(qz[idx], MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                          OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                          CACHE_FILE = CACHE_FILE, individual_templates=individual_templates, scale_flambda=True)
    
    zgrid, pz = getEazyPz(qz[idx], MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                                   OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                                   CACHE_FILE = CACHE_FILE)

    ##### start plot
    fig = plt.figure()
    
    #### Plot parameters
    plotsize=35
    alph=0.9
    
    if fnu:
        temp_sed *= (lambdaz / 5500.)**2
        fobs *= (lci/5500.)**2
        efobs *= (lci/5500.)**2
        obs_sed *= (lci/5500.)**2

    #### Full best-fit template
    axp = fig.add_subplot()

    
    ##### P(z)
    if pz is not None:            
        axp.plot(zgrid, pz, linewidth=1.0, color='orange',alpha=alph)
        axp.fill_between(zgrid,pz,np.zeros(zgrid.size),color='yellow')

        if zout.z_spec[qz[idx]] > 0:
            axp.plot(zout.z_spec[qz[idx]]*np.ones(2), np.array([0,1e6]),color='red',alpha=0.4)

        #### Set axis range and titles
        axp.set_xlim(0,np.ceil(np.max(zgrid)))
        axp.set_ylim(0,1.1*max(pz))
        axp.set_xlabel(r'$z$')
        axp.set_ylabel(r'$p(z)$')
        
    if (writePNG is not False) & (axes is None):
        if isinstance(writePNG, str):
            out=writePNG
        else:
            out='/tmp/test.pdf'
            
        fig.savefig(out,dpi=100)

    if axes is None:
        return [ax, axp]
