####Script to run inversion processes#####
#VJS 9/2016


def setup_run_inversion(home,dbpath,dbname,ncoeff,rng,sdist,Mc,smth,mdep_ffdf):
    '''
    Make the necessary matrices, invert, and save model output
    Input:
        home:           String with project home (i.e., anza)
        db:             String to database path
        dbname:         String with database path name (i.e., db2013 for pckl/db2013/)
        ncoeff:         Number of coefficients used in inversion
        rng:            Array of ranges to use in constraining the inversion
        sdist:          Number of distances to include in smoothing - there will be this many extra
                        equations added on at each range boundary
        Mc:             M squared centering term (8.5 in ASK2014)
        smth:           Smoothing
        mdep_ffdf:      Flag to add mag dependent ffdf (0/1=no/yes)
    Output: 
        inversion object:   Stored in model path, based on ranges etc. used in inversion
        
        '''
        
    import cPickle as pickle
    import cdefs as cdf
    import inversion as inv
    import numpy as np

    #Get directories for things:
    obj_dir=home+'/models/pckl/'+dbname+'/'        
                        
    #Open database object:
    dbfile=open(dbpath,'r')
    db=pickle.load(dbfile)
    dbfile.close()
        
                
    #Invert:
    #Make matrices
    G,d=inv.iinit_pga(db,ncoeff,rng,sdist,Mc,smth,mdep_ffdf)
    #Invert
    m, resid, L2norm, VR, rank, svals=inv.invert(G,d)
    
    ##
    #Save G, d, and m.....and other things...
    #Save plots:
    
    
    #Get the string for the filename, based on the ranges:
    for k in range(len(rng)):
        if k==0:
            strname=np.str(rng[k])
        else:
            strname=strname+'_'+np.str(rng[k])
        
    basename='regr_Mc'+str(Mc)+'_'+strname+'_VR_'+np.str(np.around(VR,decimals=1))
    
    # This is a normal inversion, so set stderror and tvalue to "NaN", since they do not apply:
    stderror = float('NaN')
    tvalue = float('NaN')
    
    #Put into an inversion object:
    invdat=cdf.invinfo(G,d,m,resid,L2norm,VR,rank,svals,rng,sdist,smth,stderror,tvalue)
    fname=obj_dir+basename+'.pckl'
    datobj=open(fname,'w')
    pickle.dump(invdat,datobj)
    datobj.close()
    
    #Return the model info...
    return invdat
    
    
############
def plot_data_model(home,dbpath,dbname,modelpath,coeff_file,mdep_ffdf,sdist,Mc,axlims,bmin,bmax,vref):
    '''
    Plot the data with the model, and ASK 2014
    Input:
        home:           String with path to the project directory (ie., Anza)
        dbpath:         String with path to the database object
        dbname:         String with database path name (i.e., db2013 for pckl/db2013/)
        modelpath:      String with path to the model object
        coeff_file:     String with path to the coefficient file
        mdep_ffdf:      Flag for using mag dependent ffdf (0/1=off/on)
        sdist:          Array with distances to plot in model
        Mc:             M squared centering term (8.5 in ASK 2014)
        axlims:         Axis limits for plotting [[xmin,xmax],[ymin,ymax]]
        bmin:           Minimum value for distance bins
        bmax:           Maximum value for distance bins
        vref:           Vs30 reference value for computing prediction
    Output: 
        figure:         Figure with data and model (saves to /figs directory)
    '''
    
    import gmpe as gm
    import cPickle as pickle
    from numpy import ones,str,around,isnan
    
    #Get directories for things:
    fig_dir=home+'/models/figs/'+dbname+'/'
    
    #Read in database:
    dbfile=open(dbpath,'r')
    db=pickle.load(dbfile)
    dbfile.close()
    
    #Read in model:
    mfile=open(modelpath,'r')
    model=pickle.load(mfile)
    mfile.close()    
    
    #Get the string for the figure filename, based on the ranges:
    for k in range(len(model.rng)):
        if k==0:
            strname=str(model.rng[k])
        else:
            strname=strname+'_'+str(model.rng[k])
              
    # Get basename for figure file, depending on if it's mixed or not:
    # If it's a normal inversion, it will have a rank:
    if isnan(model.rank) == False:
        basename = 'regr_Mc'+str(Mc)+'_'+strname+'_VR_'+str(around(model.VR,decimals=1))
    #If it's from mixed effects, this is set to NaN:
    elif isnan(model.rank) == True:
        basename = 'mixedregr_Mc'+str(Mc)+'_VR_'+np.str(np.around(model.VR,decimals=1))


    ############################################
    ###Compute model from inversion and ASK#####
    ############################################
    
    #Compute the magnitude/log10pga for each distance, to plot on top of data:
    mw_model,d_model=gm.compute_model_fixeddist(model.m,model.rng,sdist,mdep_ffdf)

    #Get the NGA predictions to plot on the same figure:
    #Coefficient file:
    coeff_file=coeff_file=home+'/data/coeffs/ASK14_coeffs.m'
    #Do it just for one distance for now, say R=5km.  
    Rrup=5*ones(db.r.shape)
    #Get the NGA predictions...
    freq1,M_sort,freq1_sort=gm.ask2014_pga(db.mw,Rrup,coeff_file,1,[0,0])


    ############
    ####Plot####
    ############
    
    #Plotting params...
    #Plot against data to check:
    fig1=db.plot_rpga_withmodel(bmin,bmax,mw_model,d_model,model.rng,sdist,axlims,model.VR,M_sort,freq1_sort,vref)


    #Save figure as pdf and png:
    figname=fig_dir+basename+'.png'
    figpdf=fig_dir+'pdf/'+basename+'.pdf'
    
    fig1.savefig(figname)
    fig1.savefig(figpdf)
    
    #Show and return the figure:
    fig1.show()
    return fig1


#def write_inversion_stats(home,)


##########
def run_mixedeffects(home,codehome,run_name,dbpath,dbname,Mc,vref,c):
    '''
    Run a mixed effects model for a given database, and certain parameters.
    Input:
        home:       Working home (i.e., /media/vsahakian/katmai/anza), with no slash at the end
        codehome:   Home for code (i.e., /home/vsahakian/software), with no slash at the end
        run_name:   Databse/inversion combo run name for mixed effects
        dbpath:     Full path to the database object
        dbname:     Database name (i.e., 'test2013'), for path purposes
        Mc:         Magnitude around which to center mag squared term, i.e., (8.5 - M)**2
        vref:       Reference vs30 velocity for GMPE
        c:          ffdf parameter for GMPE
    '''
    
    import cPickle as pickle
    import inversion as inv
    import cdefs as cdf
    import numpy as np
    from os import path
    
    #Open database:
    dbfile=open(dbpath,'r')
    db=pickle.load(dbfile)
    dbfile.close()
    
    # Get path names for output files:
    run_dir=path.expanduser(home+'/models/residuals/'+run_name+'/')
    
    # Get information that is needed to run inversion:
    
    pga = db.pga_pg
    mw = db.mw
    rrup = db.r
    vs30 = db.vs30
    evnum = db.evnum
    sta = db.sta
    
    # Run Model
    log, fixed, event, site = inv.mixed_effects(codehome,home,dbname,pga,mw,rrup,vs30,evnum,sta,vref,c,Mc)
    
    # Add these to an inversion object....set unused values to nan:
    d = np.log10(pga) - 0.6*np.log(vs30/vref)
    rng = [min(mw),max(mw)]
    
    G = float('NaN')
    resid = float('NaN')
    rank = float('NaN')
    svals = float('NaN')
    sdist = float('NaN')
    smth = float('NaN')

    # Now get the ones that were included:
    model = fixed[:,0]
    stderror = fixed[:,1]
    tvalue = fixed[:,2]


    # Get L2norm and VR #
    # Basic stuff:
    ffdf_factor = np.sqrt(c**2 + rrup**2)
    prediction = model[0] + (model[1]*mw) + model[2]*((Mc - mw)**2) + model[3]*np.log(ffdf_factor) + model[4]*rrup
    modelresid = d - prediction
    
    # L2norm:
    # Square these, sum them, get the square root:
    L2norm = np.sqrt(np.sum(modelresid**2))
    
    # Variance reduction:
    # Get numerator and denominator first
    VR_top = np.sum(modelresid**2)
    VR_bot = np.sum(d**2)
    VR = (1 - (VR_top/VR_bot))*100
    
    
    # Make inversion object to store later:
    invdat = cdf.invinfo(G,d,model,resid,L2norm,VR,rank,svals,rng,sdist,smth,stderror,tvalue)
    
    
    # Now also get residuals object to store later per recording:
    totalresid = modelresid
    total_mean = np.mean(modelresid)
    total_std = np.std(modelresid)
    
    # Event
    eventterm = event[:,0]
    eventmean = np.mean(eventterm)
    eventstd = np.std(eventterm)
    
    # Within Event:
    weterm = totalresid - eventterm
    wemean = np.mean(weterm)
    westd = np.std(weterm)
    
    # Site/Station term:
    siteterm = site[:,0]
    sitemean = np.mean(siteterm)
    sitestd = np.std(siteterm)
    
    # Path term plus aleatory per recording:
    pathterm = totalresid - (eventterm + siteterm)
    pathmean = np.mean(pathterm)
    pathstd = np.std(pathterm)
    
    
    ### Make the total residuals object:
    tresid = cdf.total_residuals(mw,totalresid,total_mean,total_std)
    
    ## Make all residuals object:
    mixedresid = cdf.mixed_residuals(db,totalresid,total_mean,total_std,eventterm,eventmean,eventstd,weterm,wemean,westd,siteterm,sitemean,sitestd,pathterm,pathmean,pathstd)
    
    ## Save objects ##
    basename = 'mixedregr_Mc'+str(Mc)+'_VR_'+np.str(np.around(VR,decimals=1))
    
    ## Save inversion object:
    invpath = home + '/models/pckl/' + dbname + '/' + basename + '.pckl'
    invfile = open(invpath,'w')
    pickle.dump(invdat,invfile)
    invfile.close()
    
    print 'Saved inversion object'
    
    ## Save Mixed residuals object:
    objpath = run_dir + basename + '_robj.pckl'
    objfile = open(objpath,'w')
    pickle.dump(mixedresid,objfile)
    objfile.close()
    
    print 'Saved all residuals object'
    
    
    # return mixed residuals object:
    return invdat, invpath, tresid, mixedresid