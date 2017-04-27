####Script to run inversion processes#####
#VJS 9/2016


def setup_run_inversion(home,dbpath,dbname,ncoeff,rng,sdist,Mc,smth,vref,mdep_ffdf,predictive_parameter='pga',data_correct=-0.6):
    '''
    Make the necessary matrices, invert, and save model output
    Input:
        home:                   String with project home (i.e., anza)
        dbpath:                 String to database path
        dbname:                 String with database path name (i.e., db2013 for pckl/db2013/)
        ncoeff:                 Number of coefficients used in inversion
        rng:                    Array of ranges to use in constraining the inversion
        sdist:                  Number of distances to include in smoothing - there will be this many extra
                                   equations added on at each range boundary
        Mc:                     M squared centering term (8.5 in ASK2014)
        smth:                   Smoothing
        vref:                   Reference vs30 value (like 760 m/s)
        mdep_ffdf:              Flag to add mag dependent ffdf (0/1=no/yes)
        predictive_parameter:   Default is pga.  Else, 'pgv', or...
        data_correct:           Correct data for a term?  0/1 = no/yes - vs30 term DEFAULT: vs30 correct
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
        
    ###  DEBUGGING...
    print 'ncoeff is %i, predictive parameter is %s, and data_correct is %i' % (ncoeff,predictive_parameter,data_correct)
    #####    
        
    print 'data correct is: %f' % data_correct
    
    print 'smth is %i, vref is %i' % (smth,vref)
                
    #Invert:
    #Make matrices
    G,d=inv.iinit_predparam(db,ncoeff,rng,sdist,Mc,smth,vref,mdep_ffdf,predictive_parameter=predictive_parameter,data_correct=data_correct)
    #Invert
    m, resid, L2norm, VR, rank, svals=inv.invert(G,d)

    print 'Inversion residual is: '
    print G.dot(m) - d
    
    print 'Mean of inversion residual is: '
    print np.mean(G.dot(m) - d)

    
    #Get the string for the filename, based on the ranges:
    for k in range(len(rng)):
        if k==0:
            strname=np.str(rng[k])
        else:
            strname=strname+'_'+np.str(rng[k])
        
    basename='regr_' + predictive_parameter + '_Mc'+str(Mc)+'_'+strname+'_VR_'+np.str(np.around(VR,decimals=1))
    
    # This is a normal inversion, so set stderror and tvalue to "NaN", since they do not apply:
    stderror = float('NaN')
    tvalue = float('NaN')
    
    #### DEBUGGING....
    print '\n ln the data, pga, are: \n'
    print np.log(db.pga_pg)
    print '\n the model, what goes into model.d, is: \n'
    print d
    
    
    #Put into an inversion object:
    invdat=cdf.invinfo(G,d,m,resid,L2norm,VR,rank,svals,rng,sdist,smth,stderror,tvalue)
    fname=obj_dir+basename+'.pckl'
    datobjfile=open(fname,'w')
    pickle.dump(invdat,datobjfile)
    datobjfile.close()
    
    # Now, also print the model coefficients to a file:
    model_coeff_file = obj_dir + 'coeff_' + basename + '.txt'
    coeff=open(model_coeff_file,'w')
    coeff.write('Model Coefficients for inversion' + basename + ':')
    
    for coeff_i in range(ncoeff):
        coeff.write('\n a' + str(coeff_i+1) + ' = ' + str(m[coeff_i]))
        
    coeff.close()
    
    #Return the model info...
    return invdat
    
    
############
def plot_data_model(home,dbpath,dbname,modelpath,coeff_file,mdep_ffdf,sdist,ask_dist,Mc,axlims,bmin,bmax,vref,predictive_parameter='pga',ncoeff=5,data_correct=-0.6):
    '''
    Plot the data with the model, and ASK 2014
    Input:
        home:                   String with path to the project directory (ie., Anza)
        dbpath:                 String with path to the database object
        dbname:                 String with database path name (i.e., db2013 for pckl/db2013/)
        modelpath:              String with path to the model object
        coeff_file:             String with path to the coefficient file
        mdep_ffdf:              Flag for using mag dependent ffdf (0/1=off/on)
        sdist:                  Array with distances to plot in model
        ask_dist:               One number with the distance to plot ASK
        Mc:                     M squared centering term (8.5 in ASK 2014)
        axlims:                 Axis limits for plotting [[xmin,xmax],[ymin,ymax]]
        bmin:                   Minimum value for distance bins
        bmax:                   Maximum value for distance bins
        vref:                   Vs30 reference value for computing prediction
        predictive_parameter:   Predictive parameter. 'pga', or 'pgv'.  Default: 'pga'
        ncoeff:                 Number of coefficients that were inverted for.  Default: 5
        data_correct:           Correct data?  0/1 = no/by vs30 term.  Default: 1
    Output: 
        figure:         Figure with data and model (saves to /figs directory)
    '''
    
    import gmpe as gm
    import cPickle as pickle
    from numpy import ones,str,around,isnan,linspace
    
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
        basename = 'regr_'+predictive_parameter+'_Mc'+str(Mc)+'_'+strname+'_VR_'+str(around(model.VR,decimals=1))
    #If it's from mixed effects, this is set to NaN:
    elif isnan(model.rank) == True:
        basename = 'mixedregr_'+predictive_parameter + '_ncoeff' + str(ncoeff)+'_Mc'+str(Mc)+'_VR_'+str(around(model.VR,decimals=1))

    ############################################
    ###Compute model from inversion and ASK#####
    ############################################
    
    #Compute the magnitude/log10pga for each distance, to plot on top of data:
    mw_model,d_model=gm.compute_model_fixeddist(model.m,model.rng,sdist,Mc,mdep_ffdf,ncoeff=ncoeff)

    #Get the NGA predictions to plot on the same figure:
    #Coefficient file:
    coeff_file=coeff_file=home+'/data/coeffs/ASK14_coeffs.m'

    # get the magnitude at which to compute:
    ask_mw=linspace(axlims[0][0],axlims[0][1])
    
    #Do it just for one distance for now, say R=5km.  
    Rrup=ask_dist*ones(ask_mw.shape)
    
    #Get the NGA predictions...
    freq1,M_sort,freq1_sort=gm.ask2014(ask_mw,Rrup,coeff_file,0,[0,0])


    ############
    ####Plot####
    ############
    
    #Plotting params...
    #Plot against data to check:
    fig1=db.plot_rpga_withmodel(bmin,bmax,mw_model,d_model,model.rng,sdist,ask_dist,axlims,model.VR,M_sort,freq1_sort,vref,predictive_parameter=predictive_parameter)


    #Save figure as pdf and png:
    figname=fig_dir+basename+'.png'
    figpdf=fig_dir+'pdf/'+basename+'.pdf'
    
    fig1.savefig(figname)
    fig1.savefig(figpdf)
    
    print 'Saved figure to %s, and %s' % (figname,figpdf)
    
    #Show and return the figure:
    fig1.show()
    return fig1


#def write_inversion_stats(home,)


##########
def run_mixedeffects(home,codehome,run_name,dbpath,dbname,Mc,vref,c,predictive_parameter='pga',ncoeff=5,data_correct=-0.6,a1='none',a2='none',a3='none',a4='none',a5='none',a6='none'):
    '''
    Run a mixed effects model for a given database, and certain parameters.
    Input:
        home:                       Working home (i.e., /media/vsahakian/katmai/anza), with no slash at the end
        codehome:                   Home for code (i.e., /home/vsahakian/software), with no slash at the end
        run_name:                   Databse/inversion combo run name for mixed effects
        dbpath:                     Full path to the database object
        dbname:                     Database name (i.e., 'test2013'), for path purposes
        Mc:                         Magnitude around which to center mag squared term, i.e., (8.5 - M)**2
        vref:                       Reference vs30 velocity for GMPE
        c:                          ffdf parameter for GMPE
        predictive_parameter:       Parameter to predict.  Default: 'pga'
        ncoeff:                     Number of coefficients inverting for.  Default: 5
        data_correct:               Vs30 coefficient to correct data by, if correct at all.  0/data_correct = no correction/correction by data_correct.  Default: -0.6
        a1:                         a1 coefficient, if it's being fixed.  Default:'none'
        a2:                         a2 coefficient, if it's being fixed.  Default:'none'
        a3:                         a3 coefficient, if it's being fixed.  Default:'none'
        a4:                         a4 coefficient, if it's being fixed.  Default:'none'
        a5:                         a5 coefficient, if it's being fixed.  Default:'none'
        a6:                         a6 coefficient, if it's being fixed.  Default:'none'.  If data_correct!=0, should be the same as data_correct
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
    
    print dbpath
    
    # Get path names for output files:
    run_dir=path.expanduser(home+'/models/residuals/'+run_name+'/')
    
    # Get information that is needed to run inversion:
    
    #define the predictive parameter
    if predictive_parameter=='pga':
        pred_param = db.pga_pg
    elif predictive_parameter=='pgv':
        pred_param = db.pgv
        
    # Now define the otherv alues...
    mw = db.mw
    rrup = db.r
    vs30 = db.vs30
    evnum = db.evnum
    sta = db.sta
    
    # Run Model
    me_log, fixed, event, site, d_r_prediction, d_observed, event_mean, event_std, site_mean, site_std = inv.mixed_effects(codehome,home,dbname,pred_param,mw,rrup,vs30,evnum,sta,vref,c,Mc,predictive_parameter=predictive_parameter,ncoeff=ncoeff,data_correct=data_correct,a1=a1,a2=a2,a3=a3,a4=a4,a5=a5,a6=a6)
    
    ## Add these to an inversion object....set unused values to nan:
    #if data_correct==0:
    #    d = np.log(pred_param)
    #elif data_correct!=0:
    #    d = np.log(pred_param) - data_correct*np.log(vs30/vref)
    
    ## Add the "observed" (corrected d) to the inversion object:
    d = d_observed
  
    # Get ranges for inversion object (like ranges from traditional inversion)  
    if min(mw)>0:
        rngmin=0
    else:
        rngmin=min(mw)
        
    if max(mw)<6.5:
        rngmax=6.5
    else:
        rngmax=max(mw)
    rng = [rngmin,rngmax]
    
    G = float('NaN')
    resid = float('NaN')
    rank = float('NaN')
    svals = float('NaN')
    sdist = float('NaN')
    smth = float('NaN')

    # Now get the model values.
    # For the model, include the values that were inverted for, but also what was provided - since this goes into a saved object.
    
    # First need to make lists of the variables that were prescribed or fixed, and the variables that were inverted for.
    # Concatenate all variables into a list:
    coeff_list = [a1,a2,a3,a4,a5,a6]
    
    # Set empty lists that will store the "indices" of the fixed or inverted values,
    r_fixed=[]
    r_invert=[]
    #      and also the model, stderror, and tvalues
    model = np.zeros(ncoeff)
    stderror = np.zeros(ncoeff)
    tvalue = np.zeros(ncoeff)
    #      and an empty string for naming the model:
    modelpath_list = ''
    
    print 'fixed'
    print fixed
    
    # Also set the r model index counter to 0:
    r_invert_counter=0
    
    # Loop over the coefficients in the list (a1 - a6), to get the values that are fixed, andvalues that are inverted:
    for coeffi in range(min(ncoeff,len(coeff_list))):
        
        # If this coeff was not provided, then it's inverted for:
        if coeff_list[coeffi]=='none':
            # So append the index of it to the r_invert list (probably unnecessary):
            r_invert.append(coeffi)
            
            # And get the inverted value from the model array - the index of this is the same as r_invert_counter,
            #   as we're looping through this index separatelyk, only using it when we pull out an inverted value:
            
            # If only the intercept was inverted for, fixed will be one dimensional (length of shape will be 1), so only one value of a model coeff (fixed[0]):                
            if len(np.shape(fixed))==1:
            	model[coeffi] = fixed[0]
            	stderror[coeffi]=fixed[0]
            	tvalue[coeffi]=fixed[0]
            # Otherwise if more were inverted for, len(shape(fixed)) will have two values, so len 2:
            elif len(np.shape(fixed))==2:
                model[coeffi]=fixed[r_invert_counter,0]
                stderror[coeffi]=fixed[r_invert_counter,1]
                tvalue[coeffi]=fixed[r_invert_counter,2]

                
            # Add to counter
            r_invert_counter+=1
            
            
        # If this coeff was provided, it's fixed:
        elif coeff_list[coeffi]!='none':
            # So pull out the index for the "fixed" array, and then get its actual value for the model:
            r_fixed.append(coeffi)
            model[coeffi] = coeff_list[coeffi]
            # In addition, set stderror and tvalue to
            stderror[coeffi] = float('NaN')
            tvalue[coeffi] = float('NaN')
            
            # and also add a list of the names of fixed coefficients with their numbers:
            coeff_num = coeffi + 1
            coeff_name = 'a%i' % coeff_num
            print coeff_num
            print coeff_name
            print modelpath_list
            modelpath_list = modelpath_list + '_'+coeff_name+'_'+np.str(np.around(model[coeffi],decimals=2))
    
            print 'coeffi %i was fixed' % coeffi


    # Get L2norm and VR #
    # Basic stuff:
    
    # This is bad because the ME makes a prediction including the random effects....so this isn't the full prediction
    #R = np.sqrt(rrup**2 + c**2)
    #prediction = model[0] + (model[1]*mw) + model[2]*((Mc - mw)**2) + model[3]*np.log(R) + model[4]*rrup
    
    # Instead, use the prediction from R:
    modelresid = d - d_r_prediction
    
    
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
    
    
    ################  RESIDUALS   ############
    # Now also get residuals object to store later per recording:
    # Consider the "total" residuals to be the traditional total -
    #   meaning Gm - d.....not th emixed effects residuals, which is
    #   Gm + An - d   where An represents the site and event effects.
        
    # First get event and site terms:    
    # Event
    eventterm = event[:,0]
    event_stderr = event[:,1]

    #eventmean = np.mean(eventterm)
    #eventstd = np.std(eventterm)
    
    # Site/Station term:
    siteterm = site[:,0]
    site_stderr = site[:,1]
    #sitemean = np.mean(siteterm)
    #sitestd = np.std(siteterm)
    
    # Now get the total residual - it's the aleatory residual, from the observed
    # minus r prediction, plus the event and site terms from r:
    totalresid = (d - d_r_prediction)   # + eventterm + siteterm ..... # used to add this, but norm says this is a bad represenation.
    total_mean = np.mean(totalresid)
    total_std = np.std(totalresid)
    
    # Within Event:
    #weterm = totalresid - eventterm
    weterm = totalresid + siteterm   # Now that totalresid is just the path term, this is the same as the path + site
    wemean = np.mean(weterm)
    westd = np.std(weterm)
    
    
    
    # Path term plus aleatory per recording:
    #pathterm = totalresid - (eventterm + siteterm)   # Note this is the same as d - d_r_prediction!
    pathterm = totalresid 
    pathmean = np.mean(pathterm)
    pathstd = np.std(pathterm)
    
    
    ### Make the total residuals object:
    tresid = cdf.total_residuals(mw,totalresid,total_mean,total_std)
    
    ## Make all residuals object:
    mixedresid = cdf.mixed_residuals(db,totalresid,total_mean,total_std,eventterm,event_stderr,event_mean,event_std,weterm,wemean,westd,siteterm,site_stderr,site_mean,site_std,pathterm,pathmean,pathstd)
    
    ## Save objects ##
    basename = 'mixedregr_' + dbname + '_' + predictive_parameter + '_' + '_ncoeff' + str(ncoeff) + '_Mc_'+str(Mc)+'_VR_'+np.str(np.around(VR,decimals=1))+modelpath_list
    
    ## Save inversion object:
    invpath = home + '/models/pckl/' + dbname + '/' + basename + '.pckl'
    invfile = open(invpath,'w')
    pickle.dump(invdat,invfile)
    invfile.close()
    
    print 'Saved inversion object to %s' % invpath
    
    
    ## Pritn out the coefficients to a file:
    basename = 'mixedcoeff_' + dbname + '_' + predictive_parameter + '_' + '_ncoeff' + str(ncoeff) + '_Mc_'+str(Mc)+'_VR_'+np.str(np.around(VR,decimals=1))+modelpath_list
    coeffpath = home + '/models/pckl/' + dbname + '/' + basename + '.txt'
    headerstring = 'coefficient \t std error \t tvalue '
    outdata = np.c_[model,stderror,tvalue]
    dataformat = '%8.4f\t%9.5f\t%9.5f'
    np.savetxt(coeffpath,outdata,fmt=dataformat,header=headerstring)
    
    print 'Saved coefficients to text file %s' % coeffpath
    
    
    ## Save Mixed residuals object:
    objpath = run_dir + basename + '_robj.pckl'
    objfile = open(objpath,'w')
    pickle.dump(mixedresid,objfile)
    objfile.close()
    
    print 'Saved all residuals object'
    
    
    # return mixed residuals object:
    return invdat, invpath, tresid, mixedresid, d_r_prediction, objpath