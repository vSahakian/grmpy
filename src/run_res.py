####Run the residual analysis####
'''
Run the Residual analysis for a given parameter file
VJS 7/2016
'''


#Initialize the folders
def init(home,run_name):
    '''
    Initialize the folders for output for a new run (new database/model combo)
    Input:
        home:           String to path of the home environment
        run_name:       String with the name for this db/model combo
    Output: 
        Makes directories...
    '''
    
    from shutil import rmtree,copy
    from os import path,makedirs,environ
    
    clob='y'
    run_dir=path.expanduser(home+run_name+'/')
    if path.exists(run_dir):    #If the path exists...
        print 'Run directory is: '+run_dir
        clob=raw_input('The run directory already exists, are you sure you want to clobber it? (y/n)')
        if clob is 'y' or clob is 'Y':
            clob=raw_input('This deletes everything in '+run_dir+', are you sure you want to do that? (y/n)')
            if clob is 'y' or clob is 'Y':
                rmtree(run_dir)
            else:
                #Don't do anything
                print 'Aborting mission...not clobbering...'
        else:
            #AGain, don't do anything
            print 'Aborting mission...not clobbering...'
    if clob is 'y' or clob is 'Y':
        #Make the main directory
        makedirs(run_dir)
        #And the subdirectories:
        makedirs(run_dir+'event_objs/')
        makedirs(run_dir+'sta_objs/')
        makedirs(run_dir+'E_resids/')
        makedirs(run_dir+'W_resids/')
        makedirs(run_dir+'site_resids/')
        makedirs(run_dir+'figs/')
        
        
def get_total_res(home,run_name,dbpath,modelpath,ffdf_flag,resaxlim):
    '''
    Get the total residuals for a database.
    Input:
        home:            STring with home path
        run_name:        String with name of the run
        dbpath:          String with path to a pickle database object
        modelpath:       String with path to a gmpe model pickle object
        ffdf_flag:       Flag for mag dependent ffdf.  0=off, 1=on
        resaxlim:        Array with axis limits for the resid plot: [[xmin,xmax],[ymin,ymax]]
        
    Output:
        total_residuals:    Array with total residuals
        mw:                 Array, same dim as total_residuals, with matching mw
        mean_residual:      Mean residual
        std_dev_resid:      Standard deviation of the residuals
    '''
    
    import pickle
    import gmpe as gm
    import rescomp as rcomp
    from numpy import where
    from os import path
    import cdefs as cdf
    from matplotlib.pyplot import savefig
    
    run_dir=path.expanduser(home+run_name+'/')
    
    
    ##########Open the database object:###################
    #Filename:
    dname=dbpath
    datobj=open(dname,'r')
    db=pickle.load(datobj)
    datobj.close()
    
    ##Residual Computation####
    #Load in the model to use:
    mname=modelpath
    datobj=open(mname,'r')
    model=pickle.load(datobj)
    datobj.close()
    
    #Overall residual, 
    #In some places, vs30 is 0.  Set these to vref.
    vref=760
    mdep_ffdf=ffdf_flag
    #Where are they 0?
    vs30_0ind=where(db.vs30==0)[0]
    #Get vs30 from database...
    vs30=db.vs30
    #Set 0 entries to vref:
    vs30[vs30_0ind]=vref
    
    
    #####Get Total Residuals######
    
    #Now compute the predicted value of PGA...
    d_predicted=gm.compute_model(model.m,model.rng,db.mw,db.r,db.ffdf,vs30,vref,mdep_ffdf)
    
    #Get residuals:
    total_residuals,mean_residual,std_dev=rcomp.total_residual(db,d_predicted)
    
    #Plot residuals and save plots:
    allresid=cdf.total_residuals(db.mw,total_residuals,mean_residual,std_dev)
    f1,f2=allresid.plt_resids(run_name,resaxlim)
    
    #Figure names:
    f1name=run_dir+'figs/'+run_name+'_total.png'
    f2name=run_dir+'figs/'+run_name+'_total_hist.png'
    #Save:
    f1.savefig(f1name)
    f2.savefig(f2name)
    
    return allresid
    
    
def get_E_W_residuals():
    '''
    '''