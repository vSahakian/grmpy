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
    
    #Put into an inversion object:
    invdat=cdf.invinfo(G,d,m,resid,L2norm,VR,rank,svals,rng,sdist,smth)
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
    from numpy import ones,str,around
    
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
            
    #Get basenae for figure file:        
    basename='regr_Mc'+str(Mc)+'_'+strname+'_VR_'+str(around(model.VR,decimals=1))


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


##########
def run_mixedeffects(home,codehome,dbpath,dbname,modelpath,Mc):
    '''
    '''
    