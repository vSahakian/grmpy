######Class definition######
#VJS 6/2016

'''
Generic classes that are useful for this package...
'''


class db:
    '''
    This class describes a set of events
    '''
    
    def __init__(self,event,sta,N,ml,mw,DA,DV,r,vs30,elat,elon,edepth,stlat,stlon,stelv,source_i,receiver_i,):
        '''
        Initiate the class by giving the event number (event), 
        station name (sta), station number (N), local mag (ml), moment mag (mw), 
        PGA (DA), PGV (DV), and source to site distance (r)
        '''
        import numpy as np
        
        #"Fictitions depth" or "Finite fault dimension factor"
        c=4.5
        
        #Save pga + pgv in m/s/s, not nm/s/s
        DAm=DA*1e-9
        DVm=DV*1e-9
        
        #Get percent g:
        pga_pg=DAm/9.81
        
        #Get magnitude-dependent ffdf:
        #ASK2014 c4 coefficient:
        c4=4.5
        #Find the indices for each range:
        cr1_ind=np.where(mw>5)
        cr2_ind=np.where((mw<=5) & (mw>4))
        cr3_ind=np.where(mw<=4)
        
        #Zero out the c array:
        c=np.zeros(mw.shape)
        c[cr1_ind]=c4
        c[cr2_ind]=c4-((c4-1)*(5-mw[cr2_ind]))
        c[cr3_ind]=1
        md_ffdf=np.sqrt(r**2 + c**2)
        
        #Give these values to the db:
        self.evnum=event
        self.sta=sta
        self.stnum=N
        self.ml=ml
        self.mw=mw
        self.pga=DAm
        self.pgv=DVm
        self.pga_pg=pga_pg
        self.r=r
        self.vs30=vs30
        self.ffdf=np.sqrt(self.r**2 + c4**2)
        self.md_ffdf=md_ffdf
        self.elat=elat
        self.elon=elon
        self.edepth=edepth
        self.stlat=stlat
        self.stlon=stlon
        self.stelv=stelv
        self.source_i=source_i
        self.receiver_i=receiver_i
        
    def plot_allpga(self):
        '''
        Plots log10 of all PGA, regardless of M/r
        '''
        
        from matplotlib import pyplot as plt
        import numpy as np
        
        #Open figure
        plt.figure()
        
        #Plot...
        plt.scatter(self.mw,np.log10(self.pga_pg),marker='o')
        plt.xlabel(r"$\mathbf{M}$")
        plt.ylabel(r"$\log_{10} PGA$")
        plt.title(r"PGA vs. $\mathbf{M}$ for all distances")
        
        plt.show()
        
    def plot_rpga(self,bmin,bmax,step):
        '''
        Plots log10 PGA, for various distance ranges specified by bmin, bmax,
        and step.
        Input:
            bmin:       Min value for bins
            bmax:       Max balue for bins
            step:       Step interval for bins
        '''
        
        from matplotlib import pyplot as plt
        import numpy as np
        
        #Get bins:
        bins=np.arange(bmin,bmax,step)
        #Get the bin index for each recording:
        dinds=np.digitize(np.floor(self.r),bins)
        
        #Sort data to plot:
        
        #Open figure
        f=plt.figure()
        
        #Define color scale -
        #colormap needs floats that go from 0 to 1, so bins must be normalized.
        #bins is not a float though, and in order to divide by the scalar it 
        #must first be converted to a float. 
        colors=plt.cm.rainbow(bins.astype(float)/bins.max())
        #FIGURE OUT THE COLORBAR PROBLEM!!
        #f.colorbar(colors)
        
        #Plot a different color for each distance bin:
        for i in bins:
            #Find which data points are in this bin:
            bind=np.where(dinds==i+1)[0]
            #Make an array of size len(bind),4 for the colors, so that these can
            #be plotted in scatter as x,y,z (mw, pga, color).  In this bin, all
            #the colors should be the same, so tile the color for this bin i and 
            #multiply it by an array of ones. (maybe I don't even need to do this?)
            clrs=np.ones((len(bind),4))*np.tile(colors[i,:],(len(bind),1))
            #plot
            f=plt.scatter(self.mw[bind],np.log10(self.pga_pg[bind]),edgecolors=clrs,facecolors='none',lw=0.5)

        #Label the plot - Mbold on the x, log10PGA on the y, 
        plt.xlabel(r"$\mathbf{M}$")
        plt.ylabel(r"$\log_{10} PGA$")
        plt.title(r"PGA vs. $\mathbf{M}$, binned by distance")
        
        plt.show()
        
    def plot_rpga_withmodel(self,bmin,bmax,mw,d,rng,sdist,axlims,VR,nga_mw=True,nga_pred=True,vref=True):
        '''
        Plots log10 PGA, for various distance ranges specified by bmin, bmax,
        and step.
        Input:
            bmin:       Min value for bins for data
            bmax:       Max balue for bins for data
            mw:         Mw array from gmpe.compute_model_fixeddist
            d:          d array from compute_model_fixeddist
            rng:        Magnitude ranges, same array used for inversion
            sdist:      Distances array used for inversion
            axlims:     Array with lims: [[xmin,xmax],[ymin,ymax]]
            VR:         Variance Reduction from inversion
            vref:       Reference vs30 value (Default: 760 m/s)
            nga_mw:     mw range for NGA plotting, if provided
            nga_pred:   prediction array for NGA plotting, if provided
        '''
        
        from matplotlib import pyplot as plt
        import matplotlib.colors as colors
        import matplotlib.cm as cm
        import numpy as np
        
        #Vs30 reference:
        if vref==None:
            vref=760
            
        ##Get bins:
        #bins=np.arange(bmin,bmax,step)
        ##Get the bin index for each recording:
        #dinds=np.digitize(np.floor(self.r),bins)
        #
        
        #Get colormap
        mymap='jet'
        #Make colormap:
        colormap_radius=plt.get_cmap(mymap)
        #Make a normalized colorscale
        cNorm=colors.Normalize(vmin=bmin, vmax=bmax)
        #Apply normalization to colormap:
        scalarMap=cm.ScalarMappable(norm=cNorm, cmap=colormap_radius)
        
        #Make a fake contour plot for the colorbar:
        Z=[[0,0],[0,0]]
        levels=np.arange(bmin,bmax,0.01)
        c=plt.contourf(Z, levels, cmap=colormap_radius)
        
        #Open figure
        f=plt.figure()
        
        #get colorvalue to plot
        colorVal=scalarMap.to_rgba(self.r)
        
        plt.scatter(self.mw,np.log10(self.pga_pg),edgecolors=colorVal,facecolors='none',lw=0.5)
        
        #Add colorbar:
        cb=plt.colorbar(c)
        cb.set_label('Distance (km)')

        #Label the plot - Mbold on the x, log10PGA on the y, 
        plt.xlabel(r"$\mathbf{M}$")
        plt.ylabel(r"$\log_{10} PGA$")
        plt.title(r"PGA vs. $\mathbf{M}$, binned by distance" + "\n" + \
            "M Ranges: " + np.str(rng)+ ", Var Red="+np.str(np.around(VR,decimals=1)))
        
        colors_gmpe=plt.cm.rainbow(sdist.astype(float)/sdist.max())
        
        for j in range(len(sdist)):
            #Label for plot:
            lab="R="+np.str(sdist[j])+"km"
            
            mw_dist=mw[:,j]
            d_dist=d[:,j]
            
            #Plot
            plt.plot(mw_dist,d_dist,color=colors_gmpe[j],linewidth=2,label=lab)
            
        ##Plot dummy object
        #xdumdum=array([1e99,1e99])
        #ydumdum=array([1e99,1e99])
        #zdumdum=array([1e99,1e99])
        #dummy=plt.scatter(xdumdum,ydumdum,c=zdumdum,cmap=colors)
        #print colors
        #plt.colorbar(colors)
        
        #Limits:
        xlims=axlims[0]
        ylims=axlims[1]
        plt.xlim(xlims)
        plt.ylim(ylims)
        

        
        #If NGA data are not provided, end and return the figure.
        if nga_mw==None and nga_pred==None:
            #Add legend:
            plt.legend(loc=4)
        
            plt.show(f)
            return f
        #Otherwise, plot the NGA data:
        else:
            plt.plot(nga_mw,nga_pred,linestyle='--',linewidth=2,color='b',label='ASK2014')
            
            #Add legend:
            plt.legend(loc=4)  
            
            plt.show(f)
            return f          
            
        
        
        
class invinfo:
    '''
    Save paramters from an inversion.
       G:       Left hand side matrix
       d:       Data vector
       m:       Resulting model vector
       resid:   Residuals from inversion
       rank:    rank from inversion
       svals:   Singular values from inversion
       rng:     Magnitude ranges used in inversion
       sdist:   Distances used in smoothing for inversion
       smth:    Smoothing value used in inversion 
    '''
    def __init__(self,G,d,m,resid,L2norm,VR,rank,svals,rng,sdist,smth):
        self.G=G
        self.d=d
        self.m=m
        self.resid=resid
        self.norm=L2norm
        self.VR=VR
        self.rank=rank
        self.svals=svals
        self.rng=rng
        self.sdist=sdist
        self.smth=smth
   
class total_residuals:
    '''
    Save all total residual data
    '''
    def __init__(self,mw,total_residuals,mean_residual,std_dev_residuals):
        self.mw=mw
        self.total_residuals=total_residuals
        self.mean_residual=mean_residual
        self.std_dev=std_dev_residuals
        
    def plt_resids(self,run_name,axlims):
        '''
        Plot all residuals
        Input:
            run_name:     Name of the run
            axlims:     Array with axislimits: [[xmin,xmax],[ymin,ymax]]
        Output: 
            f1:         Plot with residuals
            f2:         Plot with histogram
        '''
            
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MultipleLocator
        from numpy import str,array,ones,around
        
        #Color:
        rgb=array([111,168,163])/255.0
        rgb.astype(float)
        color=rgb*ones((len(self.mw),3))
        
        ##
        #Set up plot to have histogram adjacent to scatter...
        ##
        
        #definitions for axes
        fudge_factor=0.02
        left, width=0.1, 0.70
        bottom, height=0.1,0.81
        left_h=left+width+fudge_factor
        width_h=0.15
        
        #define axis limits for scatter, and histogram:
        rect_scatter=[left,bottom,width,height]
        rect_histy=[left_h,bottom,width_h,height]
        
        #define axis tick locations for histogram:
        hist_xLocator=MultipleLocator(500) 
        
        #Start figure:
        f1=plt.figure()
        
        #define axes:
        axScatter=plt.axes(rect_scatter)
        axHisty=plt.axes(rect_histy)
        
        #Scatter:
        axScatter.scatter(self.mw,self.total_residuals,edgecolors=color,facecolors='none',lw=0.9)
        
        #Histogram:
        #want 4x as many bins as main plot y-axis limit units:
        nbins=(axlims[1][1]-axlims[1][0])*4
        
        #set the number of bins, adn the range to be the x axis limits (same as y axis, ln residuals):
        axHisty.hist(self.total_residuals,bins=nbins,range=[axlims[1][0],axlims[1][1]],orientation='horizontal',color=rgb)
        
        #Also plot a dashed line at 0:
        axScatter.plot(axlims[0],[0,0],linestyle='--',color='0.75')
        
        #set axis limits:
        #scatter
        axScatter.set_xlim(axlims[0])
        axScatter.set_ylim(axlims[1])
        #histogram
        axHisty.set_ylim(axlims[1])
        #set axis ticks:
        axHisty.xaxis.set_major_locator(hist_xLocator)
        #set no labels on the y axis:
        axHisty.yaxis.set_ticklabels('')
        
        #Labels
        axScatter.set_xlabel(r"$\mathbf{M}$")
        axScatter.set_ylabel('ln Residuals')
        axScatter.set_title('Total Residuals for run '+run_name+'\n'+'Mean: '+str(around(self.mean_residual,decimals=2))+' Std Dev: '+str(around(self.std_dev,decimals=2)))

        #Show
        f1.show()        
        
        return f1
             
        
class event:
    '''
    Save all data for one event, to use in residual computation
    '''
    
    def __init__(self,evnum,sta,stnum,ml,mw,pga,pgv,pga_pg,r,vs30,ffdf,md_ffdf,elat,elon,edepth,stlat,stlon,stelv,source_i,receiver_i):
        self.evnum=evnum
        self.sta=sta
        self.stnum=stnum
        self.ml=ml
        self.mw=mw
        self.pga=pga
        self.pgv=pgv
        self.pga_pg=pga_pg
        self.r=r
        self.vs30=vs30
        self.ffdf=ffdf
        self.md_ffdf=md_ffdf
        self.elat=elat
        self.elon=elon
        self.edepth=edepth
        self.stlat=stlat
        self.stlon=stlon
        self.stelv=stelv
        self.source_i=source_i
        self.receiver_i=receiver_i
    
    def add_total_resid(self,total_residuals):
        self.total_residual=total_residuals
            
    def add_E_resid(self,E_residual,E_std):
        self.E_residual=E_residual
        self.E_std=E_std
        
    def add_W_resids(self,W_residuals,W_mean,W_std):
        self.W_residuals=W_residuals
        self.W_mean=W_mean
        self.W_std=W_std
        
class station:
    '''
    Save all data for one station
    '''
    
    def __init__(self,sta,stnum,vs30,evnum,ml,mw,pga_pg,pga,pgv,ffdf,md_ffdf,elat,elon,edepth,stlat,stlon,stelv,source_i,receiver_i,total_residual,E_residual,W_residual):
        self.sta=sta
        self.stnum=stnum
        self.vs30=vs30
        self.evnum=evnum
        self.ml=ml
        self.mw=mw
        self.pga_pg=pga_pg
        self.pga=pga
        self.pgv=pgv
        self.ffdf=ffdf
        self.md_ffdf=md_ffdf
        self.elat=elat
        self.elon=elon
        self.edepth=edepth
        self.stlat=stlat
        self.stlon=stlon
        self.stelv=stelv
        self.source_i=source_i
        self.receiver_i=receiver_i
        self.total_residual=total_residual
        self.E_residual=E_residual
        self.W_residual=W_residual
        
    def get_site_resid(self):
        from numpy import mean
        
        W_residual=self.W_residual
        site_resid=mean(W_residual)
        
        self.site_resid=site_resid

    def plot_site_WE(self,sta_axis,colors,axis_lims,xlabel_toggle,ylabel_toggle):
        '''
        Plot all WE residuals for this site
        Input:
            sta_axis:       Axes to plot on. Can use this to plot iteratively
            colors:         RGB values to plot arr([r,g,b])
            axlims:         Array with axis limits: [[xmin, xmax],[ymin,ymax]]
            xlabel_toggle:  String with instructions to add x label: 'on'/'off'
            ylabel_toggle:  String with instructions to add y label: 'on'/'off'
        '''
        
        import matplotlib.pyplot
        from numpy import array,str,round
        
        #Station name?
        sta_name=self.sta
        #Station stats? (Site term)
        site_term=round(self.site_resid,2)
        
        #Get the plotting title:
        plttitle1=str(sta_name)
        plttitle2=str(site_term)
        
        #Plot WE residuals:
        sta_axis.scatter(self.mw,self.W_residual,edgecolors=colors,facecolors='none')
        
        #Get and plot the site term as a dashed line:
        site_line_x=array([axis_lims[0][0],axis_lims[0][1]])
        site_line_y=array([site_term,site_term])
        #plot, with dashed line and mid gray color:
        sta_axis.hold(True)
        sta_axis.plot(site_line_x,site_line_y,linestyle='--',color='0.75')
        
        #Set axis limits:
        sta_axis.set_xlim((axis_lims[0][0],axis_lims[0][1]))
        sta_axis.set_ylim((axis_lims[1][0],axis_lims[1][1]))
        
        #Set labels:
        #Set site name, etc.
        sta_axis.text((axis_lims[0][1]-0.8),(axis_lims[1][1]-1.3),plttitle1)
        sta_axis.text((axis_lims[0][1]-0.8),(axis_lims[1][1]-2.5),plttitle2)
       
        if xlabel_toggle=='on':
            sta_axis.set_xlabel('Moment Magnitude')
        
        if ylabel_toggle=='on':
            sta_axis.set_ylabel('ln residual')
        
        
        
        
        
#####
class residuals:
    '''
    Save database info plus residuals into one object for analysis
    '''
    
    def __init__(self,dbpath,event_list_path,station_list_path):
        '''
        Initialize database - pull necessary information and save to the object
        Input:
            dbpath:             Path to the database
            event_list_path:    Path to hte object holding the list of event objects
            station_list_path:  Path to the object holding the list of station objects
        Output:
            residual:           Object holding all data and residuals for a database
        '''
        
        import cPickle as pickle
        import dread
        from numpy import where
        
        #Load in database object:
        dname=open(dbpath,'r')
        db=pickle.load(dname)
        dname.close()
        
        #First, save the database info per recording:
        self.evnum=db.evnum
        self.elat=db.elat
        self.elon=db.elon
        self.edepth=db.edepth
        self.sta=db.sta
        self.stnum=db.stnum
        self.ml=db.ml
        self.mw=db.mw
        self.pga=db.pga
        self.pgv=db.pgv
        self.pga_pg=db.pga_pg
        self.r=db.r
        self.vs30=db.vs30
        self.ffdf=db.ffdf
        self.md_ffdf=db.md_ffdf
        self.stlat=db.stlat
        self.stlon=db.stlon
        self.stelv=db.stelv
        self.source_i=db.source_i
        self.receiver_i=db.receiver_i
        
        ###
        #Load in list of event objects:
        eobjs=dread.read_obj_list(event_list_path)
        
        #Load in list of station objects:
        sobjs=dread.read_obj_list(station_list_path)
        
        #Initialize the W residual, site term, and path term arrays: 
        total_residual=[]
        E_residual=[]
        E_std=[]
        W_residual=[]
        W_mean=[]
        W_std=[]
        site_terms=[]
        path_terms=[]
        
        #Loop through recordings to extract residuals...
        for record_i in range(len(self.evnum)):
            #Get the recorded pga, event, and station info:
            record_evnum_i=self.evnum[record_i]
            record_stnum_i=self.stnum[record_i]
            
            #Find the event object and station object that corresponds to this 
            #recording, and save the corresponding residuals...
            for event_i in range(len(eobjs)):
                #Get the event object for this event_i index: 
                event=eobjs[event_i]
                
                #Get the information from this event: evnumber, station list,mw,
                #event residuals, etc.
                evnum_i=event.evnum[0]
                
                #If this event is the same as the recording in question, continue:
                if evnum_i==record_evnum_i:
                    event_total_residual_i=event.total_residual
                    event_E_i=event.E_residual
                    event_Estd_i=event.E_std
                    event_Wmean_i=event.W_mean
                    event_Wstd_i=event.W_std
                    
                    #Save the values that correspond to this recording, which will
                    #be stored in the residuals object:
                    record_total_residual_i=event_total_residual_i
                    record_E_i=event_E_i
                    record_Estd_i=event_Estd_i
                    record_Wmean_i=event_Wmean_i
                    record_Wstd_i=event_Wstd_i
                    
                    #Append to the event term and std lists for the object:
                    total_residual.append(record_total_residual_i)
                    E_residual.append(record_E_i)
                    E_std.append(record_Estd_i)
                    W_mean.append(record_Wmean_i)
                    W_std.append(record_Wstd_i)
            
                    ########
                    #Get the station information from this event: within-event 
                    #residuals, site term, etc.:
                    for station_i in range(len(sobjs)):
                        #Get the station object for this station_i index:
                        station=sobjs[station_i]
                        
                        #Get the site number for this station:
                        station_stnum_i=station.stnum
                        
                        #Does this station correspond to the current recording?
                        #If so, store the information:
                        if station_stnum_i==record_stnum_i:
                        #Which event recorded in this station corresponds to the 
                        #current event?
                            station_evnum_ind=where(station.evnum==evnum_i)[0]
                            
                            #Take this index, and save the info from it:
                            record_W_i=station.W_residual[station_evnum_ind][0][0]
                            
                            #Also get the site term:
                            record_site_term_i=station.site_resid
                            
                            #Get the path term - it's the remainder of the within-event
                            #residual after removing the site term:
                            record_path_term_i=record_W_i-record_site_term_i
                            
                            #Save these to the recording...
                            W_residual.append(record_W_i)
                            site_terms.append(record_site_term_i)
                            path_terms.append(record_path_term_i)
                            
                            
                        #If the station doesn't match the recording, then carry on...    
                        else:
                            continue
                            
                #Close event loop, if this even tis not the same as the recording:
                else:
                    continue        
                    
            
        #Save new 
        self.total_residual=total_residual
        self.E_residual=E_residual
        self.E_std=E_std    
        self.W_residual=W_residual
        self.W_mean=W_mean
        self.W_std=W_std
        self.site_terms=site_terms
        self.path_terms=path_terms
        
    
##########
    def plot_path_term_mw(self,run_name,axlims):
        '''
        VJS 9/2016
        Plot the path terms vs. Mw
        Input:  
            run_name:           String with the database/inversion combo run name
            axlims:             Axis limits [[xmin,xmax],[ymin,ymax]]
        Output:
            fig1:               Figure with path terms 
        '''
        
        from matplotlib import pyplot as plt
        from matplotlib.ticker import MultipleLocator
        from numpy import str,array,ones,around,mean,std
        
        #Stats:
        mean_pterm=mean(self.path_terms)
        std_pterm=std(self.path_terms)
        
        #Color:
        rgb=array([111,168,163])/255.0
        rgb.astype(float)
        color=rgb*ones((len(self.mw),3))
        
        #definitions for axes
        fudge_factor=0.02
        left, width=0.1, 0.70
        bottom, height=0.1,0.81
        left_h=left+width+fudge_factor
        width_h=0.15
        
        #define axis limits for scatter, and histogram:
        rect_scatter=[left,bottom,width,height]
        rect_histy=[left_h,bottom,width_h,height]
        
        #define axis tick locations for histogram:
        hist_xLocator=MultipleLocator(500) 
        
        #Start figure:
        f1=plt.figure()
        
        #define axes:
        axScatter=plt.axes(rect_scatter)
        axHisty=plt.axes(rect_histy)
        
        #Scatter:
        axScatter.scatter(self.mw,self.path_terms,edgecolors=color,facecolors='none',lw=0.9)
        
        #Histogram:
        #want 4x as many bins as main plot y-axis limit units:
        nbins=(axlims[1][1]-axlims[1][0])*4
        
        #set the number of bins, adn the range to be the x axis limits (same as y axis, ln residuals):
        axHisty.hist(self.path_terms,bins=nbins,range=[axlims[1][0],axlims[1][1]],orientation='horizontal',color=rgb)
        
        #Also plot a dashed line at 0:
        axScatter.plot(axlims[0],[0,0],linestyle='--',color='0.75')
        
        #set axis limits:
        #scatter
        axScatter.set_xlim(axlims[0])
        axScatter.set_ylim(axlims[1])
        #histogram
        axHisty.set_ylim(axlims[1])
        #set axis ticks:
        axHisty.xaxis.set_major_locator(hist_xLocator)
        #set no labels on the y axis:
        axHisty.yaxis.set_ticklabels('')
        
        #Labels
        axScatter.set_xlabel(r"$\mathbf{M}$")
        axScatter.set_ylabel('ln Residuals')
        axScatter.set_title('Path Terms for run '+run_name+'\n'+'Mean: '+str(around(mean_pterm,decimals=2))+' Std Dev: '+str(around(std_pterm,decimals=2)))

        #Show
        f1.show()        
        
        return f1

##########
    def plot_path_term_r(self,run_name,axlims):
        '''
        VJS 9/2016
        Plot the path terms vs. distances
        Input:  
            run_name:           String with the database/inversion combo run name
            axlims:             Axis limits [[xmin,xmax],[ymin,ymax]]
        Output:
            fig1:               Figure with path terms 
        '''
        
        from matplotlib import pyplot as plt
        from matplotlib.ticker import MultipleLocator
        from numpy import str,array,ones,around,mean,std
        
        #Stats:
        mean_pterm=mean(self.path_terms)
        std_pterm=std(self.path_terms)
        
        #Color:
        rgb=array([111,168,163])/255.0
        rgb.astype(float)
        color=rgb*ones((len(self.r),3))
        
        #definitions for axes
        fudge_factor=0.02
        left, width=0.1, 0.70
        bottom, height=0.1,0.81
        left_h=left+width+fudge_factor
        width_h=0.15
        
        #define axis limits for scatter, and histogram:
        rect_scatter=[left,bottom,width,height]
        rect_histy=[left_h,bottom,width_h,height]
        
        #define axis tick locations for histogram:
        hist_xLocator=MultipleLocator(500) 
        
        #Start figure:
        f1=plt.figure()
        
        #define axes:
        axScatter=plt.axes(rect_scatter)
        axHisty=plt.axes(rect_histy)
        
        #Scatter:
        axScatter.scatter(self.r,self.path_terms,edgecolors=color,facecolors='none',lw=0.9)
        
        #Histogram:
        #want 4x as many bins as main plot y-axis limit units:
        nbins=(axlims[1][1]-axlims[1][0])*4
        
        #set the number of bins, adn the range to be the x axis limits (same as y axis, ln residuals):
        axHisty.hist(self.path_terms,bins=nbins,range=[axlims[1][0],axlims[1][1]],orientation='horizontal',color=rgb)
        
        #Also plot a dashed line at 0:
        axScatter.plot(axlims[0],[0,0],linestyle='--',color='0.75')
        
        #set axis limits:
        #scatter
        axScatter.set_xlim(axlims[0])
        axScatter.set_ylim(axlims[1])
        #histogram
        axHisty.set_ylim(axlims[1])
        #set axis ticks:
        axHisty.xaxis.set_major_locator(hist_xLocator)
        #set no labels on the y axis:
        axHisty.yaxis.set_ticklabels('')
        
        #Labels
        axScatter.set_xlabel('Distance (km)')
        axScatter.set_ylabel('ln Residuals')
        axScatter.set_title('Path Terms for run '+run_name+'\n'+'Mean: '+str(around(mean_pterm,decimals=2))+' Std Dev: '+str(around(std_pterm,decimals=2)))

        #Show
        f1.show()        
        
        return f1
            
    
    ######
    def add_vp_paths(self,ray_depth,ray_lat,ray_lon):
        '''
        Add Vp raypath locations to the residuals object
        Input:
            ray_depth:      List of arrays with the depth in km for each vp ray
            ray_lat:        List of arrays with the lat in deg for each vp ray
            ray_lon:        List of arrays with the lon in deg for each vp ray
        '''
        
        self.vp_depth=ray_depth
        self.vp_lat=ray_lat
        self.vp_lon=ray_lon
        
    ######
    def add_vs_paths(self,ray_depth,ray_lat,ray_lon):
        '''
        Add Vs raypath locations to the residuals object
        Input:
            ray_depth:      List of arrays with the depth in km for each vs ray
            ray_lat:        List of arrays with the lat in deg for each vs ray
            ray_lon:        List of arrays with the lon in deg for each vs ray
        '''
        
        self.vs_depth=ray_depth
        self.vs_lat=ray_lat
        self.vs_lon=ray_lon
                
    #######
    def plot_raypaths(self,veltype,view,axlims,stations,events,by_path,mymap,faultfile):
        '''
        Plot the path terms
        Input:
            veltype:                 Velocity type to plot (vp/vs = 1/2)
            view:                    View of the plot (map=0, lat vs depth=1, lon vs depth=2)  
            axlims:                  Axis limits [[xmin, xmax], [ymin, ymax]]
            stations:                Plot stations on figure? no=0, yes=1
            events:                  Plot events on figure?  no=0, yes=1             
            by_path:                 Color black/by path term (0/1) 
            mymap:                   String with python colormap (i.e. 'jet')
            faultfile:               String to path of the pckl file with fault segments
        '''
    
        
        import matplotlib.pyplot as plt
        from matplotlib import ticker
        from numpy import zeros,unique,where,array,mean,std,c_,arange
        from matplotlib.collections import LineCollection
        import matplotlib.colors as colors
        import matplotlib.cm as cm
        import dread
        
        #Read fault file and store data - list of arrays, with each array being a segment of the fault (lon, lat):
        fault_segments=dread.read_obj_list(faultfile)
        
        #Which velocity data is being plotted, Vp or Vs?
        #Depending on what it is, specify the depth, lat and lon separately
        if veltype==1:
            depth=self.vp_depth
            lat=self.vp_lat
            lon=self.vp_lon
            
            #Plot titles:
            ptitle='Plot of raypaths for Vp'
            
        elif veltype==2:
            depth=self.vs_depth
            lat=self.vs_lat
            lon=self.vs_lon
            
            #Plot titles:
            ptitle='Plot of raypaths for Vs'
        
        #Get color for plotting:
        if by_path==0:
            pcolor='k'
        elif by_path==1:
            #Get mean and std of dataset:
            mean_pterm=mean(self.path_terms)
            std_pterm=std(self.path_terms)
            #Set hte colorscale to cover 97% of the data:
            cmin=-3*std_pterm
            cmax=3*std_pterm 
        
        #Get unique event indices for plotting events:
        unique_events=unique(self.evnum)
        #Zero out the lat, lon, and depth arrays:
        uedepth=[]
        uelat=[]
        uelon=[]
        
        for event_c in range(len(unique_events)):
            #Get the event number for each event:
            evnum_i=unique_events[event_c]
            #Get the index of the first occurrence of this event:
            unique_event_ind=where(self.evnum==evnum_i)[0][0]
            #Pull out the info from here, as it should all be the same for all instances:
            uedepth_i=self.edepth[unique_event_ind]
            uelat_i=self.elat[unique_event_ind]
            uelon_i=self.elon[unique_event_ind]
            
            #Append to arrays:
            uedepth.append(uedepth_i)
            uelat.append(uelat_i)
            uelon.append(uelon_i) 
        
        #Make them arrays:
        uedepth=array(uedepth)
        uelat=array(uelat)
        uelon=array(uelon)  
            
        
        #Define the x and y to plot based on the view:
        #Map view:   
        if view==0:
            x=lon
            y=lat
            #Stations:
            stx=self.stlon
            sty=self.stlat
            #Events:
            evx=uelon
            evy=uelat
            
            #Labels:
            xlab='Longitude (degrees)'
            ylab='Latitude (degrees)'
            
        #cross section with latitude and depth:
        elif view==1:
            x=lat
            y=depth
            #Stations:
            stx=self.stlat
            sty=self.stelv
            #Events:
            evx=uelat
            evy=-1*uedepth
            
            #Labels:
            xlab='Latitude (degrees)'
            ylab='Depth (km)'
            
        #cross section with longitude and depth
        elif view==2:
            x=lon
            y=depth
            #Stations:
            stx=self.stlon
            sty=self.stelv
            #Events:
            evx=uelon
            evy=-1*uedepth
            
            #Labels:
            xlab='Longitude (deg)'
            ylab='Depth (km)'
          
        
        ##Plot:
        #Get colormap
        #Make colormap:
        colormap_pterm=plt.get_cmap(mymap)
        #Make a normalized colorscale
        cNorm=colors.Normalize(vmin=cmin, vmax=cmax)
        #Apply normalization to colormap:
        scalarMap=cm.ScalarMappable(norm=cNorm, cmap=colormap_pterm)
        
        #Make a fake contour plot for the colorbar:
        Z=[[0,0],[0,0]]
        levels=arange(cmin,cmax,0.01)
        c=plt.contourf(Z, levels, cmap=colormap_pterm)
 
        
        #Initiate plot
        figure=plt.figure()
        #Set axis format:
        x_formatter=ticker.ScalarFormatter(useOffset=False)
        
        #Plot the dem
        
        #Plot the raypaths 
        for path_i in range(len(depth)):
            #Assign color to path term:
            colorVal = scalarMap.to_rgba(self.path_terms[path_i])
            #Get x and y
            x_i=x[path_i]
            y_i=y[path_i]
            
            plt.plot(x_i,y_i,color=colorVal)
            
        #Add colorbar:
        cb=plt.colorbar(c)
        cb.set_label('Path term (ln residual)')
            
        
        #If stations are to be plotted:    
        if stations==1:
            #Hold on:
            #plt.hold(True)
            #Scatter:
            plt.scatter(stx,sty,color='black',s=100,marker='^',zorder=len(self.mw)+5)
            
        if events==1:
            #Hold on
            #plt.hold(True)
            #Scatter events:
            plt.scatter(evx,evy,edgecolors='g',facecolors='none',s=15,linewidths=2,zorder=len(self.mw)+7)
            
        #Plot faults, if it's map view:
        if view==0:
            for segment_i in range(len(fault_segments)):
                fault=fault_segments[segment_i]
                plt.plot(fault[:,0],fault[:,1],color='k',zorder=len(self.mw)+9)
            
                #Axis limits:
                plt.xlim(axlims[0])
                plt.ylim(axlims[1])
        
        #Axis labels, etc.:
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(ptitle)
        
        #Set format of axis:
        ax=plt.gca()
        ax.xaxis.set_major_formatter(x_formatter)
        
        #Show plot:
        plt.show()
        
        #Return :
        
        return figure
        
    def plot_raypaths_cutoffval(self,veltype,view,axlims,stations,events,mymap,faultfile,cutoff_val):
        '''
        Plot the path terms above a certain cutoff value.  Plot the rest as gray.
        Input:
            veltype:                 Velocity type to plot (vp/vs = 1/2)
            view:                    View of the plot (map=0, lat vs depth=1, lon vs depth=2)  
            axlims:                  Axis limits [[xmin, xmax], [ymin, ymax]]
            stations:                Plot stations on figure? no=0, yes=1
            events:                  Plot events on figure?  no=0, yes=1 
            mymap:                   String with python colormap (i.e. 'jet')
            fautlfile:               String to path of the pckl file with fault segments
            cutoff_val:              Cutoff value to plot (i.e., only plot path term if abs(path term) >= cutoff_val) 
        '''
    
        
        import matplotlib.pyplot as plt
        from matplotlib import ticker
        from numpy import zeros,unique,where,array,std,arange
        import matplotlib.colors as colors
        import matplotlib.cm as cm
        import dread
        
        #Read fault file and store data - list of arrays, with each array being a segment of the fault (lon, lat):
        fault_segments=dread.read_obj_list(faultfile)
        
        #Which velocity data is being plotted, Vp or Vs?
        #Depending on what it is, specify the depth, lat and lon separately
        if veltype==1:
            depth=self.vp_depth
            lat=self.vp_lat
            lon=self.vp_lon
            
            #Plot titles:
            ptitle='Plot of raypaths for Vp'
            
        elif veltype==2:
            depth=self.vs_depth
            lat=self.vs_lat
            lon=self.vs_lon
            
            #Plot titles:
            ptitle='Plot of raypaths for Vs'
        
        #Get color for plotting:
        #Get mean and std of dataset:
        std_pterm=std(self.path_terms)
        #Set hte colorscale to cover 97% of the data:
        cmin=-3*std_pterm
        cmax=3*std_pterm 
        
        #Get unique event indices for plotting events:
        unique_events=unique(self.evnum)
        #Zero out the lat, lon, and depth arrays:
        uedepth=[]
        uelat=[]
        uelon=[]
        
        for event_c in range(len(unique_events)):
            #Get the event number for each event:
            evnum_i=unique_events[event_c]
            #Get the index of the first occurrence of this event:
            unique_event_ind=where(self.evnum==evnum_i)[0][0]
            #Pull out the info from here, as it should all be the same for all instances:
            uedepth_i=self.edepth[unique_event_ind]
            uelat_i=self.elat[unique_event_ind]
            uelon_i=self.elon[unique_event_ind]
            
            #Append to arrays:
            uedepth.append(uedepth_i)
            uelat.append(uelat_i)
            uelon.append(uelon_i) 
        
        #Make them arrays:
        uedepth=array(uedepth)
        uelat=array(uelat)
        uelon=array(uelon)  
            
        
        #Define the x and y to plot based on the view:
        #Map view:   
        if view==0:
            x=lon
            y=lat
            #Stations:
            stx=self.stlon
            sty=self.stlat
            #Events:
            evx=uelon
            evy=uelat
            
            #Labels:
            xlab='Longitude (degrees)'
            ylab='Latitude (degrees)'
            
        #cross section with latitude and depth:
        elif view==1:
            x=lat
            y=depth
            #Stations:
            stx=self.stlat
            sty=self.stelv
            #Events:
            evx=uelat
            evy=-1*uedepth
            
            #Labels:
            xlab='Latitude (degrees)'
            ylab='Depth (km)'
            
        #cross section with longitude and depth
        elif view==2:
            x=lon
            y=depth
            #Stations:
            stx=self.stlon
            sty=self.stelv
            #Events:
            evx=uelon
            evy=-1*uedepth
            
            #Labels:
            xlab='Longitude (deg)'
            ylab='Depth (km)'
          
        
        ##Plot:
        #Get colormap
        #Make colormap:
        colormap_pterm=plt.get_cmap(mymap)
        #Make a normalized colorscale
        cNorm=colors.Normalize(vmin=cmin, vmax=cmax)
        #Apply normalization to colormap:
        scalarMap=cm.ScalarMappable(norm=cNorm, cmap=colormap_pterm)
        
        #Make a fake contour plot for the colorbar:
        Z=[[0,0],[0,0]]
        levels=arange(cmin,cmax,0.01)
        c=plt.contourf(Z, levels, cmap=colormap_pterm)
 
        
        #Initiate plot
        figure=plt.figure()
        #Set axis format:
        x_formatter=ticker.ScalarFormatter(useOffset=False)
        
        #Plot the raypaths 
        for path_i in range(len(depth)):
            #Assign color to path term:
            #If the absolute value of the path term is below the cutoff value, 
            #color it gray:
            if abs(self.path_terms[path_i])<cutoff_val:
                colorVal=scalarMap.to_rgba(self.path_terms[path_i])
                #Make the gray tuple rgb value, completely opaque (255 at end):
                colorVal=tuple(array([184,186,186,255])/255.)
                
                #Get x and y
                x_i=x[path_i]
                y_i=y[path_i]
            
                plt.plot(x_i,y_i,color=colorVal)
            
        #Plot the raypaths above the cutoff value: 
        for path_i in range(len(depth)):
            #Assign color to path term:
            #If the path term is above/below the cutoff value, color it based on
            #the colorscale made above:
            if abs(self.path_terms[path_i])>=cutoff_val:
                colorVal=scalarMap.to_rgba(self.path_terms[path_i])
                
                #Get x and y
                x_i=x[path_i]
                y_i=y[path_i]
            
                plt.plot(x_i,y_i,color=colorVal)
            
        #Add colorbar:
        cb=plt.colorbar(c)
        cb.set_label('Path term (ln residual)')
        
        #If stations are to be plotted:    
        if stations==1:
            #Hold on:
            #plt.hold(True)
            #Scatter:
            plt.scatter(stx,sty,color='black',s=100,marker='^',zorder=len(self.mw)+5)
            
        if events==1:
            #Hold on
            #plt.hold(True)
            #Scatter events:
            plt.scatter(evx,evy,edgecolors='g',facecolors='none',s=15,linewidths=2,zorder=len(self.mw)+7)
            
        #Plot faults, if it's map view:
        if view==0:
            for segment_i in range(len(fault_segments)):
                fault=fault_segments[segment_i]
                plt.plot(fault[:,0],fault[:,1],color='k',zorder=len(self.mw)+9)
            
        #Axis limits:
        plt.xlim(axlims[0])
        plt.ylim(axlims[1])
            
        #Axis labels, etc.:
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(ptitle)
        
        #Set format of axis:
        ax=plt.gca()
        ax.xaxis.set_major_formatter(x_formatter)
        
        #Show plot:
        plt.show()
        
        #Return :
        
        return figure
        
        
    ###############
    def plot_raypaths_3d(self,veltype,stations,events,axlims,mymap,faultfile):
        '''
        Plot the raypaths in 3d
        VJS 9/2016
        Input:
            veltype:        Velocity type (1/2, Vp/Vs)
            stations:       Plot stations?  0/1 = no/yes
            events?         Plot events?  0/1 = no/yes
            axlims:         Axis limits: [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
            mymap:          STring with the colormap to plot
            faultfile:      String with path to the faultfile to plot
        Output:
            figure:         Figure with 3D raypaths
        '''
        
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import ticker
        from numpy import zeros,unique,where,array,std,arange
        import matplotlib.colors as colors
        import matplotlib.cm as cm
        import dread
        
        #Read fault file and store data - list of arrays, with each array being a segment of the fault (lon, lat):
        fault_segments=dread.read_obj_list(faultfile)
        
        #Get data to plot based on specified velocity type:
        if veltype==1:
            #Then it's vp
            ray_x=self.vp_lon
            ray_y=self.vp_lat
            ray_z=self.vp_depth
        elif veltype==2:
            ray_x=self.vs_lon
            ray_y=self.vs_lat
            ray_z=self.vs_depth
        
        ##Plot:
        #get min and max for colormap:
        cmin=-3*std(self.path_terms)
        cmax=3*std(self.path_terms)
        
        #Get colormap
        #Make colormap:
        colormap_pterm=plt.get_cmap(mymap)
        #Make a normalized colorscale
        cNorm=colors.Normalize(vmin=cmin, vmax=cmax)
        #Apply normalization to colormap:
        scalarMap=cm.ScalarMappable(norm=cNorm, cmap=colormap_pterm)
        
        #Make a fake contour plot for the colorbar:
        Z=[[0,0],[0,0]]
        levels=arange(cmin,cmax,0.01)
        c=plt.contourf(Z, levels, cmap=colormap_pterm)
        
        #Initiate figure
        f3d=plt.figure()
        
        #and axes with projection:
        ax=f3d.gca(projection='3d')
        #Set axis format:
        x_formatter=ticker.ScalarFormatter(useOffset=False)
        
        #Plot raypaths in 3d:
        for ray_i in range(len(ray_x)):
            x_i=ray_x[ray_i]
            y_i=ray_y[ray_i]
            z_i=ray_z[ray_i]
                
            #get color to plot:
            colorVal=scalarMap.to_rgba(self.path_terms[ray_i])
            
            ax.plot(x_i,y_i,z_i,color=colorVal)
            
        #Plot stations?
        if stations==1:
            sta_x=self.stlon
            sta_y=self.stlat
            sta_z=zeros(len(self.stlat))
            
            #Plot:
            ax.scatter(sta_x,sta_y,sta_z,s=150,marker='^',color='k',zorder=(len(ray_x)+10))
            
        if events==1:
            ev_x=self.elon
            ev_y=self.elat
            ev_z=self.edepth
            
            #Make event depths negative:
            ev_z=-1*ev_z
            
            #Plot:
            ax.scatter(ev_x,ev_y,ev_z,s=20,edgecolors='g',facecolors='none',linewidths=1,zorder=(len(ray_x)+15))
        
        #Plot faults:
        for segment_i in range(len(fault_segments)):
            fault=fault_segments[segment_i]
            fault_z=zeros(len(fault))
            ax.plot(fault[:,0],fault[:,1],fault_z,color='k',zorder=len(self.mw)+17)
    
        #Add colorbar:
        cb=plt.colorbar(c)
        cb.set_label('Path term (ln residual)')
        
        #Set labels:
        #Set velocity handle:
        if veltype==1:
            vtype='Vp'
        elif veltype==2:
            vtype='Vs'
            
        ax.set_xlabel('Longitude (deg)')
        ax.set_ylabel('Latitude (deg)')
        ax.set_zlabel('Depth (km)')
        ax.set_title('Raypaths for '+vtype)
        
        #Set limits:
        ax.set_xlim(axlims[0])
        ax.set_ylim(axlims[1])
        ax.set_zlim(axlims[2])
        
        #Set format of axis:
        ax=plt.gca()
        ax.xaxis.set_major_formatter(x_formatter)        
             
        #Return:
        return ax