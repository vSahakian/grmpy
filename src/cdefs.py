# -*- coding: utf-8 -*-
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
        Initiate the class by giving database information.
        Input:
            event:          Array with event number (event)
            sta:            Array/list with station name (sta)
            N:              Array with station number (N)
            ml:             Array with local mag (ml)
            mw:             Array with moment mag (mw)
            DA:             Array with PGA in m/s/s
            DV:             Array with PGV in m/s
            r:              Array with Rrup, source to site distance (r)
            vs30:           Array with vs30 (in m/s)
            elat:           Array with event latitude
            elon:           Array with event longitude
            edepth:         Array with event depth (km), positive
            stlat:          Array with station latitude
            stlon:          Array with station longitude
            stelv:          Array with statin elevation (km), positive
            source_i:       Array with source index for raytracing
            receiver_i:     Array with receiver index for raytracing
        '''
        import numpy as np
        
        #"Fictitions depth" or "Finite fault dimension factor"
        c=4.5
        
        ##Save pga + pgv in m/s/s, not nm/s/s
        #DAm=DA*1e-9
        #DVm=DV*1e-9
        
        ##Get percent g:
        #pga_pg=DAm/9.81

        #Get percent g:
        pga_pg=DA/9.81 
                       
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
        self.pga=DA
        self.pgv=DV
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
        
        
    def plot_mpga(self,bmin,bmax,axlims):
        '''
        Plots log10 PGA, for various magnitude ranges specified by bmin, bmax,
        and step.
        Input:
            bmin:       Min value for bins
            bmax:       Max balue for bins
            axlims:     [[xmin,xmax],[ymin,ymax]]
        '''
        
        from matplotlib import pyplot as plt
        import numpy as np
        import matplotlib.colors as colors
        import matplotlib.cm as cm
                
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
        
        #Open figure:
        f1=plt.figure()
        #get colorvalue to plot
        colorVal=scalarMap.to_rgba(self.mw)
        
       #Plot...
        plt.scatter(self.r,np.log10(self.pga_pg),edgecolors=colorVal,facecolors='none',lw=0.5)

        #Add colorbar:
        cb=plt.colorbar(c)
        cb.set_label(r"$\mathbf{M}$")

        #Label the plot - Mbold on the x, log10PGA on the y, 
        plt.xlabel('Distance (km)')
        plt.ylabel(r"$\log_{10} PGA$")
        plt.title(r"PGA vs. Distance, binned by $\mathbf{M}$")
        
        #Axis labels:
        plt.xlim(axlims[0][0],axlims[0][1])
        plt.ylim(axlims[1][0],axlims[1][1])
        
        plt.show()
        
        return f1
        
    def plot_rpga(self,bmin,bmax,axlims):
        '''
        Plots log10 PGA, for various distance ranges specified by bmin, bmax,
        and step.
        Input:
            bmin:       Min value for colorscale
            bmax:       Max balue for colorscale
            axlims:     [[xmin,xmax],[ymin,ymax]]
        '''
        
        from matplotlib import pyplot as plt
        import numpy as np
        import matplotlib.colors as colors
        import matplotlib.cm as cm
                
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
        
        #Open figure:
        f1=plt.figure()
        #get colorvalue to plot
        colorVal=scalarMap.to_rgba(self.r)
        
       #Plot...
        plt.scatter(self.mw,np.log10(self.pga_pg),edgecolors=colorVal,facecolors='none',lw=0.5)
        
        #Add colorbar:
        cb=plt.colorbar(c)
        cb.set_label('Distance (km)')

        #Label the plot - Mbold on the x, log10PGA on the y, 
        plt.xlabel(r"$\mathbf{M}$")
        plt.ylabel(r"$\log_{10} PGA$")
        plt.title(r"PGA vs. $\mathbf{M}$, binned by distance")

        #Label the plot - Mbold on the x, log10PGA on the y, 
        plt.xlabel(r"$\mathbf{M}$")
        plt.ylabel(r"$\log_{10} PGA$")
        plt.title(r"PGA vs. $\mathbf{M}$, binned by distance")
        
        #Axis labels:
        plt.xlim(axlims[0][0],axlims[0][1])
        plt.ylim(axlims[1][0],axlims[1][1])
        
        plt.show()
        
        return f1
        
    def plot_rpga_withmodel(self,bmin,bmax,mw,d_ln,rng,sdist,ask_dist,axlims,VR,nga_mw=True,nga_pred=True,vref=True,annotate_mag=4,rotation_angle=35,predictive_parameter='pga'):
        '''
        Plots log10 PGA, for various distance ranges specified by bmin, bmax,
        and step.
        Input:
            bmin:           Min value for bins for data
            bmax:           Max balue for bins for data
            mw:             Mw array from gmpe.compute_model_fixeddist
            d_ln:           d array from compute_model_fixeddist - IN LN PGA or PGV!!
            rng:            Magnitude ranges, same array used for inversion
            sdist:          Distances array used for inversion
            ask_dist:       Distance at which ASK is plotted
            axlims:         Array with lims: [[xmin,xmax],[ymin,ymax]]
            VR:             Variance Reduction from inversion
            nga_mw:         mw range for NGA plotting, if provided
            nga_pred:       prediction array for NGA plotting, if provided
            vref:           Reference vs30 value (Default: 760 m/s)
            annotate_mag:   Magnitue at which to annotate.  Default is M 4.
            rotation_angle: Rotation angle for text.  Default is 35 degrees 
            predictive_parameter:   Predictive parameter: 'pga', or 'pgv'.  Default: 'pga'
        '''
        
        from matplotlib import pyplot as plt
        import matplotlib.colors as colors
        import matplotlib.cm as cm
        import numpy as np
        
        #Vs30 reference:
        if vref==None:
            vref=760
            
        # Convert d_ln to log10 for plotting with data:
        d = np.log10(np.exp(d_ln))
        
        # What's the predictive parameter to plot?  Whatever it is, convert to log10 for plotting:
        if predictive_parameter=='pga':
            log10predparam=np.log10(self.pga_pg)
        elif predictive_parameter=='pgv':
            log10predparam=np.log10(self.pgv)
        
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
        
        plt.scatter(self.mw,log10predparam,edgecolors=colorVal,facecolors='none',lw=0.5)
        
        #Add colorbar:
        cb=plt.colorbar(c)
        cb.set_label('Distance (km)')

        #Label the plot - Mbold on the x, log10PGA on the y, 
        plt.xlabel(r"$\mathbf{M}$")
        
        if predictive_parameter=='pga':
            plt.ylabel(r"$\log_{10} PGA$")
        elif predictive_parameter=='pgv':
            plt.ylabel(r"$\log_{10} PGV$")
            
        plt.title(r"PGA vs. $\mathbf{M}$, binned by distance" + "\n" + \
            "M Ranges: " + np.str(rng)+ ", Var Red="+np.str(np.around(VR,decimals=1)))
        
        ##Colors for GMPE:
        #colors_gmpe=plt.cm.rainbow(sdist.astype(float)/sdist.max())
        
        # Plot GMPE as dashed lines with annotations
        for j in range(len(sdist)):
            #Label for plot:
            #lab="R="+np.str(sdist[j])+"km"
            lab=np.str(sdist[j])+" km"

            
            mw_dist=mw[:,j]
            d_dist=d[:,j]
            
            #Plot
            #plt.plot(mw_dist,d_dist,color=colors_gmpe[j],linewidth=2,label=lab)
            plt.plot(mw_dist,d_dist,'--',c='#313332',linewidth=1.7)
            
            #Now annotate the lines.  
            #First find what is the value of PGA for my annotation magnitude, and this plotting distance:
            mag_where = np.argmin(abs(mw_dist - annotate_mag))
            
            #Make the x,y coordinates of where it annotates
            x_annotate = annotate_mag
            y_annotate = d_dist[mag_where]
            
            #Make annotation text:
            text_annotate = lab
            
            #plot it:
            #plt.annotate(text_annotate,xy=(x_annotate,y_annotate),rotation=rotation_angle)
            plt.text(x_annotate,y_annotate,text_annotate,rotation=rotation_angle,va='center',ha='center',backgroundcolor='w',fontsize=8,color='k')
        
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
            plt.plot(nga_mw,nga_pred,linestyle='-.',linewidth=2.5,color='#C0C0C0',label='ASK2014')
            
            #Now annotate the lines.  
            #First find what is the value of PGA for my annotation magnitude, and this plotting distance:
            # Add 0.5 to annotate_ma so its 0.5 M away from the normal line annotations and doesn't impinge on them
            mag_where = np.argmin(abs(nga_mw - (annotate_mag + 0.6)))
            
            #Make the x,y coordinates of where it annotates
            # Again , ad 0.5 to annotate_mag in x_annotate so the ASK label is not on top of the normal labels
            x_annotate = annotate_mag + 0.6
            y_annotate = nga_pred[mag_where]
            
            #Make annotation text:
            text_annotate = 'ASK %.0fkm' % ask_dist            
            
            # Plot it
            plt.text(x_annotate,y_annotate,text_annotate,rotation=rotation_angle-5,va='center',ha='center',fontsize=8,color='k')
            
            plt.show(f)            
            return f  
            
            
    def plot_histogram(self,axlims,nbins,type_flag):
        '''
        Plot a histogram of distance or magnitude for this database
        Input:
            axlims:         Axis limits for plot [[xmin,xmax],[ymin,ymax]] 
            nbins:          Number of bins for the histogram
            type_flag:      Plot distance or M?  0=Rrup/1=M
        Output:           
            f1:             Plot with histogram
        '''
        
        import matplotlib.pyplot as plt
        # Get quantity for histogram:        
        if type_flag==0:            
            histquantity = self.r            
            histname = 'Rrup'        
        elif type_flag==1:            
            histquantity = self.mw            
            histname = 'M'
        
        
        # Initiate figure:
        fhist = plt.figure()
        plt.hist(histquantity,nbins)
        
        plt.xlim(axlims[0])
        plt.ylim(axlims[1])
        
        plt.xlabel(histname)
        plt.ylabel('# Counts')
        plt.title('Histogram of %s' % histname)
        
        return fhist
        
class invinfo:
    '''
    Save paramters from an inversion.
    G:        Left hand side matrix
    d:        Data vector
    m:        Resulting model vector
    resid:    Residuals from inversion
    rank:     rank from inversion
    svals:    Singular values from inversion
    rng:      Magnitude ranges used in inversion
    sdist:    Distances used in smoothing for inversion
    smth:     Smoothing value used in inversion 
    stderror: Standard error if running mixed effects
    tvalue:   T value if running mixed effects
    '''
        
    def __init__(self,G,d,m,resid,L2norm,VR,rank,svals,rng,sdist,smth,stderror,tvalue):
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
        self.stderror=stderror
        self.tvalue=tvalue
            
            
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
        
    def add_E_resid(self,E_residual,E_stderr,E_std=float('NaN')):
        # E_std is the standard deviation of ALL events
        # E_stderr is the standard error/standard deviation for each event
        self.E_residual=E_residual
        self.E_std=E_std
        self.E_stderr=E_stderr
         
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
        from numpy import mean,std
        
        W_residual=self.W_residual
        site_resid=mean(W_residual)
        site_stderr=std(W_residual)
        
        self.site_resid=site_resid
        self.site_stderr=site_stderr
        
    def add_site_resid(self,site_term,site_stderr=float('NaN')):
        self.site_resid=site_term
        self.site_stderr=site_stderr
        
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
        site_line_0=array([0,0])
        #plot, with dashed line and mid gray color:
        sta_axis.hold(True)
        sta_axis.plot(site_line_x,site_line_y,linestyle='--',color='0.75')
        sta_axis.plot(site_line_x,site_line_0,linestyle='-',color='0.75')
        
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
    
    def __init__(self,dbpath,event_list_path,station_list_path,init_style='basic',
                    evnum=None,elat=None,elon=None,edepth=None,sta=None,stnum=None,ml=None,mw=None,
                    pga=None,pgv=None,pga_pg=None,r=None,vs30=None,ffdf=None,md_ffdf=None,stlat=None,
                    stlon=None,stelv=None,source_i=None,receiver_i=None,total_residual=None,E_residual=None,
                    E_mean=None,E_std=None,W_residual=None,W_mean=None,W_std=None,site_terms=None,site_mean=None,site_stderr=None,site_std=None,
                    path_terms=None,path_mean=None,path_std=None):
                   
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
        from numpy import where,mean,std,array,unique
        
        if init_style != 'basic':
            ##
            self.evnum=evnum
            self.elat=elat
            self.elon=elon
            self.edepth=edepth
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
            self.stlat=stlat
            self.stlon=stlon
            self.stelv=stelv
            self.source_i=source_i
            self.receiver_i=receiver_i
            
            self.total_residual=total_residual
            self.E_residual=E_residual
            self.E_mean=E_mean
            self.E_std=E_std   
                
            self.W_residual=W_residual
            self.W_mean=W_mean
            self.W_std=W_std
            
            self.site_terms=site_terms
            self.site_stderr=site_stderr
            self.site_mean=site_mean
            self.site_std=site_std
            
            self.path_terms=path_terms
            self.path_mean=path_mean
            self.path_std=path_std
            
                        
        # Otherwise, if it's used in the basic sense of making it for the first time, then do as follows...
        elif init_style == 'basic':
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
        E_stderr=[]
        W_residual=[]
        W_mean=[]
        W_std=[]
        site_terms=[]
        site_stderr=[]
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
                    event_E_i=event.E_residual
                    event_Estderr_i=event.E_stderr
                    event_Wmean_i=event.W_mean
                    event_Wstd_i=event.W_std
                    
                    #Save the values that correspond to this recording, which will
                    #be stored in the residuals object:
                    record_E_i=event_E_i
                    record_Estderr_i=event_Estderr_i
                    
                    record_Wmean_i=event_Wmean_i
                    record_Wstd_i=event_Wstd_i
                    
                    #Append to the event term and std lists for the object:
                    
                    E_residual.append(record_E_i)
                    E_stderr.append(record_Estderr_i)
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
                            record_total_residual_i=station.total_residual[station_evnum_ind][0]
                            record_W_i=station.W_residual[station_evnum_ind][0][0]
                            
                            #Also get the site term:
                            record_site_term_i=station.site_resid
                            record_site_stderr_i=station.site_stderr
                            
                            #Get the path term - it's the remainder of the within-event
                            #residual after removing the site term:
                            record_path_term_i=record_W_i-record_site_term_i
                            
                            #Save these to the recording...
                            total_residual.append(record_total_residual_i)
                            
                            W_residual.append(record_W_i)
                            
                            site_terms.append(record_site_term_i)
                            site_stderr.append(record_site_stderr_i)
                            
                            path_terms.append(record_path_term_i)
                            
                                                        
                        #If the station doesn't match the recording, then carry on...    
                        else:
                            continue
                            
                #Close event loop, if this even tis not the same as the recording:
                else:
                    continue        
                    
        total_residual = array(total_residual)
        
        #Save new 
        self.total_residual=total_residual
    
        self.E_residual=E_residual
        self.E_stderr=E_stderr   
        
        self.W_residual=W_residual
        self.W_mean=W_mean
        self.W_std=W_std
        
        self.site_terms=site_terms
        self.site_stderr=site_stderr
        
        self.path_terms=path_terms
        self.path_mean=mean(path_terms)
        self.path_std=std(path_terms)
        
        # Get the mean and std for Event and patht erms - has to be for unique events/stations:
        unique_ev,uev_ind = unique(self.evnum,return_index=True)
        unique_sta,usta_ind = unique(self.sta,return_index=True)
        
        self.E_mean = mean(array(self.E_residual)[uev_ind])
        self.E_std = std(array(self.E_residual)[uev_ind])
        
        self.site_mean = mean(array(self.site_terms)[usta_ind])
        self.site_std = std(array(self.site_terms)[usta_ind])
            
            
            
        
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
    def add_material_values(self,value_list,value_flag,ray_type):
        '''
        Add information about the material values to a residuals object
        Input:
            value_list:         List of arrays of values of a particular material model.
                                    The length of each array in the list should be the same as
                                    the number of points along the correspoinding ray.
            value_flag:         Flag for what type of value is being added
                                    0=Vp, 1=Vs, 2=Vp/Vs, 3=Qp, 4=Qs, 5=Qp/Qs
            ray_type:           0=p-wave, 1=s-wave
        '''
        
        #Save for p-waves?
        if ray_type==0:
            print 'Adding values for p-rays into object'
            if value_flag==0:
                self.rayval_p_vp=value_list
            elif value_flag==1:
                self.rayval_p_vs=value_list
            elif value_flag==2:
                self.rayval_p_vpvs=value_list
            elif value_flag==3:
                self.rayval_p_qp=value_list
            elif value_flag==4:
                self.rayval_p_qs=value_list
            elif value_flag==5:
                self.rayval_p_qpqs=value_list
            
                    
                            
        #S- rays?
        elif ray_type==1:
            print 'Adding values for s-rays into object'
            if value_flag==0:
                self.rayval_s_vp=value_list
            elif value_flag==1:
                self.rayval_s_vs=value_list
            elif value_flag==2:
                self.rayval_s_vpvs=value_list
            elif value_flag==3:
                self.rayval_s_qp=value_list
            elif value_flag==4:
                self.rayval_s_qs=value_list
            elif value_flag==5:
                self.rayval_s_qpqs=value_list  
            
    #######
    def add_indices(self,indices,indextype,ray_type,value_flag):
        '''
        Add indices to the residuals object.
        Input:
            indices:            Array with index values for every ray
            indextype:          Flag for type of index: 
                                    0=path integral,1=normalized path integral, 
                                    2=gradient path integral
            ray_type:           0=P-wave, 1=S-wave
            value_flag:         Type of material model.  
                                    0=Vp, 1=Vs, 2=Vp/Vs, 3=Qp, 4=Qs, 5=Qp/Qs
        '''
        
        #P-rays?
        if ray_type==0:
            #material model?
            #Vp
            if value_flag==0:
                if indextype==0:
                    self.ind_p_vp_pathint=indices
                elif indextype==1:
                    self.ind_p_vp_normpathint=indices
                elif indextype==2:
                    self.ind_p_vp_gradpathint=indices
            #Vs       
            elif value_flag==1:
                if indextype==0:
                    self.ind_p_vs_pathint=indices
                elif indextype==1:
                    self.ind_p_vs_normpathint=indices
                elif indextype==2:
                    self.ind_p_vs_gradpathint=indices      
            #Vp/Vs        
            elif value_flag==2:
                if indextype==0:
                    self.ind_p_vpvs_pathint=indices
                elif indextype==1:
                    self.ind_p_vpvs_normpathint=indices
                elif indextype==2:
                    self.ind_p_vpvs_gradpathint=indices  
            #Qp        
            elif value_flag==3:
                if indextype==0:
                    self.ind_p_qp_pathint=indices
                elif indextype==1:
                    self.ind_p_qp_normpathint=indices
                elif indextype==2:
                    self.ind_p_qp_gradpathint=indices         
            #Qs        
            elif value_flag==4:
                if indextype==0:
                    self.ind_p_qs_pathint=indices
                elif indextype==1:
                    self.ind_p_qs_normpathint=indices
                elif indextype==2:
                    self.ind_p_qs_gradpathint=indices   
            #Qp /Qs       
            elif value_flag==5:
                if indextype==0:
                    self.ind_p_qpqs_pathint=indices
                elif indextype==1:
                    self.ind_p_qpqs_normpathint=indices
                elif indextype==2:
                    self.ind_p_qpqs_gradpathint=indices      
                
        #S-rays?
        if ray_type==1:
            #material model?
            #Vp
            if value_flag==0:
                if indextype==0:
                    self.ind_s_vp_pathint=indices
                elif indextype==1:
                    self.ind_s_vp_normpathint=indices
                elif indextype==2:
                    self.ind_s_vp_gradpathint=indices
            #Vs       
            elif value_flag==1:
                if indextype==0:
                    self.ind_s_vs_pathint=indices
                elif indextype==1:
                    self.ind_s_vs_normpathint=indices
                elif indextype==2:
                    self.ind_s_vs_gradpathint=indices      
            #Vp/Vs        
            elif value_flag==2:
                if indextype==0:
                    self.ind_s_vpvs_pathint=indices
                elif indextype==1:
                    self.ind_s_vpvs_normpathint=indices
                elif indextype==2:
                    self.ind_s_vpvs_gradpathint=indices  
            #Qp        
            elif value_flag==3:
                if indextype==0:
                    self.ind_s_qp_pathint=indices
                elif indextype==1:
                    self.ind_s_qp_normpathint=indices
                elif indextype==2:
                    self.ind_s_qp_gradpathint=indices         
            #Qs        
            elif value_flag==4:
                if indextype==0:
                    self.ind_s_qs_pathint=indices
                elif indextype==1:
                    self.ind_s_qs_normpathint=indices
                elif indextype==2:
                    self.ind_s_qs_gradpathint=indices   
            #Qp /Qs       
            elif value_flag==5:
                if indextype==0:
                    self.ind_s_qpqs_pathint=indices
                elif indextype==1:
                    self.ind_s_qpqs_normpathint=indices
                elif indextype==2:
                    self.ind_s_qpqs_gradpathint=indices        
                    
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
            
            print 'Raypaths plotted'
            
                    
            #If stations are to be plotted:             
            if events==1:
                #Hold on
                #plt.hold(True)
                #Scatter events:
                plt.scatter(evx,evy,edgecolors='g',facecolors='none',s=15,linewidths=1.5,zorder=len(self.mw)+5)
                
                print 'Events plotted'
                
            if stations==1:
                #Hold on:
                #plt.hold(True)
                #Scatter:
                plt.scatter(stx,sty,color='black',s=100,marker='^',zorder=len(self.mw)+7)   
                
                print 'Stations plotted'
                
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
        
        print 'Raypaths plotted'
        
        #If stations are to be plotted:    
        if events==1:
            #Hold on
            #plt.hold(True)
            #Scatter events:
            plt.scatter(evx,evy,edgecolors='g',facecolors='none',s=15,linewidths=1.5,zorder=len(self.mw)+5)
            
            print 'Events plotted'
            
        if stations==1:
            #Hold on:
            #plt.hold(True)
            #Scatter:
            plt.scatter(stx,sty,color='black',s=100,marker='^',zorder=len(self.mw)+7)
            
            print 'Stations plotted'
            
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
            sta_z=self.stelv
            
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
        #Set velocity handle:ᰜ
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
        
                
                
###########
##########
## Mixed effects residuals class - different from above only in that the residuals are already computed, 
##   don't need event and station terms to make the initial object...

class mixed_residuals:
    '''
    Save database info plus residuals into one object for analysis
    '''
    
    def __init__(self,db,total_resid,tresidmean,tresidstd,evresid,ev_stderr,evmean,evstd,weresid,wemean,westd,siteresid,site_stderr,sitemean,sitestd,pathresid,pathmean,pathstd):
        '''
        Initialize database - pull necessary information and save to the object
        Input:
            db:                 Database object
            total_resid:        Array of total residuals
            tresidmean:         Mean of total residuals
            tresidstd:          Standard deviation of total residuals
            evresid:            Array of event terms
            ev_stderr:          Array of standard error per event term
            evmean:             Mean of event terms
            evstd:              Standard deviation of event terms
            weresid:            Array of within-event terms
            wemean:             Mean of within-event terms
            westd:              Standard deviation of within-event terms
            siteresid:          Array of site terms
            site_stderr:        Array of standard error per site term
            sitemean:           Mean of site terms
            sitestd:            Standard deviation of site terms
            pathresid:          Array of path terms
            pathmean:           Mean of path terms
            pathstd:            Standard deviation of path terms
        Output:
            residual:           Object holding all data and residuals for a database＜
        '''
        
                
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
        
        self.total_residual=total_resid
        
        self.E_residual=evresid
        self.E_stderr=ev_stderr
        self.E_mean=evmean
        self.E_std=evstd   
         
        self.W_residual=weresid
        self.W_mean=wemean
        self.W_std=westd
        
        self.site_terms = siteresid
        self.site_stderr = site_stderr
        self.site_mean = sitemean
        self.site_std = sitestd
        
        self.path_terms = pathresid
        self.path_mean = pathmean
        self.path_std = pathstd
        
            
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
        Add Vp raypath locations to the residuals object
        Input:
            ray_depth:      List of arrays with the depth in km for each vs ray
            ray_lat:        List of arrays with the lat in deg for each vs ray
            ray_lon:        List of arrays with the lon in deg for each vs ray
        '''
        
        self.vs_depth=ray_depth
        self.vs_lat=ray_lat
        self.vs_lon=ray_lon
    
    
    ############
    def add_material_values(self,value_list,value_flag,ray_type):
        '''
        Add information about the material values to a residuals object
        Input:
            value_list:         List of arrays of values of a particular material model.
                                The length of each array in the list should be the same as
                                the number of points along the correspoinding ray.
            value_flag:         Flag for what type of value is being added:
                                0=Vp, 1=Vs, 2=Vp/Vs, 3=Qp, 4=Qs, 5=Qp/Qs
            ray_type:           0=p-wave, 1=s-wave
        '''
        
        #Save for p-waves?
        if ray_type==0:
            print 'Adding values for p-rays into object'
            if value_flag==0:
                self.rayval_p_vp=value_list
            elif value_flag==1:
                self.rayval_p_vs=value_list
            elif value_flag==2:
                self.rayval_p_vpvs=value_list
            elif value_flag==3:
                self.rayval_p_qp=value_list
            elif value_flag==4:
                self.rayval_p_qs=value_list
            elif value_flag==5:
                self.rayval_p_qpqs=value_list            
        
        #S- rays?
        elif ray_type==1:
            print 'Adding values for s-rays into object'
            if value_flag==0:
                self.rayval_s_vp=value_list
            elif value_flag==1:
                self.rayval_s_vs=value_list
            elif value_flag==2:
                self.rayval_s_vpvs=value_list
            elif value_flag==3:
                self.rayval_s_qp=value_list
            elif value_flag==4:
                self.rayval_s_qs=value_list
            elif value_flag==5:
                self.rayval_s_qpqs=value_list
    
    
             
       
        
    #######
    def add_indices(self,indices,indextype,ray_type,value_flag):
        '''
        Add indices to the residuals object.
        Input:
            indices:            Array with index values for every ray
            indextype:          Flag for type of index: 
                                    0=path integral,1=normalized path integral, 
                                    2=gradient path integral
            ray_type:           0=P-wave, 1=S-wave
            value_flag:         Type of material model.  
                                    0=Vp, 1=Vs, 2=Vp/Vs, 3=Qp, 4=Qs, 5=Qp/Qs
        '''
        
        #P-rays?
        if ray_type==0:
            #material model?
            #Vp
            if value_flag==0:
                if indextype==0:
                    self.ind_p_vp_pathint=indices
                elif indextype==1:
                    self.ind_p_vp_normpathint=indices
                elif indextype==2:
                    self.ind_p_vp_gradpathint=indices
             #Vs       
            elif value_flag==1:
                if indextype==0:
                    self.ind_p_vs_pathint=indices
                elif indextype==1:
                    self.ind_p_vs_normpathint=indices
                elif indextype==2:
                    self.ind_p_vs_gradpathint=indices      
            #Vp/Vs        
            elif value_flag==2:
                if indextype==0:
                    self.ind_p_vpvs_pathint=indices
                elif indextype==1:
                    self.ind_p_vpvs_normpathint=indices
                elif indextype==2:
                    self.ind_p_vpvs_gradpathint=indices  
            #Qp        
            elif value_flag==3:
                if indextype==0:
                    self.ind_p_qp_pathint=indices
                elif indextype==1:
                    self.ind_p_qp_normpathint=indices
                elif indextype==2:
                    self.ind_p_qp_gradpathint=indices         
            #Qs        
            elif value_flag==4:
                if indextype==0:
                    self.ind_p_qs_pathint=indices
                elif indextype==1:
                    self.ind_p_qs_normpathint=indices
                elif indextype==2:
                    self.ind_p_qs_gradpathint=indices   
            #Qp /Qs       
            elif value_flag==5:
                if indextype==0:
                    self.ind_p_qpqs_pathint=indices
                elif indextype==1:
                    self.ind_p_qpqs_normpathint=indices
                elif indextype==2:
                    self.ind_p_qpqs_gradpathint=indices      
                    
        #S-rays?
        if ray_type==1:
            #material model?
            #Vp
            if value_flag==0:
                if indextype==0:
                    self.ind_s_vp_pathint=indices
                elif indextype==1:
                    self.ind_s_vp_normpathint=indices
                elif indextype==2:
                    self.ind_s_vp_gradpathint=indices
             #Vs       
            elif value_flag==1:
                if indextype==0:
                    self.ind_s_vs_pathint=indices
                elif indextype==1:
                    self.ind_s_vs_normpathint=indices
                elif indextype==2:
                    self.ind_s_vs_gradpathint=indices      
            #Vp/Vs        
            elif value_flag==2:
                if indextype==0:
                    self.ind_s_vpvs_pathint=indices
                elif indextype==1:
                    self.ind_s_vpvs_normpathint=indices
                elif indextype==2:
                    self.ind_s_vpvs_gradpathint=indices  
            #Qp        
            elif value_flag==3:
                if indextype==0:
                    self.ind_s_qp_pathint=indices
                elif indextype==1:
                    self.ind_s_qp_normpathint=indices
                elif indextype==2:
                    self.ind_s_qp_gradpathint=indices         
            #Qs        
            elif value_flag==4:
                if indextype==0:
                    self.ind_s_qs_pathint=indices
                elif indextype==1:
                    self.ind_s_qs_normpathint=indices
                elif indextype==2:
                    self.ind_s_qs_gradpathint=indices   
            #Qp /Qs       
            elif value_flag==5:
                if indextype==0:
                    self.ind_s_qpqs_pathint=indices
                elif indextype==1:
                    self.ind_s_qpqs_normpathint=indices
                elif indextype==2:
                    self.ind_s_qpqs_gradpathint=indices
    
    
    
    
    
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
            sta_z=self.stelv
            
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




###########
##########

##Material Model object:

class material_model:
    '''
    A class for storing and manipulating material models 
    '''
    
    def __init__(self,x,y,z,nx,ny,nz,material_model):
        '''
        Make a material model object to store the information relating to 
        a velocity or attenuation model, etc.
        Input:
            x:                  Array with values of x nodes
            y:                  Array with values of y nodes
            z:                  Array with values of z nodes
            nx:                 Number of x nodes
            ny:                 Number of y nodes
            nz:                 Number of z nodes
            material_model:     Three-dimensional array with dimensions being (z,x,y)
        Output:
            material_object:    An object with the above info
        '''
        
        self.x=x
        self.y=y
        self.z=z
        self.nx=nx
        self.ny=ny
        self.nz=nz
        self.materials=material_model
        
    def plot_zslice(self,z_val,colormap,climits,xlab,ylab,modellab):
        '''
        Plot a z slice
        Input:
            zval:           Depth at which to plot
            colormap:       String with colormap to use
            climits:        Array with limits for color: [cmin, cmax]
            xlab:           String with xlabel
            ylab:           String with ylabel
            modellab:       String with label for the colorbar (i.e., 'Vs')
        '''
        
        import matplotlib.pyplot as plt
        from numpy import argmin
        
        #Find z value in model closest to the input slice depth:
        z_dist=abs(self.z-z_val)
        min_zdist_i=argmin(z_dist)
        
        #Print which distance:
        print 'Closest z node to requested value is '+str(self.z[min_zdist_i])
        
        #Get array to plot:
        slice_array=self.materials[min_zdist_i]
        #Get X and Y to plot:
        X=self.x
        Y=self.y
        
        #Initiate plot:
        #f1=plt.figure()
        #plt.pcolormesh(X,Y,slice_array,cmap=colormap,vmin=climits[0],vmax=climits[1])
        #plt.axis([X.min(),X.max(),Y.min(),Y.max()])
        #plt.colorbar()
        
        f1=plt.figure()
        plt.imshow(slice_array,cmap=colormap,vmin=climits[0],vmax=climits[1],extent=[X.min(),X.max(),Y.min(),Y.max()],interpolation='spline36',origin='lower')
        cbar=plt.colorbar()
        cbar.set_label(modellab+' (km/s)')
    
        ###Title:
        #depth being plotted:
        z_plot=self.z[min_zdist_i]
        ptitle='Depth slice at '+str(z_plot)+' km'
        
        plt.title(ptitle)
        plt.xlabel(xlab)
        plt.ylabel(ylab)

        return f1
        
    def plot_yslice(self,y_val,colormap,climits,xlab,zlab,modellab,aspectr):
        '''
        Plot a z slice
        Input:
            yval:           y node at which to plot
            colormap:       String with colormap to use
            climits:        Array with limits for color: [cmin, cmax]
            xlab:           String with xlabel
            zlab:           String with zlabel
            modellab:       String with label for the colorbar (i.e., 'Vs')
            aspectr:        Aspect ratio for plot
        '''
        
        import matplotlib.pyplot as plt
        from numpy import argmin
        
        #Find z value in model closest to the input slice depth:
        y_dist=abs(self.y-y_val)
        min_ydist_i=argmin(y_dist)
        
        #Print which distance:
        print 'Closest y node to requested value is '+str(self.y[min_ydist_i])
        
        #Get array to plot:
        slice_array=self.materials[:,:,min_ydist_i]
        #Get X and Y to plot:
        X=self.x
        Y=-1*self.z
        
        #Initiate plot:
        #f1=plt.figure()
        #plt.pcolormesh(X,Y,slice_array,cmap=colormap,vmin=climits[0],vmax=climits[1])
        #plt.axis([X.min(),X.max(),Y.min(),Y.max()])
        #plt.colorbar()
        
        f1=plt.figure()
        plt.imshow(slice_array,cmap=colormap,vmin=climits[0],vmax=climits[1],extent=[X.min(),X.max(),Y.min(),Y.max()],interpolation='spline36',origin='upper',aspect=aspectr)
        cbar=plt.colorbar()
        cbar.set_label(modellab+' (km/s)')
        
        ###Title:
        #depth being plotted:
        y_plot=self.y[min_ydist_i]
        ptitle='Latitude slice at '+str(y_plot)+' degrees'
        
        plt.title(ptitle)
        plt.xlabel(xlab)
        plt.ylabel(zlab)

        return f1
        
    def plot_xslice(self,x_val,colormap,climits,ylab,zlab,modellab,aspectr):
        '''
        Plot a z slice
        Input:
            xval:           x node at which to plot
            colormap:       String with colormap to use
            climits:        Array with limits for color: [cmin, cmax]
            ylab:           String with ylabel
            zlab:           String with zlabel
            modellab:       String with label for the colorbar (i.e., 'Vs')
            aspectr:        Aspect ratio for plot
        '''
        
        import matplotlib.pyplot as plt
        from numpy import argmin
        
        #Find z value in model closest to the input slice depth:
        x_dist=abs(self.x-x_val)
        min_xdist_i=argmin(x_dist)
        
        #Print which distance:
        print 'Closest x node to requested value is '+str(self.x[min_xdist_i])
        
        #Get array to plot:
        slice_array=self.materials[:,:,min_xdist_i]
        #Get X and Y to plot:
        X=self.y
        Y=-1*self.z
        
        #Initiate plot:
        #f1=plt.figure()
        #plt.pcolormesh(X,Y,slice_array,cmap=colormap,vmin=climits[0],vmax=climits[1])
        #plt.axis([X.min(),X.max(),Y.min(),Y.max()])
        #plt.colorbar()
        
        f1=plt.figure()
        plt.imshow(slice_array,cmap=colormap,vmin=climits[0],vmax=climits[1],extent=[X.min(),X.max(),Y.min(),Y.max()],interpolation='spline36',origin='upper',aspect=aspectr)
        cbar=plt.colorbar()
        cbar.set_label(modellab+' (km/s)')
        
        ###Title:
        #depth being plotted:
        x_plot=self.x[min_xdist_i]
        ptitle='Longitude slice at '+str(x_plot)+' degrees'
        
        plt.title(ptitle)
        plt.xlabel(ylab)
        plt.ylabel(zlab)

        return f1        
        
        
        
#########PLACEHOLDER########
##Object to hold the raypath path term grid
class pterm_3dgrid:
    def __init__(self,statistic,binedges,binnumber):
        '''
        Initiate the path term grid object.
        Input: 
            statistic:      Array (dims (nx, ny, nz)) with the path term statistic
            binedges:       Array with three arrays, containing bin edges (lon, lat, depth)
                            dims: (nx+1, ny+1, nz+1)
            binnumber:      Array with indices of bin number
        '''
        self.statistic=statistic
        self.binedges=binedges
        self.binnumber=binnumber
        
    def plot_slice(self,sliceaxis,slicecoord,coordtype,aspectr,cmap,climits,axlims):
        '''
        Plot a slice of the path term grid model.
        Input:
            sliceaxis:          Axes handle to plot on
            slicecoord:         Coordinate (lon, lat in degree, depth in km) of the slice to plot.
                                If depth, it must be positive.
            coordtype:          Type of coordinate: 'lon', 'lat', 'depth'
            aspectr:            Aspect ratio to plot
            cmap:               String with colormap, i.e., 'jet'
            climits:            Value for colorscale [cmin,cmax]
            axlims:             Axis limits for plot [[xmin,xmax],[ymin,ymax]]
        Output:
            pterm_ax:           Axis with path t            plt.show()erm grid slice plotted
        '''
        
        from numpy import argmin
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        
       
        #Get the slice from the 3D array that corresponds to this lat/lon:
        if coordtype=='lon':
            binind=argmin(abs(self.binedges[0]+slicecoord))
            #If it's the right hand side or last edge of tha tdimension, use the bin to the inside:
            if binind==len(self.binedges[0]):
                binind=binind=1
            #Get array to plot:
            statistic=self.statistic[binind,:,:]
            
            #Get axes and extent information:
            xmin=min(self.binedges[1])
            xmax=max(self.binedges[1])
            ymin=min(self.binedges[2])
            ymax=max(self.binedges[2])
            
            #Labels:
            xlabel='Latitude'
            ylabel='Depth (km)'
            ptitle='Path term slice at Longitude '+str(slicecoord)
            
            #Plot:
            sliceaxis.imshow(statistic.T,origin='lower',aspect=aspectr,extent=[xmin,xmax,ymin,ymax],interpolation='spline36',vmin=climits[0],vmax=climits[1])
            cbar=plt.colorbar()
            cbar.set_label('ln Residual')
            
        elif coordtype=='lat':
            binind=argmin(abs(self.binedges[1]-slicecoord))
            #If it's the right hand side or last edge of tha tdimension, use the bin to the inside:
            if binind==len(self.binedges[0]):
                binind=binind=1
            #Get array to plot:
            statistic=self.statistic[:,binind,:]
            
            #Get axes and extent information:
            xmin=min(self.binedges[0])
            xmax=max(self.binedges[0])
            ymin=min(self.binedges[2])
            ymax=max(self.binedges[2])
            
            xlabel='Longitude'
            ylabel='Depth (km)'
            ptitle='Path term slice at Latitude '+str(slicecoord)
            
            #Plot:
            sliceaxis.imshow(statistic.T,origin='lower',aspect=aspectr,extent=[xmin,xmax,ymin,ymax],interpolation='spline36',vmin=climits[0],vmax=climits[1])
            divider=make_axes_locatable(sliceaxis)
            caxis=divider.append_axes('right',size='25%',pad=0.05)
            cbar=plt.colorbar(sliceaxis,cax=caxis)
            cbar.set_label('ln Residual')     
                     
        elif coordtype=='depth':
            binind=argmin(abs(self.binedges[2]+slicecoord))
            #If it's the right hand side or last edge of tha tdimension, use the bin to the inside:
            if binind==len(self.binedges[0]):
                binind=binind=1
            #Get array to plot:
            statistic=self.statistic[:,:,binind]
            
            #Get axes and extent information:
            xmin=min(self.binedges[0])
            xmax=max(self.binedges[0])
            ymin=min(self.binedges[1])
            ymax=max(self.binedges[1])
            
            xlabel='Longitude'
            ylabel='Latitude'
            ptitle='Path term slice at Depth '+str(slicecoord)
            
            #Plot:
            sliceaxis.imshow(statistic.T,origin='lower',aspect=aspectr,extent=[xmin,xmax,ymin,ymax],interpolation='spline36',vmin=climits[0],vmax=climits[1])
            cbar=plt.colorbar()
            cbar.set_label('ln Residual')
            

        
        #Titles and limits:
        sliceaxis.set_xlabel(xlabel)
        sliceaxis.set_ylabel(ylabel)
        sliceaxis.set_title(ptitle)
        
        sliceaxis.set_xlim(axlims[0][0],axlims[0][1])
        sliceaxis.set_ylim(axlims[1][0],axlims[1][1])
        
        #Return:
        return sliceaxis
        
        
        
        
        
        