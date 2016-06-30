######Class definition######
#VJS 6/2016

'''
Generic classes that are useful for this package...
'''


class db:
    '''
    This class describes a set of events
    '''
    
    def __init__(self,event,sta,N,ml,mw,DA,DV,r,vs30,lat,lon,depth):
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
        self.lat=lat
        self.lon=lon
        self.depth=depth
        
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
        
    def plot_rpga_withmodel(self,bmin,bmax,step,mw,d,rng,sdist,axlims,VR,vref=True):
        '''
        Plots log10 PGA, for various distance ranges specified by bmin, bmax,
        and step.
        Input:
            bmin:       Min value for bins for data
            bmax:       Max balue for bins for data
            step:       Step interval for bins for data
            mw:         Mw array from gmpe.compute_model_fixeddist
            d:          d array from compute_model_fixeddist
            rng:        Magnitude ranges, same array used for inversion
            sdist:      Distances array used for inversion
            axlims:     Array with lims: [[xmin,xmax],[ymin,ymax]]
            VR:         Variance Reduction from inversion
            vref:       Reference vs30 value (Default: 760 m/s)
        '''
        
        from matplotlib import pyplot as plt
        import numpy as np
        
        #Vs30 reference:
        if vref==None:
            vref=760
            
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
        plt.title(r"PGA vs. $\mathbf{M}$, binned by distance" + "\n" + \
            "M Ranges: " + np.str(rng)+ ", Var Red="+np.str(np.around(VR,decimals=1)))
        
        colors_gmpe=plt.cm.rainbow(sdist.astype(float)/sdist.max())
        
        for j in range(len(sdist)):
            #Label for plot:
            lab="R="+np.str(sdist[j])+"km"
            
            mw_dist=mw[:,j]
            d_dist=d[:,j]
            
            #Plot
            f=plt.plot(mw_dist,d_dist,color=colors_gmpe[j],label=lab)

        #Add legend:
        plt.legend(loc=4)
            
        #Limits:
        xlims=axlims[0]
        ylims=axlims[1]
        plt.xlim(xlims)
        plt.ylim(ylims)
        
        plt.show(f)
        
        
        
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
        
        
class event:
    '''
    Save all data for one event, to use in residual computation
    '''
    
    def __init__(self,evnum,sta,stnum,ml,mw,pga,pgv,pga_pg,r,vs30,ffdf,md_ffdf,lat,lon,depth):
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
        self.lat=lat
        self.lon=lon
        self.depth=depth
        

        
    
#class hdr:
#    '''
#    This class describes the header info for a db (like AB's for now)
#    '''
#    
#    def __init__(self,dmin,dmax,mag,dist,az,year,day,hour,minu,sec,msec,nevid,idep):
#        '''
#        Has header information for events...
#        '''
#        
#        #Give header values:
#        self.dmin=dmin
#        self.dmax=dmax
#        self.mag=mag
#        self.r=dist
#        self.


