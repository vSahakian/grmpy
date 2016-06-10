######Class definition######
#VJS 6/2016

'''
Generic classes that are useful for this package...
'''


class db:
    '''
    This class describes a set of events
    '''
    
    def __init__(self,event,sta,N,ml,mw,DA,DV,dist,vs30):
        '''
        Initiate the class by giving the event number (event), 
        station name (sta), station number (N), local mag (ml), moment mag (mw), 
        PGA (DA), PGV (DV), and source to site distance (dist)
        '''
        import numpy as np
        
        #"Fictitions depth" or "Finite fault dimension factor"
        c=4.5
        
        #Save pga + pgv in m/s/s, not nm/s/s
        DAm=DA*1e-9
        DVm=DV*1e-9
        
        #Get percent g:
        pga_pg=DAm/9.81
        
        #Give these values to the db:
        self.evnum=event
        self.sta=sta
        self.stnum=N
        self.ml=ml
        self.mw=mw
        self.pga=DAm
        self.pgv=DVm
        self.pga_pg=pga_pg
        self.r=dist
        self.vs30=vs30
        self.ffdf=np.sqrt(self.r**2 + c**2)
        
    def plot_apga(self):
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
        
    def plot_rpga_withmodel(self,bmin,bmax,step,m,rng,sdist,axlims,resid,vref=True):
        '''
        Plots log10 PGA, for various distance ranges specified by bmin, bmax,
        and step.
        Input:
            bmin:       Min value for bins
            bmax:       Max balue for bins
            step:       Step interval for bins
            m:          Model vector from inversion
            rng:        Magnitude ranges, same array used for inversion
            sdist:      Distances array used for inversion
            axlims:     Array with lims: [[xmin,xmax],[ymin,ymax]]
            resid:      Residual from inversion
            vref:       Reference vs30 value (Default: 760 m/s)
        '''
        
        from matplotlib import pyplot as plt
        import numpy as np
        
        #Vs30 reference:
        if vref==None:
            vref=760
            
        #Fictitious depth:
        c=4.5
        
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
        plt.title(r"PGA vs. $\mathbf{M}$, binned by distance; Residual="+np.str(resid))
        
        colors_gmpe=plt.cm.rainbow(sdist.astype(float)/sdist.max())
        
        for j in range(len(sdist)):
            ffdf=np.sqrt(sdist[j]**2 + c**2)
            
            #Label for plot:
            lab="R="+np.str(sdist[j])+"km"
            
            for i in range(len(rng)-1):
                #Get the magnitudes to plot against:
                mw=np.linspace(rng[i],rng[i+1],100)
                
                #Get the coefficients for this range:
                a1=m[i*5]
                a2=m[(i*5)+1]
                a3=m[(i*5)+2]
                a4=m[(i*5)+3]
                a5=m[(i*5)+4]
                
                #GMPE:
                d=a1+a2*mw + a3*(8.5-mw)**2 + a4*np.log(ffdf) + \
                    a5*sdist[j] 
                    # Don't add this yet...I think it should only go with the 
                    #data for residuals... 
                    #+ 0.6*np.log(self.vs30/vref)
                
                #Plot
                if i==0:
                    f=plt.plot(mw,d,color=colors_gmpe[j],label=lab)
                else:
                    f=plt.plot(mw,d,color=colors_gmpe[j],label=None)
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
    def __init__(self,G,d,m,resid,rank,svals,rng,sdist,smth):
        self.G=G
        self.d=d
        self.m=m
        self.resid=resid
        self.rank=rank
        self.svals=svals
        self.rng=rng
        self.sdist=sdist
        self.smth=smth
        
    
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


#class event:
#    '''
#    This class describes one event
#    '''
#    
#    def __init__(self,event,sta,N,ml,mw,DA,DV)