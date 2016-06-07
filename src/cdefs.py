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
        
        #Give these values to the db:
        self.evnum=event
        self.sta=sta
        self.stnum=N
        self.ml=ml
        self.mw=mw
        self.pga=DA
        self.pgv=DV
        self.r=dist
        self.vs30=vs30
        
    def plot_apga(self):
        '''
        Plots all PGA, regardless of M/r
        '''
        
        from matplotlib import pyplot as plt
        import numpy as np
        
        #Open figure
        plt.figure()
        
        #Plot...
        plt.scatter(self.mw,np.log10(self.pga),marker='o')
        plt.xlabel(r"$\mathbf{M}$")
        plt.ylabel(r"$\log_{10} PGA$")
        plt.title(r"PGA vs. $\mathbf{M}$ for all distances")
        
        plt.show()
        
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
#    def __init__(self,event,sta,N,ml,mw,DA,DV,