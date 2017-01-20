    
    #definitions for axes
    fudge_factor=0.02
    left, width=0.1, 0.70
    bottom, height=0.1,0.70
    left_h=left+width+fudge_factor
    bottom_h=bottom+height+fudge_factor
    width_h=0.15
    height_h=0.15
    
        
    #define axis limits for scatter, and histogram:
    rect_scatter=[left,bottom_h,width,height]
    rect_histy=[left_h,bottom_h,width_h,height]
    rect_histx=[left,bottom,width, height_h]
        
    #define axis tick locations for histogram:
    hist_xLocator=MultipleLocator(5000)
    hist_yLocator=MultipleLocator(5000) 
        
    #Start figure:
    eventfig = plt.figure()
        
    #define axes:
    axScatter=plt.axes(rect_scatter)
    axHisty=plt.axes(rect_histy)
    axHistx-plt.axes(rect_histx)
        
    # Scatter:
    axScatter.scatter(tradevent,mixedevent,marker='o',s=symbol_size[0],color=colorVal)
    axScatter.plot(xe,ye,color='gray',linewidth=1.5)        
        

        
    #Histograms:
    #want 4x as many bins as main plot y-axis limit units:
    nbins_y=(evaxlims[1][1]-evaxlims[1][0])*4
    nbins_x=(evaxlims[0][1]-evaxlims[0][0])*4
        
    #set the number of bins, adn the range to be the x axis limits (same as y axis, ln residuals):
    axHisty.hist(mixedevent,bins=nbins_y,range=[evaxlims[1][0],evaxlims[1][1]],orientation='horizontal',color=rgb)
        
    #set it for x:
    axHistx.hist(tradevent,bins=nbins_x,range=[evaxlims[0][0],evaxlims[0][1]],orientation='vertical',color=rgb)
        
        #set axis limits:
        #scatter
        axScatter.set_xlim(evaxlims[0])
        axScatter.set_ylim(evaxlims[1])
        #histogram
        axHisty.set_ylim(evaxlims[1])
        axHistx.set_xlim(evaxlims[0])
        #set axis ticks:
        axHisty.xaxis.set_major_locator(hist_yLocator)
        axHistx.yaxis.set_major_locator(hist_xLocator)
        #set no labels on the y axis:
        axHisty.yaxis.set_ticklabels('')
        axHistx.xaxis.set_ticklabels('')
        
        #Labels
        axScatter.set_xlabel(r"$\mathbf{M}$")
        axScatter.set_ylabel('ln Residuals')
        axScatter.set_title('Total Residuals for run '+run_name+'\n'+'Mean: '+str(around(self.mean_residual,decimals=2))+' Std Dev: '+str(around(self.std_dev,decimals=2)))
    

    
    
    
    

    
    # Colrobar:
    cb = axScatter.colorbar(c)
    cb.set_label('Number of stations per event')