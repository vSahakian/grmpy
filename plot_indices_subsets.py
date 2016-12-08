## Plot indices for subsets of data
# VJS 12/2016

import res_analysis as ra
from os import path
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from numpy import array,arange
import matplotlib.colors as colors
import matplotlib.cm as cm

index=['ind_s_vs_gradpathint']
color_by=['r','mw']
axlims=[[[-3,6],[100,8650]],[[-3,6],[0.66,0.85]],[[-3,6],[0,6.5]]]
cmap='jet'
cvals=[[0,140],[1.2,3]]

#for indexi in range(len(index)):
#    for color_by_i in range(len(color_by)):
#        ra.plot_pathterms_colored(home,run_name,robj,index[indexi],axlims[indexi],color_by[color_by_i],cvals[color_by_i],cmap)

        
                        
        
#What to plot?
x=getattr(robj,term)
y=getattr(robj,index)
    
#get correlation coefficient:
pcoeff,tails=pearsonr(x,y)
pcoeff=round(pcoeff,2)

#Title:
ptitle='Plot of '+termtitle+' vs. '+indname+'\n Pearson coefficient: '+str(pcoeff)

#Get colormap
#Make colormap:
colormap_mwdist=plt.get_cmap(mymap)
#Make a normalized colorscale
cmin=cvals[0]
cmax=cvals[1]
cNorm=colors.Normalize(vmin=cmin, vmax=cmax)
#Apply normalization to colormap:
scalarMap=cm.ScalarMappable(norm=cNorm, cmap=colormap_mwdist)

#Make a fake contour plot for the colorbar:
Z=[[0,0],[0,0]]
levels=arange(cmin,cmax,0.01)
c=plt.contourf(Z, levels, cmap=colormap_mwdist)   

#Assign values to colormap
colorVal = scalarMap.to_rgba(getattr(robj,color_by))

#Plot:
f1=plt.figure()
plt.scatter(x,y,facecolors='none',edgecolors=colorVal,lw=0.5)

#Add colorbar:
cb=plt.colorbar(c)
if color_by=='r':
    cbarlabel='Distance (km)'
elif color_by=='mw':
    cbarlabel='M'
cb.set_label(cbarlabel)

#Set axis limits and labels:
plt.xlim(axlims[0])
plt.ylim(axlims[1])

plt.xlabel(termname)
plt.ylabel(indname)
plt.title(ptitle)

#Save figure:
pngname=fig_dir+index+'_'+term+'_'+color_by+'.png'
pdfname=pdf_dir+index+'_'+term+'_'+color_by+'.pdf'
plt.savefig(pngname)
plt.savefig(pdfname)