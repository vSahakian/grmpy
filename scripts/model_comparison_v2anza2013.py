# Compare event, site, and path terms between models
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

#Change this parameter depending on where you run:
#0=desktop
#1=mac

what_home=1

if what_home==0:
    #Desktop:
    HOME='/media/vsahakian/katmai'
    codehome='/home/vsahakian'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'
    codehome='/Users/vsahakian'

# Set homes:
home=HOME+'/anza/models/residuals/'

# Model 1 - Inversion of Anza data
#   Residuals object:
residpath_mod1=home+'v2anza2013_Mc8.5_pgrid_5sta/v2anza2013_Mc8.5_pgrid_5sta_robj.pckl'

# Model 2 - Mixed Effects inversion of Anza data
residpath_mod2=home+'mixedregr_v2anza2013_Mc_8.5/mixedregr_Mc8.5_VR_99.9_robj.pckl'

# Figure directory:
fig_dir=home+'v2anza2013_modelcompare/figs/'

# Plotting colormap:
mymap='gnuplot'

evaxlim = [[-6,6],[-6,6]]
paxlim = [[-6,6],[-6,6]]
staxlim = [[-6,6],[-6,6]]

##############################################################################

# Import models:
res1file=open(residpath_mod1,'r')
res1 = pickle.load(res1file)
res1file.close()

res2file=open(residpath_mod2,'r')
res2 = pickle.load(res2file)
res2file.close()

# Get info to plot:
evnums,evind = np.unique(res1.evnum,return_index=True)
r1event = np.array(res1.E_residual)[evind]
r2event = np.array(res2.E_residual)[evind]

r1path = res2.path_terms
r2path = res2.path_terms

sta,stind = np.unique(res1.sta,return_index=True)
stnum=res1.stnum[stind]
r1site = np.array(res1.site_terms)[stind]
r2site = np.array(res2.site_terms)[stind]


#################################
######  Auxilliary info    ######
#################$###############







#### How many stations record each unique event? ###
# Get an array with the number of stations recording each event:
event_numstas = []
unique_events,unevind=np.unique(res1.evnum,return_index=True)

for eventi in range(len(unique_events)):
    evwhere = np.where(res1.evnum==unique_events[eventi])[0]
    numstas = len(evwhere)
    # append the number of statoin srecording this event to the list:
    event_numstas.append(numstas)
    
# At the end, now turn into an array:
event_numstas = np.array(event_numstas)

################
#### How many events are recorded at each station? ###
station_numevs = []
unique_stas,unstind=np.unique(res1.sta,return_index=True)

for stai in range(len(unique_stas)):
    stwhere = np.where(res1.sta==unique_stas[stai])[0]
    numevs = len(stwhere)
    # append the number of events recorded on this station to the list:
    station_numevs.append(numevs)
    
# Turn into an array:
station_numevs = np.array(station_numevs)
  

###########################################################################################
###########################
#######  Plotting    ######
##################$########
#
#print 'Plotting event terms...'
#
### First event terms ##
## Make a straight line for relationship:
#xe=np.linspace(evaxlim[0][0],evaxlim[0][1],2)
#ye=np.linspace(evaxlim[1][0],evaxlim[1][1],2)
#
## Plotting colored by number of stations recording each event:
#cmin = min(event_numstas)
#cmax = max(event_numstas)-4
#
###Plot:
##Get colormap
##Make colormap:
#colormap_numstas=plt.get_cmap(mymap)
##Make a normalized colorscale
#cNorm=colors.Normalize(vmin=cmin, vmax=cmax)
##Apply normalization to colormap:
#scalarMap=cm.ScalarMappable(norm=cNorm, cmap=colormap_numstas)
#
##Make a fake contour plot for the colorbar:
#Z=[[0,0],[0,0]]
#levels=np.arange(cmin,cmax,0.01)
#c=plt.contourf(Z, levels, cmap=colormap_numstas)
#
##Assign values to colormap
#colorVal = scalarMap.to_rgba(event_numstas)
#
#eventfig = plt.figure()
##plt.scatter(r1event,r2event,marker='o',s=7,color='#333333')
#plt.scatter(r1event,r2event,marker='o',s=7,color=colorVal)
#plt.plot(xe,ye,color='gray',linewidth=1.5)
#
## Colrobar:
#cb = plt.colorbar(c)
#cb.set_label('Number of stations per event')
#
## Axis limits:
#plt.xlim(evaxlim[0][0],evaxlim[0][1])
#plt.ylim(evaxlim[1][0],evaxlim[1][1])
#
#plt.xlabel('Simple inversion event term (ln residual)')
#plt.ylabel('Mixed effects inversion event term (ln residual)')
#plt.title('Plot of Simple vs. Mixed effects event terms')
#
##Show the figures
#eventfig.show()
#
#evpngfile = fig_dir+'event_comp.png'
#evpdffile = fig_dir+'pdfs/event_comp.pdf'
#
#print 'Saving event figures'
#eventfig.savefig(evpngfile)
#eventfig.savefig(evpdffile)
#
#
#
################################################
### Then make path terms ##
#xp=np.linspace(paxlim[0][0],paxlim[0][1],2)
#yp=np.linspace(paxlim[1][0],paxlim[1][1],2)
#
#print 'Plotting path terms'
#pathfig = plt.figure()
#plt.scatter(r1path,r2path,marker='o',s=11,color='#333333')
#plt.plot(xp,yp,color='#9C9C9C',linewidth=1.5)
#
## Axis limits:
#plt.xlim(paxlim[0][0],paxlim[0][1])
#plt.ylim(paxlim[1][0],paxlim[1][1])
#
#plt.xlabel('Simple inversion path term (ln residual)')
#plt.ylabel('Mixed effects inversion path term (ln residual)')
#plt.title('Plot of Simple vs. Mixed effects path terms')
#
##Show the figure
#pathfig.show()
#
##Save the figure
#ppngfile = fig_dir+'path_comp.png'
#ppdffile = fig_dir+'pdfs/path_comp.pdf'
#
#print 'Saving path figures'
#pathfig.savefig(ppngfile)
#pathfig.savefig(ppdffile)
#
#
################################################
#
### Then site terms: ##
#xs=np.linspace(staxlim[0][0],staxlim[0][1],2)
#ys=np.linspace(staxlim[1][0],staxlim[1][1],2)
#
## Plotting colored by number of events recorded at each station:
#cmin = min(station_numevs)
#cmax = max(station_numevs)
#
#print 'Plotting the site figure'
#
#print 'Making colormap'
###Plot:
##Get colormap
##Make colormap:
#colormap_numevs=plt.get_cmap(mymap)
##Make a normalized colorscale
#cNorm=colors.Normalize(vmin=cmin, vmax=cmax)
##Apply normalization to colormap:
#scalarMap=cm.ScalarMappable(norm=cNorm, cmap=colormap_numevs)
#
#print 'Making fake figure for colormap'
##Make a fake contour plot for the colorbar:
#Z=[[0,0],[0,0]]
#levels=np.arange(cmin,cmax,10)
#c=plt.contourf(Z, levels, cmap=colormap_numevs)
#
#print 'Assigning values to colormap'
##Assign values to colormap
#colorVal = scalarMap.to_rgba(station_numevs)
#
#print 'starting figure, plotting...'
#stafig = plt.figure()
##plt.scatter(r1event,r2event,marker='o',s=7,color='#333333')
#plt.scatter(r1site,r2site,marker='o',s=11,color=colorVal)
#plt.plot(xs,ys,color='gray',linewidth=1.5)
#
#print 'Making colorbar'
## Colrobar:
#cb = plt.colorbar(c)
#cb.set_label('Number of events per station')
#
#print 'Axis limits...'
#
## Axis limits:
#plt.xlim(staxlim[0][0],staxlim[0][1])
#plt.ylim(staxlim[1][0],staxlim[1][1])
#
#print 'Labels'
#plt.xlabel('Simple inversion site term (ln residual)')
#plt.ylabel('Mixed effects inversion site term (ln residual)')
#plt.title('Plot of Simple vs. Mixed effects site terms')
#
#stafig.show()
#
#stpngfile = fig_dir+'station_comp.png'
#stpdffile = fig_dir+'pdfs/station_comp.pdf'
#
#print 'Saving station figures'
#stafig.savefig(stpngfile)
#stafig.savefig(stpdffile)