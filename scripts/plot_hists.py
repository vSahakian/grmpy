# Plot histograms of R and M for a particular dataset
# VJS 2/2017


import cPickle as pickle
import matplotlib.pyplot as plt

what_home=1

if what_home==0:
    #Desktop:
    HOME='/media/vsahakian/katmai'
    codehome='/home/vsahakian'
elif what_home==1:
    #Mac:
    HOME='/Users/vsahakian'
    codehome='/Users/vsahakian'

#Location to store figures:
home=HOME+'/anza'
invrun='v2anza2013'
fig_dir=home+'/models/figs/'+invrun+'/'
obj_dir=home+'/models/pckl/'+invrun+'/'
model_dir=home+'/models/pckl/'+invrun+'/'


# Set database name
dbpath = home + '/data/databases/'+ 'v2anza2013/v2anza2013_pgrid_5sta.pckl'
dbname = 'v2anza2013'

# Set names for output files:
rpdf = fig_dir + 'pdf/' + dbname + '_histogram_rrup.pdf'
rpng = fig_dir + dbname + '_histogram_rrup.png'

mwpdf = fig_dir + 'pdf/' + dbname + '_histogram_mw.pdf'
mwpng = fig_dir + dbname + '_histogram_mw.png'

mrpdf = fig_dir + 'pdf/' + dbname + '_rrup_mw.pdf'
mrpng = fig_dir + dbname + '_rrup_mw.png'

# Open database
dbfile=open(dbpath,'r')
db = pickle.load(dbfile)
dbfile.close()

# Plotting parameters:
nbins = 1000

###########################
# Plot histogram of Rrup:
flag_type = 0
axlims = [[0,190],[0,1300]]

# Plot
rhist = db.plot_histogram(axlims,nbins,flag_type)

# Save:
rhist.savefig(rpdf)
rhist.savefig(rpng)

##########################
# Plot histogram of M:
flag_type = 1
axlims = [[0,5],[0,1500]]

#Plot
mwhist = db.plot_histogram(axlims,nbins,flag_type)

#Save
mwhist.savefig(mwpdf)
mwhist.savefig(mwpng)


##  Now also scatter the mw vs. rrup...

axlims = [[0,5],[0,200]]


mrfig = plt.figure()
plt.scatter(db.mw,db.r,edgecolors='#626263',facecolors='none',s=9)

plt.xlim(axlims[0])
plt.ylim(axlims[1])

plt.xlabel('M')
plt.ylabel('Rrup')
plt.title('Rrup vs. M for database')

mrfig.savefig(mrpdf)
mrfig.savefig(mrpng)
