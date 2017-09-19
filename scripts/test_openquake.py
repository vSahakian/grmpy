## Test OpenQuake


from openquake.hazardlib.gsim.abrahamson_2014 import AbrahamsonEtAl2014
from openquake.hazardlib import imt, const
from openquake.hazardlib.gsim.base import RuptureContext
from openquake.hazardlib.gsim.base import DistancesContext
from openquake.hazardlib.gsim.base import SitesContext
import numpy as np
import gmpe as gm
import matplotlib.pyplot as plt

fig_dir = '/Users/vsahakian/anza/models/statistics/misc/oq_vs_matlab/'


## This all works..... ##

#ASK14 = AbrahamsonEtAl2014()
#
#IMT = imt.PGA()
#rctx = RuptureContext()
#dctx = DistancesContext()
#sctx = SitesContext()
#sctx_rock = SitesContext()
#
#rctx.rake = 0.0
#rctx.dip = 90.0
#rctx.ztor = 7.13
#rctx.mag = 3.0
#rctx.mag = np.linspace(0.1,5.)
#rctx.width = 10.0
#rctx.hypo_depth = 8.0
#
##dctx.rrup = np.logspace(1,np.log10(200),100)
#dctx.rrup = np.logspace(np.log10(10),np.log10(10.0),1)
#
#
## Assuming average ztor, get rjb:
#dctx.rjb = np.sqrt(dctx.rrup**2 - rctx.ztor**2)
#dctx.rhypo = dctx.rrup
#dctx.rx = dctx.rjb
#dctx.ry0 = dctx.rx
#
#sctx.vs30 = np.ones_like(dctx.rrup) * 760.0
#sctx.vs30measured = np.full_like(dctx.rrup, False, dtype='bool')
#sctx.z1pt0 = np.ones_like(dctx.rrup) * 0.05
#
#lmean_ask14, sd_ask14 = ASK14.get_mean_and_stddevs(
#    sctx, rctx, dctx, IMT, [const.StdDev.TOTAL])





###############################################################################
#
### Now try using mien which uses this... ##
#
## Read in annemarie's ask:
#ask_data_ab = np.genfromtxt('/Users/vsahakian/anza/data/coeffs/ASK_for_Valerie.txt',skip_header=1)
#
#M = ask_data_ab[:,0]
##Rrup = np.arange(0, np.log10(200),0.1)
##Rrup = np.array([np.log10(10),np.log10(80)])
#Rrup = np.array([10,80])
#vs30 = np.array([200,400,600,760,1100])
#ztor = np.array([0,5,6,7.13,10])
#z1pt0 = np.array([0.005,0.05,0.10,0.20])
#
#L2 = []
#L1 = []
#descriptions = []
#
#for ivs30 in range(len(vs30)):
#    for jztor in range(len(ztor)):
#        for kz1pt0 in range(len(z1pt0)):
#
#            mediangm,sdgm = gm.oq_ask2014(M,Rrup,vs30=vs30[ivs30],ztor=ztor[jztor],z1pt0=z1pt0[kz1pt0])
#            
#            mediangm = np.log10(np.exp(mediangm))
#            
#            # Difference, L2:
#            x2 = (np.sum((mediangm[:,0] - ask_data_ab[:,1])**2))/len(mediangm[:,0])
#            L2.append(x2)
#            
#            L1.append(np.sum(np.abs(mediangm[:,0] - ask_data_ab[:,1]))/len(mediangm[:,0]))
#            
#            descriptions.append('vs30'+np.str(vs30[ivs30]) + '_ztor' + np.str(ztor[jztor]) + '_z1pt0' + np.str(z1pt0[kz1pt0]))
#            
#            
#            # Plot together:
#            plt.figure()
#            plt.plot(M,mediangm[:,0],linewidth=2,linestyle='--',color='black',label='OpenQuake 10km')
#            plt.plot(M,np.log10(ask_data_ab[:,1]),linewidth=2,linestyle='-.',color='black',label='Matlab 10km')
#            
#            plt.plot(M,mediangm[:,1],linewidth=2,linestyle='--', color='red',label='OpenQuake 80km')
#            plt.plot(M,np.log10(ask_data_ab[:,2]),linestyle = '-.',linewidth=2,color='red',label='Matlab 80km')
#            
#            plt.legend(loc=4)
#            plt.xlabel('Magnitude')
#            plt.ylabel('log10 PGA')
#            plt.title('Vs30 = ' + np.str(vs30[ivs30]) + ', ztor = ' + np.str(ztor[jztor]) + ', z1pt0 = ' + np.str(z1pt0[kz1pt0]))
#
#            figname = 'Vs30_' + np.str(vs30[ivs30]) + '_ztor_' + np.str(ztor[jztor]) + '_z1pt0_' + np.str(z1pt0[kz1pt0]) + '.png'
#            plt.savefig(fig_dir + figname)
#
#L2 = np.array(L2)
#L1 = np.array(L1)
#descriptions = np.array(descriptions)
#
#bestl2fit_ind = np.argmin(L2)
#print 'Best L2 fit is ' + descriptions[bestl2fit_ind]
#
#bestl1fit_ind = np.argmin(L1)
#print 'Best L1 fit is ' + descriptions[bestl1fit_ind]



##############################################################################

## Now try using mien which uses this... ##

# Read in annemarie's ask:
ask_data_ab = np.genfromtxt('/Users/vsahakian/anza/data/coeffs/ASK_for_Valerie.txt',skip_header=1)

M = ask_data_ab[:,0]
#Rrup = np.arange(0, np.log10(200),0.1)
#Rrup = np.array([np.log10(10),np.log10(80)])
Rrup = np.array([np.log10(5),np.log10(10),np.log10(80),np.log10(150)])
vs30 = 1000
ztor = 0
z1pt0 = 0.05


mediangm,sdgm = gm.oq_ask2014(M,Rrup,vs30=vs30,ztor=ztor,z1pt0=z1pt0)

mediangm = np.log10(np.exp(mediangm))



# Plot together:
plt.figure()
plt.plot(M,mediangm[:,0],linewidth=2,linestyle='--',color='red',label='OpenQuake 5km')

plt.plot(M,mediangm[:,1],linewidth=2,linestyle='--', color='black',label='OpenQuake 10km')

plt.plot(M,mediangm[:,2],linewidth=2,linestyle='--', color='green',label='OpenQuake 80km')

plt.plot(M,mediangm[:,3],linewidth=2,linestyle='--', color='cyan',label='OpenQuake 150km')


plt.legend(loc=4)
plt.xlabel('Magnitude')
plt.ylabel('log10 PGA')
plt.title('Vs30 = ' + np.str(vs30) + ', ztor = ' + np.str(ztor) + ', z1pt0 = ' + np.str(z1pt0))

figname = 'logR_ztor0.png'
plt.savefig(fig_dir + figname)

