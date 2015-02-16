import numpy as np
import astropy.io.fits as apif
import matplotlib.pyplot as plt
import numpy.random as npr
from matplotlib import cm
from scipy.spatial import KDTree
import LRmodule as LR
from numpy.random import shuffle
import LRfiles as LRf

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica Neue']})
rc('font',**{'family':'serif'})
rc('xtick',labelsize=15)
rc('ytick',labelsize=15)
rc('axes',labelsize=15)
#rc('ytick.major',pad=5)
#rc('text',usetex=True)


'''This code mirrors fig 2 in Wang'13, arxiv:1312.0552, shows the radial distribution of background input and output sources
    
    Specify the input and output catalogues. Converts the ra and decs into x and y. Uses a kd-tree method to quickly calculate the distances, then scales the x axis to be in arcminutes. Y axis number if matches per extracted source. Data is also binned by flux difference. '''


#catfilein ="/Users/cc401/datascripts/fits/laceyfilt/lacey_07012015_MillGas.ALLVOLS_RAND_cat_PSW_FLS.fits"
#catfileout="/Users/cc401/datascripts/fits/laceyfilt/FLS_RAND_PSW_hipe_SXT.fits"
#associatedmapfile="/Users/cc401/datascripts/fits/laceyfilt/fls_itermap_lacey_07012015_RAND_simulated_observation_w_noise_PSW_hipe_filtered.fits"

#catfilein ="/Users/cc401/datascripts/fits/laceyfilt/uds_catalog_clustering_bias_1.0_cut_1000mJy_0_PLW.fits"
#catfileout="/Users/cc401/datascripts/fits/laceyfilt/uds_itermap_simulation_clustering_bias_1.0_cut_1000mJy_w_noise_0_PLW_SXT.fits"
#associatedmapfile ="/Users/cc401/datascripts/fits/laceyfilt/uds_itermap_simulation_clustering_bias_1.0_cut_1000mJy_w_noise_0_PLW.fits"

catfilein=LRf.catfilein
catfileout=LRf.catfileout
associatedmapfile=LRf.associatedmapfile
optRad=LRf.radius

catin  = apif.open(catfilein)[1].data
catout = apif.open(catfileout)[1].data
hdr    = apif.open(associatedmapfile)[1].header
map    = apif.open(associatedmapfile)[1].data

rai,dei,fluxi=( catin["ra"], catin["dec"], catin["flux"]*LRf.inputfluxfactor)
rao,deo,fluxo=(catout["ra"],catout["dec"],catout["flux"])

#what is size of a pixel?
pixsize=3600*np.abs(hdr["cd1_1"]) #in arcsecs

#change from ra and dec to x and y
xi,yi=LR.adxy(hdr,rai,dei)
xo,yo=LR.adxy(hdr,rao,deo)

fluxbinedges  =np.arange(LRf.fluxminedge,LRf.fluxmaxedge,LRf.fluxbinsize)
fluxbincentres=LR.linBinCentres(fluxbinedges)
qhisttotreal  =np.zeros(len(fluxbinedges)-1)
qhisttotrand  =np.zeros(len(fluxbinedges)-1)
qhisttot      =np.zeros(len(fluxbinedges)-1)


#randomise the xo positions
xorand,yorand=LR.randomXYsInMap(map,xo,yo)

#q(dq)
#matching between real out and real in
Tree=KDTree(np.array([xi,yi]).T,leafsize=100)

for i in range(len(xo)):
    if i%100==0: print i," out of ",len(xo)
    #need to pull out the sources within the optimal search radius
    ind=LR.findSourcesInRadius(xo[i],yo[i],xi,yi,optRad/pixsize,kdTree=Tree)
    hist,junk=np.histogram(fluxi[ind]-fluxo[i],bins=fluxbinedges)
    qhisttotreal+=hist

for i in range(len(xorand)):
    if i%100==0: print i," out of ",len(xorand)
    #need to pull out the sources within the optimal search radius
    ind=LR.findSourcesInRadius(xorand[i],yorand[i],xi,yi,optRad/pixsize,kdTree=Tree)
    hist,junk=np.histogram(fluxi[ind]-fluxo[i],bins=fluxbinedges)
    qhisttotrand+=hist

#must normalised histogram to the average expected number of counterparts per source 
#AND the counts per bin
qhisttot=qhisttotreal-qhisttotrand
qhisttotreal/=(np.float(len(xorand)) * (fluxbinedges[1:]-fluxbinedges[:-1]))
qhisttotrand/=(np.float(len(xorand)) *(fluxbinedges[1:]-fluxbinedges[:-1]))
qhisttot/= (np.float(len(xorand)) * (fluxbinedges[1:]-fluxbinedges[:-1]))
qhisttoterr=np.sqrt(qhisttotreal**2+qhisttotrand**2)/np.sqrt(len(xo)) #error is poission count

  
#loop over the histograms
hists=[qhisttotreal,qhisttotrand,qhisttot]

##plotting 


plt.plot(fluxbincentres,hists[0],color='b',drawstyle='steps-mid',label="real output matched with real input")
plt.plot(fluxbincentres,hists[1],color='r',drawstyle='steps-mid',label="random output matched with real input")
plt.errorbar(fluxbincentres,hists[2],yerr=qhisttoterr,color='k',drawstyle='steps-mid',label="residual showing true matches")
plt.plot(fluxbincentres,np.zeros(len(fluxbincentres)),'g:')
plt.legend(loc="best")

plt.xlabel(r"fluxin-fluxout (mJy)", fontsize=15)
plt.ylabel(r"counts",fontsize=15)
#axes[1].yaxis.labelpad=10
#axes[2].xaxis.labelpad=12
#
#cbar_ax=fig.add_axes([0.75,0.15,0.05,0.7])
#cbar=fig.colorbar(im, cax=cbar_ax)
#cbar.set_label(r"Counts", fontsize=15)
plt.savefig(LRf.qdsplotfile)
np.savetxt(LRf.qdsfile, zip(fluxbincentres,qhisttot))
plt.show()




