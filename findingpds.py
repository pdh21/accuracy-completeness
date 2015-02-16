import numpy as np
import astropy.io.fits as apif
import matplotlib.pyplot as plt
import LRmodule as LR
import LRfiles as LRf

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica Neue']})
rc('font',**{'family':'serif'})
rc('xtick',labelsize=15)
rc('ytick',labelsize=15)
rc('axes',labelsize=15)
#rc('ytick.major',pad=5)
#rc('text',usetex=True)


'''finding the surface density of background sources as a function of dS.
    This is quite calculation intensive as its a histogram of lots of things so it's saved out.
    
    First,
    rho(S, c), the background number counts, are determined by
    binning the entire matching catalogue (Section 2.2) as a
    function of the two flux densities and colour, and dividing
    by the survey area - from Chapin'''


#catfilein ="/Users/cc401/datascripts/fits/laceyfilt/lacey_07012015_MillGas.ALLVOLS_RAND_cat_PSW_FLS.fits"
#catfileout="/Users/cc401/datascripts/fits/laceyfilt/FLS_RAND_PSW_hipe_SXT.fits"
#associatedmapfile="/Users/cc401/datascripts/fits/laceyfilt/fls_itermap_lacey_07012015_RAND_simulated_observation_w_noise_PSW_hipe_filtered.fits"


#catfilein ="/Users/cc401/datascripts/fits/laceyfilt/uds_catalog_clustering_bias_1.0_cut_1000mJy_0_PLW.fits"
#catfileout="/Users/cc401/datascripts/fits/laceyfilt/uds_itermap_simulation_clustering_bias_1.0_cut_1000mJy_w_noise_0_PLW_SXT.fits"
#associatedmapfile ="/Users/cc401/datascripts/fits/laceyfilt/uds_itermap_simulation_clustering_bias_1.0_cut_1000mJy_w_noise_0_PLW.fits"

catfilein=LRf.catfilein
catfileout=LRf.catfileout
associatedmapfile=LRf.associatedmapfile

catin  = apif.open(catfilein)[1].data
catout = apif.open(catfileout)[1].data
hdr    = apif.open(associatedmapfile)[1].header
map    = apif.open(associatedmapfile)[1].data


fluxi=( catin["flux"]*LRf.inputfluxfactor)
fluxo=(catout["flux"])
pixsize=3600*np.abs(hdr["cd1_1"]) #in arcsecs

minfluxbin=np.floor(np.min(fluxi)-np.max(fluxo))
maxfluxbin=np.ceil(np.max(fluxi)-np.min(fluxo))

print 'looking at fluxes between ', minfluxbin,maxfluxbin

fluxbinedges  =np.arange(minfluxbin,maxfluxbin,1)
fluxbincentres=LR.linBinCentres(fluxbinedges)
phisttot      =np.zeros(len(fluxbinedges)-1)

for i in range(len(fluxo)):
    if i%100 == 0: print i,' out of ', len(fluxo)
    hist,junk=np.histogram(fluxi-fluxo[i],bins=fluxbinedges)
    phisttot+=hist

#calculate map area
indx,indy=np.where(~np.isnan(map))
maparea=len(indx)*pixsize**2

phisttot/=np.float(maparea*len(fluxo)) #have to divide by number of output sources and by the surface area in arcsec

np.savetxt(LRf.pdsfile, zip(fluxbincentres,phisttot), header="dS mJy/tsurfacedensity per arcsec", comments="#")

plt.plot(fluxbincentres,phisttot,'k', drawstyle='steps-mid',label='q(ds)')
plt.xlabel('ds (mJy)')
plt.ylabel('surface density')
plt.savefig(LRf.pdsplotfile)
plt.show()


