import numpy as np
import astropy.io.fits as apif
import matplotlib.pyplot as plt
import LRmodule as LR
import LRfiles as LRf



'''This code mirrors fig 3 in Wang'13, arxiv:1312.0552, shows the radial distribution of positional offsets between extracted sources and input sources per extracted sources. The radial distribution of background input and output follows a Rayleigh distribution.
    
    Also calculates the optimal search radius. NB not theoretical
    
    Specify the input and output catalogues. Converts the ra and decs into x and y. Uses a kd-tree method to quickly calculate the distances, then scales the x axis to be in arcminutes. Y axis number if matches per extracted source???? Poisson errors are used'''



#catfilein ="/Users/cc401/datascripts/fits/laceyfilt/lacey_07012015_MillGas.ALLVOLS_RAND_cat_PSW_FLS.fits"
#catfileout="/Users/cc401/datascripts/fits/laceyfilt/FLS_RAND_PSW_hipe_SXT.fits"
#associatedmapfile="/Users/cc401/datascripts/fits/laceyfilt/fls_itermap_lacey_07012015_RAND_simulated_observation_w_noise_PSW_hipe_filtered.fits"


#catfilein ="/Users/cc401/datascripts/fits/laceyfilt/lacey_07012015_MillGas.ALLVOLS_cat_PMW_COSMOS.fits"
#catfileout="/Users/cc401/datascripts/fits/laceyfilt/COSMOS_PMW_hipe_SXT.fits"
#associatedmapfile="/Users/cc401/datascripts/fits/laceyfilt/cosmos_itermap_lacey_07012015_simulated_observation_w_noise_PMW_hipe.fits.gz"


#catfilein ="/Users/cc401/datascripts/fits/laceyfilt/uds_catalog_clustering_bias_1.0_cut_1000mJy_0_PSW.fits"
#catfileout="/Users/cc401/datascripts/fits/laceyfilt/uds_itermap_simulation_clustering_bias_1.0_cut_1000mJy_w_noise_0_PSW_SXT.fits"
#associatedmapfile ="/Users/cc401/datascripts/fits/laceyfilt/uds_itermap_simulation_clustering_bias_1.0_cut_1000mJy_w_noise_0_PSW.fits"

catfilein=LRf.catfilein
catfileout=LRf.catfileout
associatedmapfile=LRf.associatedmapfile


catin  = apif.open(catfilein)[1].data
catout = apif.open(catfileout)[1].data
hdr    = apif.open(associatedmapfile)[1].header

rai,dei=( catin["ra"], catin["dec"])
rao,deo=(catout["ra"],catout["dec"])

#what is size of a pixel?
pixsize=3600*np.abs(hdr["cd1_1"]) #in arcsecs

#change from ra and dec to x and y
xi,yi=LR.adxy(hdr,rai,dei)
xo,yo=LR.adxy(hdr,rao,deo)

#define binedges
radbinedges=np.linspace(0,LRf.radiusmax,endpoint=True)

#get the num input objects matched to output source as fn of distance
histtot,histtoterr = LR.histDistance(xi,yi,xo,yo,radbinedges/pixsize)

#get binCentres
bincentres=LR.linBinCentres(radbinedges)

plt.errorbar(bincentres,histtot,yerr=histtoterr,color='r',linestyle='steps',drawstyle="steps-mid",ecolor='r',label="data matches")
plt.xlabel('Radius r (arcsec)')
plt.ylabel('N(r)')

#fit the data to find sigma and the bg scaling*r counts amount
sigma,scaling,scaling2,pcov=LR.findSigmaR(bincentres,histtot, error=histtoterr)

#plotting the stuff up :2 
plt.plot(bincentres, LR.rayleighPlus(bincentres,sigma,scaling,scaling2),'k-',label="fit")
plt.plot(bincentres, histtot-LR.rayleighPlus(bincentres,sigma, scaling,scaling2),color='b',linestyle='steps',label="residual of fit",drawstyle="steps-mid")
plt.plot(bincentres, LR.rayleighScaled(bincentres,sigma,scaling), linestyle='--', color='k',label="real matches")
plt.plot(bincentres, scaling2*bincentres,'k:',label="bg counts")
plt.legend(loc="best")
plt.show()


optRad,SNRmax,indRad=LR.findOptimalRadius(bincentres, sigma, scaling,scaling2, tolerance=1e-4,plot=True)






