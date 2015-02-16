import numpy as np
import astropy.io.fits as apif
import matplotlib.pyplot as plt
import LRmodule as LR

'''just checking that kdTree stuff gives real distances and not approximations, that way I can speed things up'''

catfilein ="/Users/cc401/datascripts/fits/laceyfilt/lacey_07012015_MillGas.ALLVOLS_RAND_cat_PSW_FLS.fits"
catfileout="/Users/cc401/datascripts/fits/laceyfilt/FLS_RAND_PSW_hipe_SXT.fits"
associatedmapfile="/Users/cc401/datascripts/fits/laceyfilt/fls_itermap_lacey_07012015_RAND_simulated_observation_w_noise_PSW_hipe_filtered.fits"


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
radbinedges80=np.arange(0,51)
bincentres80=LR.linBinCentres(radbinedges80)
radbinedges50=np.arange(0,51)
bincentres50=LR.linBinCentres(radbinedges50)


#get the num input objects matched to output source as fn of distance
histtot80,histtoterr = LR.histDistance(xi,yi,xo,yo,radbinedges80/pixsize,distupper=80/pixsize)
histtot50,histtoterr = LR.histDistance(xi,yi,xo,yo,radbinedges50/pixsize,distupper=50/pixsize)

print histtot80-histtot50

#THIS RETURNS ALL ZEROS. GOOD. YAY! 
