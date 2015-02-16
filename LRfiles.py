
'''put in all file names here! carries through the rest of the code :)
        '''

##################
#catfilein=
#catfileout=
#associatedmapfile=
#
#qdsfile=
#pdsfile=
#pdsplotfile=
#qdsplotfile=
#################


catfilein="/Users/cc401/datascripts/fits/laceyfilt/cosmos_catalog_clustering_bias_1.0_cut_500mJy_0_PSW.fits"
catfileout="/Users/cc401/datascripts/fits/laceyfilt/cosmos_itermap_simulation_clustering_bias_1.0_cut_500mJy_0_PSW_SXT.fits"
associatedmapfile="/Users/cc401/datascripts/fits/laceyfilt/cosmos_itermap_simulation_clustering_bias_1.0_cut_500mJy_0_PSW.fits"

qdsfile="/Users/cc401/Dropbox/thesis/counts/helms/cosmosPSWqds.dat"
pdsfile="/Users/cc401/Dropbox/thesis/counts/helms/udscosmosPSWpds.dat"
pdsplotfile="/Users/cc401/Dropbox/thesis/counts/helms/cosmosPSWpds.eps"
qdsplotfile="/Users/cc401/Dropbox/thesis/counts/helms/cosmosPSWqdsthreeplots.eps"

catfilein="/Users/cc401/datascripts/fits/laceyfilt/cosmos_catalog_clustering_bias_1.0_cut_500mJy_0_PLW.fits"
catfileout="/Users/cc401/datascripts/fits/laceyfilt/cosmos_itermap_simulation_clustering_bias_1.0_cut_500mJy_0_PLW_SXT.fits"
associatedmapfile="/Users/cc401/datascripts/fits/laceyfilt/cosmos_itermap_simulation_clustering_bias_1.0_cut_500mJy_0_PLW.fits"

qdsfile="/Users/cc401/Dropbox/thesis/counts/helms/cosmosPLWqds.dat"
pdsfile="/Users/cc401/Dropbox/thesis/counts/helms/udscosmosPLWpds.dat"
pdsplotfile="/Users/cc401/Dropbox/thesis/counts/helms/cosmosPLWpds.eps"
qdsplotfile="/Users/cc401/Dropbox/thesis/counts/helms/cosmosPLWqdsthreeplots.eps"


inputfluxfactor=1000. #in case things in Jy rather than mJy

sigma=5.19249936131#PSW#7.62732607131#PMW#12.1017095101#PLW
sigma=12.1017
radius=23. #PLW
LRThresh=1.04567237079#1.86#PLW

#bins
#radius
radiusmax=60. #max binedge in arcmin default 50
#fluxes
fluxminedge=-42 #default -42
fluxmaxedge=22 #default 22
fluxbinsize=2 #default 2 #MAKE IT A MUTLIPLE OF 1 TROLOLOLOL