import numpy as np
import astropy.io.fits as apif
import matplotlib.pyplot as plt
import astropy.wcs as apw
from scipy.spatial import KDTree
from scipy.optimize import curve_fit
import numpy.random as npr
from matplotlib import cm

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica Neue']})
rc('font',**{'family':'serif'})
rc('xtick',labelsize=15)
rc('ytick',labelsize=15)
rc('axes',labelsize=15)
#rc('ytick.major',pad=5)
rc('text',usetex=True)


'''This code mirrors fig 2 in Wang'13, arxiv:1312.0552, shows the radial distribution of background input and output sources. 
    
    Specify the input and output catalogues. Converts the ra and decs into x and y. Uses a kd-tree method to quickly calculate the distances, then scales the x axis to be in arcminutes. Y axis number if matches per extracted source. Data is also binned by flux difference. '''


catfilein ="/Users/cc401/datascripts/fits/laceyfilt/lacey_07012015_MillGas.ALLVOLS_RAND_cat_PSW_FLS.fits"
catfileout="/Users/cc401/datascripts/fits/laceyfilt/FLS_RAND_PSW_hipe_SXT.fits"
associatedmapfile="/Users/cc401/datascripts/fits/laceyfilt/fls_itermap_lacey_07012015_RAND_simulated_observation_w_noise_PSW_hipe_filtered.fits"


catin  = apif.open(catfilein)[1].data
catout = apif.open(catfileout)[1].data
hdr    = apif.open(associatedmapfile)[1].header
map    = apif.open(associatedmapfile)[1].data

rai,dei,fluxi=( catin["ra"], catin["dec"], catin["s250"])
rao,deo,fluxo=(catout["ra"],catout["dec"],catout["flux"])



#what is size of a pixel?
pixsize=3600*np.abs(hdr["cd1_1"]) #in arcsecs

#change from ra and dec to x and y
WCSmap=apw.WCS(hdr)
xi,yi=WCSmap.wcs_world2pix(rai,dei,0)
xo,yo=WCSmap.wcs_world2pix(rao,deo,0)

radbinedges=np.arange(0,51)
fluxbinedges=np.arange(-42,24,2)
histtot=np.zeros((len(radbinedges)-1,len(fluxbinedges)-1))

#make KDTree
Tree=KDTree(np.array([xi,yi]).T,leafsize=100)

#real output posns matching
for i in range(len(xo)):
    if i%100==0: print i," out of ",len(xo)
    d,ind=Tree.query(np.array([xo[i],yo[i]]),k=None,distance_upper_bound=80/pixsize)
    hist,junk,junk2=np.histogram2d(d,fluxi[ind]-fluxo[i],bins=[radbinedges/pixsize,fluxbinedges])
    #if i%100==0:
    #plt.subplot(2,1,1)
    #plt.plot((xi[ind]-xo[i])*pixsize, (yi[ind]-yo[i])*pixsize, 'bx')
    #plt.subplot(2,1,2)
    #plt.plot(radbinedges[:-1]+0.5*(radbinedges[1:]-radbinedges[:-1]),hist)
    #plt.show()
    histtot+=hist

#randomise the xo positions
stilltorandom=np.ones(len(xo))
xorand=xo.copy()
yorand=yo.copy()
mapsize=map.shape
#loop until all assigned
while np.sum(stilltorandom) != 0:
    indtorand,=np.where(stilltorandom == 1)
    numtorand=len(indtorand)
    print 'left to random: ',numtorand
    tempxo=npr.random(numtorand)*mapsize[1]
    tempyo=npr.random(numtorand)*mapsize[0]
    gdind,=np.where(~np.isnan(map[tempyo.astype(int),tempxo.astype(int)]))
    indyes=indtorand[gdind]
    xorand[indyes]=tempxo[gdind]
    yorand[indyes]=tempyo[gdind]
    stilltorandom[indyes]=0


#rand output posns matching
histtotrand=np.zeros((len(radbinedges)-1,len(fluxbinedges)-1))

for i in range(len(xorand)):
    if i%100==0: print i," out of ",len(xo)
    d,ind=Tree.query(np.array([xorand[i],yorand[i]]),k=None,distance_upper_bound=80/pixsize)
    hist,junk,junk2=np.histogram2d(d,fluxi[ind]-fluxo[i],bins=[radbinedges/pixsize,fluxbinedges])
    #if i%100==0:
    #plt.subplot(2,1,1)
    #plt.plot((xi[ind]-xo[i])*pixsize, (yi[ind]-yo[i])*pixsize, 'bx')
    #plt.subplot(2,1,2)
    #plt.plot(radbinedges[:-1]+0.5*(radbinedges[1:]-radbinedges[:-1]),hist)
    #plt.show()
    histtotrand+=hist

#creating residual array
histtotresi=histtot-histtotrand    

#creating a vmin_vmax for the colorbar
hists=[histtot,histtotrand,histtotresi]
allhists=np.array(hists).flatten()
maxhist=np.max(allhists)
minhist=np.min(allhists)


#plotting 
aspectratio = 1.0*(fluxbinedges[-1] - fluxbinedges[0])/(1.0*radbinedges[-1] - radbinedges[0])

fig,axes = plt.subplots(nrows=3,ncols=1,figsize=(5,11))

for i,ax in enumerate(axes.flat):
    im=ax.imshow(hists[i],extent=[fluxbinedges[0],fluxbinedges[-1],radbinedges[0],radbinedges[-1]],aspect=aspectratio,vmin=minhist,vmax=maxhist,cmap=cm.hot, origin='lower')
fig.subplots_adjust(right=0.7)

axes[2].set_xlabel(r"fluxin-fluxout (mJy)", fontsize=15)
axes[1].set_ylabel(r"distance (arcsec)",fontsize=15)
axes[1].yaxis.labelpad=10
axes[2].xaxis.labelpad=12

cbar_ax=fig.add_axes([0.75,0.15,0.05,0.7])
cbar=fig.colorbar(im, cax=cbar_ax)
cbar.set_label(r"Counts", fontsize=15)
plt.savefig("/Users/cc401/Dropbox/thesis/counts/helms/fluxpositionmatches.eps")
plt.show()




