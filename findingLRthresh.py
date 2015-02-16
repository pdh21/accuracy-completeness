import numpy as np
import astropy.io.fits as apif
import matplotlib.pyplot as plt
import astropy.wcs as apw
import LRmodule as LR
from scipy.interpolate import interp1d
import LRfiles as LRf



'''This code mirrors the LR stuff in Wang'13, arxiv:1312.0552
    '''


#catfilein ="/Users/cc401/datascripts/fits/laceyfilt/lacey_07012015_MillGas.ALLVOLS_RAND_cat_PSW_FLS.fits"
#catfileout="/Users/cc401/datascripts/fits/laceyfilt/FLS_RAND_PSW_hipe_SXT.fits"
#associatedmapfile="/Users/cc401/datascripts/fits/laceyfilt/fls_itermap_lacey_07012015_RAND_simulated_observation_w_noise_PSW_hipe_filtered.fits"


#catfilein ="/Users/cc401/datascripts/fits/laceyfilt/uds_catalog_clustering_bias_1.0_cut_1000mJy_0_PLW.fits"
#catfileout="/Users/cc401/datascripts/fits/laceyfilt/uds_itermap_simulation_clustering_bias_1.0_cut_1000mJy_w_noise_0_PLW_SXT.fits"
#associatedmapfile ="/Users/cc401/datascripts/fits/laceyfilt/uds_itermap_simulation_clustering_bias_1.0_cut_1000mJy_w_noise_0_PLW.fits"


#qdsfile="/Users/cc401/Dropbox/thesis/counts/helms/UDSPLWqds.dat"
#pdsfile="/Users/cc401/Dropbox/thesis/counts/helms/udsPLWpds.dat"
#optradius=49. #23. #49.
#sigma=12.055 #5.546 #12.055

catfilein=LRf.catfilein
catfileout=LRf.catfileout
associatedmapfile=LRf.associatedmapfile
qdsfile=LRf.qdsfile
pdsfile=LRf.pdsfile
optradius=LRf.radius
sigma=LRf.sigma

print "loading files"
qdsBinCentres, qds=np.loadtxt(qdsfile, unpack=True)
pdsBinCentres, pds=np.loadtxt(pdsfile,unpack=True)

catin  = apif.open(catfilein)[1].data
catout = apif.open(catfileout)[1].data
hdr    = apif.open(associatedmapfile)[1].header
map    = apif.open(associatedmapfile)[1].data

rai,dei,fluxi=( catin["ra"], catin["dec"], catin["flux"]*LRf.inputfluxfactor)
rao,deo,fluxo=(catout["ra"],catout["dec"],catout["flux"])

#what is size of a pixel?
pixsize=3600*np.abs(hdr["cd1_1"]) #in arcsecs

#change from ra and dec to x and y
"changing to euclidean coordinates"
xi,yi=LR.adxy(hdr,rai,dei)
xo,yo=LR.adxy(hdr,rao,deo)
#making randomised catalogue
xorand,yorand=LR.randomXYsInMap2(map,xo,yo)

print "defining bins and rebinning"
#define binedges
radbinedges  =np.linspace(0,LRf.radiusmax,endpoint=True)
fluxbinedges =np.arange(LRf.fluxminedge,LRf.fluxmaxedge,LRf.fluxbinsize)
#get binCentres
radBinCentres=LR.linBinCentres(radbinedges)
fluxBinCentres=LR.linBinCentres(fluxbinedges)

#rebin the pds thing 
pdsBinEdges=LR.linBinEdges(pdsBinCentres)
pds=LR.rebinDataBins(pds,pdsBinCentres, pdsBinEdges, fluxbinedges)

print "making kdTree and doing interpolations"
#make kdTree
Tree=LR.getkdTree(xi,yi)
#do interpolation of pds and qds
pdsInterp =interp1d(fluxBinCentres,pds)
qdsInterp =interp1d(fluxBinCentres,qds)
frInterp  =interp1d(np.insert(radBinCentres,0,0.0),np.insert(LR.rayleigh(radBinCentres,sigma),0,0.0)) #adding 0,0 as bincentres for the interpolation too as we know it goes to 0 there

print "calculating LRs"
allLRs=[]
for i in range(len(xo)):
    if i%100 == 0: print "source ",i," out of ",len(xo)
    #find the sources that are within the optimal search range only
    r,rInd=Tree.query(np.array([xo[i],yo[i]]),k=None,distance_upper_bound=optradius/pixsize)
    rInd=np.array(rInd)
    if len(rInd) == 0: continue
    r=np.array(r*pixsize)
    #find sources in the flux range:
    dFlux=fluxi[rInd]-fluxo[i]
    fInd,=np.where((dFlux <= fluxBinCentres[-1]) & (dFlux >= fluxBinCentres[0]))
    if len(fInd) > 0:
        gdxi=xi[rInd[fInd]]
        gdyi=yi[rInd[fInd]]
        dFlux=dFlux[fInd]
        r=r[fInd]
    
        
        LikeRatios=(frInterp(r)*qdsInterp(dFlux))/(pdsInterp(dFlux)*np.pi*2*r)

        if np.any(LikeRatios <0):
            print "LR:",LikeRatios,"\nr: ",r,"\ndS: ",dFlux,"\nfr: ",frInterp(r)," \nqds: ",qdsInterp(dFlux)," \npds: ", pdsInterp(dFlux)
        allLRs.append(list(LikeRatios))

allLRs=np.array([LR1 for LRs in allLRs for LR1 in LRs])
print np.min(allLRs),np.max(allLRs)
LRhist,LRBinEdges=np.histogram(allLRs,bins=50)
plt.plot(LR.linBinCentres(LRBinEdges),LRhist,label='real posns')


print "calculating rand LRs"
allrandLRs=[]
for i in range(len(xorand)):
    if i%100 == 0: print "source ",i," out of ",len(xo)
    #find the sources that are within the optimal search range only
    r,rInd=Tree.query(np.array([xorand[i],yorand[i]]),k=None,distance_upper_bound=optradius/pixsize)
    rInd=np.array(rInd)
    if len(rInd) == 0: continue
    r=np.array(r*pixsize)
    ind=LR.findSourcesInRadius(xorand[i],yorand[i],xi,yi,optradius/pixsize,kdTree=Tree)
    #find sources in the flux range:
    dFlux=fluxi[rInd]-fluxo[i]
    fInd,=np.where((dFlux <= fluxBinCentres[-1]) & (dFlux >= fluxBinCentres[0]))
    if len(fInd) > 0:
        gdxi=xi[rInd[fInd]]
        gdyi=yi[rInd[fInd]]
        dFlux=dFlux[fInd]
        r=r[fInd]
        
        LikeRatios=(frInterp(r)*qdsInterp(dFlux))/(pdsInterp(dFlux)*np.pi*2*r)
        
        #if np.any(LikeRatios <0):
        #print "LR:",LikeRatios,"fr: ",frInterp(r)," qds: ",qdsInterp(dFlux)," pds", pdsInterp(dFlux)
        allrandLRs.append(list(LikeRatios))

allrandLRs=np.array([LR1 for LRs in allrandLRs for LR1 in LRs])
LRBinCentres = LR.linBinCentres(LRBinEdges)
print np.min(allrandLRs),np.max(allrandLRs)
LRrandhist,LRBinEdges=np.histogram(allrandLRs,bins=LRBinEdges)
plt.plot(LRBinCentres,LRrandhist,label='rand posns')
plt.legend(loc='best')
plt.show()

#10% false id rate, all LRs, cumulative fraction above certain LR
# the find when 10% of sources have been matched at that LR
LRhistcum=np.array([np.sum(LRhist[i:]) for i in range(len(LRhist))])
LRhistrandcum=np.array([np.sum(LRrandhist[i:]) for i in range(len(LRrandhist))])
LRratio=LRhistrandcum/np.float(LRhistrandcum[0])

plt.plot(LRBinCentres,LRhistcum,label='real posns')
plt.plot(LRBinCentres,LRhistrandcum,label='rand posns')
plt.legend(loc='best')
plt.show()

#false identification rate
plt.plot(LRBinCentres,LRratio)

LRinterp=interp1d(LRratio[::-1],LRBinCentres[::-1]) 
LRlimit=LRinterp(0.1)
   
print LRlimit
    
plt.plot([LRlimit,LRlimit],[0,0.1],'r-')
plt.plot([LRlimit],[0.1],'r*')
plt.show()
