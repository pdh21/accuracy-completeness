import numpy as np
import astropy.io.fits as apif
import matplotlib.pyplot as plt
import LRmodule as LR
from scipy.interpolate import interp1d
import LRfiles as LRf

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica Neue']})
rc('font',**{'family':'serif'})
rc('xtick',labelsize=15)
rc('ytick',labelsize=15)
rc('axes',labelsize=15)
#rc('ytick.major',pad=5)
#rc('text',usetex=True)

'''Given all the derived info return a list of matches of output source to input sources as a .dat file, then this can be used when read in with the two catalogues to match them up to then go on and calculate things~~'''

#catfilein ="/Users/cc401/datascripts/fits/laceyfilt/uds_catalog_clustering_bias_1.0_cut_1000mJy_0_PLW.fits"
#catfileout="/Users/cc401/datascripts/fits/laceyfilt/uds_itermap_simulation_clustering_bias_1.0_cut_1000mJy_w_noise_0_PLW_SXT.fits"
#associatedmapfile ="/Users/cc401/datascripts/fits/laceyfilt/uds_itermap_simulation_clustering_bias_1.0_cut_1000mJy_w_noise_0_PLW.fits"


#qdsfile="/Users/cc401/Dropbox/thesis/counts/helms/UDSPLWqds.dat"
#pdsfile="/Users/cc401/Dropbox/thesis/counts/helms/udsPLWpds.dat"
#optradius=49.#23.
#sigma=12.055#5.4601
#LRThresh=2.93#3.53

catfilein=LRf.catfilein
catfileout=LRf.catfileout
associatedmapfile=LRf.associatedmapfile
qdsfile=LRf.qdsfile
pdsfile=LRf.pdsfile
optradius=LRf.radius
sigma=LRf.sigma
LRThresh=LRf.LRThresh

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
frInterp  =interp1d(np.insert(radBinCentres,0,0.0),np.insert(LR.rayleigh(radBinCentres,sigma),0,0.0000000)) #adding 0,0001 as bincentres for the interpolation too as we know it goes to 0 there



print "calculating LRs"
allLRs=[]
allind=""
allrs=""
allfs=""
for i in range(len(xo)):
    if i%100 == 0: print "source ",i," out of ",len(xo)
    #find the sources that are within the optimal search range only
    r,rInd=Tree.query(np.array([xo[i],yo[i]]),k=None,distance_upper_bound=optradius/pixsize)
    rInd=np.array(rInd)
    if len(rInd) > 0:
        r=np.array(r*pixsize)

        #ind=LR.findSourcesInRadius(xo[i],yo[i],xi,yi,optradius/pixsize,kdTree=Tree)
        #find sources in the flux range:
        dFlux=fluxi[rInd]-fluxo[i]
        fInd,=np.where((dFlux <= fluxBinCentres[-1]) & (dFlux >= fluxBinCentres[0]))
        if len(fInd) > 0:
            gdxi=xi[rInd[fInd]]
            gdyi=yi[rInd[fInd]]
            dFlux=dFlux[fInd]
            r=r[fInd]
            
            
            LikeRatios=(frInterp(r)*qdsInterp(dFlux))/(pdsInterp(dFlux)*np.pi*2*r)
            ind,=np.where(LikeRatios >= LRThresh)
            if len(ind) == 0:
                allLRs.append(',')
                allind=allind+str(i)+"\t,\n"
                allrs=allrs+str(i)+"\t,\n"
                allfs=allfs+str(i)+"\t,\n"
            else:
                allLRs.append(list(LikeRatios[ind]))
                allind=allind+str(i)+'\t'+str(list(rInd[fInd[ind]])).strip('[]')+'\n'
                #print np.array([qdsInterp(dFlux[ind]), pdsInterp(dFlux[ind]), frInterp(r[ind]), r[ind]]).T
                allrs=allrs+str(i)+'\t'+str(list(r[ind])).strip('[]')+'\n'
                allfs=allfs+str(i)+'\t'+str(list(dFlux[ind])).strip('[]')+'\n'
        else:
            allLRs.append(',')
            allind=allind+str(i)+"\t,\n"
            allrs=allrs+str(i)+"\t,\n"
            allfs=allfs+str(i)+"\t,\n"
    else:
        allLRs.append(',')
        allind=allind+str(i)+"\t,\n"
        allrs=allrs+str(i)+"\t,\n"
        allfs=allfs+str(i)+"\t,\n"

#print str(allind).replace('\'','').strip('[]')
buffOut=open("/Users/cc401/Dropbox/thesis/counts/helms/cosmosmatchesPLW.dat","w")
delme=open("/Users/cc401/Dropbox/thesis/counts/helms/delme.dat","w")
delmef=open("/Users/cc401/Dropbox/thesis/counts/helms/delmef.dat","w")
#np.savetxt("/Users/cc401/Dropbox/thesis/counts/helms/UDSmatchesPLW.dat",allind)

buffOut.write(str(allind).replace('\'','').strip('[]'))
delme.write(str(allrs).replace('\'','').strip('[]'))
delmef.write(str(allfs).replace('\'','').strip('[]'))
buffOut.close()
delme.close()
delmef.close()



