import numpy as np
import astropy.io.fits as apif
import matplotlib.pyplot as plt
import LRmodule as LR
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

''' read in catalogues and read in the matches list, re pick them out and then go through and analyse the results by finding the accuracy of the methods'''

catfilein ="/Users/cc401/datascripts/fits/laceyfilt/uds_catalog_clustering_bias_1.0_cut_1000mJy_0_PLW.fits"
catfileout="/Users/cc401/datascripts/fits/laceyfilt/uds_itermap_simulation_clustering_bias_1.0_cut_1000mJy_w_noise_0_PLW_SXT.fits"
matchfile="/Users/cc401/Dropbox/thesis/counts/helms/UDSmatchesPLW.dat"


catin  = apif.open(catfilein)[1].data
catout = apif.open(catfileout)[1].data

rai,dei,fluxi=( catin["ra"], catin["dec"], catin["flux"])
rao,deo,fluxo=(catout["ra"],catout["dec"],catout["flux"])

ids=[]
matches=[]
#read in the match file and get an id and match list
with open(matchfile) as f:
    content=f.readlines()
    for line in content:
        line=line.split('\t')
        ids.append(line[0])
        line1=line[1][:-2].split(',')
        if len(line1) > 1:
            #print [float(x) for x in line1]
            matches.append([float(x) for x in line1])
        else: matches.append([])
f.close()

ids=np.array([float(x) for x in ids])
matches=np.array(matches)


sInarr=[]
sOutarr=[]

for id in ids:
    if len(matches[id]) != 0:
        sIn=fluxi[matches[id]]
        sOut=fluxo[id]
        sInarr.append(sIn)
        sOutarr.append(np.ones(len(sIn))*sOut)
#linearise
sInarr= np.array([x1 for xs in sInarr for x1 in xs])
sOutarr=np.array([x1 for xs in sOutarr for x1 in xs])

plt.plot(sInarr, (sOutarr-sInarr)/sInarr,'r.')

#bin the sIn and sOut along sIn in log
#sBinEdges=np.linspace(np.log10(1),np.log10(500), 12, endpoint=True)
sBinEdges=10**np.linspace(np.log10(np.min(sInarr)),np.log10(np.max(sInarr)), 12, endpoint=True)
sIn=LR.logBinCentres(sBinEdges)

#bin sOut according to sIn bins, find average flux of sOut and err
sOutmean=np.zeros(len(sIn))
sOuterr=np.zeros(len(sIn))


#need to normalise the histogrammmmmmmmm??????
for i in range(len(sIn)):
    ind,=np.where((sInarr<sBinEdges[i+1]) & (sInarr>=sBinEdges[i]))
    if len(ind)>0:
        sOutmean[i]=np.mean(sOutarr[ind])
        sOuterr[i]=np.std(sOutarr[ind])
plt.errorbar(sIn,(sOutmean-sIn)/sIn,yerr=sOuterr/sIn,color='b',linestyle='--')
plt.xscale('log')
plt.yscale('linear',nonposy='clip')


def accLine(sIn, a0,a1,a2):
    return a0*sIn**a1 +a2

params,pcov=curve_fit(accLine,sIn,(sOutmean-sIn)/sIn,sigma=sOuterr/sIn)
a0,a1,a2=params
a0e,a1e,a2e=(pcov[0,0],pcov[1,1],pcov[2,2])

print 'a0: '+str(a0)+' +/- '+str(a0e)
print 'a1: '+str(a1)+' +/- '+str(a1e)
print 'a2: '+str(a2)+' +/- '+str(a2e)
plt.plot(sBinEdges, accLine(sBinEdges,a0,a1,a2),'g:')
plt.show()
