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

fluxi= catin["flux"]
fluxo=catout["flux"]

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

#find the set of indices that are in matches, not repeating.
flatmatches=np.array([x for xs in matches for x in xs])
sourcesMatched=np.array(list(set(flatmatches)),dtype='int')

sBinEdges=10**np.linspace(np.log10(np.min(fluxi)),np.log10(np.max(fluxi)), 15, endpoint=True)
sIn=LR.logBinCentres(sBinEdges)
print sIn
print np.max(fluxi)

histIn,junk=np.histogram(fluxi,bins=sBinEdges)
histMatched,junk=np.histogram(fluxi[sourcesMatched],bins=sBinEdges)

plt.plot(histMatched/histIn)
plt.show()

def compCurve(sIn, Q,M,B,N):
    return 1./((1.+Q*np.exp((-sIn-M)/B))**N)

params,pcov=curve_fit(compCurve,sIn,histMatched/histIn)
Q,M,B,N=params
plt.plot(sIn,compCurve(sIn, Q,M,B,N))
plt.show()

