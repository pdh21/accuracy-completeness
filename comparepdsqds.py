import numpy as np
import matplotlib.pyplot as plt
import LRmodule as LR
from scipy.interpolate import interp1d

qdsfile="/Users/cc401/Dropbox/thesis/counts/helms/UDSPLWqds.dat"
pdsfile="/Users/cc401/Dropbox/thesis/counts/helms/udsPLWpds.dat"

qdsBinCentres, qds=np.loadtxt(qdsfile, unpack=True)
pdsBinCentres, pds=np.loadtxt(pdsfile,unpack=True)

radbinedges   =np.arange(0,51)
fluxbinedges  =np.arange(-42,22,2)

#get binCentres
radBinCentres=LR.linBinCentres(radbinedges)
fluxBinCentres=LR.linBinCentres(fluxbinedges)


plt.figure()
plt.plot(fluxBinCentres, qds)
plt.plot(pdsBinCentres, pds*1000)
plt.show()

#rebin the pds thing 
pdsBinEdges=LR.linBinEdges(pdsBinCentres)
pds=LR.rebinDataBins(pds,pdsBinCentres, pdsBinEdges, fluxbinedges)

pdsInterp =interp1d(fluxBinCentres,pds)
qdsInterp =interp1d(fluxBinCentres,qds)


plt.figure()
plt.plot(fluxBinCentres,qds/pds)
print qdsBinCentres,fluxBinCentres
plt.figure()
plt.plot(fluxBinCentres, qds)
plt.figure()
plt.plot(fluxBinCentres, pds)

plt.show()