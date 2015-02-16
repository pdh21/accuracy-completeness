import numpy as np
import astropy.io.fits as apif
import matplotlib.pyplot as plt
import astropy.wcs as apw
import numpy.random as npr
from scipy.spatial import KDTree
from scipy.optimize import curve_fit

'''set of codes to use when calculating likelyhood ratios. Doesn't really need to be in a separate module but does make things easier in the long run if I need to seperate things'''


def adxy(header, ra, dec):
    '''give it a header and some ra and dec and it'll return the x and y without you having to think too hard'''
    WCSmap=apw.WCS(header)
    x,y=WCSmap.wcs_world2pix(ra,dec,0)
    return x,y


def linBinCentres(edges):
    '''given linear spaced bins, calculates the bin centres, (length -1 of bin edges)'''
    return edges[:-1]+0.5*(edges[1:]-edges[:-1])


def logBinCentres(edges):
    '''given linear spaced bins, calculates the bin centres, (length -1 of bin edges)'''
    return 10**(np.log10(edges[:-1])+0.5*(np.log10(edges[1:])-np.log10(edges[:-1])))


def linBinEdges(centres):
    '''given linear spaced bin centres, calculates the bin edges (length +1 of bin centres)'''
    edges=np.zeros(len(centres)+1)
    edges[0]=centres[0]-0.5*(edges[1]-edges[0])
    edges[1:]=centres[0]+0.5*(edges[1]-edges[0])
    return edges


def getkdTree(xi,yi):
    '''shortcut of having to remember the correct notation
        '''
    return KDTree(np.array([xi,yi]).T,leafsize=100)

def rayleigh(r,sigma): 
    return r*np.exp(-r**2/(2*sigma**2))/(sigma**2.)

def rayleighScaled(r,sigma,scaling): 
    return scaling*r*np.exp(-r**2/(2*sigma**2))/(sigma**2.)

def rayleighPlus(r,sigma,scaling, scaling2):
    return scaling*r*np.exp(-r**2/(2*sigma**2))/(sigma**2.) + scaling2*r

def histDistance(xi,yi,xo,yo,bins,distupper=None):
    '''calculates a kdTree and gives the normalised (by number of xo yo sources) counts of matches as a function of radius. xs,ys,bins, and distupper must all be in the same distance units.
        NB distupper is max dist the Tree query stops at, defaults to the last bin edge you give. Making it bigger than your bins won't have any effect on the output but making it smaller will.'''
    
    histtot=np.zeros(len(bins)-1)
    if distupper == None:
        distupper=bins[-1]
    Tree=KDTree(np.array([xi,yi]).T,leafsize=100)
    
    for i in range(len(xo)):
        if i%100==0: print i," out of ",len(xo)
        d,ind=Tree.query(np.array([xo[i],yo[i]]),k=None,distance_upper_bound=distupper)
        hist,junk=np.histogram(d,bins=bins)
        #if i%100==0:
        #plt.subplot(2,1,1)
        #plt.plot((xi[ind]-xo[i])*pixsize, (yi[ind]-yo[i])*pixsize, 'bx')
        #plt.subplot(2,1,2)
        #plt.plot(radbinedges[:-1]+0.5*(radbinedges[1:]-radbinedges[:-1]),hist)
        #plt.show()
        histtot+=hist
    
    #find mean
    histtot/=np.float(len(xo))
    #need to make a poisson error thing
    histtoterr=histtot/np.sqrt(len(xo)) 
    
    return histtot,histtoterr


def randomXYsInMap(map,x,y):
    '''Create a list of randomised x and y positions within a given map and header, where NaNs denote outside the map.
        '''
    
    stilltorandom=np.ones(len(x))
    xorand=x.copy()
    yorand=y.copy()
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
    return xorand, yorand

def randomXYsInMap2(map,x,y):
    '''Swaps around the x and y positions rather than randomising in the map'''
    ind=np.arange(len(x))
    npr.shuffle(ind)
    #calculate perturbation
    print (npr.random(len(x))-0.5)*10*31./8.3333
    x+=(npr.random(len(x))-0.5)*10*31./8.3333
    y+=(npr.random(len(y))-0.5)*10*31./8.3333
    return x[ind],y[ind]


def findSigmaR(bincentres,data,error=None):
    '''fits Rayleigh function + bg noise linear in r to the counts and optional relative errors. Returns the fits as sigma + scaling, returns the covariance matrix for those parameters respectively.'''

    params,pcov=curve_fit(rayleighPlus, bincentres,data,sigma=error)
    #params,pcov=curve_fit(rayleighPlus, bincentres,data)#,sigma=error)
    
    print "sigma: ",params[0]," arcsec p/m",np.sqrt(pcov[0,0])," arcsec"
    print "excess scale: ",params[1]," scaling p/m",np.sqrt(pcov[1,1])
    print "random scale: ",params[2]," scaling p/m",np.sqrt(pcov[2,2])

    return params[0], params[1], params[2], pcov


def findSourcesInRadius(xo,yo,xi,yi,radius,kdTree=None):
    '''xo yo are scalars.
        Find sources within a search radius. All distances/positions must be in the same units. Optionally can pass in a kdTree so it doesn't make one each function call.'''
    if kdTree == None:
        Tree=KDTree(np.array([xi,yi]).T,leafsize=100)
    else: Tree= kdTree
    d,ind=Tree.query(np.array([xo,yo]),k=None,distance_upper_bound=radius)

    return ind



def rebinDataBins(data,binCentresToBin, binEdgesToBin, binEdgesToUse):
    '''given bin CENTRES that you have and bin EDGES that you want, code takes data and rebins it to the new bin sizes. Works on 1D bins only and when the new binsize is a multiple of the old one, so it's really basic'''
    
    #un-normalise data:
    data*=(binCentresToBin[1]-binCentresToBin[0])
    dataout=np.zeros(len(binEdgesToUse)-1)
    for i in range(len(binEdgesToUse)-1):
        #look for where the bin centres are in the new bins
        ind, =np.where((binCentresToBin > binEdgesToUse[i]) & (binCentresToBin< binEdgesToUse[i+1]))
        dataout[i]=np.sum(data[ind])

    #re-normalise data:
    dataout/=(binEdgesToUse[1:]-binEdgesToUse[:-1])
    return dataout



def findOptimalRadius(bincentres, sigma, scaling, scaling2, tolerance=1e-4,plot=False):
    '''given data/full model information, construct where the optimum radius is.
    This is defined as the Signal to Noise Ratio, SNR, (excess of counts)/(poisson noise in excess_. Excess counts are the number of counts above what's expected from background counts alone, i.e. counts as fn of radius - model of bg counts. Poisson noise in the excess is sqrt(excess counts) as it's shot noise. 
        
        scaling must always be given, it's how the bg counts scale with r. if sigma is given as a scalar, then the model is used, if it's an array then that's assumed to be all the counts (bg + excess)
        
        tolerance is optional and is the difference between neighbouring SNR calcuatiions at which to stop searching. If set to 0 will find the max SNR value which may be a silly radius so be careful. Optimal search radius will be massively sensitive to this but not q(dS) (I hope)'''

    #checking if sigma's a scalar or an array.
    if ~hasattr(sigma, "__len__"):
        histtotmodel=rayleighPlus(bincentres,sigma,scaling,scaling2)
    else: 
        histtotmodel=sigma
    #cumulatively sum the counts
    histtotcum = np.array([np.sum(histtotmodel[0:i]) for i in range(len(histtotmodel))])

    #make the background counts and cumulatively sum them
    bgcounts = scaling2*bincentres
    bgcountscum= np.array([np.sum(bgcounts[0:i]) for i in range(len(bgcounts))])
    #SNR signalcounts/sqrt(signalcounts)
    SNR=(2*(histtotcum-bgcountscum))/np.sqrt(2*histtotcum)
    for ind in range(len(SNR)):
        diff=SNR[ind+1]-SNR[ind]
        if diff < tolerance: break
    print "Optimal radius is ", bincentres[ind], " with SNR: ", SNR[ind], " at index ", ind

    plt.plot(bincentres, histtotcum)
    plt.plot(bincentres, bgcountscum)
    plt.plot(bincentres, SNR**2)
    plt.show()
    if plot==True:
        plt.plot(bincentres, SNR, 'k-', label="SNR(r)")
        plt.plot([bincentres[ind],bincentres[ind]],[0,SNR[ind]],'k--')
        plt.plot(bincentres[ind],SNR[ind], 'r*',label="optimal radius", markersize=20)
        plt.xlabel("Radius")
        plt.ylabel("Signal to Noise Ratio")
        plt.legend(loc="best")
        plt.show()

    return bincentres[ind], SNR[ind],ind

    