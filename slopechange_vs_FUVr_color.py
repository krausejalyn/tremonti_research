import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import bootstrap
import numpy as np
import math
import statistics
import random
import scipy.stats.distributions as dist # from Ewan Cameron's paper

# fraction of galaxies where the metallicity profile changes slopes vs FUV - r color for fits file v2
def hdu():
    hdu = fits.open('zooinverse_summary_v2.fit')
    log_mass = hdu[1].data['LOG_MASS']
    slope_change = hdu[1].data['SLOPE_CHANGE'] # fraction of galaxies with slope changes as a function of stellar mass
    FUV_r = hdu[1].data['FNUGRIZ_ABSMAG'] # array of 7 colors in this order:  FUV, NUV, u, g, r, i, z

    fuvr = FUV_r[:,0] - FUV_r[:,4]
    
    nbins = 10
    binsize = np.linspace(np.nanmin(fuvr), np.nanmax(fuvr), nbins)
    bins_cens = np.zeros(nbins)

    for j in range(0, nbins-1):
        bins_cens[j] = ((binsize[j] + binsize[j+1]) / 2.0)
    
    frac = np.zeros(len(bins_cens)) # returns a new array filled with zeros
    c = np.zeros(len(bins_cens))
    err = np.zeros(len(bins_cens))

    for i in range(0, len(binsize)-1):
        in_bin = np.where((fuvr > binsize[i]) & (fuvr < binsize[i+1]))[0] # Need the '[0]' because np.where returns a tuple with an array of indices.
        total = len(in_bin) # total number of galaxies in the bin
        b = np.sum(slope_change[np.where((fuvr > binsize[i]) & (fuvr < binsize[i+1]) & (slope_change > 0))])
        d = slope_change[np.where((fuvr > binsize[i]) & (fuvr < binsize[i + 1]) & (slope_change > 0))]
        boot = bootstrap(d, bootnum=1000, samples=None, bootfunc=None) 
        boot_frac = np.zeros(1000)
        j = 0
        frac[i] =  b / total
        for row in boot:
            c = np.sum(row)
            boot_frac[j] = c / total 
            j+=1
        err[i] = statistics.stdev(boot_frac, xbar=None)  # xbar=none means that the mean is automatically calculated 
        
    # plotted figure of fraction of galaxies with slope changes vs. FUV-r color
    #plt.plot(binsize, frac, linewidth = 3)
    #plt.scatter(binsize, frac)
    plt.hist(fuvr, bins=nbins, facecolor='blue')
    plt.errorbar(bins_cens, frac, err) # divide by 760 because that is the total amount of galaxies
    plt.title('Fraction of Galaxies that have a Slope Change vs. FUV-r Color')
    plt.xlabel('FUV-r Color')
    plt.ylabel('Fraction')
    plt.show()
    #plt.savefig('/Users/jalynkrause/Documents/astro/frac_slope_FUV-r.png')
    plt.close()
    
def bootstrap(d, bootnum=1000, samples=None, bootfunc=None): #(data, bootnum=100, samples=None, bootfunc=None)
    # from http://docs.astropy.org/en/v0.3/_modules/astropy/stats/funcs.html#bootstrap
    if samples is None:
        samples = d.shape[0]
    
    if bootfunc is None:
        resultdims = (bootnum,) + (samples,) + d.shape[1:]
        boot = np.empty(resultdims)
    else:
        resultdims = (bootnum,)
        boot = np.empty(resultdims)
    
    for i in range(0, bootnum):
        bootarr = np.random.randint(low=-1, high=d.shape[0], size=samples)
        boot[i] = d[bootarr]
        
    return boot # Returns numpy.ndarray Bootstrapped data. Each row is a bootstrap resample of the data.    

'''
print(hdu[1].header)

# from Ewan Cameron's paper
p_lower = dist.beta.ppf((1-c)/2.,k+1,n-k+1)
p_upper = dist.beta.ppf(1-(1-c)/2.,k+1,n-k+1)
'''

# main method
fig = plt.figure(figsize = (20,9), facecolor = 'white')

hdu()