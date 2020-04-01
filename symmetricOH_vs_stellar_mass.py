import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import bootstrap
import numpy as np
import math
import statistics
import random
import scipy.stats.distributions as dist # from Ewan Cameron's paper
    
# newest fits file fraction of symmetric O/H map vs. stellar mass
def fits_v2():
    hdu_v2 = fits.open('zooinverse_summary_v2.fit')
    plateifu = hdu_v2[1].data['PLATEIFU']
    sym_OH = hdu_v2[1].data['SYMMETRIC_OH']
    log_mass = hdu_v2[1].data['LOG_MASS']
    data = hdu_v2[1].data

    bins = np.linspace(9, 11, 9) # sets up the edges of the mass bins
    bins_cens = np.linspace(9.125, 10.875, 8) # centers of the bins. This should be done in a smarter way if you have arbitrary bin edges

    frac = np.zeros(len(bins_cens))  # returns a new array filled with zeros
    c = np.zeros(len(bins_cens))
    
    err = np.zeros(len(bins_cens))
    
    for i in range(0, len(bins) - 1): 
        a = np.where((log_mass > bins[i]) & (log_mass < bins[i+1]))[0] #Need the '[0]' because np.where returns a tuple with an array of indices.
        total = len(a) # total number of galaxies in the bin
        b = np.sum(sym_OH[np.where((log_mass > bins[i]) & (log_mass < bins[i + 1]) & (sym_OH > 0))])
        d = sym_OH[np.where((log_mass > bins[i]) & (log_mass < bins[i + 1]) & (sym_OH > 0))]
        boot = bootstrap(d, bootnum=1000, samples=None, bootfunc=None)
        boot_frac = np.zeros(1000)
        j = 0
        frac[i] =  b / total 
        for row in boot:
            c = np.sum(row)
            boot_frac[j] = c / total 
            j += 1
        err[i] = statistics.stdev(boot_frac, xbar=None)  # xbar=none means that the mean is automatically calculated   
         
    # plotted line for fits_v2 file
    plt.plot(bins_cens, frac, label='Summary v2 Fits File', linewidth=3)
    plt.scatter(bins_cens, frac, lw = 3)
    plt.errorbar(bins_cens, frac, err) # from before (np.sqrt(c)/760)
    plt.title('Fraction of Galaxies with Symmetric O/H Maps')
    plt.ylabel('Fraction')
    plt.xlabel('$log_{10}$ Stellar Mass')
    plt.show()
    #plt.savefig('/Users/jalynkrause/Documents/astro/fracOH_stellarmass.png')
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
        bootarr = np.random.randint(low=0, high=d.shape[0], size=samples)
        boot[i] = d[bootarr]
        
    return boot # Returns numpy.ndarray Bootstrapped data. Each row is a bootstrap resample of the data.    
    
'''
##print(hdu[1].header)

# from Ewan Cameron's paper
p_lower = dist.beta.ppf((1-c)/2.,k+1,n-k+1)
p_upper = dist.beta.ppf(1-(1-c)/2.,k+1,n-k+1)
'''

# main method
fig = plt.figure(figsize = (20,9), facecolor = 'white')

fits_v2()