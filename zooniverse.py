import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import math
import scipy.stats.distributions as dist # from Ewan Cameron's paper

# fraction of galaxies that have a symmetric O/H map as a function of stellar mass
def fits_v0():
    hdu_v0 = fits.open('zooinverse_summary_v0.fit')
    plateifu = hdu_v0[1].data['PLATEIFU']
    sym_OH = hdu_v0[1].data['SYMMETRIC_OH']
    log_mass = hdu_v0[1].data['LOG_MASS']
    data = hdu_v0[1].data

    bins = np.linspace(9, 11, 5) # sets up the edges of the mass bins
    bins_cens = np.linspace(9.25, 10.75, 4) # centers of the bins. This should be done in a smarter way if you have arbitrary bin edges

    frac = np.zeros(len(bins_cens)) # returns a new array filled with zeros

    for i in range(0, len(bins) - 1):
        a = np.where((log_mass>bins[i]) & (log_mass<bins[i+1]))[0] #Need the '[0]' because np.where returns a tuple with an array of indices.
        total = len(a) # total number of galaxies in the bin
        b = np.where((log_mass>bins[i]) & (log_mass<bins[i+1]) & (sym_OH>0))[0]
        n_sym = len(b) # total number of symmetric galaxies in bin
        frac[i] = n_sym/total
        
    fits_v2() # calls other definition
    
    # plotted figure of symmetric O/H map & stellar mass
    plt.plot(bins_cens, frac, label = 'Summary v0 Fits File', linewidth = 3)
    plt.scatter(bins_cens, frac, lw = 3)
    plt.legend(loc = 'upper left')
    plt.title('Fraction of Galaxies with Symmetric O/H Maps')
    plt.ylabel('Fraction')
    plt.xlabel('$log_{10}$ Stellar Mass')
    plt.show()
    #plt.savefig('/Users/jalynkrause/Documents/astro/fracOH_stellarmass.png')
    plt.close()
    
# newest fits file fraction of symmetric O/H map vs. stellar mass
def fits_v2():
    hdu_v2 = fits.open('zooinverse_summary_v2.fit')
    plateifu = hdu_v2[1].data['PLATEIFU']
    sym_OH = hdu_v2[1].data['SYMMETRIC_OH']
    log_mass = hdu_v2[1].data['LOG_MASS']
    data = hdu_v2[1].data

    bins = np.linspace(9, 11, 5) # sets up the edges of the mass bins
    bins_cens = np.linspace(9.25, 10.75, 4) # centers of the bins. This should be done in a smarter way if you have arbitrary bin edges

    frac = np.zeros(len(bins_cens))  # returns a new array filled with zeros

    for i in range(0, len(bins) - 1):
        a = np.where((log_mass>bins[i]) & (log_mass<bins[i+1]))[0] #Need the '[0]' because np.where returns a tuple with an array of indices.
        total = len(a) # total number of galaxies in the bin
        b = np.where((log_mass>bins[i]) & (log_mass<bins[i+1]) & (sym_OH>0))[0]
        n_sym = len(b) # total number of symmetric galaxies in bin
        frac[i] = n_sym/total

    # plotted line for fits_v2 file
    plt.plot(bins_cens, frac, label = 'Summary v2 Fits File', linewidth = 3)
    plt.scatter(bins_cens, frac, lw = 3)
    
'''
##print(hdu[1].header)

# from Ewan Cameron's paper
p_lower = dist.beta.ppf((1-c)/2.,k+1,n-k+1)
p_upper = dist.beta.ppf(1-(1-c)/2.,k+1,n-k+1)
'''

# main method
fig = plt.figure(figsize = (20,9), facecolor = 'white')

fits_v0()