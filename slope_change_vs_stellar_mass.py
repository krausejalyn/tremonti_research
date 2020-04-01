import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import math
import scipy.stats.distributions as dist # from Ewan Cameron's paper

# fraction of galaxies where the metallicity profile changes slopes vs stellar mass for first fits file (v0)
def fits_v0():
    hdu_v0 = fits.open('zooinverse_summary_v0.fit')
    plateifu = hdu_v0[1].data['PLATEIFU']
    slope_change = hdu_v0[1].data['SLOPE_CHANGE']
    log_mass = hdu_v0[1].data['LOG_MASS']
    data = hdu_v0[1].data

    bins = np.linspace(9, 11, 9) # sets up the edges of the mass bins
    bins_cens = np.linspace(9.125, 10.875, 8) # centers of the bins. This should be done in a smarter way if you have arbitrary bin edges

    fract = np.zeros(len(bins_cens)) # returns a new array filled with zeros
        
    fits_v2() # calls other definition
     
    for j in range(0, len(bins) - 1):
        a = np.where((log_mass>bins[j]) & (log_mass<bins[j+1]))[0] #Need the '[0]' because np.where returns a tuple with an array of indices.
        tot = len(a) # total number of galaxies in the bin
        b = np.where((log_mass>bins[j]) & (log_mass<bins[j+1]) & (slope_change>0))[0]
        n_change = len(b) # total number of symmetric galaxies in bin
        fract[j] = n_change/tot

    # plotted figure of fraction of galaxies with slope changes & stellar mass
    plt.plot(bins_cens, fract, label='Summary v0 Fits File', linewidth = 3)
    plt.scatter(bins_cens, fract, lw = 3)
    plt.legend(loc = 'upper left')
    plt.title('Fraction of Galaxies with Slope Changes vs. Stellar Mass')
    plt.xlabel('$log_{10}$ Stellar Mass')
    plt.ylabel('Fraction')
    plt.show()
    ##plt.savefig('/Users/jalynkrause/Documents/astro/frac_slope_stellarmass.png')
    plt.close()
    
# fraction of galaxies where the metallicity profile changes slopes vs stellar mass from fits file (v2)
def fits_v2():
    hdu_v2 = fits.open('zooinverse_summary_v2.fit')
    plateifu = hdu_v2[1].data['PLATEIFU']
    sym_OH = hdu_v2[1].data['SYMMETRIC_OH']
    log_mass = hdu_v2[1].data['LOG_MASS']
    data = hdu_v2[1].data

    bins = np.linspace(9, 11, 9) # sets up the edges of the mass bins
    bins_cens = np.linspace(9.125, 10.875, 8) # centers of the bins. This should be done in a smarter way if you have arbitrary bin edges

    fract = np.zeros(len(bins_cens))  # returns a new array filled with zeros
    
    # fraction of galaxies with slope changes as a function of stellar mass
    slope_change = hdu_v2[1].data['SLOPE_CHANGE']
    
    fract = np.zeros(len(bins_cens))
    
    for j in range(0, len(bins) - 1):
        a = np.where((log_mass>bins[j]) & (log_mass<bins[j+1]))[0] #Need the '[0]' because np.where returns a tuple with an array of indices.
        tot = len(a) # total number of galaxies in the bin
        b = np.where((log_mass>bins[j]) & (log_mass<bins[j+1]) & (slope_change>0))[0]
        n_change = len(b) # total number of symmetric galaxies in bin
        fract[j] = n_change/tot

    # plotted figure of fraction of galaxies with slope changes & stellar mass
    plt.plot(bins_cens, fract, label='Summary v2 Fits File', linewidth = 3)
    plt.scatter(bins_cens, fract, lw = 3)
    
    
    
'''
##print(hdu[1].header)

# from Ewan Cameron's paper
p_lower = dist.beta.ppf((1-c)/2.,k+1,n-k+1)
p_upper = dist.beta.ppf(1-(1-c)/2.,k+1,n-k+1)
'''

# main method
fig = plt.figure(figsize = (20,9), facecolor = 'white')

fits_v0()