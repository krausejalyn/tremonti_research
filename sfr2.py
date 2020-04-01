import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM as flc
cosmo = flc(H0=70,Om0=0.3)

def main():
    # For me the directory is /Users/admin/Desktop/astro_research/SFRD_and_SFRT
    plateifu = '8454-6104'
   
    hdu = fits.open('/Users/jalynkrause/Documents/astro/SFRD_and_SFRT/' + plateifu + '_SFRD.fits')
    # Error on the SFRD
    plateifu = hdu[1].data['plateifu'][0]
    sfrd = hdu[1].data['sfrd']
    sfrt = hdu[1].data['sfrt']

    tot = str(sfrt[0]) # measured in solar masses / year
    
    shape = sfrd.shape[1]
    shapemap = [-.25*shape, .25*shape, -.25*shape, .25*shape]
    
    plt.title(r'SFR H$\alpha$')
    #plt.title('PLATEIFU = '+ plateifu + ', SFR TOTAL = ' + tot )
    plt.xlabel('Arcseconds')
    plt.ylabel('Arcseconds')
    plt.imshow(sfrd, cmap = "viridis", origin = 'lower', extent = shapemap, zorder = 1)
    cb_sfr = plt.colorbar(shrink=0.7)
    cb_sfr.set_label('', rotation=270, labelpad=25)
    plt.show()
    #plt.savefig('/Users/jalynkrause/Documents/astro/SFR/' + plateifu + '(2).png')
    plt.close()

main()
