from astropy.io import fits
import numpy as np
from astropy.cosmology import FlatLambdaCDM as flc
cosmo = flc(H0=70,Om0=0.3)

p3ddir='/Users/jalynkrause/Documents/astro/pipe3d_maps/'
drpall='/Users/jalynkrause/Documents/astro/drpall-v2_4_3.fits'


def get_redshift(plateifu):
    '''
    Get the redshift for each galaxy
    '''
    hdu = fits.open(drpall)
    indx = hdu[1].data['PLATEIFU'] == plateifu
    return hdu[1].data['NSA_Z'][indx][0]

def get_ba(plateifu):
    '''
    Get the b/a ratio for each galaxy
    '''
    hdu = fits.open(drpall)
    indx = hdu[1].data['PLATEIFU'] == plateifu
    return hdu[1].data['nsa_elpetro_ba'][indx][0]

def get_massmap(plateifu):
    #Open the file and extract the data
    hdu = fits.open(p3ddir + 'manga-' + plateifu + '.Pipe3D.cube.fits.gz')
    mass = hdu[1].data[19,:,:] #dust-corrected, log10(M*/spaxel)
    mass_err = hdu[1].data[20,:,:]
    #convert from mass/spaxel to mass/kpc^2 with an inclination correction
    z = get_redshift(plateifu)
    ba = get_ba(plateifu)
    sig_cor = 2*np.log10((0.5/cosmo.arcsec_per_kpc_proper(z).value))
    mass = mass-sig_cor + np.log10(ba)
    mass_err = mass_err-sig_cor + np.log10(ba)
    return mass, mass_err


mass,mass_error = get_massmap('7495-12704')
