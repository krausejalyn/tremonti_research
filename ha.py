import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import requests
import matplotlib.image as img

def get_filenames(url):
    file_names = np.genfromtxt(url, usecols = (0), skip_header = 1, dtype = str, delimiter = ',')
    return file_names

def get_ha_map(files):
    hdulist = get_hdu(files)
    logcube = get_logcube(files)
    plate = hdulist['PRIMARY'].header['PLATEID']
    iden = hdulist['PRIMARY'].header['PLATEIFU']
    fiber = hdulist['PRIMARY'].header['IFUDSGN']
    Ha = hdulist['EMLINE_GFLUX'].data[18,...]
    errs = (hdulist['EMLINE_GFLUX_IVAR'].data)**-0.5
    Ha_err = errs[18,...]
    bad = ((Ha/Ha_err) < 3)
    Ha[bad] = np.nan
    eff = hdulist['SPX_ELLCOO'].data[1]
    contours_i = logcube['IIMG'].data
    contours_i[bad] = np.nan
    shape = Ha.shape[1]
    shapemap = [-.25*shape, .25*shape, -.25*shape, .25*shape]
    
    get_image(plate, fiber)
    get_metalicity(iden)
    a = fig.add_subplot(1,3,2)
    
    plt.imshow(Ha, cmap = "viridis", extent = shapemap, zorder = 1)
    plt.gca().invert_yaxis()
    plt.gca().contour(eff*2,[2], extent = shapemap, origin = 'upper', colors = 'r', zorder = 3)
    plt.gca().contour(contours_i, 8, colors = 'black', alpha = 0.6, extent = shapemap, origin = 'upper', zorder = 2)
    cb = plt.colorbar(shrink = .7)
    cb.set_label('H-alpha Flux [10$^{17}$ergs/s/cm$^{2}$/pix]', rotation = 270, labelpad = 15)
    plt.xlabel('Arcseconds')
    plt.ylabel('Arcseconds')
    plt.title('H-alpha Flux')
    #plt.show()
    plt.savefig('/Users/jalynkrause/Documents/astro/halpha_map.png')
    plt.close('all')
    
def get_metalicity(file):
    # all same from prev function or taken from Celeste code
    hdulist = get_hdu(files)
    logcube = get_logcube(files)
    plate = hdulist['PRIMARY'].header['PLATEID']
    fiber = hdulist['PRIMARY'].header['IFUDSGN']
    Ha = hdulist['EMLINE_GFLUX'].data[18,...]
    Hb = hdulist['EMLINE_GFLUX'].data[1,...]
    fluxes = hdulist['EMLINE_GFLUX'].data
    errs = (hdulist['EMLINE_GFLUX_IVAR'].data)**-0.5
    Ha_err = errs[18,...]
    O3 = fluxes[13,:,:]
    O3_err = errs[13,:,:]
    H_beta = fluxes[11,:,:]
    H_beta_err = errs[11,:,:]
    N2 = fluxes[19,:,:]
    N2_err = errs[19,:,:]
    
    # from Pettini paper
    O3_HB = O3/Hb
    N2_HA = N2/Ha
    R = O3_HB/N2_HA	
    logR = np.log10(R)
    # from eq3 in Pettini
    logOH12 = 8.73 - 0.32 * np.log10(R) 
    #bad = ((logOH12err > max_err) | (Ha/Ha_err <3))
    #logOH12[bad] = np.nan
    eff = hdulist['SPX_ELLCOO'].data[1]
    contours_i = logcube['IIMG'].data
    #contours_i[bad] = np.nan
    shape = logOH12.shape[1]
    shapemap = [-.25*shape, .25*shape, -.25*shape, .25*shape]
        
    b = fig.add_subplot(1,3,3)
    
    # plots figure   
    plt.imshow(logOH12, cmap = "viridis", extent = shapemap, zorder = 1)
    plt.gca().invert_yaxis()
    plt.gca().contour(eff*2,[2], extent = shapemap, origin = 'upper', colors = 'r', zorder = 3)
    plt.gca().contour(contours_i, 8, colors = 'black', alpha = 0.6, extent = shapemap, origin = 'upper', zorder = 2)
    plt.title('Metalicity Map')
    plt.xlabel('Arcseconds')
    plt.ylabel('Arcseconds')
    cb = plt.colorbar(shrink = .7)
    cb.set_label('12+log(O/H)', rotation = 270, labelpad = 15)
    #plt.show()
    #plt.savefig('/Users/jalynkrause/Documents/astro/metalicity_map.png')
    #plt.close('all')
                
           
def get_image(plate, fiber):
    r = requests.get('https://data.sdss.org/sas/mangawork/manga/spectro/redux/v2_4_3/' + str(plate) + '/stack/images/' + str(fiber) + '.png', auth=('sdss', '2.5-meters'))

    ##Saves the image
    with open('/Users/jalynkrause/Documents/astro/' + str(plate) + '-' + str(fiber) + '.png', 'wb') as fd:
        for chunk in r.iter_content(chunk_size=128):
            fd.write(chunk)
            
    a = fig.add_subplot(1,3,1)    

    try:
        image = img.imread('/Users/jalynkrause/Documents/astro/' + str(plate) + '-' + str(fiber) + '.png')
    except ValueError:
        print("No image.")
        
    imgplot = plt.imshow(image)    
    
def get_logcube(iden):
    logcube = fits.open('/Users/jalynkrause/Documents/astro/manga-' + str(iden) + '-LOGCUBE.fits.gz')
    return logcube

def get_hdu(iden):
    hdulist = fits.open('/Users/jalynkrause/Documents/astro/manga-' + str(iden) + '-MAPS-HYB10-GAU-MILESHC.fits.gz')
    return hdulist

filename = '/Users/jalynkrause/Documents/astro/filename.txt'
files = get_filenames(filename)
print(files)

fig = plt.figure(figsize = (20,9), facecolor = 'white')

for x in range(0, 1):
    get_ha_map(files)
