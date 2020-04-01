################################################################
########### ALL CODE GENERATED FROM CAM ########################
################################################################

from astropy.io import fits
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.table as t
from astropy.table import Table, Column
import os
from astropy.visualization import quantity_support
np.seterr(divide='ignore', invalid='ignore')

# to get SRFD for whole galaxy
galaxy_filepaths = '/Users/jalynkrause/Documents/goodmaps/'

def start(fits_file):
    while True:
        try:
            filePath = galaxy_filepaths + '/' + fits_file
            holder = fits.open(filePath)
            break
        except (FileNotFoundError):
            print("No FITS file with name: " + fits_file)
    return holder

def one_start(fits_file):
    while True:
        try:
            filePath = galaxy_filepaths + '/' + fits_file
            holder = fits.open(filePath)
            break
        except (FileNotFoundError):
            print("No FITS file with name: " + fits_file)
            # to get SRFD for whole galaxy
            fits_file = 'manga-' + input("Enter (plate-ifu): ") + '-MAPS-HYB10-GAU-MILESHC.fits.gz'
    return holder, fits_file

def main(hdu):
    lines = hdu['EMLINE_GFLUX'].data
    mask = hdu['EMLINE_GFLUX_MASK'].data
    lines[mask!=0] = np.nan
    ha = lines[18,:,:]
    n2 = lines[19,:,:]
    o3 = lines[13,:,:]
    hb = lines[11,:,:]
    X = np.log10(n2/ha)
    Y = np.log10(o3/hb)
    element_list = [ha,n2,o3,hb,X,Y]
    return element_list

# element_list:   [0] is ha, [1] is n2, [2] is o3, [3] is hb, [4] is log10(n2/ha)
#                 [5] is log10(o3/hb)

def bptCut(element_list):
    cut_list=element_list
    inds=np.where(np.logical_or(cut_list[5] > 0.61/(cut_list[4] -0.05) + 1.3,cut_list[4] > 0))
    (cut_list[4])[inds]=np.nan
    (cut_list[5])[inds]=np.nan
    (cut_list[0])[inds]=np.nan
    (cut_list[1])[inds]=np.nan
    (cut_list[2])[inds]=np.nan
    (cut_list[3])[inds]=np.nan
    return cut_list

def galaxy_plots(mappable,user_choice):
    while True:
        if user_choice=='1':
            plt.figure()
            xk = np.linspace(-.73,-.15, 500)
            yk = 0.61/(xk -0.05) + 1.3
            # Kauffmann+03
            plt.plot(xk, yk, '--', color='blue', lw=2, label='Kauffmann+03' )
            plt.title('BPT with Kauffmann + 03 line')
            plt.xlabel('log([NII]/Ha)')
            plt.ylabel('log[OIII]/Hb')
            X=mappable[4]
            Y=mappable[5]
            plt.scatter(X,Y,c=hdu[2].data[0,:,:],cmap='plasma',s=2,marker='^')
            plt.colorbar()
            plt.show()
            user_choice=input('Choose a plot: 1)BPT 2)log(Ha) 3)log(N[II]/Ha) 4)exit: ')
        elif user_choice=='2':
            plt.figure()
            plt.imshow(np.log10(mappable[0]),origin='lower')
            plt.colorbar()
            plt.show()
            user_choice=input('Choose a plot: 1)BPT 2)log(Ha) 3)log(N[II]/Ha) 4)exit:  ')
        elif user_choice=='3':
            plt.figure()
            plt.imshow(mappable[4],origin='lower')
            plt.colorbar()
            plt.show()
            user_choice=input('Choose a plot: 1)BPT 2)log(Ha) 3)log(N[II]/Ha) 4)exit:  ')
        elif user_choice=='exit':
            break
        elif user_choice=='5':
            plt.figure()
            plt.imshow(mappable,origin='lower')
            plt.colorbar()
            plt.show()
            break
        else:
            print('Not a choice')
            break
            
################################################################################

def dust_correction(h_alpha,h_beta):
    ratio=h_alpha/h_beta
    ha_correct=h_alpha*(ratio/2.86)**2.36
    ha_correct[np.where(ratio>30)]=np.nan
    print(np.ndim(ha_correct))
    return ha_correct
cosmo=FlatLambdaCDM(H0=70,Om0=0.3)

def lum(ha_correct,z_value):
    dl=cosmo.luminosity_distance(z_value)
    print(dl)
    dl=dl.to('cm')/u.cm
    dl=float(dl)
    print('Luminosity Distance: ',dl)
    L_halpha=(ha_correct)*(4*np.pi)*(dl**2)*(10**-17)/10000000
    print('Luminosity: ',L_halpha[22][22])
    sfr=L_halpha/(2.16*(10**34))
    return sfr

# sfr_t=np.nansum(sfr)
# print(sfr_t)

def inclination_correction(z_value,inclination_value):
    arcsecs_per_kpc=cosmo.arcsec_per_kpc_proper(z_value)*u.kpc/u.arcsec
    arcsecs_per_kpc=float(arcsecs_per_kpc)
    inclination=inclination_value
    print('Inclination: ', inclination)
    spax_area=((.5/arcsecs_per_kpc)**2)/inclination
    print('Spaxel Area: ', spax_area)
    return spax_area

def sfr_density(sfr,spaxel_area):
    sfr=sfr/spaxel_area
    return (sfr)

def sfr_total(sfr):

    sfr[np.where(np.isinf(sfr))]=np.nan


    return np.nansum(sfr)

hdu2=fits.open('/Users/jalynkrause/Documents/astro/drpall-v2_4_3.fits')
info=hdu2[1].data
plateifus_drpall=info['plateifu']
zs_drpall=info['nsa_z']
inclinations_drpall=info['nsa_elpetro_ba']
masses_drpall=info['nsa_elpetro_mass']

plateifus_drpall_dict={}
zs_drpall_dict={}
inclinations_drpall_dict={}
mass_drpall_dict={}

for plateifu,z,inclination,mass in zip(plateifus_drpall,zs_drpall,inclinations_drpall,masses_drpall):
    zs_drpall_dict[plateifu]=z
    inclinations_drpall_dict[plateifu]=inclination
    mass_drpall_dict[plateifu]=mass


galaxyd={}
plateifus_goodmaps={}
zs_goodmaps={}
inclinations_goodmaps={}
mass_goodmaps={}

print(len(os.listdir(galaxy_filepaths)))

for galaxy in sorted(os.listdir(galaxy_filepaths)):
    if 'fits' in galaxy:
        galaxyd[galaxy]=galaxy

while True:
    view_choice = input('Work with all galaxies or one galaxy (all/one/exit)? ')

    if view_choice == 'one':
        # to get SRFD for whole galaxy
        galaxy_choice = 'manga-' + input("Enter (plate-ifu): ") + '-MAPS-HYB10-GAU-MILESHC.fits.gz'
        hdu,galaxy_choice = one_start(galaxy_choice)

        plateifu = hdu[0].header['plateifu']
        z_galaxy = zs_drpall_dict[plateifu]
        inclination_galaxy = inclinations_drpall_dict[plateifu]

        cut_element_list = bptCut(main(hdu))
        h_alpha_correct = dust_correction(cut_element_list[0],cut_element_list[3])
        print("")
        sfr = lum(h_alpha_correct,z_galaxy)
        sfrd = sfr_density(sfr,inclination_correction(z_galaxy,inclination_galaxy))
        sfr_t = sfr_total(sfr)

        print('SFR_D[22][22]: ',sfrd[22][22])
        print('SFR_T: ',sfr_t)
        print(z_galaxy)
        print("")
        
        while True:
            choice1=input('Plots with (cuts, no cuts, sfrd) or (exit): ')
            if choice1=='cuts':
                new_hdu,ignore=one_start(galaxy_choice)
                list_wo_cuts=main(new_hdu)
                elements_with_cuts=bptCut(list_wo_cuts)
                choice2=input('Choose a plot: 1)BPT 2)log(Ha) 3)log(N[II]/Ha) 4)exit: ')
                galaxy_plots(elements_with_cuts,choice2)
            elif choice1=='no cuts':
                new_hdu,ignore=one_start(galaxy_choice)
                elements_without_cuts=main(new_hdu)
                choice2=input('Choose a plot: 1)BPT 2)log(Ha) 3)log(N[II]/Ha) 4)exit: ')
                galaxy_plots(elements_without_cuts,choice2)
            elif choice1=='dustc':
                new_hdu,ignore=one_start(galaxy_choice)
                cut_element_list=bptCut(main(new_hdu))
                h_alpha_correct=dust_correction(cut_element_list[0],cut_element_list[3])
                galaxy_plots(h_alpha_correct,'2')
            elif choice1=='sfrd':
                galaxy_plots(sfrd,'5')
            elif choice1=='exit':
                break
            else:
                print('Not a choice')

    elif view_choice=='all':
        # sfrd_all=[]
        plateifu_all=[]
        sfr_t_all=[]
        sfr_t_large=[]
        mass_all=[]
        t2=Table()
        for galaxy in galaxyd:
            t=Table()
            hdu=start(galaxy)

            plateifu=hdu[0].header['plateifu']
            print("")
            print(plateifu)
            plateifu_all.append(plateifu)
            z_galaxy=zs_drpall_dict[plateifu]
            inclination_galaxy=inclinations_drpall_dict[plateifu]
            mass_galaxy=mass_drpall_dict[plateifu]/(.7**2)
            mass_all.append(mass_galaxy)
            list_wo_cuts=main(hdu)
            cut_element_list=bptCut(main(hdu))
            h_alpha_correct=dust_correction(cut_element_list[0],cut_element_list[3])
            sfr=lum(h_alpha_correct,z_galaxy)
            sfrd=sfr_density(sfr,inclination_correction(z_galaxy,inclination_galaxy))
            sfr_t=sfr_total(sfr)
            sfr_t_all.append(sfr_t)
            
            print('SFR_D: ',sfrd[22][22])
            print('SFR_T: ',sfr_t)
            print("")
            # t['SFRD']=Column(sfrd,description='SFR Densities')
            # t['SFRT']=Column(sfr_t,description='Total SFR for the galaxy',length=0)
            # t.write('/Users/admin/Desktop/astro_research/SFRD_and_SFRT/'+plateifu+':SFRD.fits')
        print(len(plateifu_all))
        t2['PLATEIFU']=Column(plateifu_all)
        t2['SFRT']=Column(sfr_t_all)
        t2['MASS70']=Column(mass_all)
        t2.write('/Users/jalynkrause/Documents/astro/MASS70_SFRT.fits')
        # sfr_t_all=np.asarray(sfr_t_all)
        # mass_all=np.asarray(mass_all)
        # t2['MASS']=Column(mass_all,description='Masses')
        # t2['SFRT']=Column(sfr_t_all,description='Total SFR for a galaxy')
        # t2.write('/Users/admin/Desktop/astro_research/Masses_and_SFRT.fits')
        # plt.figure
        # plt.scatter(np.log10(mass_all),np.log10(sfr_t_all))
        # plt.show()
        break

    elif view_choice=='exit':
        break
        
    else:
        print('Not a choice')