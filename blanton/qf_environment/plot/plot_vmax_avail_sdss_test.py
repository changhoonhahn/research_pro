import numpy as np
import pyfits as fits
import matplotlib.pyplot as plt
import matplotlib.font_manager
from matplotlib import rc
import scipy.optimize as optimize
from scipy import interpolate
import matplotlib.pylab as py
import yanny as yan
import sys

rc('text', usetex=True)
rc('font', family='serif')
sdss_run    = int(sys.argv[1])
nran        = int(sys.argv[2])
nrsk        = int(sys.argv[3])
par = yan.yanny('/home/users/hahn/qf_environment/zerod_environment_parameters.par',np=True)

rad = str(par['PARAMETERS']['cylrad'][sdss_run-1])
height = str(par['PARAMETERS']['cylheight'][sdss_run-1])
thresh = str(par['PARAMETERS']['thresh'][sdss_run-1])
nbin = str(par['PARAMETERS']['nbin'][sdss_run-1])

vmax_avail_dir  = '/global/data/scr/chh327/primus/data/vmax_avail/'
ransack_dir     = '/global/data/scr/chh327/primus/data/ransack/'
envcount_dir    = '/global/data/scr/chh327/primus/data/envcount/'

random_vmax_name = 'vmax_avail_cylr'+rad+'h'+height+'_thresh'+thresh+'_nbin'+nbin+'_sdss_ran'+str(nran)+'_rsk'+str(nrsk)+'.fits'
random_vmax_data = fits.open(vmax_avail_dir+random_vmax_name)
random_vmax_z           = random_vmax_data[1].data.field('z')
random_vmax_vmax_avail  = random_vmax_data[1].data.field('v_max_avail')

fig1 = py.figure(1) 
ax10 = fig1.add_subplot(111)
ax10.plot(random_vmax_z, random_vmax_vmax_avail, color='k',label='Random')

ransack = np.loadtxt(ransack_dir+'ransack_sdss_'+str(nrsk)+'.dat',skiprows=1)
fig2 = py.figure(2,figsize=(15,10))
ax20 = fig2.add_subplot(111)
ax20.scatter(ransack[:,0], ransack[:,1], color='r', label='Ransack') 

sfq = ['active','quiescent']
sfq_color = ['b','r']
sfq_label = ['Active', 'Quiescent']
for i in range(len(sfq)): 
    envcount_vmax_name = 'vmax_avail_cylr'+rad+'h'+height+'_thresh'+thresh+'_nbin'+nbin+'_sdss_ran'+str(nran)+'_rsk'+str(nrsk)+'_'+sfq[i]+'_sdsstest.fits'
    envcount_vmax_data = fits.open(vmax_avail_dir+envcount_vmax_name)
    envcount_vmax_z           = envcount_vmax_data[1].data.field('z')
    envcount_vmax_vmax_avail  = envcount_vmax_data[1].data.field('v_max_avail')
    ax10.plot(envcount_vmax_z, envcount_vmax_vmax_avail, color=sfq_color[i], label=sfq_label[i]+' Galaxies')

    envcount_name = 'envcount_cylr'+rad+'h'+height+'_thresh'+thresh+'_nbin'+nbin+'_sdss_'+sfq[i]+'_EDP-sdss-z00375_0145-numden.fits'
    envcount_data = fits.open(envcount_dir+envcount_name) 
    envcount_ra     = envcount_data[1].data.field('ra')
    envcount_dec    = envcount_data[1].data.field('dec')
    if i==0: 
        data_ra     = envcount_ra
        data_dec    = envcount_dec 
    else: 
        data_ra     = np.concatenate([data_ra, envcount_ra])
        data_dec    = np.concatenate([data_dec, envcount_dec])
ax10.set_xlabel('Redshift ($z$)',fontsize=20)
ax10.set_ylabel(r'$V_{\rm{max},\rm{avail}}$',fontsize=20)
ax10.legend(loc='upper left')
ax10.text(0.1,random_vmax_vmax_avail[0], r'Radius '+rad+' Height '+height, fontsize=18)

mfdata_dir  = '/global/data/scr/chh327/primus/mfdata/mfdata_chh/new/'
mfdata_name = 'mfdata_all_supergrid01_sdss.fits.gz'
mfdata_data = fits.open(mfdata_dir+mfdata_name)
mfdata_ra   = mfdata_data[1].data.field('ra')
mfdata_dec  = mfdata_data[1].data.field('dec')

ax20.scatter(mfdata_ra, mfdata_dec, color='k',label='Parent Data')
ax20.scatter(data_ra, data_dec, color='g',label='Envcount Data')
ax20.legend(loc='upper left')

fig1.text(0.5,0.95,r"SDSS $V_{\rm{max},\rm{avail}}$ Comparison", ha='center', va='center', fontsize=18)
fig2.text(0.5,0.95,r"SDSS Ransack, Environment Count Data, and Parent Data $RA$ and $Dec$ Comparison", ha='center', va='center', fontsize=18)
fig1.savefig('/home/users/hahn/figures/primus/vmax_avail/sdss_random_data_vmax_avail_comparison_cylr'+rad+'h'+height+'.png')
fig2.savefig('/home/users/hahn/figures/primus/vmax_avail/sdss_ransack_envcount_mfdata_ra_dec_comparison.png')
py.show()
