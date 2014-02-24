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
############################################################################
rc('text', usetex=True)
rc('font', family='serif')

fig=plt.figure(1,figsize=(10,10))
ax = fig.add_subplot(111)

mfdata_name = '/home/users/hahn/qf_environment/dump/sdss_masslimit_smftest/mfdata_smf.fits'
target_name = '/home/users/hahn/qf_environment/dump/sdss_masslimit_smftest/target_smf.fits'
target_novmax_name = '/home/users/hahn/qf_environment/dump/sdss_masslimit_smftest/target_noedgecut_smf.fits'
target_incorrect_name = '/home/users/hahn/qf_environment/dump/sdss_masslimit_smftest/target_incorrect_smf.fits'

mfdata_smfdata = fits.open(mfdata_name)
target_smfdata = fits.open(target_name)
target_novmax_smfdata = fits.open(target_novmax_name)
target_incorrect_smfdata = fits.open(target_incorrect_name)

mfdata_mass = mfdata_smfdata[1].data.field('mass')
mfdata_phi  = mfdata_smfdata[1].data.field('phi')

target_mass = target_smfdata[1].data.field('mass')
target_phi  = target_smfdata[1].data.field('phi')

target_novmax_mass = target_novmax_smfdata[1].data.field('mass')
target_novmax_phi  = target_novmax_smfdata[1].data.field('phi')

target_incorrect_mass = target_incorrect_smfdata[1].data.field('mass')
target_incorrect_phi  = target_incorrect_smfdata[1].data.field('phi')

ax.plot(mfdata_mass[0],np.log10(mfdata_phi[0]),c='k',linewidth=3,label='Entire SDSS sample')
ax.plot(target_mass[0],np.log10(target_phi[0]),c='r',linewidth=3,label=r"$V_{\rm{max}}$ SDSS Targets")
ax.plot(target_novmax_mass[0],np.log10(target_novmax_phi[0]),c='b',linewidth=3,label=r"$V_{\rm{max},\rm{avail}}$ SDSS Targets")
ax.plot(target_incorrect_mass[0],np.log10(target_incorrect_phi[0]),c='g',linewidth=3,label=r"Uncorrected $V_{\rm{max}}$ SDSS Targets")
ax.vlines(9.0,-6,0,color='r',label="9.0")
ax.vlines(9.5,-6,0,color='m',label="9.5")
ax.vlines(10.0,-6,0,color='b',label="10.0")

ax.set_xlim([7.5, 11.80])
#ax.set_ylim([-5.75,-1.75])
ax.set_xlabel('Mass',fontsize=16)
ax.set_ylabel('$\Phi$',fontsize=16)
ax.legend(loc='best')

mfdata_file = '/global/data/scr/chh327/primus/mfdata/mfdata_all_supergrid01_sdss_lit.fits.gz' 
target_file = '/global/data/scr/chh327/primus/data/target/target_cylr2h50_thresh75_nbin5_all_sdss_lit.fits' 

mfdata = fits.open(mfdata_file)
target = fits.open(target_file)

mfmass = mfdata[1].data.field('mass')
mfz = mfdata[1].data.field('z')
mfvmax_evol = mfdata[1].data.field('vmax_evol')

targz = target[1].data.field('redshift')
targmass = target[1].data.field('mass')
targzmax_evol = target[1].data.field('zmax_evol')
targzmin_evol = target[1].data.field('zmin_evol')
targvmax_evol = target[1].data.field('vmax_evol')
targvmax = target[1].data.field('vmax')

fig2 = plt.figure(2,figsize=(10,10))
ax2 = fig2.add_subplot(111)
ax2.scatter(mfz,mfmass,color='k')
ax2.scatter(targz,targmass,color='r')
ax2.set_xlim([0.0,0.2])
ax2.set_xlabel('Redshift',fontsize=20)
ax2.set_ylabel(r'$\rm{log} \: \mathcal{M}/\mathcal{M}_\odot$',fontsize=20)

fig3 = plt.figure(3,figsize=(10,10))
ax3 = fig3.add_subplot(111)
ax3.scatter(mfmass,mfvmax_evol,color='g',label=r'$V_{\rm{max}}$ Parent')
ax3.scatter(targmass,targvmax_evol,color='k', label=r'Not corrected $V_{\rm{max}}$ SDSS Targets')
ax3.scatter(targmass,targvmax,color='r', label=r'Corrected $V_{\rm{max}}$ SDSS Targets')
ax3.set_xlabel('Mass',fontsize=20)
ax3.set_ylabel(r'$V_{\rm{max}}$',fontsize=20)
ax3.legend(loc='best')

fig4 = plt.figure(4,figsize=(10,10))
ax4 = fig4.add_subplot(111)
ax4.scatter(targmass,targvmax/targvmax_evol,color='k')
ax4.set_xlabel('Mass',fontsize=20)
ax4.set_ylabel(r'$V_{\rm{max},\rm{corrected}}/V_{\rm{max},\rm{not corrected}}$',fontsize=20)

fig5 = plt.figure(5,figsize=(10,10))
ax5 = fig5.add_subplot(111)
ax5.scatter(targmass,targzmax_evol,color='k')
ax5.axhline(y=0.145,color='r')
ax5.set_xlabel('Mass',fontsize=20)
ax5.set_ylabel(r'$z_{\rm{max,evol}}$',fontsize=20)

fig6 = plt.figure(6,figsize=(10,10))
ax6 = fig6.add_subplot(111)
ax6.scatter(targmass,targzmin_evol,color='k')
ax6.set_xlabel('Mass',fontsize=20)
ax6.set_ylabel(r'$z_{\rm{min,evol}}$',fontsize=20)
ax6.axhline(y=0.06,color='r')
plt.show()
