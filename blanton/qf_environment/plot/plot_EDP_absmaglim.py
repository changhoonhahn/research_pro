import numpy as np
import pylab as py
import matplotlib.pyplot as plt
from matplotlib import rc
import pyfits as fits

rc('text',usetex=True)
rc('font',**{'family':'monospace','monospace':['Courier']})

envdir     = '/global/data/scr/chh327/primus/data/EDP/' 
sdss    = np.loadtxt(envdir+'env_zbin0.100000_lim.dat')
lowz    = np.loadtxt(envdir+'env_zbin0.300000_lim.dat')
midz    = np.loadtxt(envdir+'env_zbin0.500000_lim.dat')  
highz   = np.loadtxt(envdir+'env_zbin0.700000_lim.dat') 
higherz = np.loadtxt(envdir+'env_zbin0.900000_lim.dat') 

targetdir       = '/global/data/scr/chh327/primus/data/target/'

target_sdss     = fits.open(targetdir+'target_cylr2h50_thresh75_nbin5_all_sdss_lit.fits')
target_primus   = fits.open(targetdir+'target_cylr2h50_thresh75_nbin20_all_lit.fits') 

target_sdss_mg  = target_sdss[1].data.field('MG_01') 
target_sdss_z   = target_sdss[1].data.field('REDSHIFT') 
target_sdss_mass        = target_sdss[1].data.field('mass') 
target_sdss_masslimit   = target_sdss[1].data.field('masslimit') 
target_primus_mg= target_primus[1].data.field('MG_01')
target_primus_z = target_primus[1].data.field('REDSHIFT')
target_primus_mass        = target_primus[1].data.field('mass') 
target_primus_masslimit   = target_primus[1].data.field('masslimit') 

fig1 = py.figure(1,figsize=(17,6))

target_sdss_z_lim   = target_sdss_z[ target_sdss_mass > target_sdss_masslimit ] 
target_sdss_mg_lim  = target_sdss_mg[ target_sdss_mass > target_sdss_masslimit ]

ax01 = fig1.add_subplot(151)
ax01.text(0.045,-25.,'SDSS: $z=0.0375-0.145$',fontsize=15)
ax01.scatter(target_sdss_z_lim,target_sdss_mg_lim,color='k',s=1,marker='o')
ax01.scatter(sdss[:,1], sdss[:,0],color='r',s=1,marker='o')
ax01.set_xlim([0.0375, 0.145])
ax01.set_ylim([-25.7, -15.0])
ax01.set_xticks([0.05,0.075,0.1,0.125])
#ax01.locator_params(axis='x',nbins=6)
ax01.set_ylabel(r'$M_{\rm{g}}$', fontsize=24)
plt.gca().invert_yaxis()

target_primus_mg_lim0 = target_primus_mg[ (target_primus_z > 0.2) & (target_primus_z < 0.4) & (target_primus_mass > target_primus_masslimit)] 
target_primus_z_lim0  = target_primus_z[ (target_primus_z > 0.2) & (target_primus_z < 0.4) & (target_primus_mass > target_primus_masslimit)] 

ax02 = fig1.add_subplot(152)
ax02.text(0.21,-25.,'PRIMUS: $z=0.2-0.4$',fontsize=15)
ax02.scatter(target_primus_z_lim0,target_primus_mg_lim0,color='k',s=1,marker='o')
ax02.scatter(lowz[:,1], lowz[:,0],color='r',s=1,marker='o')
ax02.set_xlim([0.2, 0.4])
ax02.set_ylim([-25.7, -15.0])
ax02.get_yaxis().set_ticklabels([])
plt.gca().invert_yaxis()

target_primus_mg_lim1 = target_primus_mg[ (target_primus_z > 0.4) & (target_primus_z < 0.6) & (target_primus_mass > target_primus_masslimit)] 
target_primus_z_lim1  = target_primus_z[ (target_primus_z > 0.4) & (target_primus_z < 0.6) & (target_primus_mass > target_primus_masslimit)] 

ax03 = fig1.add_subplot(153)
ax03.text(0.41,-25.,'PRIMUS: $z=0.4-0.6$',fontsize=15)
ax03.scatter(target_primus_z_lim1,target_primus_mg_lim1,color='k',s=1,marker='o')
ax03.scatter(midz[:,1], midz[:,0],color='r',s=1,marker='o')
ax03.set_xlim([0.4, 0.6])
ax03.set_ylim([-25.7, -15.0])
ax03.set_xlabel(r'{\bf Redshift}', fontsize=20)
ax03.get_yaxis().set_ticklabels([])
plt.gca().invert_yaxis()

target_primus_mg_lim2 = target_primus_mg[ (target_primus_z > 0.6) & (target_primus_z < 0.8) & (target_primus_mass > target_primus_masslimit)] 
target_primus_z_lim2  = target_primus_z[ (target_primus_z > 0.6) & (target_primus_z < 0.8) & (target_primus_mass > target_primus_masslimit)] 

ax04 = fig1.add_subplot(154)
ax04.text(0.61,-25.,'PRIMUS: $z=0.6-0.8$',fontsize=15)
ax04.scatter(target_primus_z_lim2,target_primus_mg_lim2,color='k',s=1,marker='o')
ax04.scatter(highz[:,1], highz[:,0],color='r',s=1,marker='o')
ax04.set_xlim([0.6,0.8])
ax04.set_ylim([-25.7, -15.0])
ax04.get_yaxis().set_ticklabels([])
plt.gca().invert_yaxis()

target_primus_mg_lim3 = target_primus_mg[ (target_primus_z > 0.8) & (target_primus_z < 1.0) & (target_primus_mass > target_primus_masslimit)] 
target_primus_z_lim3  = target_primus_z[ (target_primus_z > 0.8) & (target_primus_z < 1.0) & (target_primus_mass > target_primus_masslimit)] 

ax05 = fig1.add_subplot(1,5,5)
ax05.text(0.81,-25.0,'PRIMUS: $z=0.8-1.0$',fontsize=15)
ax05.scatter(target_primus_z_lim3,target_primus_mg_lim3,color='k',s=1,marker='o')
ax05.scatter(higherz[:,1], higherz[:,0],color='r',s=1,marker='o')
ax05.set_xlim([0.8,1.0])
ax05.set_ylim([-25.7, -15.0])
ax05.get_yaxis().set_ticklabels([])
plt.gca().invert_yaxis()

fig1.subplots_adjust(left=0.10,right=0.90,wspace=0.0,hspace=0.0)

bins = [-25.75,-25.5,-25.25,-25.0,-24.75,-24.50,-24.25,-24.0,-23.75,-23.5,-23.25,-23.0,-22.75,-22.5,-22.25,-22.0,-21.75,-21.5,-21.25,-21.0,-20.75,-20.50,-20.25,-20.0,-19.75,-19.5,-19.25,-19,-18.75,-18.5,-18.25,-18.0,-17.75,-17.50,-17.25,-17.0]

fig2 = py.figure(2,figsize=(20,7))
ax21 = fig2.add_subplot(151)
#py.text(0.05, 0.85,'Low Z 90%: -18.51',ha='left', va='top', transform = ax.transAxes)
ax21.hist(sdss[:,0], bins, normed=1, histtype='bar')
ax21.set_xlabel(r'$M_{\rm{g}}$', fontsize=15)

ax22 = fig2.add_subplot(152)
#py.text(0.05, 0.85,'Low Z 90%: -18.51',ha='left', va='top', transform = ax.transAxes)
ax22.hist(lowz[:,0], bins, normed=1, histtype='bar')
ax22.set_xlabel(r'$M_{g}$', fontsize=15)

ax23 = fig2.add_subplot(153)
#py.text(0.05, 0.85,'Mid Z 90%: -19.65',ha='left', va='top', transform = ax.transAxes)
ax23.hist(midz[:,0], bins, normed=1)
ax23.set_xlabel(r'$M_{g}$', fontsize=15)

ax24 = fig2.add_subplot(154)
#py.text(0.05, 0.85,'High Z 90%: -20.37',ha='left', va='top', transform = ax.transAxes)
ax24.hist(highz[:,0],bins, normed=1, histtype='bar')
ax24.set_xlabel(r'$M_{g}$', fontsize=15)

ax25 = fig2.add_subplot(1,5,5)
#py.text(0.05, 0.85,'High Z 90%: -20.37',ha='left', va='top', transform = ax.transAxes)
ax25.hist(higherz[:,0],bins, normed=1, histtype='bar')
ax25.set_xlabel(r'$M_{g}$', fontsize=15)

fig3 = py.figure(3,figsize=(10,10))
ax31 = fig3.add_subplot(111)                               
ax31.scatter(sdss[:,1], sdss[:,0],s=1,marker='o')
z1 = 0.145*np.ones(500)
z2 = 0.035*np.ones(500)
z3 = 0.04*np.ones(500)
z4 = 0.05*np.ones(500)
y = -26.0+0.1*np.arange(500)
ax31.plot(z1,y,'r',linewidth=2)
ax31.plot(z2,y,'r',linewidth=2)
ax31.plot(z3,y,'b',linewidth=2)
ax31.plot(z4,y,'g',linewidth=2)
ax31.set_xlim([0.0, 0.2])
ax31.set_ylim([-25.7, -17.0])
ax31.set_ylabel(r'Absolute Magnitude $M_{g}$', fontsize=15)
ax31.set_xlabel(r'Redshift $(z)$', fontsize=15)
plt.gca().invert_yaxis()

py.show()
