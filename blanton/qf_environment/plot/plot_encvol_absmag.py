import numpy as np
import pylab as py
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text',usetex=True)
rc('font',family='serif')

#lowz = np.loadtxt('env_encvol_absmag_lowz.dat')
#midz = np.loadtxt('env_encvol_absmag_midz.dat')
#highz= np.loadtxt('env_encvol_absmag_highz.dat')

#path = '/global/data/scr/chh327/primus/completeness/'
#lowz = np.loadtxt(path+'primus_parent_completeness__encvol_v_absmag_lowz.dat')
#midz = np.loadtxt(path+'primus_parent_completeness__encvol_v_absmag_midz.dat')
#highz= np.loadtxt(path+'primus_parent_completeness__encvol_v_absmag_highz.dat')
#sdss = np.loadtxt(path+'sdss_parent_completeness_encvol-v-absmag.dat')

dir     = '/global/data/scr/chh327/primus/data/EDP/' 
sdss    = np.loadtxt(dir+'env_zbin0.100000.dat')
lowz    = np.loadtxt(dir+'env_zbin0.300000.dat')
midz    = np.loadtxt(dir+'env_zbin0.500000.dat')  
highz   = np.loadtxt(dir+'env_zbin0.700000.dat') 
higherz = np.loadtxt(dir+'env_zbin0.900000.dat') 

fig1 = py.figure(1,figsize=(10,7))
ax01 = fig1.add_subplot(151)
ax01.scatter(sdss[:,1], sdss[:,0],s=1,marker='o')
ax01.set_xlim([0.0, 0.2])
ax01.set_ylim([-25.7, -17.0])
ax01.set_ylabel(r'Absolute Magnitude $M_{r}$', fontsize=15)
plt.gca().invert_yaxis()

ax02 = fig1.add_subplot(152)
ax02.scatter(lowz[:,1], lowz[:,0],s=1,marker='o')
ax02.set_xlim([0.2, 0.4])
ax02.set_ylim([-25.7, -17.0])
ax02.get_yaxis().set_ticklabels([])
plt.gca().invert_yaxis()

ax03 = fig1.add_subplot(153)
ax03.scatter(midz[:,1], midz[:,0],s=1,marker='o')
ax03.set_xlim([0.4, 0.6])
ax03.set_ylim([-25.7, -17.0])
ax03.set_xlabel('Redshift (z)', fontsize=15)
ax03.get_yaxis().set_ticklabels([])
plt.gca().invert_yaxis()

ax04 = fig1.add_subplot(154)
ax04.scatter(highz[:,1], highz[:,0],s=1,marker='o')
ax04.set_xlim([0.6,0.8])
ax04.set_ylim([-25.7, -17.0])
ax04.get_yaxis().set_ticklabels([])
plt.gca().invert_yaxis()

ax05 = fig1.add_subplot(155)
ax05.scatter(higherz[:,1], higherz[:,0],s=1,marker='o')
ax05.set_xlim([0.8,1.0])
ax05.set_ylim([-25.7, -17.0])
ax05.get_yaxis().set_ticklabels([])
plt.gca().invert_yaxis()

bins = [-25.75,-25.5,-25.25,-25.0,-24.75,-24.50,-24.25,-24.0,-23.75,-23.5,-23.25,-23.0,-22.75,-22.5,-22.25,-22.0,-21.75,-21.5,-21.25,-21.0,-20.75,-20.50,-20.25,-20.0,-19.75,-19.5,-19.25,-19,-18.75,-18.5,-18.25,-18.0,-17.75,-17.50,-17.25,-17.0]

fig2 = py.figure(2,figsize=(20,7))
ax21 = fig2.add_subplot(151)
#py.text(0.05, 0.85,'Low Z 90%: -18.51',ha='left', va='top', transform = ax.transAxes)
ax21.hist(sdss[:,0], bins, normed=1, histtype='bar')
ax21.set_xlabel(r'$M_{r}$', fontsize=15)

ax22 = fig2.add_subplot(152)
#py.text(0.05, 0.85,'Low Z 90%: -18.51',ha='left', va='top', transform = ax.transAxes)
ax22.hist(lowz[:,0], bins, normed=1, histtype='bar')
ax22.set_xlabel(r'$M_{r}$', fontsize=15)

ax23 = fig2.add_subplot(153)
#py.text(0.05, 0.85,'Mid Z 90%: -19.65',ha='left', va='top', transform = ax.transAxes)
ax23.hist(midz[:,0], bins, normed=1)
ax23.set_xlabel(r'$M_{r}$', fontsize=15)

ax24 = fig2.add_subplot(154)
#py.text(0.05, 0.85,'High Z 90%: -20.37',ha='left', va='top', transform = ax.transAxes)
Mg = -20.75*np.ones(6)
y = 0.1*np.arange(6)
ax24.plot(Mg,y,'r',linewidth=3)
ax24.hist(highz[:,0],bins, normed=1, histtype='bar')
ax24.set_xlabel(r'$M_{r}$', fontsize=15)

ax25 = fig2.add_subplot(155)
#py.text(0.05, 0.85,'High Z 90%: -20.37',ha='left', va='top', transform = ax.transAxes)
ax25.hist(higherz[:,0],bins, normed=1, histtype='bar')
Mg = -20.75*np.ones(7)
y = 0.1*np.arange(7)
ax25.plot(Mg,y,'r',linewidth=3)
ax25.set_xlabel(r'$M_{r}$', fontsize=15)

py.show()
