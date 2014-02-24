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
smfdir = '/global/data/scr/chh327/primus/data/smf/'

fiducial_mass = 10.5

if (sys.argv[1]=='lit'):
    zbin = ['0810z','0608z','0406z','0204z','nobinz']
    redshift = ["$0.8<z<1.0$","$0.6<z<0.8$","$0.4<z<0.6$","$0.2<z<0.4$","$0.01<z<0.2$"]
    color = ['k','b','m','r','g']
    litsuffix='_lit'

if (sys.argv[1]=='nolit'):
    zbin = ['0810z','06508z','05065z','0405z','0304z','0203z','nobinz']
    redshift = ["$0.8<z<1.0$","$065<z<0.8$","$0.5<z<0.65$","$0.4<z<0.5$","$0.3<z<0.4$","$0.2<z<0.3$","$0.01<z<0.2$"]
    color = ['k','b','c','m','r','y','g']
    litsuffix=''
if (sys.argv[2]=='john'): 
    johnsuffix='_john'
else: 
    johnsuffix=''

sfq = ['active', 'quiescent']
sfqtitle = ['Star Forming', 'Quiescent', 'Quiescent Fraction']
markers = ['o', 's', '^']
smf_mass = []
sf = []
q = []
sferr = []
qerr = []

rc('text', usetex=True)
rc('font', family='serif')

mfdir = '/global/data/scr/chh327/primus/mfdata/'

fig=plt.figure(1,figsize=(10,10))
fig2=plt.figure(2,figsize=(10,10))
for i in range(len(sfq)):
    ax = fig.add_subplot(2,1,i+1)
    ax2= fig2.add_subplot(2,1,i+1)

    john_file   = 'mf_'+sfq[i]+'_supergrid01_bin0.10_evol'+litsuffix+'.fits.gz'
    john_data   = fits.open(mfdir+john_file) 
    john_mass   = john_data[1].data.field('mass')
    john_phi    = john_data[1].data.field('phi')
    john_limit  = john_data[1].data.field('limit')
    john_number = john_data[1].data.field('number')
    john_zbins  = len(john_mass)
    
    good_john_mass = []
    good_john_phi  = []
    for jj in range(len(john_mass)):
        good_mass   = (john_mass[jj])[(john_limit[jj]==1) & (john_number[jj]>3)]
        good_phi    = (john_phi[jj])[(john_limit[jj]==1) & (john_number[jj]>3)]
        ax.plot(good_mass,np.log10(good_phi),color[john_zbins-jj-1]+'--')
        good_john_mass.append(john_mass[jj])
        good_john_phi.append(john_phi[jj])

    sdss_john_file   = 'mf_'+sfq[i]+'_supergrid01_bin0.10_evol_sdss.fits.gz'
    sdss_john_data  = fits.open(mfdir+sdss_john_file) 
    sdss_john_mass  = sdss_john_data[1].data.field('mass')
    sdss_john_phi   = sdss_john_data[1].data.field('phi')
    sdss_john_limit = sdss_john_data[1].data.field('limit')
    sdss_john_number= sdss_john_data[1].data.field('number')

    sdss_good_mass   = (sdss_john_mass[0])[(sdss_john_limit[0]==1) & (sdss_john_number[0]>3)]
    sdss_good_phi    = (sdss_john_phi[0])[(sdss_john_limit[0]==1) & (sdss_john_number[0]>3)]
    ax.plot(sdss_good_mass,np.log10(sdss_good_phi),'g--',linewidth=2)
    
    sdss_john_mass = sdss_john_mass[0]
    sdss_john_phi = sdss_john_phi[0] 

    for ii in range(len(zbin)):
        sample = ''
        if (zbin[ii]=='nobinz'): 
            sample = '_sdss'
            fname = 'smf_mfdata_new_'+sfq[i]+'_supergrid01'+sample+'_'+zbin[ii]+johnsuffix+'.fits'
        else:
            fname = 'smf_mfdata_new_'+sfq[i]+'_supergrid01'+sample+litsuffix+'_'+zbin[ii]+johnsuffix+'.fits'
        print fname
        smfdata = fits.open(smfdir+fname)
        mass = smfdata[1].data.field('mass')
        phi  = smfdata[1].data.field('phi')
        limit   = smfdata[1].data.field('limit')
        number  = smfdata[1].data.field('number')
        phierr_poisson  = smfdata[1].data.field('phierr_poisson')
        phierr_lower    = smfdata[1].data.field('phierr_lower')
        phierr_upper    = smfdata[1].data.field('phierr_upper')
        phierr_cv       = smfdata[1].data.field('phierr_cv')
        phierr          = smfdata[1].data.field('phierr')
        phi_err         = phierr_poisson
        error = 0.434*(phi_err[0]/phi[0])

        good_mass   = (mass[0])[(limit[0]==1) & (number[0] > 3)]
        good_phi    = (phi[0])[(limit[0]==1) & (number[0] > 3)]
        ax.scatter(good_mass,np.log10(good_phi),color=color[ii],label=redshift[ii])
#        ax.errorbar(mass[0],phi[0],yerr=error,fmt=color[ii]+'o')
        if (zbin[ii]=='nobinz'):
            good_phi_ratio = (phi[0]/sdss_john_phi)[(limit[0]==1) & (number[0] > 3)]
            ax2.plot(good_mass, good_phi_ratio, 'g', linewidth=2)
        else: 
            good_phi_ratio = (phi[0]/good_john_phi[john_zbins-ii-1])[(limit[0]==1) & (number[0] > 3)]
            ax2.plot(good_mass, good_phi_ratio, color[john_zbins-ii-1], linewidth=2)

        if i==0: 
            sf.append(phi[0])
            sferr.append(phierr[0])
            smf_mass.append(mass[0])
        else: 
            q.append(phi[0])
            qerr.append(phierr[0])

        if i == 0:
            leg_prop = matplotlib.font_manager.FontProperties(size=15)
            lg=ax.legend(loc='lower left', prop=leg_prop)
            lg.get_frame().set_linewidth(0)

        ax.set_ylabel(r"$\bf{log(\Phi / Mpc^{-3} dex^{-1})}$", fontsize=15)
        ax2.set_ylabel(r"$\bf{\Phi_{\rm{Hahn}}/\Phi_{\rm{Moustakas}}}$", fontsize=20)
        ax.tick_params(labelsize=15)
        axx = ax.twinx()
        axx.set_ylabel(sfqtitle[i], fontsize=15)
        axx.get_yaxis().set_ticklabels([])
        ax.set_xlim([8.20, 11.80])
        ax.set_ylim([-5.75,-1.75])
        ax2.set_xlim([8.20, 11.80])
        ax2x = ax2.twinx()
        ax2x.set_ylabel(sfqtitle[i], fontsize=15)
        ax2x.get_yaxis().set_ticklabels([])
fig.subplots_adjust(left=0.1,right=0.85,bottom=0.05,top=0.95,wspace=0.0,hspace=0.0)
fig2.subplots_adjust(left=0.1,right=0.85,bottom=0.05,top=0.95,wspace=0.0,hspace=0.0)
if (sys.argv[2]=='john'):
    file_out = '/home/users/hahn/figures/primus/smf/smf_comparisontoJohnFig11'+litsuffix+'.png'
    file_ratio_out = '/home/users/hahn/figures/primus/smf/smf_comparisontoJohnFig11'+litsuffix+'_smfratio.png'
else: 
    file_out = '/home/users/hahn/figures/primus/smf/smf_comparisontoJohnFig11'+litsuffix+'_mfdata_chh.png'
    file_ratio_out = '/home/users/hahn/figures/primus/smf/smf_comparisontoJohnFig11'+litsuffix+'_mfdata_chh_smfratio.png'
fig.savefig(file_out)
fig2.savefig(file_ratio_out)
plt.show()
