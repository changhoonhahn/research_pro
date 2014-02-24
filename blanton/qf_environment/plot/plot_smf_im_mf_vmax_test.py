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

fiducial_mass = 11.0
binsize = '0.10'

litsuffix = '_lit'
johnsuffix = ''

#zbin = ['0810z','0608z','0406z','0204z']
#zbinval = [0.9, 0.7, 0.5, 0.3]
#redshift = ["$0.8<z<1.0$","$0.6<z<0.8$","$0.4<z<0.6$","$0.2<z<0.4$"]
#color = ['k','b','m','r']

zbin = ['0810z','0608z','0406z','0204z','0002z']
zbinval = [0.9, 0.7, 0.5, 0.3, 0.1025]
redshift = ["$0.8<z<1.0$","$0.6<z<0.8$","$0.4<z<0.6$","$0.2<z<0.4$","$0.0375<z<0.145$"]
color = ['k','b','m','r','g']

sfq = ['active', 'quiescent']
sfqtitle = ['Star Forming', 'Quiescent','Quiescent Fraction']
markers = ['o', 's', '^']
smf_mass = []
sf = []
q = []
sferr = []
qerr = []

rc('text', usetex=True)
rc('font', family='serif')
primus_run = int(sys.argv[1])
par = yan.yanny('/home/users/hahn/research/pro/blanton/qf_environment/zerod_environment_parameters.par',np=True)

rad = str(par['PARAMETERS']['cylrad'][primus_run-1])
height = str(par['PARAMETERS']['cylheight'][primus_run-1])
thresh = str(par['PARAMETERS']['thresh'][primus_run-1])

sdss_run = int(sys.argv[2])
par = yan.yanny('/home/users/hahn/research/pro/blanton/qf_environment/zerod_environment_parameters.par',np=True)

sdss_rad = str(par['PARAMETERS']['cylrad'][sdss_run-1])
sdss_height = str(par['PARAMETERS']['cylheight'][sdss_run-1])
sdss_thresh = str(par['PARAMETERS']['thresh'][sdss_run-1])

fig=plt.figure(1,figsize=(8,10))
fig2 = plt.figure(2,figsize=(8,10))

mfdir = '/global/data/scr/chh327/primus/mfdata/'
for i in range(len(sfq)):
    ax = fig.add_subplot(2,1,i+1)
    ax2 = fig2.add_subplot(2,1,i+1)
    
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
#    sdss_john_file   = 'mf_'+sfq[i]+'_supergrid01_bin0.10_evol_sdss.fits.gz' 
#    sdss_john_data  = fits.open(mfdir+sdss_john_file)    
    sdss_john_data  = fits.open('/global/data/scr/chh327/primus/data/smf/smf_mfdata_'+sfq[i]+'_supergrid01_sdss_nobinz_envzrange.fits')
    sdss_john_mass  = sdss_john_data[1].data.field('mass')
    sdss_john_phi   = sdss_john_data[1].data.field('phi')
    sdss_john_limit = sdss_john_data[1].data.field('limit')
    sdss_john_number= sdss_john_data[1].data.field('number')

    sdss_good_mass   = (sdss_john_mass[0])[(sdss_john_limit[0]==1) & (sdss_john_number[0]>3)]
    sdss_good_phi    = (sdss_john_phi[0])[(sdss_john_limit[0]==1) & (sdss_john_number[0]>3)]
    ax.plot(sdss_good_mass,np.log10(sdss_good_phi),'g--',linewidth=2)

    sdss_john_mass = sdss_john_mass[0]
    sdss_john_phi = sdss_john_phi[0]

    ax = fig.add_subplot(2,1,i+1)
    ax2 = fig2.add_subplot(2,1,i+1)
    for iii in range(len(zbin)):
        fname = 'smf_primus_cylr'+rad+'h'+height+'_thresh'+thresh+'_lit_'+sfq[i]+'_'+zbin[iii]+'_bin'+binsize+'_noenv.fits'
        if zbin[iii]=='0002z': 
            fname = 'smf_sdss_cylr'+sdss_rad+'h'+sdss_height+'_thresh'+sdss_thresh+'_'+sfq[i]+'_nobinz_bin'+binsize+'_noenv.fits'
        print fname
        smfdata = fits.open(smfdir+fname)
        mass    = smfdata[1].data.field('mass')
        phi     = smfdata[1].data.field('phi')
        limit   = smfdata[1].data.field('limit')
        number  = smfdata[1].data.field('number')
        phierr_poisson  = smfdata[1].data.field('phierr_poisson')
        phierr_lower    = smfdata[1].data.field('phierr_lower')
        phierr_upper    = smfdata[1].data.field('phierr_upper')
        phierr_cv       = smfdata[1].data.field('phierr_cv')
        phierr          = smfdata[1].data.field('phierr')
        phi_err         = phierr_cv
        error = 0.434*(phi_err[0]/phi[0])

        good_mass   = (mass[0])[(limit[0]==1) & (number[0] > 3)]
        good_phi    = (phi[0])[(limit[0]==1) & (number[0] > 3)]
        good_error  = (error)[(limit[0]==1) & (number[0] > 3)]
        ax.scatter(good_mass,np.log10(good_phi),color=color[iii],label=redshift[iii])
#            ax.errorbar(good_mass,np.log10(good_phi),yerr=good_error,fmt=color[iii]+'o')

        if (zbin[iii]=='0002z'): 
            good_mass = sdss_john_mass[(limit[0]==1) & (number[0] > 3)]
            print sdss_john_phi
            good_phi_ratio = (phi[0]/sdss_john_phi)[(limit[0]==1) & (number[0] > 3)]
            print good_phi_ratio
        else: 
            good_phi_ratio = (phi[0]/good_john_phi[john_zbins-iii-1])[(limit[0]==1) & (number[0] > 3)]
        ax2.plot(good_mass, good_phi_ratio, color[iii], linewidth=2)
        
#        print np.log10((phi[0])[(limit[0]==1) & (number[0] > 3)])
#        print np.log10((good_john_phi[john_zbins-iii-1])[(limit[0]==1) & (number[0] > 3)])
        if i==0: 
            sf.append(good_phi)
            sferr.append(good_error)
            smf_mass.append(good_mass)
        else: 
            q.append(good_phi)
            qerr.append(good_error)

    if i == 0:
        leg_prop = matplotlib.font_manager.FontProperties(size=12)
        lg=ax.legend(loc='lower left', prop=leg_prop)
        lg.get_frame().set_linewidth(0)

    ax.tick_params(labelsize=15)

    axx = ax.twinx()
    axx.set_ylabel(sfqtitle[i], fontsize=15)
    axx.get_yaxis().set_ticklabels([])

    ax.set_xlim([8.20, 11.80])
    ax.set_ylim([-5.75,-1.75])
    ax.get_xaxis().set_ticklabels([])
    ax.set_xlabel(r"$\bf{log(M/M_{\odot})}$",fontsize=18)
   
    if i == 0: 
        ax2.get_xaxis().set_ticklabels([])
    ax2.set_xlim([8.20, 11.80])
    ax2.set_xlabel(r"$\bf{log(M/M_{\odot})}$",fontsize=18)
    ax2x = ax2.twinx()
    ax2x.set_ylabel(sfqtitle[i], fontsize=15)
    ax2x.get_yaxis().set_ticklabels([])
fig.subplots_adjust(left=0.10,right=0.90,bottom=0.05,top=0.95,wspace=0.0,hspace=0.0)
fig2.subplots_adjust(left=0.10,right=0.90,bottom=0.05,top=0.95,wspace=0.0,hspace=0.0)

fig.text(0.5, 0.97, r'Radius '+str(rad)+' Height '+str(height), ha='center', va='center', fontsize=18)
fig2.text(0.5, 0.97, r'Radius '+str(rad)+' Height '+str(height), ha='center', va='center', fontsize=18)

fig.text(0.03,0.5,r"$\bf{log(\Phi / Mpc^{-3} dex^{-1})}$", ha='center', va='center', rotation='vertical', fontsize=18)
fig2.text(0.03,0.5,r"$\bf{\Phi_{V_{\rm{max,avail}}}/\Phi}$", ha='center', va='center', rotation='vertical', fontsize=18)

#fig.savefig('/home/users/hahn/figures/primus/smf/smf_mfdata_r'+str(rad)+'h'+str(height)+'edgecut_vmaxavail_comparison.png')
#fig2.savefig('/home/users/hahn/figures/primus/smf/smf_mfdata_r'+str(rad)+'h'+str(height)+'edgecut_vmaxavail_comparison_ratio.png')
plt.show()

##################################################################################
## QUIESCENT FRACTION
##################################################################################
#ax1=fig.add_subplot(3,1,3)
#for jj in range(len(zbin))[::-1]:
#    qf_fname = 'qf_primus_cylr'+rad+'h'+height+'_thresh'+thresh+'_nbin'+nbin+'_lit_'+zbin[jj]+'_bin'+binsize+'_test.fits'
#    if zbin[jj]=='0002z':
#        qf_fname = 'qf_sdss_cylr'+sdss_rad+'h'+sdss_height+'_thresh'+sdss_thresh+'_nbin'+sdss_nbin+'_nobinz_bin'+binsize+'_test.fits'
#        print qf_fname
#    qfdata = fits.open(smfdir+qf_fname)
#    qflimit = qfdata[1].data.field('limit')
#    qfmass = qfdata[1].data.field('mass')
#    qfqf = qfdata[1].data.field('qf')
#    qferr = qfdata[1].data.field('err')
#    qferrcv  = qfdata[1].data.field('err_cv')
#
#    qfzero = (qflimit==1) & (qfqf > 0.0) & (qfqf < 1.0)# & (qfmass > 10.5) & (qfmass < 11.5)
#    qfno0 = qfqf[qfzero]
#    qferrno0 = qferr[qfzero]
#    massno0 = qfmass[qfzero]
#    ax1.plot(massno0,qfno0,color[jj],label=redshift[jj])
#    ax1.errorbar(massno0,qfno0,yerr=qferrno0,fmt=color[jj]+'o', capsize=3, elinewidth=1)
#    ax1.set_xlim([8.20, 11.80])
#    ax1.set_ylim([0.0,1.0])
#    ax1.set_xlabel(r"$\bf{log(M/M_{\odot})}$",fontsize=18)
