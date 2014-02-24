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
ebin = ['hienv','midenv','lowenv']
etitle = ['High Environment','Mid Environment','Low Environment']
zbin = ['0810z','0608z','0406z','0204z']
zbinval = [0.9, 0.7, 0.5, 0.3]
sfq = ['active', 'quiescent']
sfqtitle = ['Star Forming', 'Quiescent','Quiescent Fraction']
color = ['k','b','m','r']
markers = ['o', 's', '^']
redshift = ["$0.8<z<1.0$","$0.6<z<0.8$","$0.4<z<0.6$","$0.2<z<0.4$"]
smf_mass = []
sf = []
q = []
sferr = []
qerr = []

rc('text', usetex=True)
rc('font', family='serif')
run = int(sys.argv[1])
par = yan.yanny('/home/users/hahn/qf_environment/zerod_environment_parameters.par',np=True)

rad = str(par['PARAMETERS']['cylrad'][run-1])
height = str(par['PARAMETERS']['cylheight'][run-1])
thresh = str(par['PARAMETERS']['thresh'][run-1])
nbin = str(par['PARAMETERS']['nbin'][run-1])

fig=plt.figure(1,figsize=(10,10))
for i in range(len(sfq)):
    for ii in range(len(ebin)):
        ax = fig.add_subplot(3, 3, 3*i + ii + 1)
        for iii in range(len(zbin)):
            fname = 'smf_primus_cylr'+rad+'h'+height+'_thresh'+thresh+'_nbin'+nbin+'_lit_'+sfq[i]+'_'+ebin[ii]+'_'+zbin[iii]+'.fits'
            print fname
            smfdata = fits.open(smfdir+fname)
            mass = smfdata[1].data.field('mass')
            phi  = smfdata[1].data.field('phi')
            phierr_poisson  = smfdata[1].data.field('phierr_poisson')
            phierr_lower    = smfdata[1].data.field('phierr_lower')
            phierr_upper    = smfdata[1].data.field('phierr_upper')
            phierr_cv       = smfdata[1].data.field('phierr_cv')
            phierr          = phierr_cv
            error = 0.434*(phierr[0]/phi[0])
            ax.plot(mass[0],np.log10(phi[0]),c=color[iii],label=redshift[iii])
            ax.errorbar(mass[0],np.log10(phi[0]),yerr=error,fmt=color[iii]+'o')
            if i==0: 
                sf.append(phi[0])
                sferr.append(phierr[0])
                smf_mass.append(mass[0])
            else: 
                q.append(phi[0])
                qerr.append(phierr[0])

        if i == 0:
            ax.set_title(etitle[ii])
            if ii== 0:
                leg_prop = matplotlib.font_manager.FontProperties(size=12)
                lg=ax.legend(loc='lower left', prop=leg_prop)
                lg.get_frame().set_linewidth(0)

        if ii == 0:
            ax.set_ylabel(r"$\bf{log(\Phi / Mpc^{-3} dex^{-1})}$", fontsize=15)
            ax.tick_params(labelsize=15)
        else:
            ax.get_yaxis().set_ticklabels([])

        if ii == 2:
            axx = ax.twinx()
            axx.set_ylabel(sfqtitle[i], fontsize=15)
            axx.get_yaxis().set_ticklabels([])
        ax.set_xlim([8.20, 11.80])
        ax.set_ylim([-5.75,-1.75])
        ax.get_xaxis().set_ticklabels([])
#################################################################################
# QUIESCENT FRACTION
#################################################################################
qffitfig = py.figure(2)
qflinfitfig = py.figure(3,figsize=(10,5))

offsets = np.zeros([len(ebin),len(zbin)])

for j in range(len(ebin))[::-1]:
    ax1=fig.add_subplot(3,3,7+j)
    qflinfitax = qflinfitfig.add_subplot(1,3,j+1)
    for jj in range(len(zbin))[::-1]:
        sfdata = sf[len(zbin)*j+jj]
        sferrdata = sferr[len(zbin)*j+jj]

        qdata = q[len(zbin)*j+jj]
        qerrdata = qerr[len(zbin)*j+jj]
        massdata = smf_mass[len(zbin)*j+jj]

        qf = qdata/(sfdata+qdata)
        #qferr = np.sqrt((sfdata/(sfdata+qdata)**2*qerrdata)**2+(qdata/(sfdata+qdata)**2*sferrdata)**2)

        qf_fname = 'qf_primus_cylr'+rad+'h'+height+'_thresh'+thresh+'_nbin'+nbin+'_lit_'+ebin[j]+'_'+zbin[jj]+'.fits'
        print qf_fname
        qfdata = fits.open(smfdir+qf_fname)
        qfmass = qfdata[1].data.field('mass')
        qfqf = qfdata[1].data.field('qf')
        qferr = qfdata[1].data.field('err')
        qferrcv  = qfdata[1].data.field('err_cv')

        qfzero = qfqf > 0
        qfno0 = qfqf[qfzero]
        qferrno0 = qferr[qfzero]
        massno0 = qfmass[qfzero]
        ax1.plot(massno0,qfno0,color[jj]+'--')
                # +'o-',label=redshift[jj])
        ax1.errorbar(massno0,qfno0,yerr=qferrno0,fmt=color[jj]+'o', capsize=5, elinewidth=2)

        massout = massno0
        qfout = qfno0
        qferrout = qferrno0

        qffitax = qffitfig.add_subplot(111)

        def linfit(x, a, b):
            x = x.astype(np.float)
            a = a.astype(np.float)
            b = b.astype(np.float)
            return a*x+b

        def yintfit(x,b):
            x = x.astype(np.float)
            b = b.astype(np.float)
            return slope*x+b

        if (j==len(ebin)-1 and jj==len(zbin)-1):
            guess = [0.3,-2.7]
            slope = (optimize.curve_fit(linfit, massout, qfout, guess, qferrout)[0])[0]
            yint = (optimize.curve_fit(linfit, massout, qfout, guess, qferrout)[0])[1]
            offsets[j][jj] = slope*(10.5)+yint
            print 'ebin'+ebin[j]+'zbin'+zbin[jj]+'slope=', slope, 'offset at 10.5=', offsets[j][jj]
            ax1.plot(massout,linfit(massout,slope,yint),color[jj]+'--',linewidth=2.5,label=redshift[jj])
            qflinfitax.plot(massout,linfit(massout,slope,yint),color[jj]+'--',linewidth=2.5,label=redshift[jj])
        else:
            guess = [-4.0]
            #yint = (optimize.curve_fit(yintfit,massout,qfout,guess,qferrout)[0])[0]
            yint = (optimize.curve_fit(yintfit,massout,qfout,guess)[0])[0]
            offsets[j][jj] = slope*(10.5)+yint
            print 'ebin'+ebin[j]+'zbin'+zbin[jj]+'slope=', slope, 'offset at 10.5=', offsets[j][jj]
            massouty = [ slope*massout[k]+yint for k in range(len(massout)) ]
            qflinfitax.plot(massout,massouty,color[jj]+'--',linewidth=2.5)

    if (j==0):
        ax1.set_ylabel(r'Quiescent Fraction',fontsize=15)
        ax1.tick_params(labelsize=14)
        qflinfitax.set_ylabel(r'Quiescent Fraction',fontsize=15)
        qflinfitax.tick_params(labelsize=14)
    else:
        ax1.get_yaxis().set_ticklabels([])
        ax1.tick_params(labelsize=14)
        qflinfitax.get_yaxis().set_ticklabels([])
        qflinfitax.tick_params(labelsize=14)
    ax1.set_xlim([8.20, 11.80])
    ax1.set_ylim([0.0,1.0])
    ax1.set_xlabel(r"$\bf{log(M/M_{\odot})}$",fontsize=18)

    qflinfitax.set_xlim([8.20, 11.80])
    qflinfitax.set_ylim([0.0,1.0])
    qflinfitax.set_xlabel(r"$\bf{log(M/M_{\odot})}$",fontsize=18)

for k in range(len(ebin)):
    qffitax.plot( zbinval, offsets[k,:], 'k'+markers[k]+'-', label=etitle[k])

qffitax.set_xlim([0.2, 1.0])
qffitax.set_xlabel('Redshift (z)',fontsize=18)
qffitax.set_ylabel(r'Linear fit value at $log(M/M_{sun})=10.5$',fontsize=18)
qffitax.legend(loc='best')

fig.subplots_adjust(left=0.10,right=0.90,bottom=0.05,top=0.95,wspace=0.0,hspace=0.0)
qflinfitfig.subplots_adjust(left=0.10,right=0.90,wspace=0.0,hspace=0.0)
plt.show()
