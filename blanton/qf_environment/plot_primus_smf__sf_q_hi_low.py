import numpy as np
import math as m
import matplotlib.pyplot as plt
import matplotlib.font_manager
from matplotlib import rc
import scipy.optimize as optimize
from scipy import interpolate
import matplotlib.pylab as py
import yanny as yan
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
mass = []
sf = []
q = []
hist_sf = []
hist_q = []
hist_sferr = []
hist_qerr = []
sferr = []
qerr = []

rc('text', usetex=True)
rc('font', family='serif')

par = yan.yanny('/home/users/hahn/primus/pro/spec0d/envcount/zerod_environment_parameters.par',np=True)
run = int(raw_input('Run: '))

rad = str(par['PARAMETERS']['cylrad'][run-1])
height = str(par['PARAMETERS']['cylheight'][run-1])
thresh = str(par['PARAMETERS']['thresh'][run-1])
nbin = str(par['PARAMETERS']['nbin'][run-1])

fig=plt.figure(1,figsize=(10,10))
for i in range(len(sfq)):
    for ii in range(len(ebin)):
        ax = fig.add_subplot(3, 3, 3*i + ii + 1)
        for iii in range(len(zbin)):
            fname = 'smf_primus_cylr'+rad+'h'+height+'_thresh'+thresh+'_nbin'+nbin+'_lit_'+sfq[i]+'_'+ebin[ii]+'_'+zbin[iii]+'.dat'
            print fname
            data=np.loadtxt(smfdir+fname)
            hist_data=np.loadtxt(smfdir+'hist_'+fname)

            ax.plot(data[:,0],data[:,1],c=color[iii],label=redshift[iii])
            ax.errorbar(data[:,0],data[:,1],yerr=data[:,2],fmt=color[iii]+'o')
            if i==0:
                sf.append(data[:,1])
                sferr.append(data[:,2])
                mass.append(data[:,0])

                hist_sf.append(hist_data[:,1])
                hist_sferr.append(hist_data[:,2])
            else:
                q.append(data[:,1])
                qerr.append(data[:,2])

                hist_q.append(hist_data[:,1])
                hist_qerr.append(hist_data[:,2])

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
        histsf = hist_sf[len(zbin)*j+jj]
        histsf_err = hist_sferr[len(zbin)*j+jj]

        qdata = q[len(zbin)*j+jj]
        qerrdata = qerr[len(zbin)*j+jj]
        histq = hist_q[len(zbin)*j+jj]
        histq_err = hist_qerr[len(zbin)*j+jj]

        massdata = mass[len(zbin)*j+jj]

        qf = 1.0/(10.0**(sfdata-qdata)+1.0)
        qf_hist = histq/(histsf+histq)
        qferr = qf_hist*np.sqrt( (histsf_err/histsf)**2 + (np.sqrt(histsf_err**2+histq_err**2)/(histsf+histq))**2 )

        qfzero = qf > 0
        qfno0 = qf[qfzero]
        qferrno0 = qferr[qfzero]
        massno0 = massdata[qfzero]
        ax1.plot(massno0,qfno0,color[jj]+'o-',label=redshift[jj])
        ax1.errorbar(massno0,qfno0,yerr=qferrno0,fmt=color[jj]+'o', capsize=5, elinewidth=2)

        if (j==2):
            maxm = max(massno0)
        muplim = massno0 < maxm
        qfout = qfno0[muplim]
        qferrout = qferrno0[muplim]
        massout = massno0[muplim]

        qffitax = qffitfig.add_subplot(111)
        if (j==len(ebin)-1 and jj==len(zbin)-1):
            def linfit(x, a, b):
                return a*x+b
            guess = [ 0.3, -2.7 ]
            #slope = (optimize.curve_fit(linfit, massout, qfout, guess, qferrout)[0])[0]
            #yint = (optimize.curve_fit(linfit, massout, qfout, guess, qferrout)[0])[1]
            slope = (optimize.curve_fit(linfit, massout, qfout, guess)[0])[0]
            yint = (optimize.curve_fit(linfit, massout, qfout, guess)[0])[1]
            offsets[j][jj] = slope*(10.5)+yint
            print 'ebin'+ebin[j]+'zbin'+zbin[jj]+'slope=', slope, 'offset at 10.5=', offsets[j][jj]
            ax1.plot(massout,linfit(massout,slope,yint),color[j]+'--',linewidth=2.5,label=redshift[jj])
            qflinfitax.plot(massout,linfit(massout,slope,yint),color[jj]+'--',linewidth=2.5,label=redshift[jj])
        else:
            def yintfit(x,b):
                return slope*x+b
            guess = [-2.7]
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
