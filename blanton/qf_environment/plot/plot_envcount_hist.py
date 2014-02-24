import numpy as np
import pylab as py
import pyfits as fits
from matplotlib import rc
import yanny as yan
import sys

rc('text', usetex=True)
rc('font', family='serif')

run = int(sys.argv[1])
par = yan.yanny('/home/users/hahn/research/pro/blanton/qf_environment/zerod_environment_parameters.par',np=True)

rad             = str(par['PARAMETERS']['cylrad'][run-1])
height          = str(par['PARAMETERS']['cylheight'][run-1])
thresh          = str(par['PARAMETERS']['thresh'][run-1])
highenvthresh   = par['PARAMETERS']['highenvthresh'][run-1] 
lowenvthresh    = par['PARAMETERS']['lowenvthresh'][run-1]

sdss_run = int(sys.argv[2])

sdss_envfile = 'EDP-sdss-z00375_0145-numden.fits'
envfile = 'EDP-primus-z0210-numden.fits'
if (sys.argv[3] == 'lit'): 
    litsuffix = '_lit'
else: 
    litsuffix = ''
zbin = [0.1, 0.3, 0.5, 0.7, 0.9]
zbin_labels = ['SDSS $0.0375 - 0.145$','PRIMUS $0.2 - 0.4$', 'PRIMUS $0.4 - 0.6$', 'PRIMUS $0.6 - 0.8$', 'PRIMUS $0.8 - 1.0$']
sf_q = ['active', 'quiescent'] 

all_envcount = []
all_redshift = []
envcount_dir    = '/global/data/scr/chh327/primus/data/envcount/'
for i in range(len(sf_q)):
    envcount_fname  = 'envcount_cylr'+rad+'h'+height+'_thresh'+thresh+'_sdss_'+sf_q[i]+'_'+sdss_envfile
    envcount_data   = fits.open(envcount_dir+envcount_fname)
    envcount = envcount_data[1].data.field('envcount')
    redshift = envcount_data[1].data.field('redshift')

    if (sf_q[i] == 'active'): 
        sdss_active_envcount = envcount
        sdss_active_redshift = redshift 
    if (sf_q[i] == 'quiescent'): 
        sdss_quiescent_envcount = envcount
        sdss_quiescent_redshift = redshift 

    all_envcount.extend(envcount)
    all_redshift.extend(redshift)

for i in range(len(sf_q)):
    envcount_fname  = 'envcount_cylr'+rad+'h'+height+'_thresh'+thresh+'_'+sf_q[i]+litsuffix+'_'+envfile
    envcount_data   = fits.open(envcount_dir+envcount_fname)
    envcount = envcount_data[1].data.field('envcount')
    redshift = envcount_data[1].data.field('redshift')

    if (sf_q[i] == 'active'): 
        active_envcount = envcount
        active_redshift = redshift 
    if (sf_q[i] == 'quiescent'): 
        quiescent_envcount = envcount
        quiescent_redshift = redshift 

    all_envcount.extend(envcount)
    all_redshift.extend(redshift)

all_envcount = np.array(all_envcount)
all_redshift = np.array(all_redshift)

fig1 = py.figure(1)
envhist = fig1.add_subplot(111)
env_bin  = 1.0*np.array(range(int(int(5.0*np.ceil(np.max(all_envcount)/5.0)))+1))
env_hist = envhist.hist(np.array(all_envcount), env_bin, normed=True) 

fig2 = py.figure(2)
envhist_log = fig2.add_subplot(111)
env_x   = [ (env_hist[1][i]+env_hist[1][i+1])/2.0 for i in range(len(env_hist[1])-1) ] 
envhist_log.plot(env_x, env_hist[0], linewidth=3, label=r'Environment Counts')
envhist_log.set_xlabel('Environment')
envhist_log.legend(loc='best')

fig3 = py.figure(3,figsize=(12,10))
fig4 = py.figure(4,figsize=(12,10))
fig5 = py.figure(5)
dump = fig5.add_subplot(111)
for i in range(len(zbin)): 
    envhist_z   = fig3.add_subplot(len(zbin),1,i+1)
    if (zbin[i]==0.1): 
        envcount_z  = all_envcount[ (all_redshift < zbin[i]+0.1)] 
        print np.max(all_redshift[ (all_redshift < zbin[i]+0.1)])
        sf_envcount = sdss_active_envcount
        q_envcount  = sdss_quiescent_envcount
    else: 
        envcount_z  = all_envcount[ (all_redshift < zbin[i]+0.1) & (all_redshift >= zbin[i]-0.1) ] 
        sf_envcount = active_envcount[ (active_redshift < zbin[i]+0.1) & (active_redshift >= zbin[i]-0.1) ] 
        q_envcount  = quiescent_envcount[ (quiescent_redshift < zbin[i]+0.1) & (quiescent_redshift >= zbin[i]-0.1) ] 
    env_hist_z  = envhist_z.hist(envcount_z, env_bin, color='k', normed=True, alpha=0.5, label='All')
   
    print len(envcount_z)
    print len(sf_envcount)+len(q_envcount)
    sf_hist_z   = dump.hist(sf_envcount, env_bin, color='b',histtype='step') 
    envhist_z.plot((sf_hist_z[1])[0:-1], sf_hist_z[0]/float(len(envcount_z)), color='b', linewidth=3, drawstyle='steps', label='active')
    q_hist_z    = dump.hist(q_envcount, env_bin,  color='r',histtype='step') 
    envhist_z.plot((q_hist_z[1])[0:-1], q_hist_z[0]/float(len(envcount_z)), color='r', linewidth=3, linestyle='--', drawstyle='steps', label='quiescent')
   
    envhist_z.set_xlim([0.0,10.0])
    envhist_z.set_ylim([0.0,0.5])
    envhist_z.vlines(highenvthresh,0.0,1.0,linewidth=2)
    envhist_z.vlines(lowenvthresh,0.0,1.0,linewidth=2)
    envhist_z.get_yaxis().set_ticklabels([0.0,0.1,0.2,0.3,0.4])
  
    env_hist_line   = fig4.add_subplot(len(zbin),1,i+1)
    env_x_z         = [ (env_hist_z[1][j]+env_hist_z[1][j+1])/2.0 for j in range(len(env_hist_z[1])-1) ] 
    env_hist_line.plot(env_x_z, env_hist_z[0], linewidth=3, label=zbin_labels[i])
    if (i!=len(zbin)-1): 
        envhist_z.get_xaxis().set_ticklabels([])
        env_hist_line.get_xaxis().set_ticklabels([])
    if (i==0): 
        envhist_z.legend(loc='lower right')
    env_hist_line.vlines(highenvthresh,0.0,1.0,linewidth=2)
    env_hist_line.vlines(lowenvthresh,0.0,1.0,linewidth=2)
    env_hist_line.set_xlim([0.0,10.0])
    env_hist_line.set_ylim([0.0,1.0])
    env_hist_line.get_yaxis().set_ticklabels([0.0,0.2,0.4,0.6,0.8])

    lowenv_frac     = float(len(envcount_z[(envcount_z < lowenvthresh)]))/float(len(envcount_z))*100.0
    highenv_frac    = float(len(envcount_z[(envcount_z >= highenvthresh)]))/float(len(envcount_z))*100.0
    midenv_frac     = 100.0-lowenv_frac-highenv_frac
    print lowenv_frac, midenv_frac, highenv_frac
    envhist_z.text(0.15,0.2, r''+str(round(lowenv_frac))+'$\%$', fontsize=20)
    envhist_z.text(lowenvthresh+1.0,0.2, r''+str(round(midenv_frac))+'$\%$', fontsize=20)
    envhist_z.text(highenvthresh+2.0,0.2, r''+str(round(highenv_frac))+'$\%$', fontsize=20)
    envhist_z.text(8.0,0.4, zbin_labels[i], fontsize=15)


fig3.text(0.5, 0.02, r'$n_{\rm{env}}$', ha='center', va='center', fontsize=20)
fig3.text(0.5, 0.98, r'Radius '+rad+'Height '+height, ha='center', va='center', fontsize=20)
fig3.text(0.03, 0.5, r'', ha='center', va='center', rotation='vertical',fontsize=20)

fig4.text(0.5, 0.02, r'$n_{\rm{env}}$', ha='center', va='center', fontsize=20)
fig4.text(0.5, 0.98, r'Radius '+rad+' Height '+height, ha='center', va='center', fontsize=20)
fig4.text(0.03, 0.5, r'', ha='center', va='center', rotation='vertical',fontsize=20)

fig3.subplots_adjust(left=0.10,right=0.90,bottom=0.05,top=0.95,wspace=0.0,hspace=0.0)
fig4.subplots_adjust(left=0.10,right=0.90,bottom=0.05,top=0.95,wspace=0.0,hspace=0.0)
py.show()
