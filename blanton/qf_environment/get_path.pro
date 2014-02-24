function get_path, chh327=chh327, primus=primus, envt=envt, envcount=envcount, $
	envdist=envdist, smf=smf, powercode=powercode, chichipio=chichipio,$
	FFTdat=FFTdat, nbar=nbar, power=power, mfdata=mfdata, ransack=ransack, $
	repo=repo, parent=parent, isedfit=isedfit, target=target, vmaxavail=vmaxavail, $
        qfenvpara=qfenvpara, envsmfdata=envsmfdata, qfenvcount=qfenvcount
	
    chh327dir = '/global/data/scr/chh327/'
    homedir = '/home/users/hahn/'
    chidir = '/mount/chichipio2/hahn/'
    smfdir = '/mount/moon1/ioannis/research/projects/primus/mf/2165/mfs_v20/'
    repodir = '/home/users/hahn/research/pro/blanton/qf_environment/'
    mfdir = '/global/data/scr/chh327/primus/science/mf/2165/'
    mfarchive = '/mount/moon1/ioannis/archive/primus/mf/'
    if keyword_set(chh327) then return, chh327dir
    if keyword_set(primus) then return, chh327dir+'primus/'
    if keyword_set(envt) then return, chh327dir+'primus/data/EDP/'
    if keyword_set(powercode) then return, chh327dir+'powercode/'
    if keyword_set(smf) then return, chh327dir+'primus/data/smf/'
    if keyword_set(qfenvcount) then return, chh327dir+'primus/data/qf/'
    if keyword_set(ransack) then return, chh327dir+'primus/data/ransack/'
    if keyword_set(envcount) then return, chh327dir+'primus/data/envcount/'
    if keyword_set(envdist) then return, chh327dir+'primus/envdist/'
    if keyword_set(FFTdat) then return, chidir+'FFT/'
    if keyword_set(nbar) then return, chidir+'nbar/'
    if keyword_set(power) then return, chidir+'power/'
    if keyword_set(repo) then return, repodir
    if keyword_set(qfenvpara) then return, repodir
    if keyword_set(isedfit) then return, mfdir+'parent_v20/'  
    if keyword_set(parent) then return, mfdir+'parent_v20/'
    if keyword_set(mfdata) then return, chh327dir+'primus/mfdata/'
    if keyword_set(target) then return, chh327dir+'primus/data/target/'
    if keyword_set(vmaxavail) then return, chh327dir+'primus/data/vmax_avail/'
    if keyword_set(envsmfdata) then return, chh327dir+'primus/data/smfdata/'
end
