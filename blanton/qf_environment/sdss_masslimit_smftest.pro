pro sdss_masslimit_smftest
    mfdata_dir = "/global/data/scr/chh327/primus/mfdata/"
    target_dir = "/global/data/scr/chh327/primus/data/target/"
    
    mfdata_fname = "mfdata_all_supergrid01_sdss_lit.fits.gz"
    target_fname = "target_cylr2h50_thresh75_nbin5_all_sdss_lit.fits"

    mfdata = mrdfits(mfdata_dir+mfdata_fname,1)
    target = mrdfits(target_dir+target_fname,1)

    binsize = mf_binsize(bin=0.1,minmass=minmass,maxmass=maxmass)

    mfdata_smf = im_mf_vmax(mfdata.mass,mfdata.weight/mfdata.vmax_evol,binsize=binsize,minmass=minmass,maxmass=maxmass,masslimit=mfdata.masslimit)

    target_smf = im_mf_vmax(target.mass,target.weight/target.vmax,binsize=binsize,minmass=minmass,maxmass=maxmass,masslimit=target.masslimit)
    target_incorrect_smf = im_mf_vmax(target.mass,target.weight/target.vmax_evol,binsize=binsize,minmass=minmass,maxmass=maxmass,masslimit=target.masslimit)
    target_noedgecut_smf = im_mf_vmax(target.mass,target.weight/target.vmaxavail,binsize=binsize,minmass=minmass,maxmass=maxmass,masslimit=target.masslimit)

    for i=0L,n_elements(mfdata_smf.mass)-1L do begin 
        print, (mfdata_smf.mass)[i], (mfdata_smf.phi)[i], (target_smf.phi)[i], (target_incorrect_smf.phi)[i]
    endfor

    mwrfits,mfdata_smf,"/home/users/hahn/qf_environment/dump/sdss_masslimit_smftest/mfdata_smf.fits",/create
    mwrfits,target_smf,"/home/users/hahn/qf_environment/dump/sdss_masslimit_smftest/target_smf.fits",/create
    mwrfits,target_incorrect_smf,"/home/users/hahn/qf_environment/dump/sdss_masslimit_smftest/target_incorrect_smf.fits",/create
    mwrfits,target_noedgecut_smf,"/home/users/hahn/qf_environment/dump/sdss_masslimit_smftest/target_noedgecut_smf.fits",/create
end
