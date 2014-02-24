pro mf_sdss_kcorrect_edp, sdss=sdss, clobber=clobber, $
  just_ugrizkhk_01=just_ugrizjhk_01
; compute K-corrections with and without the zeropoint
; corrections 

    mfpath = '/global/data/scr/chh327/primus/science/mf/kcorrect/edp_kcorrect/'

    field = ['sdss'] 
    nfield = n_elements(field)

    prefix = '_kcorr'
    suffix = ''
    
    kcorrfile = mfpath+field+prefix+suffix+'.fits'
    if im_file_test(kcorrfile+'.gz',clobber=clobber) then return
       
    maggies = get_mf_sdss_photometry_edp(field,zobj=zobj,filters=filters,$
        ivarmaggies=ivarmaggies,nozptoffsets=keyword_set(sdss))

    junk = mf_maglimit(field,targ_filter=targ_filter)
    kcorr = mf_do_sdss_kcorrect_edp(zobj,maggies,ivarmaggies,filterlist=filters,$
        targ_filter=targ_filter,just_ugrizkhk_01=just_ugrizjhk_01,$
        maxiter=maxiter,sdss=sdss)

    im_mwrfits, kcorr, kcorrfile, clobber=clobber
return
end
