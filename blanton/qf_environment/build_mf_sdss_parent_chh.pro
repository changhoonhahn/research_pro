pro build_mf_sdss_parent_chh, sdss=sdss, clobber=clobber
    mfpath = '/global/data/scr/chh327/primus/science/mf/2165/parent_v20/edp_parent/'
    field = 'sdss'
    realfield = 'sdss'
    nfield = n_elements(field)

    zbins = mf_zbins(sdss=sdss,zmin=zmin,zmax=zmax)
    
; loop on each field    
    t0 = systime(1)
    outfile = mfpath+'mf_parent_'+field+'.fits'

    uber = read_mf_ubersample_chh(realfield,photo=photo,xray=xray)
    kcorr = read_mf_kcorrect_chh(realfield)

    filters = get_mf_filters_chh(field)
    faint = mf_maglimit(field,bright=bright,targ_filter=targ_filter,$
        ch1faint=ch1faint,ch1bright=ch1bright,nice_targ_filter=nice)
    nuber = n_elements(uber)

; select the parent sample          
    parent = mf_select_parent(uber,field=field,realfield=realfield,$
        photo=photo,logfile=logfile,sdss=sdss)
    ngal = n_elements(parent)

    sample = uber[parent]
    sample = struct_addtags(replicate({uber_position: 0L},ngal),temporary(sample))
    sample.uber_position = parent
    sample = struct_addtags(temporary(sample),kcorr[parent])

; compute the statistical weights (the statistical weight for the SDSS
; sample is computed in BUILD_MF_UBERSAMPLE); also change
; .Z-->.ZPRIMUS for the SDSS sample
    alltags = tag_names(sample) & newtags = alltags
    newtags[where(strlowcase(strtrim(alltags,2)) eq 'z')] = 'zprimus'
    sample = im_struct_trimtags(temporary(sample),select=alltags,newtags=newtags)
    bad = where((sample.final_weight le 0) or (finite(sample.final_weight) eq 0))
    if bad[0] ne -1 then stop

    ; compute ZMIN and ZMAX with and without evolution
    moretags = {$
        zmin_opt_noevol: -999.0, zmax_opt_noevol: -999.0, zmin_opt_evol: -999.0, zmax_opt_evol: -999.0, $
        zmin_ch1_noevol: -999.0, zmax_ch1_noevol: -999.0, zmin_ch1_evol: -999.0, zmax_ch1_evol: -999.0, $
        zmin_noevol: -999.0, zmax_noevol: -999.0, zmin_evol: -999.0, zmax_evol: -999.0}
    moretags = replicate(moretags,ngal)
    sample = struct_addtags(temporary(sample),moretags)
    
    evol_opt = im_zminzmax(sample.zprimus,sample.targ_mag,sample.k_coeffs,$
        offset=evol_opt_offset,vname=vname,h100=h100,bright=bright,faint=faint,$
        filter=targ_filter,q0=q0_opt,qz0=qz0,debug=debug,sample_zmin=zmin,sample_zmax=zmax)
    sample.zmin_opt_evol = evol_opt.zmin
    sample.zmax_opt_evol = evol_opt.zmax

    noevol_opt = im_zminzmax(sample.zprimus,sample.targ_mag,sample.k_coeffs,$
        offset=noevol_opt_offset,vname=vname,h100=h100,bright=bright,faint=faint,$
        filter=targ_filter,q0=0.0,qz0=qz0,debug=debug,sample_zmin=zmin,sample_zmax=zmax)
    sample.zmin_opt_noevol = noevol_opt.zmin
    sample.zmax_opt_noevol = noevol_opt.zmax

    if (ch1faint gt -900.0) then begin
        splog, 'Computing [ch1] ZMIN, ZMAX'
        noevol_ch1 = im_zminzmax(sample.zprimus,sample.ch1mag,sample.k_coeffs,$
            offset=noevol_ch1_offset,vname=vname,h100=h100,bright=ch1bright,faint=ch1faint,$
            filter='spitzer_irac_ch1.par',q0=0.0,qz0=qz0,debug=debug,sample_zmin=zmin,sample_zmax=zmax)
        sample.zmin_ch1_noevol = noevol_ch1.zmin
        sample.zmax_ch1_noevol = noevol_ch1.zmax

        evol_ch1 = im_zminzmax(sample.zprimus,sample.ch1mag,sample.k_coeffs,$
            offset=evol_ch1_offset,vname=vname,h100=h100,bright=ch1bright,faint=ch1faint,$
            filter='spitzer_irac_ch1.par',q0=q0_ch1,qz0=qz0,debug=debug,sample_zmin=zmin,sample_zmax=zmax)
        sample.zmin_ch1_evol = evol_ch1.zmin
        sample.zmax_ch1_evol = evol_ch1.zmax

        sample.zmin_noevol = sample.zmin_opt_noevol>sample.zmin_ch1_noevol
        sample.zmax_noevol = sample.zmax_opt_noevol<sample.zmax_ch1_noevol
        sample.zmin_evol = sample.zmin_opt_evol>sample.zmin_ch1_evol
        sample.zmax_evol = sample.zmax_opt_evol<sample.zmax_ch1_evol
    endif else begin
        sample.zmin_noevol = sample.zmin_opt_noevol
        sample.zmax_noevol = sample.zmax_opt_noevol
        sample.zmin_evol = sample.zmin_opt_evol
        sample.zmax_evol = sample.zmax_opt_evol
    endelse

    im_mwrfits, sample, outfile, clobber=clobber

   print, 'Total time = ', (systime(1)-t0)/60.0
return
end
