pro build_mf_sdss_ubersample_edp, sdss=sdss, clobber=clobber
; jm10aug25ucsd - build the "uber" sample of galaxies (for each field)
; that will be fitted by isedfit in PRIMUS_MF_ISEDFIT; select primary
; *galaxies* (i.e., ALLTARGETS=0) in the formal window function
; (NOWINDOW=0); we fit this uber-sample so that we can play around
; with different parent sample selection criteria (see
; BUILD_MF_PARENT)

; the output of this routine can be read with READ_MF_UBERSAMPLE() 
    mfpath = '/global/data/scr/chh327/primus/science/mf/2165/ubersample/edp_ubersample/'

; ##################################################    
; SDSS - build the VAGC-GALEX ubersample
    vagc_sample = get_mf_vagc_chh(sample=sample,$
        letter=letter,poststr=poststr,zminmax=zminmax)
    vagc_windowfile = '/global/data/scr/chh327/primus/science/mf/2165/'+$
        'dr72bsafe0.ply'
    print, zminmax
    post = mrdfits('/global/data/sdss/lss/dr72/bsafe/0/post_catalog.dr72bsafe0.fits',1)
    keep = where((post.z ge zminmax[0]) and (post.z le zminmax[1]))
    post = post[keep]
    ngal = n_elements(post)

    sdssphot = mrdfits('/global/data/scr/chh327/primus/science/mf/2165/ubersample/'+$
        'object_sdss_imaging.fits', row=post.object_position,1)
    galex = mrdfits('/global/data/scr/chh327/primus/science/mf/2165/ubersample/'+$
        'sdss_ubersample.fits.gz',1)
    sdsstwomass = mrdfits('/global/data/scr/chh327/primus/science/mf/2165/ubersample/'+$
        'object_twomass.fits', row=post.object_position,1)
    sdsswise = mrdfits('/global/data/scr/chh327/primus/science/mf/2165/ubersample/'+$
        'object_wise.fits', row=post.object_position,1)

; build the final catalog by keeping just the tags we want
    sample = struct_trimtags(post,select=['object_position','ra','dec','z','absm'])
    sample = struct_addtags(temporary(sample),struct_trimtags(sdssphot,$
        select=['modelflux*','petroflux*','fracpsf','devflux*',$
        'expflux*','extinction','petror*']))             


    match, galex.object_position, post.object_position, mgalex,mpost
    galexx = struct_trimtags(galex, select=['tile','fov_radius','nuv_*','fuv_*'])
    
    galex_tags = clear_struct(galexx)
    galex_vagc = replicate(galex_tags[0],n_elements(post))
    galex_vagc[mpost] = galexx[mgalex]

    sample = struct_addtags(temporary(sample),$
        struct_trimtags(galex_vagc,select=['tile','fov_radius','nuv_*','fuv_*']))
    sample = struct_addtags(temporary(sample),$
        struct_trimtags(sdsstwomass,select=['*_ext']))
    sample = struct_addtags(temporary(sample),$
        struct_trimtags(sdsswise,select=['*mpro*']))

; add the targeting mag
    sample = struct_addtags(temporary(sample),replicate({targ_mag: 0.0, final_weight: 1.0},ngal))
    sdss_to_maggies, calib=sample, mm, flux='petro'
    sample.targ_mag = reform(-2.5*alog10(mm[2,*])-(primus_sdss_zpoffset())[2]) ; targeting mag


; apply the vagc window function 
; carrying through the spectroscopic incompleteness (fgot), and the
; redshift limits
    in_window = im_is_in_window(vagc_windowfile,$
        ra=sample.ra,dec=sample.dec,/silent,polyid=allpolyid)
    final = where(in_window,ngal)
    out = sample[final]
    polyid = allpolyid[final]

    read_mangle_polygons, vagc_windowfile, win ; from the window function
    out.final_weight = 1.0/win[polyid].weight            ; fgot, from the VAGC

    im_mwrfits, out, mfpath+'sdss_ubersample.fits', /clobber
return
end
    
