pro build_edgecuttest_sdss_envcount, run_sdss
; Import parameter file:
    parpath = get_path(/repo)
    parameters = yanny_readone(parpath+'zerod_environment_parameters.par', hdr=hdr)
    run_index   = where(parameters.run eq run_sdss)
    param       = parameters[run_index]
    cylrad      = param.cylrad
    cylheight   = param.cylheight
    threshold   = param.thresh
    thrshld     = float(param.thresh)/100.0

    radius_string       = strtrim(string(cylrad),2)
    height_string       = strtrim(string(cylheight),2)
    threshold_string    = strtrim(string(threshold),2)

    zbins = mf_zbins_chh(nzbins, zmin=zmin, zmax=zmax, /sdss)

    sfq = ['active','quiescent']
    sdss_ra = []
    sdss_dec = []
    sdss_redshift = []
    for i=0L,n_elements(sfq)-1L do begin
; Import Environment Count files for both active and quiescent
        envcount_dir = '/global/data/scr/chh327/primus/data/envcount/'
        envcount_file = 'EDGECUTTEST_envcount_cylr'+radius_string+'h'+height_string+$
            '_thresh'+threshold_string+'_sdss_'+sfq[i]+'_EDP-sdss-z00375_0145-numden.fits'

        envcount_data = mrdfits(envcount_dir+envcount_file, 1)
        print, envcount_dir+envcount_file
        sdss_ra = [sdss_ra, envcount_data.ra]
        sdss_dec = [sdss_dec, envcount_data.dec]
        sdss_redshift = [sdss_redshift, envcount_data.redshift]
    endfor
    sdss_combined = replicate({ra:0., dec:0., redshift:0.}, n_elements(sdss_ra))
    sdss_combined.ra = sdss_ra
    sdss_combined.dec = sdss_dec
    sdss_combined.redshift = sdss_redshift

    smallwindow = where(sdss_combined.dec LT 40.0 AND sdss_combined.dec GT 35.0 $
        AND sdss_combined.ra LT 150.0 AND sdss_combined.ra GT 100.0) 
    sdss_sub = sdss_combined[smallwindow]
    print, n_elements(sdss_sub)

    ransackfile     = 'ransack_sdss_edpfield_1000000.fits'
    ransack_data    = mrdfits(get_path(/ransack)+ransackfile,1)
    print, get_path(/ransack)+ransackfile
    ran_ra          = ransack_data.ra
    ran_dec         = ransack_data.dec
    edgecut = get_edgecut(ran_ra, ran_dec, sdss_sub ,zmin=zmin,zmax=zmax,$
        rad=cylrad,thresh=0.75,/sdss)
    output  = struct_addtags(sdss_sub, edgecut)

    fname   = get_path(/envcount)+'EDGECUTTEST_envcount_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+$
        '_sdss.fits'
    print, fname
    mwrfits, output, fname, /create
    return
end
