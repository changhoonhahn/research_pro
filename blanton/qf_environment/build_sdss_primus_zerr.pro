pro build_sdss_primus_zerr, run,seed1,seed2
; Generated aritifical redshift deviations based on PRIMUS redshift errors. 

    parpath = get_path(/repo)
    parameters = yanny_readone(parpath+'zerod_environment_parameters.par', hdr=hdr)
    para = where(parameters.run eq run)
    param = parameters[para]
    cylrad      = param.cylrad
    cylheight   = param.cylheight
    nbin        = param.nbin
    threshold   = param.thresh

    radius_string = strtrim(string(cylrad),2)
    height_string = strtrim(string(cylheight),2)
    threshold_string = strtrim(string(threshold),2)
    nbin_string = strtrim(string(nbin),2)

    envfname = 'EDP-sdss-z00375_0145-numden.fits'

    sf_q = ['active','quiescent']

    sigma_z = 0.005
     
    for i=0L,n_elements(sf_q)-1L do begin 
        envcount_file = get_path(/envcount)+'envcount_cylr'+radius_string+'h'+height_string+$
            '_thresh'+threshold_string+'_nbin'+nbin_string+'_sdss'+'_'+sf_q[i]+'_'+envfname
        envcount_data = mrdfits(envcount_file,1)
        redshift = envcount_data.redshift

        for j=0L,n_elements(redshift)-1L do begin 
            dz = (-3.0+randomu(seed1)*6.0)*sigma_z
            random2 = randomu(seed2)
            pofz = exp(-0.5*(dz)^2/sigma_z^2)
            seed1 = seed1+1
            seed2 = seed2+1
            while (pofz LE random2) do begin 
                dz = (-3.0+randomu(seed1)*6.0)*sigma_z
                random2 = randomu(seed2)
                pofz = exp(-0.5*(dz)^2/sigma_z^2)
                seed1 = seed1+1
                seed2 = seed2+1
            endwhile 
            redshift[j] = redshift[j]+dz
        endfor 
        envcount_data.redshift = redshift
        sdss_primuszerr_file = get_path(/envcount)+'envcount_cylr'+radius_string+'h'+height_string+$
                '_thresh'+threshold_string+'_nbin'+nbin_string+'_sdss'+'_'+sf_q[i]+'_primuszerr_'+$
                envfname
        mwrfits, envcount_data, sdss_primuszerr_file, /create
    endfor
end
