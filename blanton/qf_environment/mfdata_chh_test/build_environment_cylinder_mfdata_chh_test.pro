pro build_environment_cylinder_mfdata_chh_test,run,Nransack,primus=primus, sdss=sdss, literature=literature
    sf_q = ['all', 'active', 'quiescent']

    t0 = systime(1)
    print, t0
    
    parpath = get_path(/repo)
    parameters = yanny_readone(parpath+'zerod_environment_parameters.par', hdr=hdr)
    para = where(parameters.run eq run)
    param = parameters[para]
    cylrad      = param.cylrad
    cylheight   = param.cylheight
    nbin        = param.nbin
    threshold   = param.thresh
    thrshld     = float(param.thresh)/100.0
    
    radius_string = strtrim(string(cylrad),2)
    height_string = strtrim(string(cylheight),2)
    threshold_string = strtrim(string(threshold),2)
    nbin_string = strtrim(string(nbin),2)
    
    if keyword_set(primus) then begin 
        samples = ''
        fields = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire']
        if keyword_set(literature) then litsuffix = '_lit' $
            else litsuffix = ''
        if keyword_set(literature) then zbins = mf_zbins(nzbins, zmin=zmin, zmax=zmax, /literature) $
            else zbins = mf_zbins(nzbins, zmin=zmin, zmax=zmax) 
        envfname = 'EDP-primus-z0210-numden.fits' 
    endif 
    if keyword_set(sdss) then begin
        samples = '_sdss'
        fields = ['sdss']
        litsuffix = ''
        zbins = mf_zbins_chh(nzbins, zmin=zmin, zmax=zmax, /sdss)
        envfname = 'EDP-sdss-z00375_0145-numden.fits'
    endif 
    
    print, 'zmin=',zmin,' zmax=',zmax,' cylrad=',cylrad,' cylheight=',cylheight,' nbin=',nbin,' threshold=',thrshld
    
    envfile = get_path(/envt)+envfname
    print, envfile
   
    survey=''
    for i=0L,n_elements(fields)-1L do survey = survey+fields[i]+'_'
    ransackfile     = 'ransack_'+survey+strtrim(string(Nransack),2)+'.fits'
    print, ransackfile
    ransack_data    = mrdfits(get_path(/ransack)+ransackfile,1) 
    ran_ra          = ransack_data.ra
    ran_dec         = ransack_data.dec

    for i = 0L,n_elements(sf_q)-1L do begin
        mffile = get_path(/target)+'mfdata_chh_test/target_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+$
            '_nbin'+nbin_string+'_'+sf_q[i]+samples+litsuffix+'.fits'
            print, mffile
        mfdata = mrdfits(mffile, 1)
        ngal = n_elements(mfdata)
        target = replicate({ra:0.D, dec:0.D, zprimus:0.}, ngal)

        target.ra       = mfdata.ra
        target.dec      = mfdata.dec
        target.zprimus  = mfdata.redshift 

        nmatch = replicate({envcount:0.}, n_elements(target))
        nmatch.envcount = get_environment_cylinder(envfile=envfile, targ=target, rad=cylrad, h=cylheight)
            
        if keyword_set(primus) then edgecut = get_edgecut(ran_ra, ran_dec, target ,zmin=zmin,zmax=zmax,$
            rad=cylrad,nbin=nbin,thresh=thrshld,/primus)
        if keyword_set(sdss) then edgecut = get_edgecut(ran_ra, ran_dec, target ,zmin=zmin,zmax=zmax,$
            rad=cylrad,nbin=nbin,thresh=thrshld,/sdss)
            
        print, 'time=', (systime(1)-t0)/60.0
       
        output  = struct_addtags(struct_addtags(mfdata, nmatch), edgecut)
        fname   = get_path(/envcount)+'mfdata_chh_test/envcount_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+'_nbin'$
            +nbin_string+samples+'_'+sf_q[i]+litsuffix+'_'+envfname
        print, fname
        mwrfits, output, fname, /create
    endfor
    print, 'total time=', (systime(1)-t0)/60.0
end 
