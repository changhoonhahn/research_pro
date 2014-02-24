pro build_vmax_avail,run,n_rsk,n_ran,primus=primus,sdss=sdss
; Generates tables for V_max,avail interpolation later used in build_target_sample.pro
    parpath = get_path(/repo)
    parameters = yanny_readone(parpath+'zerod_environment_parameters.par', hdr=hdr)
    para = where(parameters.run eq run)
    param = parameters[para]
    zmin = param.zmin 
    zmax = param.zmax
    cylrad = param.cylrad
    cylheight = param.cylheight
    threshold = param.thresh
    thrshld = float(param.thresh)/100.0
    nz = 50 
    
    radius_string       = strtrim(string(cylrad),2)
    height_string       = strtrim(string(cylheight),2)
    threshold_string    = strtrim(string(threshold),2)
    ran_string          = strtrim(string(n_ran),1)
    rsk_string          = strtrim(string(n_rsk),1)

; set log file to record    
    logfile = '/home/users/hahn/research/pro/blanton/qf_environment/envcount_run.log'
    splog, file=logfile, /append
    splog, '-------------------------', systime(0)
    t0 = systime(1)

; set survey and field
    if keyword_set(primus) then fields = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire']
    if keyword_set(sdss) then fields = ['sdss']
    survey = ''
    for i=0L,n_elements(fields)-1L do begin 
        survey = survey+fields[i]+'_'
    endfor
    splog, 'survey: '+survey
    splog, 'R = '+radius_string+' H = '+height_string+' Threshold = '+threshold_string

; import ransack and random files
    rsk_dir = get_path(/ransack)
    if keyword_set(sdss) then begin
        ransack_fname   = rsk_dir+'ransack_'+survey+'edpfield_'+strtrim(string(n_rsk),1)+'.fits'
        random_fname    = rsk_dir+'random_'+survey+'targetfield_'+strtrim(string(n_ran),1)+'.fits'
    endif else begin  
        ransack_fname   = rsk_dir+'ransack_'+survey+strtrim(string(n_rsk),1)+'.fits'
        random_fname    = rsk_dir+'random_'+survey+strtrim(string(n_ran),1)+'.fits'
    endelse 
    splog, 'Ransack and Random file used:'
    splog, ransack_fname
    splog, random_fname
  
    ransack = mrdfits(ransack_fname,1)
    random  = mrdfits(random_fname,1)
    print, ransack_fname
    print, random_fname

    ran_vol = randomn(rseed,n_ran,/uniform) 
    ran_z   = dblarr(n_ran)

    vmax_vmin   = lf_comvol(zmax)-lf_comvol(zmin)
    vmin        = lf_comvol(zmin)
    for i=0L,n_ran-1L do begin  
        ran_vol[i]  = ran_vol[i]*(vmax_vmin)+vmin
        ran_z[i]    = lf_vtoz(ran_vol[i])
    endfor    
    randomdata = replicate({ra:0.,dec:0.,redshift:0.},n_ran) 
    randomdata.ra       = random.ra
    randomdata.dec      = random.dec
    randomdata.redshift  = ran_z

    rsk_ra = ransack.ra
    rsk_dec= ransack.dec

    edgecut = get_edgecut(rsk_ra,rsk_dec,randomdata,zmin=zmin,zmax=zmax,$
        rad=cylrad,thresh=thrshld,primus=primus,sdss=sdss)
    ran_data = struct_addtags(randomdata,edgecut)

    zbin     = dblarr(nz)
    eff_vmax = dblarr(nz)
    for i=0L,nz-1L do begin 
        zbin[i] = zmin+float(i+1L)*(zmax-zmin)/float(nz)
        z_indx  = ran_data.redshift LE zbin[i] AND ran_data.redshift GE zmin
        edgecut_indx = ran_data.edgecut eq 1 ; edgecut = 1 then galaxy is within the edge

        edgecut_pts     = ran_data[where(z_indx AND edgecut_indx, edgecut_count)] 
        unedgecut_pts   = ran_data[where(z_indx, unedgecut_count)]
        eff_vmax[i] = float(edgecut_count)/float(unedgecut_count)*lf_comvol(zbin[i])
        print, 'edgecut_count: ', edgecut_count, '  un-edgecut_count: ', unedgecut_count
        print, zbin[i],eff_vmax[i]/lf_comvol(zbin[i]),float(edgecut_count)/float(unedgecut_count)
    endfor 
    out = replicate({z:0.,v_max_avail:0.},nz)
    out.z = zbin
    out.v_max_avail = eff_vmax

    outname = get_path(/vmaxavail)+'vmax_avail_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+$
        '_'+survey+'ran'+ran_string+'_rsk'+rsk_string+'.fits'
    splog, 'Output: ', outname
    splog, 'Total time = ', (systime(1)-t0)/60.0
    splog, ' ' 
    mwrfits,out,outname,/create
    splog, /close
end 
