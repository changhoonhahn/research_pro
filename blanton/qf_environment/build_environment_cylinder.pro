pro build_environment_cylinder,run,Nransack,primus=primus, sdss=sdss, literature=literature
    parpath = get_path(/repo)
    parameters = yanny_readone(parpath+'zerod_environment_parameters.par', hdr=hdr)
    para = where(parameters.run eq run)
    param = parameters[para]
    cylrad      = param.cylrad
    cylheight   = param.cylheight
    threshold   = param.thresh
    thrshld     = float(param.thresh)/100.0
    
    radius_string = strtrim(string(cylrad),2)
    height_string = strtrim(string(cylheight),2)
    threshold_string = strtrim(string(threshold),2)

; set log file to record
    logfile = '/home/users/hahn/research/pro/blanton/qf_environment/envcount_run.log'
    splog, file=logfile, /append
    splog, '-------------------------', systime(0)
    t00 = systime(1) 
    
    sf_q = ['active', 'quiescent']
    
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
    splog, 'Fields ', fields
    
    splog, 'zmin=',zmin,' zmax=',zmax,' cylrad=',cylrad,' cylheight=',cylheight,' threshold=',thrshld
    
    envfile = get_path(/envt)+envfname
    splog, 'EDP file'+envfile
   
    survey=''
    for i=0L,n_elements(fields)-1L do survey = survey+fields[i]+'_'
    if keyword_set(sdss) then ransackfile = 'ransack_'+survey+'edpfield_'+strtrim(string(Nransack),2)+'.fits'
    if keyword_set(primus) then ransackfile = 'ransack_'+survey+strtrim(string(Nransack),2)+'.fits'
    ransack_data    = mrdfits(get_path(/ransack)+ransackfile,1) 
    splog, 'Ransack file'+get_path(/ransack)+ransackfile
    ran_ra          = ransack_data.ra
    ran_dec         = ransack_data.dec

    for i = 0L,n_elements(sf_q)-1L do begin
        t0 = systime(1)
        targetfile = get_path(/target)+'target_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+$
            '_'+sf_q[i]+samples+litsuffix+'.fits'
        splog, sf_q[i]
        splog, 'Target file '+targetfile
        targetdata = mrdfits(targetfile, 1)
        ngal = n_elements(targetdata)

        nmatch = replicate({envcount:0.}, n_elements(targetdata))
        nmatch.envcount = get_environment_cylinder(envfile=envfile, targ=targetdata, rad=cylrad, h=cylheight)
        splog, 'Envcount mean=',mean(nmatch.envcount),' median=',median(nmatch.envcount)$
            ,' min=',min(nmatch.envcount),' max=',max(nmatch.envcount)
            
        if keyword_set(primus) then edgecut = get_edgecut(ran_ra, ran_dec, targetdata ,zmin=zmin,zmax=zmax,$
            rad=cylrad,thresh=thrshld,/primus)
        if keyword_set(sdss) then edgecut = get_edgecut(ran_ra, ran_dec, targetdata ,zmin=zmin,zmax=zmax,$
            rad=cylrad,thresh=thrshld,/sdss)
            
        output  = struct_addtags(struct_addtags(targetdata, nmatch), edgecut)
        splog, 'Edgecut mean=',mean(output.edgecut),' median=',median(output.edgecut)$
            ,' min=',min(output.edgecut),' max=',max(output.edgecut)
        if (mean(output.edgecut) EQ 0.0) then STOP 

        fname   = get_path(/envcount)+'envcount_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+$
            samples+'_'+sf_q[i]+litsuffix+'_'+envfname
        mwrfits, output, fname, /create
        splog, 'Output: ', fname 
        splog, 'Edgecut/Environment Time = ', (systime(1)-t0)/60.0
        splog, ' '
    endfor
    splog, 'Total Time=', (systime(1)-t00)/60.0
end 
