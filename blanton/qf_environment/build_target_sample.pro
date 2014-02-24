function get_vmax_avail, z
    common vmaxavail_shared, vz
    vmaxavail = interpol(vz.v_max_avail,vz.z,z,/spline)
    return, vmaxavail
end 

pro build_target_sample,run,Nransack,Nrandom,literature=literature, primus=primus, sdss=sdss
    common vmaxavail_shared, vz
    mfdatapath  = get_path(/mfdata)
    parpath = get_path(/repo)
    parameters = yanny_readone(parpath+'zerod_environment_parameters.par', hdr=hdr)
    para = where(parameters.run eq run)
    param= parameters[para]
    cylrad      = param.cylrad
    cylheight   = param.cylheight
    threshold   = param.thresh
    
    radius_string   = strtrim(string(cylrad),2)
    height_string   = strtrim(string(cylheight),2)
    threshold_string= strtrim(string(threshold),2)
    rsk_string      = strtrim(string(Nransack),1)
    ran_string      = strtrim(string(Nrandom),1)
   
    h100 = mf_h100()

; set log file to record
    logfile = '/home/users/hahn/research/pro/blanton/qf_environment/envcount_run.log'
    splog, file=logfile, /append
    splog, '-------------------------', systime(0)
    t0 = systime(1)

    if keyword_set(primus) then begin
        if keyword_set(literature) then litsuffix = '_lit' $
            else litsuffix = '' 
        if keyword_set(literature) then zbins = mf_zbins(nzbins, zmin=zmin, zmax=zmax, /literature) $
            else zbins = mf_zbins(nzbins, zmin=zmin, zmax=zmax)
        jm_zbins = zbins
        jm_zmin = zmin 
        jm_zmax = zmax
        area = get_poly_area(/primus,/sr)
        samples = ''
        fields = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire'] 
    endif
    if keyword_set(sdss) then begin
        litsuffix = ''
        zbins = mf_zbins_chh(nzbins, zmin=zmin, zmax=zmax, /sdss) 
        jm_zbins = mf_zbins(jm_nzbins,zmin=jm_zmin,zmax=jm_zmax, /sdss)
        area = get_poly_area(/sdss,/target,/sr)
        samples = '_sdss'
        fields = ['sdss']
    endif 
    splog, 'Fields ', fields

    struct_print, zbins
    struct_print, jm_zbins

    sf_q   = ['all', 'active', 'quiescent']

    survey=''
    for i=0L,n_elements(fields)-1L do begin
        survey=survey+fields[i]+'_'
    endfor

;Vmax,avail file generated from build_vmax_avail.pro
    vmaxavail_name=get_path(/vmaxavail)+'vmax_avail_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+$
        '_'+survey+'ran'+ran_string+'_rsk'+rsk_string+'.fits'
    splog, 'Vmaxavail file '+vmaxavail_name
    print, vmaxavail_name
    vz = mrdfits(vmaxavail_name,1)

    for i =0L, n_elements(sf_q)-1L do begin 
        mfdatafile = mfdatapath+'mfdata_chh/new/env_zrange/mfdata_'+sf_q[i]+'_supergrid01'+samples+litsuffix+'.fits.gz'
        splog, 'mfdata file'+mfdatafile
        print, mfdatafile
        mfdata = mrdfits(mfdatafile,1)

        zlim = mfdata.z GE zmin AND mfdata.z LT zmax
        data = mfdata[where(zlim)]
        ngal = n_elements(data)	
        data = struct_addtags(data,replicate({vmax:0.,vmaxavail:0.},ngal))

        for j=0L,nzbins-1L do begin 
            zindx = where(data.z GT zbins[j].zlo AND data.z LE zbins[j].zup,zindx_count)
            data_zbin = data[zindx]

            jm_zmin_evol = data_zbin.zmin_evol>jm_zbins[j].zlo
            jm_zmax_evol = data_zbin.zmax_evol<jm_zbins[j].zup
            zmin_evol = data[zindx].zmin_evol>zbins[j].zlo
            zmax_evol = data[zindx].zmax_evol<zbins[j].zup
            vmax_evol = data[zindx].vmax_evol
;                data[zindx].vmax        = (vmax_evol/(lf_comvol(jm_zmax_evol)-lf_comvol(jm_zmin_evol)))$
;                    *(lf_comvol(zmax_evol)-lf_comvol(zmin_evol))
;            data[zindx].vmaxavail   = (vmax_evol/(lf_comvol(jm_zmax_evol)-lf_comvol(jm_zmin_evol))) $
;                *(get_vmax_avail(zmax_evol,vz)-get_vmax_avail(zmin_evol,vz))

            data[zindx].vmax = area/3.0*(lf_comvol(zmax_evol)-lf_comvol(zmin_evol))*(1/h100)^3.0
            data[zindx].vmaxavail = (area/3.0)*(get_vmax_avail(zmax_evol)-get_vmax_avail(zmin_evol))*(1/h100)^3.0
            if (min(data[zindx].vmaxavail) LT 0.0) then begin
                print, "ERROR Negative Vmax,avail"
                print, "Try using fewer points for Vmax,avail interpol"
                STOP
            endif 
        endfor 
           
        target = replicate({ra:0.D,dec:0.D,redshift:0.,mr_01:0.,mg_01:0.,mb_00:0.,mass:0.,$
            masslimit:0.,age:0.,SFR:0.,zmin_evol:0.,zmax_evol:0., $ 
            field:' ',weight:0.,vmax_evol:0.,vmax:0.,vmaxavail:0.},ngal)
        if keyword_set(sdss) then target.field = strtrim('sdss',2) $
            else target.field = strtrim(data.field,2)
        target.ra       = data.ra
        target.dec      = data.dec
        target.redshift = data.z
        target.mr_01    = data.mr_01
        target.mg_01    = data.mr_01+data.gmr_01
        target.mb_00    = data.mb_00
        target.mass     = data.mass
        target.masslimit= data.masslimit
        target.age      = data.age
        target.SFR      = data.SFR
        target.weight   = data.weight
        target.zmin_evol= data.zmin_evol
        target.zmax_evol= data.zmax_evol
        target.vmax_evol= data.vmax_evol
        target.vmax     = data.vmax
        target.vmaxavail= data.vmaxavail
        
        dir = get_path(/target)
        fname = dir+'target_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+$
            '_'+sf_q[i]+samples+litsuffix+'.fits'
        mwrfits, target, fname, /create
        splog, 'Output: ', fname
        splog, 'Total time = ', (systime(1)-t0)/60.0
        splog, ' '
    endfor 
    splog, /close
end
