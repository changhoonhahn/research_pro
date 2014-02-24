function get_vmax_avail,z,r=r,h=h,thr=thr,nbin=nbin,primus=primus, sdss=sdss
    path=get_path(/ransack)
    radius_string = strtrim(string(r),2)
    height_string = strtrim(string(h),2)
    threshold_string = strtrim(string(thr),2)
    nbin_string = strtrim(string(nbin),2)
    ran_string = strtrim(string(10000),1)
    rsk_string = strtrim(string(10000),1)
    if keyword_set(primus) then fields = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire']
    if keyword_set(sdss) then fields = ['sdss']

    survey=''
    for i=0L,n_elements(fields)-1L do begin
        survey=survey+fields[i]+'_'
    endfor

    if keyword_set(primus) then begin 
        fname='vmax_avail_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+$
            '_nbin'+nbin_string+'_'+survey+'ran'+ran_string+'_rsk'+rsk_string+'.fits'
    endif 
    if keyword_set(sdss) then begin 
        fname='sdss_v_max_avail.fits'
    endif
    
    vz=mrdfits(path+'vmax_avail/'+fname,1)
    vmaxavail=interpol(vz.v_max_avail,vz.z, z,/spline)

    return, vmaxavail
end
