function get_poly_area_old, smf=smf, sdss=sdss, sr=sr, field=field

    if keyword_set(smf) then begin 
        dir ='/global/data/scr/chh327/primus/survey_regions/fields/'
        read_mangle_polygons, dir+'cdfs_field_galex_window_2mask.ply', cdfs
        read_mangle_polygons, dir+'cfhtls_xmm_field_galex_window_2mask.ply', cfhtls
        read_mangle_polygons, dir+'cosmos_field_galex_window_2mask.ply', cosmos
        read_mangle_polygons, dir+'es1_field_galex_window_2mask.ply', es1
        read_mangle_polygons, dir+'xmm_swire_field_galex_window_2mask.ply', xmmswire

        area_cdfs = total(cdfs.str, /double)
        area_cfhtls = total(cfhtls.str, /double)
        area_cosmos = total(cosmos.str, /double)
        area_es1 = total(es1.str, /double)
        area_xmmswire = total(xmmswire.str, /double)

        area = area_cdfs+area_cosmos+area_es1+area_cfhtls+area_xmmswire
    
        if keyword_set(sr) then begin
            return, area
        endif else begin
            return, area*(180.0/!PI)^2
        endelse 
    endif

    if keyword_set(sdss) then begin
        dir = '/mount/moon1/ioannis/research/projects/primus/mf/2165/'
        read_mangle_polygons, dir+'dr72bsafe0_galex_final.ply', sdss
        area = total(sdss.str, /double)
        if keyword_set(sr) then begin
            return, area
        endif else begin
            return, area*(180.0/!PI)^2
        endelse 
    endif

    if n_elements(field) ne 0 then begin
        dir ='/global/data/primus/survey_regions/fields/'
        read_mangle_polygons, dir+'cdfs_field_galex_window_2mask.ply', cdfs
        read_mangle_polygons, dir+'cfhtls_xmm_field_galex_window_2mask.ply', cfhtls
        read_mangle_polygons, dir+'cosmos_field_galex_window_2mask.ply', cosmos
        read_mangle_polygons, dir+'es1_field_galex_window_2mask.ply', es1
        read_mangle_polygons, dir+'xmm_swire_field_galex_window_2mask.ply', xmmswire
        
        fields = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire']
        areas = dblarr(n_elements(fields))
        areas[0] = total(es1.str, /double)
        areas[1] = total(cosmos.str, /double)
        areas[2] = total(cfhtls.str, /double)
        areas[3] = total(cdfs.str, /double)
        areas[4] = total(xmmswire.str, /double)

        indx = where(fields eq field)
        if keyword_set(sr) then begin 
            return, areas[indx]       
        endif else begin
            return, areas[indx]*(180.0/!PI)^2
        endelse
    endif 
end
