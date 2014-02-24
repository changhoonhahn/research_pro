function get_poly_area, primus=primus, sdss=sdss, sr=sr, edp=edp, target=target, field=field, $
    silent=silent
; Get the area of PRIMUS or SDSS polygons
; if /sr is specific it returns answer in steradians 
    if keyword_set(primus) then begin  
        dir ='/global/data/scr/chh327/primus/survey_regions/fields/'
        fields = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire']
        if keyword_set(field) then fields = [field]
        areas=dblarr(n_elements(fields))

        for i=0L,n_elements(fields)-1L do begin
            read_mangle_polygons, dir+fields[i]+'_field_galex_window_2mask.ply',polygon
            areas[i]=total(polygon.str,/double)
        endfor

        totarea=total(areas)
        if keyword_set(sr) then return, totarea $
            else return, (totarea)*(180.0/!PI)^2
    endif

    if keyword_set(sdss) then begin
;        dir = '/mount/moon1/ioannis/research/projects/primus/mf/2165/'
        dir = '/global/data/scr/chh327/primus/science/mf/2165/'
        if keyword_set(edp) then read_mangle_polygons, dir+'dr72bsafe0.ply', sdss
        if keyword_set(target) then read_mangle_polygons, dir+'dr72bsafe0_galex_final.ply', sdss
        area = total(sdss.str, /double)
        if keyword_set(sr) then return, area $
            else return, area*(180.0/!PI)^2
    endif 
end

