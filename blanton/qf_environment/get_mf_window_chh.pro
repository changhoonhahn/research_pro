function get_mf_window_chh, field, area=area, poly=poly, sr=sr, $
  irac=irac, just2mask=just2mask
; jm11jul28ucsd - return the window function

    nfield = n_elements(field)
    if nfield gt 1 then begin
       for ii = 0, nfield-1 do begin
          windowfile1 = get_mf_window_chh(field[ii],area=area1,$
            poly=poly1,sr=sr,irac=irac,just2mask=just2mask)
          if ii eq 0 then begin
             windowfile = windowfile1
             area = area1
          endif else begin
             windowfile = [windowfile,windowfile1]
             area = [area,area1]
          endelse
       endfor
       return, windowfile
    endif
    
    if field eq 'sdss' then begin
       vagc = get_mf_vagc()
       windowfile = '/global/data/scr/chh327/primus/science/mf/2165/'+vagc+'_galex_final.ply'
    endif else begin
       if keyword_set(just2mask) then begin ; ignore the galex, irac window
          suffix = ''
       endif else begin
          suffix = '_galex'
          if keyword_set(irac) or field eq 'test_cfhtls_xmm' then suffix = suffix+'_irac'
       endelse
       windowfile = '/global/data/scr/chh327/primus/survey_regions/fields/'+$
         strtrim(field,2)+'_field'+suffix+'_window_2mask.ply'
    endelse
    if (file_test(windowfile) eq 0) then $
      message, 'No GALEX window for field '+field+'!'+windowfile

    if arg_present(area) then begin
       read_mangle_polygons, windowfile, win
       area = total(win.str,/double)
       if (keyword_set(sr) eq 0) then area *= (180.0/!pi)^2 ; [deg^2]
    endif

return, windowfile
end
