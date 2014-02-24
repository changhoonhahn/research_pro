function get_zbin, z, moustakas=moustakas, sdss=sdss
    if keyword_set(moustakas) then zbin = [0.01,0.2,0.3,0.4,0.5,0.65,0.8]
    if keyword_Set(primus) then zbin = [0.2, 0.4, 0.6, 0.8, 1.0]
    if keyword_Set(sdss) then zbin = [0.06, 0.145] 

    bin = fltarr(2)
    if (z lt zbin[0] or z gt zbin[n_elements(zbin)-1L]) then begin
        print, "Error! Out of z range!"
        bin[0] = zbin[0]
        bin[1] = zbin[0]
    endif else begin
        for i=0L,n_elements(zbin)-2L do begin 
            if (z ge zbin[i] and z lt zbin[i+1]) then begin
                bin[0] = zbin[i]
                bin[1] = zbin[i+1]
            endif
        endfor
    endelse
return, bin
end

