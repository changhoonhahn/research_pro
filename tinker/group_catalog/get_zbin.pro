function get_zbin, z, zlo=zlo, zup=zup,$
    primus=primus
; Finds the z bins that the z belongs to
    zbins = z_bins(nzbins,minz=minz,maxz=maxz,primus=primus)
    ;struct_print, zbins
    
    i_match = where((zbins.zlo LT z) and (zbins.zup GT z),n_match)
    if (n_match NE 1) then begin 
        print, "OUT OF RANGE"
        STOP     ; Sanity check to ensure that there's not more than one match
    endif

    zlo = zbins[i_match].zlo
    zup = zbins[i_match].zup
    return, [zlo, zup]
end 
