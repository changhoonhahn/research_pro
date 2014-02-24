function z_bins, nzbin, minz=minz, maxz=maxz, primus=primus 
; Sets the z bins to be used
    if keyword_set(primus) then begin 
        zlo = [ 0.2, 0.4, 0.6, 0.8 ] 
        zup = [ 0.4, 0.6, 0.8, 1.0 ]
    endif 
    nzbin = n_elements(zlo)
    zbins = replicate({zbin:0.0, zlo:0.0 , zup: 0.0, zsigma:0.0}, nzbin)

    minz = min(zlo)
    maxz = max(zup)
    
    zbins.zlo = zlo
    zbins.zup = zup
    zbins.zbin = (zup-zlo)/2.0+zlo
    zbins.zsigma = (zup-zlo)/2.0
return, zbins
end
