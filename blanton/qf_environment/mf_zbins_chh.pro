function mf_zbins_chh, nzbins, sdss=sdss, zmin=zmin, zmax=zmax, $
  literature=literature
; jm10aug25ucsd - preselected redshift bins
; redshift invertals correspond to ~1 Gyr of cosmic time
    if keyword_set(sdss) then begin
       junk = get_mf_vagc_chh(zminmax=zminmax)
       zzmin = zminmax[0]
       zzmax = zminmax[1]
;      zzmin = [zminmax[0],0.1]
;      zzmax = [0.1,zminmax[1]]
    endif else begin
;      zzmin = [0.15,0.25,0.35,0.45,0.6,0.8]
;      zzmax = [0.25,0.35,0.45,0.6,0.8,1.0]  
;      zzmin = [0.2,0.3,0.4,0.55,0.7,0.9]
;      zzmax = [0.3,0.4,0.55,0.7,0.9,1.1]

       if keyword_set(literature) then begin
          zzmin = [0.2,0.4,0.6,0.8]
          zzmax = [0.4,0.6,0.8,1.0]
       endif else begin
          zzmin = [0.2,0.3,0.4,0.5,0.65,0.8]
          zzmax = [0.3,0.4,0.5,0.65,0.8,1.0]
       endelse
    endelse
    zmin = min(zzmin)
    zmax = max(zzmax)
    nzbins = n_elements(zzmin)
    zbins = replicate({zbin: 0.0, zlo: 0.0, zup: 0.0, zsigma: 0.0},nzbins)
    zbins.zlo = zzmin
    zbins.zup = zzmax
    if keyword_set(sdss) then begin
       zbins.zbin = 0.1 ; median 
       zbins.zsigma = 0.1
    endif else begin
       zbins.zbin = (zzmax-zzmin)/2.0+zzmin
       zbins.zsigma = (zbins.zup-zbins.zlo)/2.0
    endelse
return, zbins
end
