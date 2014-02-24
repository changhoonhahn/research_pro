function get_umr_bimodality_coeff_chh, z, mrpivot=mrpivot, sdss_redpeak=sdss_redpeak, $
  sdss_greenvalley=sdss_greenvalley
; read the output from mf_umr_bimodality

    mfpath = '/mount/moon1/ioannis/research/projects/primus/mf/2165/'
    final = mrdfits(mfpath+'umr_evolution.fits.gz',1,/silent)

; get the intercept of the minimum in the bimodality at redshift Z
    ngal = n_elements(z)
    mrpivot = final.mrpivot
    umrmin = poly(z-final.zpivot,final.umrmin_coeff)
    coeff = transpose([[umrmin],[replicate(final.slope,ngal)]])

; for comparison return my measurement of the SDSS green valley and
; the Wyder+07 coefficients 
    sdss = mrdfits(mfpath+'sdss_umr_redpeak.fits.gz',1,/silent)
    sdss_greenvalley = sdss.greenvalley_coeff
    sdss_redpeak = sdss.redpeak_coeff
;   wyder_redpeak = [1.9+mrpivot*(-0.175),-0.175]
    
return, coeff
end

function mf_select_umr_red_chh, umr, mr, z=z, umr_cut=umr_cut, $
  coeff=coeff, sdss_redpeak=sdss_redpeak, sdss_greenvalley=sdss_greenvalley, $
  mrpivot=mrpivot, blue=blue
; jm11sep21ucsd - select quiescent and bluely star-forming galaxies
;   using the evolving 0.1(NUV-r) color cut (see fit_umr_bimodality)

; get the coefficients of the bimodality minimum for every redshift
    nzz = n_elements(z)
    coeff = get_umr_bimodality_coeff_chh(z,sdss_redpeak=sdss_redpeak, $
      sdss_greenvalley=sdss_greenvalley,mrpivot=mrpivot)
    
    ngal = n_elements(umr)
    if (ngal eq 0L) then return, -1 ; just return the coefficients
    if (n_elements(mr) ne ngal) then begin
       splog, 'Dimensions of UMR and MR do not agree!'
       return, -1
    endif

; evolving sloping color-cut of the form (u-z) = (u-z)_0 + a*(Mass+21) + b*(z-0.1)
    if (nzz eq 1) then umr_cut = poly(mr-mrpivot,coeff) else begin
       umr_cut = mr*0.0
       for ii = 0L, ngal-1 do umr_cut[ii] = poly(mr[ii]-mrpivot,coeff[*,ii])
    endelse
    
    iquiescent = where(umr gt umr_cut,comp=iblue)
    if keyword_set(blue) then indx = iblue else $
      indx = iquiescent
    
return, indx
end
