function get_gmr_bimodality_coeff, z, mrpivot=mrpivot, sdss_redpeak=sdss_redpeak, $
  sdss_greenvalley=sdss_greenvalley
; read the output from mf_gmr_bimodality

    mfpath = mf_path()
    final = mrdfits(mfpath+'gmr_evolution.fits.gz',1,/silent)

; get the intercept of the minimum in the bimodality at redshift Z
    ngal = n_elements(z)
    mrpivot = final.mrpivot
    gmrmin = poly(z-final.zpivot,final.gmrmin_coeff)
    coeff = transpose([[gmrmin],[replicate(final.slope,ngal)]])

; for comparison return my measurement of the SDSS green valley and
; the Wyder+07 coefficients 
    sdss = mrdfits(mfpath+'sdss_gmr_redpeak.fits.gz',1,/silent)
    sdss_greenvalley = sdss.greenvalley_coeff
    sdss_redpeak = sdss.redpeak_coeff
;   wyder_redpeak = [1.9+mrpivot*(-0.175),-0.175]
    
return, coeff
end

function mf_select_gmr_red, gmr, mr, z=z, gmr_cut=gmr_cut, $
  coeff=coeff, sdss_redpeak=sdss_redpeak, sdss_greenvalley=sdss_greenvalley, $
  mrpivot=mrpivot, blue=blue
; jm12jan05ucsd - select red and blue galaxies using the evolving
;   0.1(g-r) color cut (see fit_gmr_bimodality)

; get the coefficients of the bimodality minimum for every redshift
    nzz = n_elements(z)
    coeff = get_gmr_bimodality_coeff(z,sdss_redpeak=sdss_redpeak, $
      sdss_greenvalley=sdss_greenvalley,mrpivot=mrpivot)
    
    ngal = n_elements(gmr)
    if (ngal eq 0L) then return, -1 ; just return the coefficients
    if (n_elements(mr) ne ngal) then begin
       splog, 'Dimensions of GMR and MR do not agree!'
       return, -1
    endif

; evolving sloping color-cut of the form (u-z) = (u-z)_0 + a*(Mass+21) + b*(z-0.1)
    if (nzz eq 1) then gmr_cut = poly(mr-mrpivot,coeff) else begin
       gmr_cut = mr*0.0
       for ii = 0L, ngal-1 do gmr_cut[ii] = poly(mr[ii]-mrpivot,coeff[*,ii])
    endelse
    
    iquiescent = where(gmr gt gmr_cut,comp=iblue)
    if keyword_set(blue) then indx = iblue else $
      indx = iquiescent
    
return, indx
end
