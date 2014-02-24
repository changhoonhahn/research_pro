function get_nuvmr_bimodality_coeff, z, masspivot=masspivot, sdss_redpeak=sdss_redpeak, $
  sdss_greenvalley=sdss_greenvalley
; read the output from mf_nuvmr_bimodality

    mfpath = mf_path()
    final = mrdfits(mfpath+'nuvmr_evolution.fits.gz',1,/silent)

; get the intercept of the minimum in the bimodality at redshift Z
    ngal = n_elements(z)
    masspivot = final.masspivot
    nuvmrmin = poly(z-final.zpivot,final.nuvmrmin_coeff)
    coeff = transpose([[nuvmrmin],[replicate(final.slope,ngal)]])

; for comparison return my measurement of the SDSS green valley and
; the Wyder+07 coefficients 
    sdss = mrdfits(mfpath+'sdss_nuvmr_redpeak.fits.gz',1,/silent)
    sdss_greenvalley = sdss.greenvalley_coeff
    sdss_redpeak = sdss.redpeak_coeff
;   wyder_redpeak = [1.9+masspivot*(-0.175),-0.175]
    
return, coeff
end

function mf_select_nuvmr_quiescent, nuvmr, mass, z=z, nuvmr_cut=nuvmr_cut, $
  coeff=coeff, sdss_redpeak=sdss_redpeak, sdss_greenvalley=sdss_greenvalley, $
  masspivot=masspivot, active=active
; jm11sep21ucsd - select quiescent and actively star-forming galaxies
;   using the evolving 0.1(NUV-r) color cut (see fit_nuvmr_bimodality)

; get the coefficients of the bimodality minimum for every redshift
    nzz = n_elements(z)
    coeff = get_nuvmr_bimodality_coeff(z,sdss_redpeak=sdss_redpeak, $
      sdss_greenvalley=sdss_greenvalley,masspivot=masspivot)
    
    ngal = n_elements(nuvmr)
    if (ngal eq 0L) then return, -1 ; just return the coefficients
    if (n_elements(mass) ne ngal) then begin
       splog, 'Dimensions of NUVMR and MASS do not agree!'
       return, -1
    endif

; evolving sloping color-cut of the form (u-z) = (u-z)_0 + a*(Mass+21) + b*(z-0.1)
    if (nzz eq 1) then nuvmr_cut = poly(mass-masspivot,coeff) else begin
       nuvmr_cut = mass*0.0
       for ii = 0L, ngal-1 do nuvmr_cut[ii] = poly(mass[ii]-masspivot,coeff[*,ii])
    endelse
    
    iquiescent = where(nuvmr gt nuvmr_cut,comp=iactive)
    if keyword_set(active) then indx = iactive else $
      indx = iquiescent
    
;; three-sided selection box a la Williams+09
;    pars = {x1: -0.5, x2: 0.9, x3: 1.7, y1: 4.5, y2: 6.0, y3: 8.0}
;    slope = (pars.y2-pars.y1)/(pars.x3-pars.x2)
;    int = pars.y1-slope*pars.x2
;
;    if (n_elements(nuvmr) gt 0) then begin
;       iquiescent = where(((zmk lt pars.x2) and (nuvmr gt pars.y1)) or $
;         ((zmk gt pars.x2) and (zmk lt pars.x3) and $
;         (nuvmr gt poly(zmk,[int,slope]))),nquiescent,comp=iactive)
;       if keyword_set(active) then indx = iactive else $
;         indx = iquiescent
;    endif else indx = -1   
    
return, indx
end
