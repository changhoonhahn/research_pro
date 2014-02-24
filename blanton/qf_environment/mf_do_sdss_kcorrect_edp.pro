function mf_do_sdss_kcorrect_edp, in_redshift, in_maggies, in_ivarmaggies, $
  filterlist=filterlist, vname=vname, just_ugrizkhk_01=just_ugrizjhk_01, $
  targ_filter=targ_filter, maxiter=maxiter, sdss=sdss
; jm10aug30ucsd - simple wrapper script to compute a general set of
; K-corrections for the PRIMUS/MF project (called by MF_ISEDFIT)
; jm11mar15ucsd - added 2MASS absolute magnitudes

    ngal = n_elements(in_redshift)
    nfilt = n_elements(filterlist)

    h100 = mf_h100()
    vname = mf_vname()
    print, filterlist
    weff = k_lambda_eff(filterlist=filterlist)
    
    kcorr = {$
      k_zobj:                          -999.0,$
      k_maggies:                fltarr(nfilt),$
      k_ivarmaggies:            fltarr(nfilt),$
      k_maggies_mask:         intarr(nfilt)+1,$ ; good=1, bad=0
      k_bestmaggies:            fltarr(nfilt),$
      k_mass:                          -999.0,$
      k_coeffs:                     fltarr(5),$
      k_chi2:                          -999.0,$

      k_fnuvugrizjhk_absmag_01:         fltarr(10)-999.0,$
      k_fnuvugrizjhk_absmag_ivar_01:    fltarr(10)-999.0,$
      k_fnuvugrizjhk_kcorrect_01:       fltarr(10)-999.0,$

      k_synth_fnuvugrizjhk_absmag_01:    fltarr(10)-999.0}
    kcorr = replicate(kcorr,ngal)

    kcorr.k_zobj = in_redshift
    kcorr.k_maggies = in_maggies
    kcorr.k_ivarmaggies = in_ivarmaggies

    out_filterlist1 = [galex_filterlist(),sdss_filterlist(),twomass_filterlist()] 
; compute k-corrections iteratively, rejecting photometric datapoints
; that contribute more than MAXCONTRIB to chi^2 and that are more than
; NSIGMA from the best-fit model; require more than NMINFILT *optical*
; datapoints (otherwise if one band gets rejected then chi^2 can
; become zero)
    nminfilt = 4
    nsigma = 3.0
    maxcontrib = 0.5

    optfilt = where((weff gt 3000.0) and (weff lt 1E4),noptfilt)

    fnuvugrizjhk_kcorrect_01 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies*$
        kcorr.k_maggies_mask,filterlist,out_filterlist1,band_shift=0.1,chi2=kchi2,$
        mass=mass,coeffs=coeffs,bestmaggies=bestmaggies,absmag=fnuvugrizjhk_absmag_01,$
        ivarabsmag=fnuvugrizjhk_absmag_ivar_01,synth_absmag=synth_fnuvugrizjhk_absmag_01,$
        clineflux=cflux,uvflux=uvflux,vname=vname,h100=h100)       ; AB, band_shift=0.1

    chi2 = total(kcorr.k_maggies_mask*in_ivarmaggies*(in_maggies-bestmaggies)^2,1,/double)
    dof = total(kcorr.k_maggies_mask*in_ivarmaggies gt 0,1)-1.0 ; degrees of freedom
    neg = where(dof le 0.0)
    if (neg[0] ne -1) then message, 'This should not happen!'
    chi2 = chi2/dof
   
; reject deviant bandpasses       
    these = where(total(kcorr.k_maggies_mask[optfilt,*]*in_ivarmaggies[optfilt,*] gt 0,1) ge nminfilt,nthese)
    if (nthese ne 0L) then begin
        optmaggies = in_maggies[optfilt,*]
        optbestmaggies = bestmaggies[optfilt,*]
        optivarmaggies = in_ivarmaggies[optfilt,*]

        optchi2 = total(kcorr[these].k_maggies_mask[optfilt]*optivarmaggies[*,these]*$
            (optmaggies[*,these]-optbestmaggies[*,these])^2,1,/double)
        chi2_perband = kcorr[these].k_maggies_mask[optfilt]*optivarmaggies[*,these]*$
            (optmaggies[*,these]-optbestmaggies[*,these])^2/$
            rebin(reform(optchi2,1,nthese),noptfilt,nthese)

        mag = maggies2mag(optmaggies[*,these],ivarmaggies=optivarmaggies[*,these],magerr=magerr)
        bestmag = maggies2mag(optbestmaggies[*,these])

        kcorr[these].k_maggies_mask[optfilt] = ((abs(mag-bestmag) gt nsigma*magerr) and $
            (abs(mag-bestmag) gt 0.2) and (chi2_perband gt maxcontrib) and $
            (optivarmaggies[*,these] gt 0)) eq 0           ; 1=good, 0=reject
    endif

    kcorr.k_bestmaggies = bestmaggies
    kcorr.k_mass = alog10(mass) ; h=0.7, Chabrier
    kcorr.k_coeffs = coeffs
    kcorr.k_chi2 = chi2
    kcorr.k_fnuvugrizjhk_absmag_01      = fnuvugrizjhk_absmag_01
    kcorr.k_fnuvugrizjhk_absmag_ivar_01 = fnuvugrizjhk_absmag_ivar_01
    kcorr.k_fnuvugrizjhk_kcorrect_01    = fnuvugrizjhk_kcorrect_01
    kcorr.k_synth_fnuvugrizjhk_absmag_01 = synth_fnuvugrizjhk_absmag_01

return, kcorr
end

