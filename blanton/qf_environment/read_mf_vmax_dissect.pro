function read_mf_vmax_dissect, field, supergrid=supergrid, binsize=binsize, avgsupergrid=avgsupergrid, $
  avgfield=avgfield, quiescent=quiescent, active=active, nuvmr_quiescent=nuvmr_quiescent, $
  nuvmr_active=nuvmr_active, red=red, blue=blue, $
  gmr_red=gmr_red, gmr_blue=gmr_blue, rho=rho, bestfit=bestfit, $
  double_bestfit=double_bestfit, noevol=noevol, silent=silent, log=log, $
  literature=literature, evolved=evolved, nolog=nolog, mffile=mffile;, final_smf=final_smf
; jm10sep16ucsd - read the output of FIT_MF_VMAX, for an arbitrary
;   combination of keywords; the default is to return the binned MF
;   (stored in the first extension) unless /BESTFIT, which is stored
;   in the second FITS extension
; jm11mar31ucsd - major update

    mfpath = mf_path(/mfs)
    subsample = 'all' ; default
    if keyword_set(quiescent) then subsample = 'quiescent'
    if keyword_set(active) then subsample = 'active'
    
    if keyword_set(noevol) then esuffix = 'noevol' else esuffix = 'evol'

    if keyword_set(avgsupergrid) then super = 'supergridavg' else $
      if (n_elements(supergrid) eq 0) then supergrid = 1
    super = 'supergrid'+string(supergrid,format='(I2.2)')

    if keyword_set(avgfield) then fld = '' else begin
       if (n_elements(field) eq 0) then begin
          splog, 'FIELD input required'
          return, -1
       endif
       fld = '_'+field
    endelse
    
    if keyword_set(literature) then litsuffix = '_lit' else litsuffix = ''

    binsize = mf_binsize(bin=binsize)
;   if keyword_set(evolved) then binsize = 0.3

    if n_elements(mffile) eq 0 then begin
       mffile = mfpath+prefix+'mf_'+subsample+'_'+super+'_bin'+string(binsize,format='(F4.2)')+$
         '_'+esuffix+fld+litsuffix+'.fits.gz'
    endif

    ext = 1
    print, mffile
    mf = mrdfits(mffile,ext,silent=silent)
    nzz = n_elements(mf)

; take the log of Phi; could modify to deal with lower/upper limits 
    if keyword_set(log) and (keyword_set(rho) eq 0) and (keyword_set(bestfit) eq 0) then begin

       nbins = mf[0].nbins
       mf = struct_addtags(mf,replicate({$
         phi_lower_stat: fltarr(nbins), phi_upper_stat: fltarr(nbins),$
         phierr_lower_stat: fltarr(nbins), phierr_upper_stat: fltarr(nbins)},nzz))

       for iz = 0, nzz-1 do begin
          gd = where(mf[iz].number gt 0,ngd)
          if (ngd ne 0) then begin

; add the CV and Poisson errors in quadrature             
             mf[iz].phierr_lower_stat = sqrt(mf[iz].phierr_lower^2 + mf[iz].phierr_cv^2)
             mf[iz].phierr_upper_stat = sqrt(mf[iz].phierr_upper^2 + mf[iz].phierr_cv^2)

; only for the FINAL structure             
             if tag_exist(mf,'phierr_lower_stat') then mf[iz].phierr_lower_stat[gd] = mf[iz].phierr_lower_stat[gd]/mf[iz].phi[gd]/alog(10)
             if tag_exist(mf,'phierr_upper_stat') then mf[iz].phierr_upper_stat[gd] = mf[iz].phierr_upper_stat[gd]/mf[iz].phi[gd]/alog(10)
             
             if tag_exist(mf,'phierr') then mf[iz].phierr[gd] = mf[iz].phierr[gd]/mf[iz].phi[gd]/alog(10)
             if tag_exist(mf,'phierr_poisson') then mf[iz].phierr_poisson[gd] = mf[iz].phierr_poisson[gd]/mf[iz].phi[gd]/alog(10)

;            splog, 'need to finish coding this!!!!!'
             
; standard set of tags             
             mf[iz].phierr_cv[gd] = mf[iz].phierr_cv[gd]/mf[iz].phi[gd]/alog(10)
             mf[iz].phierr_lower[gd] = mf[iz].phierr_lower[gd]/mf[iz].phi[gd]/alog(10)
             mf[iz].phierr_upper[gd] = mf[iz].phierr_upper[gd]/mf[iz].phi[gd]/alog(10)
             mf[iz].phi[gd] = alog10(mf[iz].phi[gd])

; additional optional tags
;            if keyword_set(final_smf) then begin ; add these only in the log!
;               mf[iz].phierr_lower_syst = mf[iz].phi-mf[iz].phi_lower_syst
;               mf[iz].phierr_upper_syst = mf[iz].phi_upper_syst-mf[iz].phi
                
             mf[iz].phi_lower_stat = mf[iz].phi-mf[iz].phierr_lower_stat
             mf[iz].phi_upper_stat = mf[iz].phi+mf[iz].phierr_upper_stat

;               mf[iz].phierr_lower_stat_syst = mf[iz].phierr_lower_stat + mf[iz].phierr_lower_syst
;               mf[iz].phierr_upper_stat_syst = mf[iz].phierr_upper_stat + mf[iz].phierr_upper_syst

;               mf[iz].phi_lower_stat_syst = mf[iz].phi-mf[iz].phierr_lower_stat_syst
;               mf[iz].phi_upper_stat_syst = mf[iz].phi+mf[iz].phierr_upper_stat_syst
;            endif

;            mf[iz].phi[gd] = 10^mf[iz].phi[gd]
;            mf[iz].phierr_cv[gd] = alog(10)*mf[iz].phi[gd]*mf[iz].phierr_cv[gd]
;            mf[iz].phierr_lower[gd] = alog(10)*mf[iz].phi[gd]*mf[iz].phierr_lower[gd]
;            mf[iz].phierr_upper[gd] = alog(10)*mf[iz].phi[gd]*mf[iz].phierr_upper[gd]
;            if tag_exist(mf,'phierr') then mf[iz].phierr[gd] = alog(10)*mf[iz].phi[gd]*mf[iz].phierr[gd]
;            if tag_exist(mf,'phierr_poisson') then mf[iz].phierr_poisson[gd] = alog(10)*mf[iz].phi[gd]*mf[iz].phierr_poisson[gd]
          endif 
       endfor 
    endif 

return, mf
end
