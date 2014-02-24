function read_mf_vmax, field, supergrid=supergrid, binsize=binsize, avgsupergrid=avgsupergrid, $
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
    if keyword_set(nuvmr_quiescent) then subsample = 'quiescent_nuvmr'
    if keyword_set(nuvmr_active) then subsample = 'active_nuvmr'
    if keyword_set(red) then subsample = 'red'
    if keyword_set(blue) then subsample = 'blue'
    if keyword_set(gmr_red) then subsample = 'red_gmr'
    if keyword_set(gmr_blue) then subsample = 'blue_gmr'

    if keyword_set(noevol) then esuffix = 'noevol' else esuffix = 'evol'
    if keyword_set(evolved) then prefix = 'evolved_' else prefix = ''

    if keyword_set(avgsupergrid) then super = 'supergridavg' else $
      if (n_elements(supergrid) eq 0) then supergrid = 1
    super = 'supergrid'+string(supergrid,format='(I2.2)')
    if keyword_set(final_smf) then super = 'final'

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
    if (keyword_set(silent) eq 0) then splog, 'Reading '+mffile    

    ext = 1
    if keyword_set(rho) then ext = 2
    if keyword_set(bestfit) then ext = 3
    if keyword_set(double_bestfit) then ext = 4
    print, mffile
;    mf = mrdfits(mffile,ext,silent=silent)
    nzz = n_elements(mf)

;; take the log of 'rho' and 'num' total
;    if keyword_set(rho) then begin
;       if tag_exist(mf,'rho') then begin
;          gd = where(mf.rho gt 0,ngd)
;          if (ngd ne 0) then begin
;             mf.rhoerr = mf.rhoerr/mf.rho/alog(10)
;             mf.rho = alog10(mf.rho)
;          endif
;          gd = where(mf.num gt 0,ngd)
;          if (ngd ne 0) then begin
;             mf.numerr = mf.numerr/mf.num/alog(10)
;             mf.num = alog10(mf.num)
;          endif
;       endif
;    endif 

; add the statistical and CV errors on RHO, NUM, and the CUMUPHI
; quantities in quadrature 
    if keyword_set(rho) then begin

; build the needed tags on the-fly
       outrho = mf
       mpprefix = ['rho','num']
       pprefix = ['rhoerr','numerr']
       ssuffix = ['95','10','105','11','95_10','10_105','105_11','11_115']
       for pp = 0, 1 do begin
          for mm = 0, n_elements(ssuffix)-1 do begin
             if tag_exist(outrho,pprefix[pp]+'_stat_'+ssuffix[mm]) eq 0 then $
               outrho = struct_addtags(outrho,replicate(create_struct(pprefix[pp]+'_stat_'+ssuffix[mm],-999.0),nzz))

             etag = tag_indx(outrho,pprefix[pp]+'_'+ssuffix[mm])
             cvtag = tag_indx(outrho,pprefix[pp]+'_cv_'+ssuffix[mm])
             outtag = tag_indx(outrho,pprefix[pp]+'_stat_'+ssuffix[mm])

             outrho.(outtag) = sqrt(outrho.(etag)^2+outrho.(cvtag)^2)

             if keyword_set(nolog) then begin
                mtag = tag_indx(outrho,mpprefix[pp]+'_'+ssuffix[mm])
                outrho.(etag) = alog(10)*outrho.(etag)*10.0^outrho.(mtag)
                outrho.(cvtag) = alog(10)*outrho.(cvtag)*10.0^outrho.(mtag)
                outrho.(outtag) = alog(10)*outrho.(outtag)*10.0^outrho.(mtag)
                outrho.(mtag) = 10.0^outrho.(mtag)
             endif
          endfor
       endfor

; now do cumuphi
       ssuffix = 'cumuphi'+['50','45','40','35','30','25']
       for mm = 0, n_elements(ssuffix)-1 do begin
          if tag_exist(outrho,'masserr_stat_'+ssuffix[mm]) eq 0 then $
            outrho = struct_addtags(outrho,replicate(create_struct('masserr_stat_'+ssuffix[mm],-999.0),nzz))

          etag = tag_indx(outrho,'masserr_'+ssuffix[mm])
          cvtag = tag_indx(outrho,'masserr_cv_'+ssuffix[mm])
          outtag = tag_indx(outrho,'masserr_stat_'+ssuffix[mm])

          outrho.(outtag) = sqrt(outrho.(etag)^2+outrho.(cvtag)^2)
       endfor 
       mf = outrho
    endif

; take the log of Phi; could modify to deal with lower/upper limits 
    if keyword_set(log) and (keyword_set(rho) eq 0) and (keyword_set(bestfit) eq 0) then begin

; add some tags on the-fly
;      if keyword_set(final_smf) then begin
       nbins = mf[0].nbins
       mf = struct_addtags(mf,replicate({$
         phi_lower_stat: fltarr(nbins), phi_upper_stat: fltarr(nbins),$
         phierr_lower_stat: fltarr(nbins), phierr_upper_stat: fltarr(nbins)},nzz))
;        phierr_lower_syst: fltarr(nbins), phierr_upper_syst: fltarr(nbins),$
;        phierr_lower_stat_syst: fltarr(nbins), phierr_upper_stat_syst: fltarr(nbins),$
;        phi_lower_stat_syst: fltarr(nbins), phi_upper_stat_syst: fltarr(nbins)},nzz))
;      endif

       for iz = 0, nzz-1 do begin
          gd = where(mf[iz].number gt 0,ngd)
          if (ngd ne 0) then begin

; add the CV and Poisson errors in quadrature             
             mf[iz].phierr_lower_stat = sqrt(mf[iz].phierr_lower^2 + mf[iz].phierr_cv^2)
             mf[iz].phierr_upper_stat = sqrt(mf[iz].phierr_upper^2 + mf[iz].phierr_cv^2)

; only for the FINAL structure             
             if tag_exist(mf,'phierr_lower_stat') then mf[iz].phierr_lower_stat[gd] = mf[iz].phierr_lower_stat[gd]/mf[iz].phi[gd]/alog(10)
             if tag_exist(mf,'phierr_upper_stat') then mf[iz].phierr_upper_stat[gd] = mf[iz].phierr_upper_stat[gd]/mf[iz].phi[gd]/alog(10)
;            if tag_exist(mf,'phierr_lower_syst') then mf[iz].phierr_lower_syst[gd] = mf[iz].phierr_lower_syst[gd]/mf[iz].phi[gd]/alog(10)
;            if tag_exist(mf,'phierr_upper_syst') then mf[iz].phierr_upper_syst[gd] = mf[iz].phierr_upper_syst[gd]/mf[iz].phi[gd]/alog(10)
;            if tag_exist(mf,'phierr_lower_priors') then mf[iz].phierr_lower_priors[gd] = mf[iz].phierr_lower_priors[gd]/mf[iz].phi[gd]/alog(10)
;            if tag_exist(mf,'phierr_upper_priors') then mf[iz].phierr_upper_priors[gd] = mf[iz].phierr_upper_priors[gd]/mf[iz].phi[gd]/alog(10)
;            if tag_exist(mf,'phierr_lower_models') then mf[iz].phierr_lower_models[gd] = mf[iz].phierr_lower_models[gd]/mf[iz].phi[gd]/alog(10)
;            if tag_exist(mf,'phierr_upper_models') then mf[iz].phierr_upper_models[gd] = mf[iz].phierr_upper_models[gd]/mf[iz].phi[gd]/alog(10)
;            if tag_exist(mf,'phierr_lower_stat_syst') then mf[iz].phierr_lower_stat_syst[gd] = mf[iz].phierr_lower_stat_syst[gd]/mf[iz].phi[gd]/alog(10)
;            if tag_exist(mf,'phierr_upper_stat_syst') then mf[iz].phierr_upper_stat_syst[gd] = mf[iz].phierr_upper_stat_syst[gd]/mf[iz].phi[gd]/alog(10)

;            if tag_exist(mf,'phi_lower_priors') then mf[iz].phi_lower_priors[gd] = alog10(mf[iz].phi_lower_priors[gd])
;            if tag_exist(mf,'phi_upper_priors') then mf[iz].phi_upper_priors[gd] = alog10(mf[iz].phi_upper_priors[gd])
;            if tag_exist(mf,'phi_lower_models') then mf[iz].phi_lower_models[gd] = alog10(mf[iz].phi_lower_models[gd])
;            if tag_exist(mf,'phi_upper_models') then mf[iz].phi_upper_models[gd] = alog10(mf[iz].phi_upper_models[gd])
;            if tag_exist(mf,'phi_lower_syst') then mf[iz].phi_lower_syst[gd] = alog10(mf[iz].phi_lower_syst[gd])
;            if tag_exist(mf,'phi_upper_syst') then mf[iz].phi_upper_syst[gd] = alog10(mf[iz].phi_upper_syst[gd])

;            if tag_exist(mf,'phi_lower_stat') then mf[iz].phi_lower_stat[gd] = alog10(mf[iz].phi_lower_stat[gd])
;            if tag_exist(mf,'phi_upper_stat') then mf[iz].phi_upper_stat[gd] = alog10(mf[iz].phi_upper_stat[gd])
;            if tag_exist(mf,'phi_lower_stat_syst') then mf[iz].phi_lower_stat_syst[gd] = alog10(mf[iz].phi_lower_stat_syst[gd])
;            if tag_exist(mf,'phi_upper_stat_syst') then mf[iz].phi_upper_stat_syst[gd] = alog10(mf[iz].phi_upper_stat_syst[gd])

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
