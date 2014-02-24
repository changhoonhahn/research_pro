pro fit_mf_vmax, noevol=noevol, sdss=sdss, fit_individual=fit_individual, $
  fit_combined=fit_combined, final_smf=final_smf, literature=literature, $
  binsize=binsize, debug=debug, clobber=clobber
; jm11aug12ucsd - fit the MF using the 1/Vmax weighted points 
    
    mfpath = mf_path(/mfs)
    if keyword_set(noevol) then esuffix = 'noevol' else esuffix = 'evol'
    
    if keyword_set(sdss) then begin
       suffix = '_sdss'
    endif else begin
       suffix = ''
    endelse

    field = get_mf_fields(sdss=sdss)
    nfield = n_elements(field)
    win = get_mf_window(field,area=area,/sr)
    totarea = total(area,/double)
    
    zbins = mf_zbins(sdss=sdss,nzbins,literature=literature)
    if keyword_set(literature) then litsuffix = '_lit' else litsuffix = ''

; select quiescent and star-forming galaxies using various criteria
    subsample = ['all','quiescent','active','red','blue']
;   subsample = ['active']
;   subsample = ['all','quiescent','active','quiescent_nuvmr','active_nuvmr','red','blue']
;   subsample = ['quiescent','active']
    if keyword_set(sdss) then subsample = [subsample,'red_gmr','blue_gmr']
;   subsample = 'quiescent'

; read the supergrid parameter file
    super = get_mf_supergrid(nsuper=nsuper,superstring=superstring)
;   splog, 'Hack!!'
;   super = super[0] & superstring = superstring[0] & nsuper = 1

; binning parameters
    binsize = mf_binsize(bin=binsize,minmass=minmass,maxmass=maxmass)
    
; just do core subsamples and supergrid01 if /literature or if the
; binsize is >0.1 dex
    if keyword_set(literature) or (binsize gt 0.1) then begin
       super = get_mf_supergrid(1,nsuper=nsuper,superstring=superstring)
       subsample = ['all','quiescent','active']
    endif

; fit the primus MFs well above the formal mass limit to avoid surface
; brightness selection effects
;   if keyword_set(sdss) then minfitmass_offset = 0.0 else minfitmass_offset = 1.0
;   if keyword_set(sdss) then masslimit_offset = 0.0 else masslimit_offset = 1.0
    
; ---------------------------------------------------------------------------
; construct and fit the SDSS+GALEX and combined PRIMUS sample
    if (keyword_set(sdss) and keyword_set(final_smf) eq 0) or keyword_set(fit_combined) or $
      keyword_set(literature) then begin
       for jj = 0, nsuper-1 do begin
          for ss = 0, n_elements(subsample)-1 do begin
             datafile = mfpath+'mfdata_'+subsample[ss]+'_supergrid'+$
               superstring[jj]+suffix+litsuffix+'.fits.gz'
             data = mrdfits(datafile,1)

;            splog, 'Hack!'
;            for iz = 2, nzbins-1 do begin
             for iz = 0, nzbins-1 do begin
                title = 'Field Average, z='+string(zbins[iz].zlo,format='(F4.2)')+'-'+$
                  string(zbins[iz].zup,format='(F4.2)')+', '+repstr(subsample[ss],'_','-')+', '+$
                  esuffix

                these = where((data.z ge zbins[iz].zlo) and (data.z lt zbins[iz].zup),ngal)
                if keyword_set(noevol) then vmax = data[these].vmax_noevol else $
                  vmax = data[these].vmax_evol

; combine the MFs from the individual fields
                mfdata1 = im_mf_vmax(data[these].mass,data[these].weight/vmax,$
                  binsize=binsize,minmass=minmass,maxmass=maxmass,$
                  masslimit=data[these].masslimit,rev=rev)
;                 masslimit=(data[these].masslimit+masslimit_offset)<11)

; get the mean SFR in each stellar mass bin
                mfdata1 = struct_addtags(mfdata1,{sfrmean: fltarr(mfdata1.nbins), sfrsigma: fltarr(mfdata1.nbins)})
                for bb = 0, mfdata1.nbins-1 do begin
                   if (rev[bb] ne rev[bb+1]) then begin ; >0 galaxies
                      sfr = (data[these].sfr)[rev[rev[bb]:rev[bb+1]-1]]
                      mfdata1.sfrmean[bb] = djs_median(sfr)
                      mfdata1.sfrsigma[bb] = djsig(sfr)
                   endif
                endfor
                
; get the number of fields contributing to each stellar mass bin                 
                mfdata1 = struct_addtags(mfdata1,{nfield: intarr(mfdata1.nbins), $
                  thesefields: strarr(5,mfdata1.nbins)})
                mfdata1.thesefields = '...'
                if keyword_set(sdss) then begin
                   mfdata1.nfield = 1
                   mfdata1.thesefields[0] = 'sdss'
;                  fitdouble = 1
;                  if ss eq 1 then fitdouble = 1 else fitdouble = 0
                endif else begin
                   for bb = 0, mfdata1.nbins-1 do begin
                      if (rev[bb] ne rev[bb+1]) then begin ; >0 galaxies
                         allfield = (data[these].field)[rev[rev[bb]:rev[bb+1]-1]]
                         thesefields = allfield[uniq(allfield,sort(allfield))]
                         mfdata1.nfield[bb] = n_elements(thesefields)
                         mfdata1.thesefields[0:mfdata1.nfield[bb]-1,bb] = thesefields
                      endif
                   endfor
                endelse
;               niceprint, mfdata1.mass, mfdata1.nfield

                parinfo = init_mf_parinfo(active=subsample[ss] eq 'active',$
                  quiescent=subsample[ss] eq 'quiescent',nuvmr_active=subsample[ss] eq 'active_nuvmr',$
                  nuvmr_quiescent=subsample[ss] eq 'quiescent_nuvmr',$
                  red=subsample[ss] eq 'red',blue=subsample[ss] eq 'blue',$
                  gmr_red=subsample[ss] eq 'red_gmr',gmr_blue=subsample[ss] eq 'blue_gmr',$
                  double_parinfo=double_parinfo)

                delvarx, uniform_weights
                if keyword_set(sdss) then begin
; based on the results of qaplot_fit_mf_vmax, when fitting the
; quiescent SDSS sample I need to assume uniform weights 
                   if subsample[ss] eq 'quiescent' then uniform_weights = 1
                endif

                mffit1 = do_vmax_mffit(mfdata1,parinfo=parinfo,debug=debug,$
                  double_parinfo=double_parinfo,mffit_double=mffit_double1,$
                  uniform_weights=uniform_weights,title=title)
                rho1 = integrate_vmax_mf(mfdata1,mffit1,mffit_double=mffit_double1)
                
; use jackknife to get the cosmic variance uncertainty
                if keyword_set(sdss) then begin
                   hist = hist_nd(transpose([[data[these].ra],[data[these].dec]]),$
                     [30D,20D],rev=rev)
                   nbins = n_elements(hist)
                   notempty = where(hist gt 1000,nfield_jack)
                   splog, 'SDSS - '+strtrim(nbins,2)+' bins, '+$
                     strtrim(nfield_jack,2)+' jackknife fields, '+$
                     strtrim(mean(hist[notempty]),2)+' galaxies on average, '+$
                     strtrim(min(hist[notempty]),2)+' min, '+strtrim(max(hist[notempty]),2)+' max'
                endif else nfield_jack = nfield

; loop on each field                
                for ii = 0, nfield_jack-1 do begin
                   if keyword_set(sdss) then begin
                      kk = notempty[ii]
                      toss = these[rev[rev[kk]:rev[kk+1]-1]]
                      keep = cmset_op(these,'and',/not2,toss)

                      reweight = 1D ; assume to be unity
                      
                   endif else begin
                      keep = where(strtrim(data[these].field,2) ne field[ii])
; reweight Vmax because we've now dropped one of the fields!
                      uarea = (data[these[keep]].area)[uniq(data[these[keep]].field)]
                      reweight = totarea/total(uarea,/double)
                   endelse

                   jack_mfdata1 = im_mf_vmax(data[these[keep]].mass,$
                     reweight*data[these[keep]].weight/vmax[keep],$
                     binsize=binsize,minmass=minmass,maxmass=maxmass,$
                     masslimit=data[these[keep]].masslimit)

                   jack_mffit1 = do_vmax_mffit(jack_mfdata1,parinfo=parinfo,debug=0,$
                     double_parinfo=double_parinfo,mffit_double=jack_mffit_double1,$
                     uniform_weights=uniform_weights,/quiet)
                   jack_rho1 = integrate_vmax_mf(jack_mfdata1,jack_mffit1,mffit_double=jack_mffit_double1)
                   
                   if (ii eq 0) then begin
                      jack_mffit = jack_mffit1
                      jack_mfdata = jack_mfdata1
                      jack_rho = jack_rho1
                   endif else begin
                      jack_mffit = [jack_mffit,jack_mffit1]
                      jack_mfdata = [jack_mfdata,jack_mfdata1]
                      jack_rho = [jack_rho,jack_rho1]
                   endelse
                endfor

                jfact = (nfield_jack-1.0)/float(nfield_jack) ; jackknife prefactor
                mffit1.phistar_cv_err = sqrt(jfact*total((jack_mffit.phistar-djs_mean(jack_mffit.phistar))^2))
                mffit1.logmstar_cv_err = sqrt(jfact*total((jack_mffit.logmstar-djs_mean(jack_mffit.logmstar))^2))
                mffit1.alpha_cv_err = sqrt(jfact*total((jack_mffit.alpha-djs_mean(jack_mffit.alpha))^2))

; get the variance just based on the fields contributing to each mass bin
                good = where(mfdata1.nfield gt 0,ngood)
                for gg = 0, ngood-1 do begin
                   if keyword_set(sdss) then begin
                      nbb = nfield_jack
                      m2 = lindgen(nfield_jack)
                   endif else begin
                      match, strtrim(mfdata1.thesefields[*,good[gg]],2), field, m1, m2
                      nbb = n_elements(m1)
                   endelse
                   if (nbb ge 3) then begin
                      mfdata1.phierr_cv[good[gg]] = sqrt((nbb-1.0)/nbb*$
                        total((jack_mfdata[m2].phi[good[gg]]-$
                        djs_mean(jack_mfdata[m2].phi[good[gg]]))^2))
                   endif else begin
                   endelse
                endfor

; for mass bins with contributions from fewer then three fields use
; the nearest variance
                good = where(mfdata1.nfield gt 0 and mfdata1.phierr_cv gt 0.0)
                fix = where(mfdata1.nfield gt 0 and mfdata1.phierr_cv le 0.0,nfix)
                if (nfix ne 0) then begin
                   get_element, good, fix, nearest
                   mfdata1.phierr_cv[fix] = mfdata1.phierr_cv[good[nearest]]
                endif
;               niceprint, mfdata1.mass, mfdata1.phi, mfdata1.phierr_cv, mfdata1.nfield, mfdata1.number

; now deal with the number and mass density
                jfact = (nfield_jack-1.0)/float(nfield_jack) ; jackknife prefactor
                rho1.masserr_cv_cumuphi50 = sqrt(jfact*total((jack_rho.mass_cumuphi50-djs_mean(jack_rho.mass_cumuphi50))^2))
                rho1.masserr_cv_cumuphi45 = sqrt(jfact*total((jack_rho.mass_cumuphi45-djs_mean(jack_rho.mass_cumuphi45))^2))
                rho1.masserr_cv_cumuphi40 = sqrt(jfact*total((jack_rho.mass_cumuphi40-djs_mean(jack_rho.mass_cumuphi40))^2))
                rho1.masserr_cv_cumuphi35 = sqrt(jfact*total((jack_rho.mass_cumuphi35-djs_mean(jack_rho.mass_cumuphi35))^2))
                rho1.masserr_cv_cumuphi30 = sqrt(jfact*total((jack_rho.mass_cumuphi30-djs_mean(jack_rho.mass_cumuphi30))^2))
                rho1.masserr_cv_cumuphi25 = sqrt(jfact*total((jack_rho.mass_cumuphi25-djs_mean(jack_rho.mass_cumuphi25))^2))
;               rho1.masserr_cv_cumuphi20 = sqrt(jfact*total((jack_rho.mass_cumuphi20-djs_mean(jack_rho.mass_cumuphi20))^2))

                pprefix = ['rho','num']
                ssuffix = ['95','10','105','11','95_10','10_105','105_11','11_115']
                for pp = 0, 1 do begin
                   for mm = 0, n_elements(ssuffix)-1 do begin
                      intag = tag_indx(jack_rho,pprefix[pp]+'_'+ssuffix[mm])
                      outtag = tag_indx(rho1,pprefix[pp]+'err_cv_'+ssuffix[mm])

                      gd = where(jack_rho.(intag) gt -900.0,ngd)
                      if (ngd ne 0) then begin
;                           if ngd ne nfield then begin
;                              djs_plot, jack_mfdata[0].mass, jack_mfdata[0].phi, $
;                                psym=-6, xr=[9.8,11.8], /ylog, yr=[1E-5,1E-2], xsty=3, ysty=3
;                              for mm=1,4 do djs_oplot, jack_mfdata[mm].mass, jack_mfdata[mm].phi, psym=-6, color='blue'
;                           endif
                         jfact = (ngd-1.0)/float(ngd) ; jackknife prefactor
                         rho1.(outtag) = sqrt(jfact*total((jack_rho[gd].(intag)-djs_mean(jack_rho[gd].(intag)))^2))
                      endif
                   endfor
                endfor 

; pack it all in                
                if (iz eq 0) then begin
                   mfdata = mfdata1
                   rho = rho1
                   mffit = mffit1
                   mffit_double = mffit_double1
                endif else begin
                   mfdata = [mfdata,mfdata1]
                   rho = [rho,rho1]
                   mffit = [mffit,mffit1]
                   mffit_double = [mffit_double,mffit_double1]
                endelse
             endfor 
             outfile = mfpath+'mf_'+subsample[ss]+'_supergrid'+superstring[jj]+$
               '_bin'+string(binsize,format='(F4.2)')+'_'+esuffix+suffix+litsuffix+'.fits'

             im_mwrfits, mfdata, outfile, clobber=clobber, /nogzip
             im_mwrfits, rho, outfile, /append, /nogzip
;            if keyword_set(fitdouble) then begin
                im_mwrfits, mffit, outfile, /append, /nogzip
                im_mwrfits, mffit_double, outfile, /append, /gzip
;            endif else im_mwrfits, mffit, outfile, /append, /gzip
          endfor 
       endfor 
    endif
    
; ---------------------------------------------------------------------------
; construct the MF in the individual PRIMUS fields - just for all
; galaxies and supergrid01
    if (keyword_set(sdss) eq 0) and keyword_set(fit_individual) then begin
       ind_binsize = 0.15
       ind_minmass = 8.5-ind_binsize/2                  ; shift by a half-bin
       ind_subsample = ['all','quiescent','active'] ; only do these three samples
       ind_super = get_mf_supergrid(1,nsuper=ind_nsuper,superstring=ind_superstring)
       for jj = 0, ind_nsuper-1 do begin
          for ss = 0, n_elements(ind_subsample)-1 do begin
             parinfo = init_mf_parinfo(active=subsample[ss] eq 'active',$
               quiescent=subsample[ss] eq 'quiescent',double_parinfo=double_parinfo)

             for ii = 0, nfield-1 do begin
                splog, field[ii]
                datafile = mfpath+'mfdata_'+ind_subsample[ss]+'_supergrid'+$
                  ind_superstring[jj]+'_'+field[ii]+litsuffix+'.fits.gz'
                data = mrdfits(datafile,1)

                for iz = 0, nzbins-1 do begin
                   title = field[ii]+', z='+string(zbins[iz].zlo,format='(F4.2)')+'-'+$
                     string(zbins[iz].zup,format='(F4.2)')+', '+repstr(ind_subsample[ss],'_','-')+', '+$
                     esuffix

                   these = where((data.z ge zbins[iz].zlo) and (data.z lt zbins[iz].zup),ngal)
                   if keyword_set(noevol) then vmax = data[these].vmax_noevol else $
                     vmax = data[these].vmax_evol

                   mfdata1 = im_mf_vmax(data[these].mass,data[these].weight/vmax,$
                     binsize=ind_binsize,minmass=ind_minmass,maxmass=maxmass,noutbins=noutbins,$
                     masslimit=data[these].masslimit)
                   mffit1 = do_vmax_mffit(mfdata1,parinfo=parinfo,double_parinfo=double_parinfo,$
                     mffit_double=mffit_double1,debug=debug,title=title)

                   if (iz eq 0) then begin
                      mffit = mffit1
                      mfdata = mfdata1
                   endif else begin
                      mffit = [mffit,mffit1]
                      mfdata = [mfdata,mfdata1]
                   endelse
                endfor

                outfile = mfpath+'mf_'+ind_subsample[ss]+'_supergrid'+ind_superstring[jj]+$
                  '_bin'+string(binsize,format='(F4.2)')+'_'+esuffix+$
                  '_'+field[ii]+litsuffix+'.fits'
                im_mwrfits, mfdata, outfile, clobber=clobber, /nogzip
                im_mwrfits, mffit, outfile, /append, /gzip
             endfor 
          endfor
       endfor
    endif 
return
end
