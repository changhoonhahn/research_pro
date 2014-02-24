function get_mf_sdss_photometry_edp, field, ivarmaggies=ivarmaggies, $
  zobj=zobj, zminmax=zminmax, filters=filters, zerod=zerod, $
  photo=photo, xray=xray, parent=parent, nozptoffsets=nozptoffsets, $
  use_kcorrect_mask=use_kcorrect_mask
; pull out the photometry we need

    zerod = mrdfits('/global/data/scr/chh327/primus/science/mf/2165/ubersample/edp_ubersample/'+$
        'sdss_ubersample.fits.gz',1)
    zerodver=''
    zobj = zerod.z
    ngal = n_elements(zobj)
; use modelmaggies for colors and scale by the cmodelmag in the r-band 
    sdss_to_maggies, modelmaggies, modelivarmaggies, calib=zerod, flux='model'
    sdss_to_maggies, cmodelmaggies, cmodelivarmaggies, calib=zerod, flux='cmodel'
    neg = where(modelmaggies[2,*] le 0)
    if (neg[0] ne -1L) then message, 'Bad!'
    factor = rebin(cmodelmaggies[2,*]/modelmaggies[2,*],5,n_elements(zerod))
    smaggies = modelmaggies*factor
    sivarmaggies = modelivarmaggies/factor^2

; get the GALEX photometry
    im_galex_to_maggies, zerod, gmaggies, givarmaggies, /allow_nuv_nondetect
;    gmaggies = fltarr(2,ngal)
;    givarmaggies = fltarr(2,ngal)

; get the WISE photometry (just channels 1 & 2)
    wise_to_maggies, zerod, wmaggies, wivarmaggies, /mpro
;    wmaggies = fltarr(4,ngal)
;    wivarmaggies = fltarr(4,ngal)

; get the 2MASS photometry; use 10-sigma upper limits on 2MASS fluxes
; where needed; see http://www.astro.virginia.edu/~mfs4n/2mass/allsky
; and Jarrett+10
    twomass_to_maggies, zerod, tmaggies, tivarmaggies
;    tmaggies = fltarr(2,ngal)
;    tivarmaggies = fltarr(2,ngal)

    maggies = [gmaggies,smaggies,tmaggies,wmaggies[0:1,*]]
    ivarmaggies = [givarmaggies,sivarmaggies,tivarmaggies,wivarmaggies[0:1,*]]
    allfilters = [galex_filterlist(),sdss_filterlist(),$
        twomass_filterlist(),(wise_filterlist())[0:1]] 
    
    filters = get_mf_filters(field)
    nfilt = n_elements(filters)
    
    match, strtrim(filters,2), allfilters, m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    maggies = maggies[m2,*]
    ivarmaggies = ivarmaggies[m2,*]

; apply the hand-masked photometry
    maskfile = '/home/users/hahn/primus/pro/science/mf/photometry_mask_'+$
       zerodver+'_'+field+'.par'
    if file_test(maskfile) then begin
        mask = yanny_readone(maskfile)
        ivarmaggies[*,mask.uber_position] = ivarmaggies[*,mask.uber_position]*mask.reject
    endif

; apply the K-correct photometric mask?
    if keyword_set(use_kcorrect_mask) then begin
        kk = read_mf_kcorrect(field,zerodver=kzerodver)
        if (n_elements(kk) ne n_elements(zobj)) or $
            (n_elements(filters) ne n_elements(kk[0].k_maggies)) or $
            (kzerodver ne zerodver) then message, 'Mismatch!!'
        ivarmaggies = ivarmaggies*kk.k_maggies_mask
    endif
return, maggies
end    
