function get_mf_filters_chh, field, nice_filters=nice_filters, noirac=noirac, $
  nogalex=nogalex, allcosmos=allcosmos
; jm10aug28ucsd - return the filters to be used, given the field;
;   FIT_FILTERINDEX specifies the filters that should be used for fitting
;   and K-corrections (see MF_ISEDFIT) 
;   splog, 'Testing which filters to use for fitting!!!'

    case field of
       'sdss': begin
          filters = [$
            primus_galex_filterlist(),$
            primus_sdss_filterlist(),$
            twomass_filterlist(),$
            (wise_filterlist())[0:1]]           
       end
       'cdfs': begin
          filters = [$
            primus_galex_filterlist(),$
            primus_cdfs_swire_filterlist_ugriz(),$
            (primus_cdfs_swire_filterlist_irac())[0:1]] ; just ch1-2
       end
       'cosmos': begin
          filters = [$
            primus_galex_filterlist(),$
            primus_cosmos_filterlist(),$
            (primus_scosmos_irac_filterlist())[0:1]]
; exclude the crummy cosmos bands (see also mf_derive_zptoffsets)
          if (keyword_set(allcosmos) eq 0) then begin
             filters = filters[where((strmatch(filters,'*ukirt*') eq 0) and $
               (strmatch(filters,'*wircam*',/fold) eq 0) and $
               (strmatch(filters,'*suprimecam_V*',/fold) eq 0) and $
               (strmatch(filters,'*suprimecam_B*') eq 0))]
          endif
       end
       'deep2_02hr': begin
          filters = [$
            primus_galex_filterlist(),$
            primus_deep2_filterlist(),$
            primus_deep2_iz_filterlist(),$
            primus_deep2_ugi_filterlist()]
       end
       'deep2_23hr': begin
          filters = [$
            primus_galex_filterlist(),$
            primus_deep2_filterlist(),$
            primus_deep2_iz_filterlist()]
       end
       'dls': begin
          filters = primus_dls_filterlist() ; BVRz
       end
       'es1': begin
          filters = [$
            primus_galex_filterlist(),$
            primus_es1_swire_filterlist_bvriz(),$
            (primus_es1_swire_filterlist_irac())[0:1]] ; just ch1-2
       end
       'xmm_swire': begin
          filters = [$
            primus_galex_filterlist(),$
            primus_cfhtls_filterlist_erben(),$
;           primus_sxds_filterlist(),$
            (primus_xmm_swire_filterlist())[0:1]] ; just ch1-2
       end
       'cfhtls_xmm': begin
          filters = [$
            primus_galex_filterlist(),$
            primus_cfhtls_filterlist_erben(),$
            (primus_xmm_swire_filterlist())[0:1]] ; just ch1-2
       end
; uber-XMM field for MF_DERIVE_ZPTOFFSETS
       'xmm': begin
          filters = [$
            primus_galex_filterlist(),$
            primus_cfhtls_filterlist_erben(),$
            (primus_xmm_swire_filterlist())[0:1]] ; just ch1-2
       end
; field to test our cdfs and es1 selection
       'test_cfhtls_xmm': begin
          filters = [$
            primus_galex_filterlist(),$
            primus_cfhtls_filterlist_erben(),$
            (primus_xmm_swire_filterlist())[0:1]] ; just ch1-2
       end
       else: message, 'Code me!'
    endcase

    if keyword_set(noirac) then filters = filters[where(strmatch(filters,'*irac*') eq 0)]
    if keyword_set(nogalex) then filters = filters[where(strmatch(filters,'*galex*') eq 0)]
    
; convert the filter name to something you can show your grandma 
    nfilt = n_elements(filters)
    nice_filters = strarr(nfilt)
    for ii = 0, nfilt-1 do nice_filters[ii] = repstr(strmid(filters[ii],$
      strpos(filters[ii],'_',/reverse_search)+1),'.par','')
;   if strmatch(field,'*cosmos*',/fold) then nice_
    
return, filters
end

