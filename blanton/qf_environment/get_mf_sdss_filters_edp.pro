function get_mf_sdss_filters_edp, field, nice_filters=nice_filters, noirac=noirac, $
  nogalex=nogalex, allcosmos=allcosmos
; jm10aug28ucsd - return the filters to be used, given the field;
;   FIT_FILTERINDEX specifies the filters that should be used for fitting
;   and K-corrections (see MF_ISEDFIT) 
;   splog, 'Testing which filters to use for fitting!!!'

    case field of
       'sdss': begin
          filters = [ primus_sdss_filterlist() ]           
       end
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

