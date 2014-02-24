function read_mf_ubersample_chh, field, photo=photo, xray=xray, zerodver=zerodver
; jm10aug23ucsd - read the desired uber-sample written out by
; BUILD_MF_UBERSAMPLE
    mfpath = '/global/data/scr/chh327/primus/science/mf/2165/ubersample/edp_ubersample/'
    if (n_elements(field) eq 0) then begin
       splog, 'FIELD input required!'
       return, -1
    endif

    if (strtrim(field,2) eq 'sdss') then $
      zerodfile = mfpath+field+'_ubersample.fits.gz' else $
        zerodfile = mfpath+field+'_zerod_'+zerodver+'.fits.gz'

    if (file_test(zerodfile) eq 0) then begin
       splog, 'Zerod file '+zerodfile+' not found!'
       return, -1
    endif
    splog, 'Reading '+zerodfile
    zerod = mrdfits(zerodfile,1)
    
    if (strtrim(field,2) eq 'sdss') then return, zerod

    if arg_present(photo) then begin
       photofile = repstr(zerodfile,'_zerod_','_photo_')
       if (file_test(photofile) eq 0) then begin
          splog, 'Photomerge file '+photofile+' not found!'
          return, zerod
       endif
       splog, 'Reading '+photofile
       photo = mrdfits(photofile,1)
    endif

    if arg_present(xray) then begin
       xrayfile = repstr(zerodfile,'_zerod_','_xray_')
       if (file_test(xrayfile) eq 0) then begin
          splog, 'X-ray file '+xrayfile+' not found!'
          return, zerod
       endif
       splog, 'Reading '+xrayfile
       xray = mrdfits(xrayfile,1)
    endif

return, zerod
end
