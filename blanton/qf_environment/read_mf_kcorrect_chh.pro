function read_mf_kcorrect_chh, field, zerodver=zerodver
    mfpath = '/global/data/scr/chh327/primus/science/mf/kcorrect/edp_kcorrect/'
    if (n_elements(field) eq 0) then begin
       splog, 'FIELD input required!'
       return, -1
    endif

    if strtrim(field,2) eq 'sdss' then suffix = '' else suffix = '_'+zerodver
    kcorrfile = mfpath+field+'_kcorr'+suffix+'.fits.gz'

    if (file_test(kcorrfile) eq 0) then begin
       splog, 'K-corrections file '+kcorrfile+' not found!'
       return, -1
    endif
    splog, 'Reading '+kcorrfile
    kcorr = mrdfits(kcorrfile,1)
    
return, kcorr
end
