function read_mf_parent_chh, field, quiescent=quiescent, active=active, indx=indx
; jm10aug25ucsd - read the desired parent sample
    
    mfpath = get_path(/parent)
    if (n_elements(field) eq 0) then begin
       splog, 'FIELD input required!'
       return, -1
    endif
    parentfile = mfpath+'mf_parent_'+field+'.fits.gz'
    if (file_test(parentfile) eq 0) then begin
       splog, 'Parent file '+parentfile+' not found!'
       return, -1
    endif
    splog, 'Reading '+parentfile
    parent = mrdfits(parentfile,1)

;; select quiescent or star-forming galaxies    
;    if keyword_set(quiescent) or keyword_set(active) then begin
;       mz = parent.k_ugrizjhk_absmag_05[4]
;       umz = parent.k_ugrizjhk_absmag_05[0]-mz
;       indx = mf_select_quiescent(umz,mz=mz,z=parent.zprimus,active=active)
;       parent = parent[indx]
;    endif else indx = lindgen(n_elements(parent))
    
return, parent
end
