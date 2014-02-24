function get_ransack_command, n=n, primus=primus, sdss=sdss
    if keyword_set(primus) then begin
;        dir ='/global/data/primus/survey_regions/fields/'
        fields = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire']
        dir = '/global/data/scr/chh327/primus/survey_regions/fields/'
        polygon=strarr(n_elements(fields))
        
        fname='ransack_'
        ransack='ransack -r '+strtrim(string(n),1)+' '
        for i=0L,n_elements(fields)-1L do begin
            polygon[i]=dir+'ransack_'+fields[i]+'.ply'
            fname=fname+fields[i]+'_'
            ransack=ransack+polygon[i]+' '
        endfor
        
        fname=fname+strtrim(string(n),1)+'.dat'
        ransack=ransack+fname
;        fname = 'ransack_cdfs_cfhtls_cosmos_es1_xmmswire_'+strtrim(string(n),1)+'.dat'
;        ransack ='ransack -r '+strtrim(string(n),1)+' '+cdfs+' '+cfhtls+' '+cosmos+' '+es1+' '+xmmswire+' '+fname
    endif
    if keyword_set(sdss) then begin
;        dir = '/mount/moon1/ioannis/research/projects/primus/mf/2165/'
        dir = '/global/data/scr/chh327/primus/science/mf/2165/' 
        sdss = dir+'dr72bsafe0_galex_final.ply'
        fname = 'ransack_sdss_'+strtrim(string(n),1)+'.dat'
        ransack = 'ransack -r '+strtrim(string(n),1)+' '+sdss+' '+fname
    endif
    return, ransack
end
