function get_min_masslimit, zbin, sfq, literature=literature
    if keyword_set(literature) then ext='_lit' else ext=''

    dir = '/global/data/scr/chh327/primus/science/mf/2165/mfs_v20/'
    masslim = mrdfits(dir+'limits_'+sfq+'_supergrid01'+ext+'.fits.gz',1)
    mlim = masslim.masslim_95

    minmass = min(mlim[zbin,*])
    
    return, minmass
end
