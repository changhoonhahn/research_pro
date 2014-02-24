pro plot_smf_mfdata
    fields = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire']
    sfq = ['quiescent','active','all']
    bsize=0.1
    zbins = mf_zbins(nzbins)
    struct_print, zbins

    uber = {redshift:0., mass:0.,vmax:0., weight:0.}
;    dir = get_path(/mfdata)
    dir = '/mount/moon1/ioannis/research/projects/primus/mf/2165/mfs_v20/'
    for i = 0L, n_elements(sfq)-1L do begin
        uberdata = mrdfits(dir+'mfdata_'+sfq[i]+'_supergrid01.fits.gz',1)
        data = replicate({redshift:0., mass:0.,vmax:0., weight:0.,masslimit:0.},n_elements(uberdata))
        data.redshift = uberdata.z
        data.mass = uberdata.mass
;        data.vmax = uberdata.vmax_evol
        data.weight = uberdata.weight/uberdata.vmax_evol
        data.masslimit = uberdata.masslimit
        z = ['2_3','3_4','4_5','5_6','6_7']

        for l=0L,nzbins-2L do begin
            zindx = data.redshift ge zbins[l].zlo and data.redshift lt zbins[l].zup
            minmass = get_min_masslimit(l,sfq[i],/primus)   
            masscomp = data.mass gt minmass

            finaldata=data[where(zindx and masscomp)]
;            hist = hist_smf(finaldata.mass, finaldata.vmax, weight=finaldata.weight, binsize=bsize)
            hist = hist_weighted(finaldata.mass,weight=finaldata.weight, xmin=8.75, binsize=bsize)
            dim=size(hist,/dimensions)

            smfpath = get_path(/smf)
            fname = 'moustakas_'+strtrim(sfq[i])+'_'+z[l]+'_hist.dat'
            print, fname
            openw, lun, smfpath+fname, /get_lun
;                for k=0L, dim[1]-1L do printf,lun, hist[0,k],alog10(hist[1,k]/bsize),$
;                          0.434*(hist[2,k]/hist[1,k]),format='(f,f,f)'
                for k=0L,dim[1]-1L do printf,lun,hist[0,k],alog10(hist[1,k]/bsize),0.0,format='(f,f,f)'
            free_lun, lun
        endfor
    endfor
end 
