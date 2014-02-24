pro build_smf_primus
; Simplified version of build_smf_mfdata using John's mfdata catalogs  
; builds SMF for PRIMUS data 
    samples = ['']
    zbins = mf_zbins_chh(nzbins, /literature) 
    zbin = ['0204', '0406', '0608', '0810']
    litsuffix = '_lit'
    struct_print, zbins

    binsize = smf_binsize(bin=0.01,minmass=minmass,maxmass=maxmass)
    minmass = 9.5 
    maxmass = 11.5-binsize/2.0
    ;    sf_q = ['active','quiescent']
    sf_q = ['all']
    
    for i=0L, n_elements(sf_q)-1L do begin
        datafile = get_path(/mfdata)+'mfdata_'+strtrim(sf_q[i])+'_supergrid01'$
        +samples[0]+litsuffix+'.fits.gz'
        print, datafile
        data = mrdfits(datafile, 1)
        print, min(data.z), max(data.z)

        for iz=0L,nzbins-1L do begin
            these = where((data.z ge zbins[iz].zlo) and (data.z lt zbins[iz].zup),ngal)
            print, zbins[iz].zlo, zbins[iz].zup
            data_zbin = data[these]
            vmax = data_zbin.vmax_evol

            mfdata1 = im_mf_vmax(data_zbin.mass,data_zbin.weight/vmax,$
                binsize=binsize,minmass=minmass,maxmass=maxmass,$
                masslimit=data_zbin.masslimit,rev=rev)
            smf_file = 'smf_mfdata_'+strtrim(sf_q[i],2)+'_supergrid01'$ 
                +samples[0]+litsuffix+'_'+zbin[iz]+'z_mbin'+strmid(strtrim(string(binsize),2),0,4)$
                +'_john.fits'
            smf_dir = '/global/data/scr/chh327/tinker/group_catalog/smf/'
            print, smf_dir+smf_file
            mwrfits, mfdata1, smf_dir+smf_file, /create
        endfor
    endfor 
end
