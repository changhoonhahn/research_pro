pro build_smf_mfdata,primus=primus, sdss=sdss, literature=literature,john=john,$
    jmzrange=jmzrange, edpzrange=edpzrange, vmaxcorr=vmaxcorr
    if keyword_set(sdss) then begin 
        samples = ['_sdss']
        if keyword_set(jmzrange) then zbins = mf_zbins(nzbins, /sdss) 
        if keyword_set(edpzrange) then zbins = mf_zbins_chh(nzbins, /sdss) 
        zbin=['nobin']
        litsuffix = ''
    endif

    if keyword_set(primus) then begin 
        samples = ['']
        if keyword_set(literature) then begin 
            zbins = mf_zbins_chh(nzbins, /literature) 
            zbin = ['0204', '0406', '0608', '0810']
        endif else begin 
            zbins = mf_zbins(nzbins)
            zbin = ['0203','0304','0405','05065','06508','0810']
        endelse 
    if keyword_set(literature) then litsuffix = '_lit'$
        else litsuffix = ''
    endif 
    struct_print, zbins

    binsize = mf_binsize(bin=binsize,minmass=minmass,maxmass=maxmass)

    sf_q = ['active','quiescent']
    
    for i=0L, n_elements(sf_q)-1L do begin
        datafile = get_path(/mfdata)+'mfdata_chh/new/mfdata_'+strtrim(sf_q[i])+'_supergrid01'$
            +samples[0]+litsuffix+'.fits.gz'
        if keyword_set(john) then $
            datafile = get_path(/mfdata)+'mfdata_'+strtrim(sf_q[i])+'_supergrid01'$
            +samples[0]+litsuffix+'.fits.gz'
        if keyword_set(sdss) then begin 
            if keyword_set(vmaxcorr) then vmaxcorrdir='env_zrange/' $
                else vmaxcorrdir=''
            datafile = get_path(/mfdata)+'mfdata_chh/new/'+vmaxcorrdir+'mfdata_'+strtrim(sf_q[i])+'_supergrid01'$
                +samples[0]+litsuffix+'.fits.gz'
        endif 
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
            fname='smf_mfdata_new_'+strtrim(sf_q[i])+'_supergrid01'$
                +samples[0]+litsuffix+'_'+zbin[iz]+'z.fits'
            if keyword_set(john) then $ 
                fname='smf_mfdata_'+strtrim(sf_q[i])+'_supergrid01'$
                +samples[0]+litsuffix+'_'+zbin[iz]+'z_john.fits'
            if keyword_set(sdss) then begin 
                if keyword_set(jmzrange) then begin 
                    zrangeflag = '_jmzrange'
                    vmaxcorrflag = ''
                endif 
                if keyword_set(edpzrange) then begin 
                    zrangeflag = '_edpzrange'
                    if keyword_set(vmaxcorr) then vmaxcorrflag = '_vmaxcorr' $
                        else vmaxcorrflag = '_novmaxcorr'
                endif 
                fname='smf_mfdata_'+strtrim(sf_q[i])+'_supergrid01'$
                    +samples[0]+litsuffix+'_'+zbin[iz]+'z'+zrangeflag+vmaxcorrflag+'.fits'
            endif
            print, fname 
            
            mwrfits, mfdata1, get_path(/smf)+fname, /create
        endfor
    endfor 
end
