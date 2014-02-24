pro build_smf_im_mf_vmax, run, nolit=nolit, twoenv=twoenv
; builds SMF for the environment count data using im_mf_vmax.pro
    parpath = get_path(/repo)
    parameters = yanny_readone(parpath+'zerod_environment_parameters.par', hdr=hdr)
    para        = where(parameters.run eq run)
    param       = parameters[para]
    survey      = strtrim(param.name,1)
    zmin        = param.zmin
    zmax        = param.zmax
    cylrad      = param.cylrad
    cylheight   = param.cylheight
    thresh      = param.thresh
    highenvthresh   = param.highenvthresh
    lowenvthresh    = param.lowenvthresh
    rad_string      = strtrim(string(cylrad),2)
    h_string        = strtrim(string(cylheight),2)
    thresh_string   = strtrim(string(thresh),2)
    
    if (survey EQ 'primus') then begin 
        fields = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire'] 
        samples = ''
        if keyword_set(nolit) then litsuffix='' $
            else litsuffix = '_lit'
        if keyword_set(nolit) then zbins = mf_zbins(nzbins) $
            else zbins = mf_zbins(nzbins, /literature) 
        if keyword_set(nolit) then zbin = ['0203','0304','0405','05065','06508','0810'] $
            else zbin = ['0204', '0406', '0608', '0810']
        envfname = 'EDP-primus-z0210-numden.fits'
    endif
    if (survey EQ 'sdss') then begin 
        fields = ['sdss']
        samples = '_sdss' 
        zbins = mf_zbins_chh(nzbins, /sdss) 
        zbin = ['nobin']
        litsuffix=''
;       envfname = 'EDP-sdss-z006_0145-numden.fits'
        envfname = 'EDP-sdss-z00375_0145-numden.fits'
    endif 
    struct_print, zbins

    binsize = mf_binsize(bin=0.25,minmass=minmass,maxmass=maxmass)

    sf_q = ['active','quiescent']
    envbin = ['hienv','midenv','lowenv']
    if keyword_set(twoenv) then envbin = ['hienv','lowenv']

    for i=0L, n_elements(sf_q)-1L do begin
; import environment count data for active/quiescent galaxies
        datafile = get_path(/envcount)+'envcount_cylr'+rad_string+'h'+h_string+$
            '_thresh'+thresh_string+samples+'_'+sf_q[i]+litsuffix+'_'+envfname
        data = mrdfits(datafile, 1)
        print, datafile 
        print, n_elements(data), sf_q[i]+' Galaxies' 
            
        edgecut = where(data.edgecut eq 1, n_edgecut)   ; impose edgecuts on the galaxies
        sfq_data = data[edgecut] 
        print, n_elements(sfq_data), sf_q[i]+' Galaxies within edgecut' 

        for iz=0L,nzbins-1L do begin
; impose redshift limits
            zindx   = where((sfq_data.redshift ge zbins[iz].zlo) and $
                (sfq_data.redshift lt zbins[iz].zup),zbin_count)
            zbin_data = sfq_data[zindx]
            print, zbin_count, ' Galaxies within redshift bin '+strmid(strtrim(string(zbins[iz].zlo),2),0,4)+'-'$
                +strmid(strtrim(string(zbins[iz].zup),2),0,4)
            for ee=0L,n_elements(envbin)-1L do begin
; impose environment limits 
                if keyword_set(twoenv) then begin 
                    if (envbin[ee] EQ 'hienv') then env_indx = where(zbin_data.envcount GT highenvthresh, n_envbin)
                    if (envbin[ee] EQ 'lowenv') then env_indx = where(zbin_data.envcount LT lowenvthresh, n_envbin) 
                endif else begin
                    if ee eq 0 then env_indx = where(zbin_data.envcount GE highenvthresh, n_envbin) 
                    if ee eq 1 then env_indx = where(zbin_data.envcount GE lowenvthresh AND $ 
                        zbin_data.envcount LT highenvthresh, n_envbin)
                    if ee eq 2 then env_indx = where(zbin_data.envcount LT lowenvthresh, n_envbin)
                endelse 
                finaldata = zbin_data[env_indx]
                print, n_envbin, ' Galaxies in '+envbin[ee]
                
                vmax = finaldata.vmaxavail

                mfdata = im_mf_vmax(finaldata.mass,finaldata.weight/vmax,$
                    binsize=binsize,minmass=minmass,maxmass=maxmass,$
                    masslimit=finaldata.masslimit,rev=rev)
                mfdata = struct_addtags(mfdata,{nfield:intarr(mfdata.nbins),thesefields:strarr(5,mfdata.nbins)})
                print, zbin[iz], sf_q[i], envbin[ee],float(n_envbin)/float(zbin_count),max(finaldata.envcount),mfdata.ngal

                ; Obtains the fields that contribute to the SMF at each mass bin: 
                if (survey EQ 'sdss') then begin
                    mfdata.nfield = 1
                    mfdata.thesefields[0] = 'sdss'
                endif else begin 
                   for bb=0L,mfdata.nbins-1L do begin 
                       if (rev[bb] ne rev[bb+1]) then begin 
                           allfield = (finaldata.field)[rev[rev[bb]:rev[bb+1]-1]]
                           thesefields = allfield[uniq(allfield,sort(allfield))]
                           mfdata.nfield[bb] = n_elements(thesefields)
                           mfdata.thesefields[0:mfdata.nfield[bb]-1,bb] = thesefields
                        endif  
                    endfor 
                endelse 
                
                ; Since SDSS only has one field, it has to divide SDSS by RA and DEC: 
                if (survey EQ 'sdss') then begin 
                    sdss_hist = hist_nd(transpose([[finaldata.ra],[finaldata.dec]]),$
                        [30D,20D],rev=sdss_rev)
                    notempty = where(sdss_hist gt 100,nfield_jack)
                endif else nfield_jack = n_elements(fields)

                for f=0L,nfield_jack-1L do begin 
                    if (survey EQ 'sdss') then begin 
                        nn = notempty[f] 
                        toss = [sdss_rev[sdss_rev[nn]:sdss_rev[nn+1]-1]]
                        keep = cmset_op(lindgen(n_elements(finaldata)),'and',/not2,toss)
                        reweight = 1D
                    endif else begin 
                        keep = where(strtrim(finaldata.field,2) ne fields[f]) 
                        uniq_field  = (finaldata[keep].field)[uniq(finaldata[keep].field)]
                        uarea       = get_poly_area(/primus,/sr,field=strtrim(uniq_field,2))
                        reweight    = get_poly_area(/primus,/sr)/uarea
                    endelse 

                    jack_hist = im_mf_vmax(finaldata[keep].mass,reweight*finaldata[keep].weight/finaldata[keep].vmaxavail,$
                        binsize=binsize,minmass=minmass,maxmass=maxmass,masslimit=finaldata[keep].masslimit)
                    if (f EQ 0L) then begin
                        jack_mfdata = jack_hist
                    endif else begin
                        jack_mfdata = [jack_mfdata, jack_hist]
                    endelse 
                endfor
                
                good = where(mfdata.nfield gt 0, ngood)
                for gg=0L,ngood-1L do begin 
                    if (survey EQ 'sdss') then begin 
                        nbb = nfield_jack
                        m2 = lindgen(nfield_jack)
                    endif else begin 
                        match, strtrim(mfdata.thesefields[*,good[gg]],2),fields ,m1,m2
                        nbb = n_elements(m1)
                    endelse 
                    if (nbb ge 3) then begin 
                        mfdata.phierr_cv[good[gg]] = sqrt((nbb-1.0)/nbb*total((jack_mfdata[m2].phi[good[gg]]$
                            -djs_mean(jack_mfdata[m2].phi[good[gg]]))^2))
                    endif 
                endfor 

                good = where(mfdata.nfield gt 0 and mfdata.phierr_cv gt 0.0) 
                fix = where(mfdata.nfield gt 0 and mfdata.phierr_cv le 0.0,nfix)
                if (nfix ne 0) then begin 
                    get_element, good, fix, nearest
                    mfdata.phierr_cv[fix] = mfdata.phierr_cv[good[nearest]]
                endif 
                
                if (survey EQ 'primus') then fname='smf_primus_cylr'+rad_string+'h'+h_string+$
                    '_thresh'+thresh_string+litsuffix+'_'+strtrim(sf_q[i])+'_'+$
                    envbin[ee]+'_'+zbin[iz]+'z_bin'+string(binsize,format='(F4.2)')+'_'+$
                    strtrim(string(n_elements(envbin)),2)+'envbin.fits'
                if (survey EQ 'sdss') then fname='smf_sdss_cylr'+rad_string+'h'+h_string+$
                    '_thresh'+thresh_string+litsuffix+'_'+strtrim(sf_q[i])+'_'+$
                    envbin[ee]+'_'+zbin[iz]+'z_bin'+string(binsize,format='(F4.2)')+'_'+$
                    strtrim(string(n_elements(envbin)),2)+'envbin.fits'
                print, fname 
                mwrfits, mfdata, get_path(/smf)+fname, /create
            endfor 
        endfor
    endfor 
end
