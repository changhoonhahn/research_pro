pro build_smf_im_mf_vmax_jackknife, run, primus=primus, sdss=sdss, literature=literature, prank=prank
    if keyword_set(primus) then fields = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire']
    if keyword_set(sdss) then fields = ['sdss']
    if keyword_set(primus) then samples = ['']
    if keyword_set(sdss) then samples = ['sdss']
    if keyword_set(literature) then litsuffix = '_lit_'$
        else litsuffix = ''

    parpath = get_path(/repo)
    parameters = yanny_readone(parpath+'zerod_environment_parameters.par', hdr=hdr)
    para        = where(parameters.run eq run)
    param       = parameters[para]
    zmin        = param.zmin
    zmax        = param.zmax
    cylrad      = param.cylrad
    cylheight   = param.cylheight
    nbin        = param.nbin
    thresh      = param.thresh

    if keyword_set(primus) then envfname = 'EDP-primus-z0210-numden.fits'
    if keyword_set(sdss) then envfname = 'EDP-sdss-z006_0145-numden.fits'

    radius_string       = strtrim(string(cylrad),2)
    height_string       = strtrim(string(cylheight),2)
    threshold_string    = strtrim(string(thresh),2)
    nbin_string         = strtrim(string(nbin),2)

    if keyword_set(sdss) then zbins = mf_zbins_chh(nzbins, /sdss) 
    if keyword_set(primus) then begin 
        if keyword_set(literature) then zbins = mf_zbins_chh(nzbins, /literature) else zbins = mf_zbins(nzbins)
    endif 
    struct_print, zbins

    binsize = mf_binsize(bin=0.3,minmass=minmass,maxmass=maxmass)

    sf_q = ['active','quiescent']
    if keyword_set(primus) then zbin = ['0204', '0406', '0608', '0810']
    if keyword_set(sdss) then zbin=['nobin']
    hilow = ['hienv','midenv','lowenv']


    for i=0L,nzbins-1L do begin
        for ii=0L,n_elements(hilow)-1L do begin 
            for iii=0L,n_elements(sf_q)-1L do begin
                if keyword_set(primus) then fname='smf_primus_cylr'+radius_string+'h'+height_string+$
                    '_thresh'+threshold_string+'_nbin'+nbin_string+litsuffix+strtrim(sf_q[iii])+'_'+$
                    hilow[ii]+'_'+zbin[i]+'z.fits'
                if keyword_set(sdss) then fname='smf_sdss_cylr'+radius_string+'h'+height_string+$
                    '_thresh'+threshold_string+'_nbin'+nbin_string+litsuffix+strtrim(sf_q[iii])+'_'+$
                    hilow[ii]+'_'+zbin[i]+'z.fits'
                if (sf_q[ii] EQ 'active') then sf_hist = mrdfits(get_path(/smf)+fname,1)
                if (sf_q[ii] EQ 'quiescent') then q_hist = mrdfits(get_path(/smf)+fname,1)
            endfor
            qf =  


                hist = struct_addtags(hist,{nfield:intarr(hist.nbins),thesefields:strarr(5,hist.nbins)})
                ; Obtains the fields that contribute to the SMF at each mass bin: 
                if keyword_set(sdss) then begin
                    hist.nfield = 1
                    hist.thesefields[0] = 'sdss'
                endif else begin 
                   for bb=0L,hist.nbins-1L do begin 
                       if (rev[bb] ne rev[bb+1]) then begin 
                           allfield = (finaldata.field)[rev[rev[bb]:rev[bb+1]-1]]
                           thesefields = allfield[uniq(allfield,sort(allfield))]
                           hist.nfield[bb] = n_elements(thesefields)
                           hist.thesefields[0:hist.nfield[bb]-1,bb] = thesefields
                        endif  
                    endfor 
                endelse 
                
                ; Since SDSS only has one field, it has to divide SDSS by RA and DEC: 
                if keyword_set(sdss) then begin 
                    sdss_hist = hist_nd(transpose([[finaldata.ra],[finaldata.dec]]),$
                        [30D,20D],rev=rev)
                    nbins = n_elements(sdss_hist)
                    notempty = where(hist gt 1000,nfield_jack)
                endif else nfield_jack = n_elements(fields)

                for f=0L,nfield_jack-1L do begin 
                    if keyword_set(sdss) then begin 
                        nn = notempty[ii] 
                        ;toss = these[rev[rev[nn]:rev[nn+1]-1]]
                        ;keep = cmset_op(these,'and',/not2,toss)
                        reweight = 1D
                    endif else begin 
                        keep = where(strtrim(finaldata.field,2) ne fields[f]) 
                        uniq_field  = (finaldata[keep].field)[uniq(finaldata[keep].field)]
                        uarea       = get_poly_area(/primus,/sr,field=uniq_field)
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

                good = where(hist.nfield gt 0, ngood)
                for gg=0L,ngood-1L do begin 
                    if keyword_set(sdss) then begin 
                        nbb = nfield_jack
                        m2 = lindgen(nfield_jack)
                    endif else begin 
                        match, strtrim(hist.thesefields[*,good[gg]],2),fields ,m1,m2
                        nbb = n_elements(m1)
                    endelse 
                    if (nbb ge 3) then begin 
                        hist.phierr_cv[good[gg]] = sqrt((nbb-1.0)/nbb*total((jack_mfdata[m2].phi[good[gg]]$
                            -djs_mean(jack_mfdata[m2].phi[good[gg]]))^2))
                    endif 
                endfor 

                good = where(hist.nfield gt 0 and hist.phierr_cv gt 0.0) 
                fix = where(hist.nfield gt 0 and hist.phierr_cv le 0.0,nfix)
                if (nfix ne 0) then begin 
                    get_element, good, fix, nearest
                    hist.phierr_cv[fix] = hist.phierr_cv[good[nearest]]
                endif 
                
                if keyword_set(primus) then fname='smf_primus_cylr'+radius_string+'h'+height_string+$
                    '_thresh'+threshold_string+'_nbin'+nbin_string+litsuffix+strtrim(sf_q[k])+'_'+$
                    hilow[kk]+'_'+zbin[iii]+'z.fits'
                if keyword_set(sdss) then fname='smf_sdss_cylr'+radius_string+'h'+height_string+$
                    '_thresh'+threshold_string+'_nbin'+nbin_string+litsuffix+strtrim(sf_q[k])+'_'+$
                    hilow[kk]+'_'+zbin[iii]+'z.fits'
                
                print, fname 
                mwrfits, hist, get_path(/smf)+fname, /create
            endfor 
        endfor
    endfor 
end
