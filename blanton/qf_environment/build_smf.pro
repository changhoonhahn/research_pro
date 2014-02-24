pro build_smf, run, primus=primus, sdss=sdss, literature=literature, prank=prank
    if keyword_set(primus) then fields = ['es1','cosmos','cfhtls_xmm','cdfs','xmm_swire']
    if keyword_set(sdss) then fields = ['sdss']
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

    fieldfactor = dblarr(n_elements(fields))
    for k=0L,n_elements(fields)-1L do begin
        if keyword_set(primus) then begin 
            fieldfactor[k]= get_poly_area(/primus,/sr)/$
                get_poly_area(/primus,/sr,field=fields[k])
        endif 
        if keyword_set(sdss) then fieldfactor[k]=1.0
    endfor

    bsize = 0.3

    sf_q = ['active','quiescent']
    if keyword_set(primus) then zbin = ['0204', '0406', '0608', '0810']
    if keyword_set(sdss) then zbin=['nobin']
    hilow = ['hienv','midenv','lowenv']

    data_struct = {redshift:0.,mass:0.,masslimit:0.,envcount:0.,vmax:0.,vmaxavail:0.,weight:0.,edgecut:0L, prank:0.,sfq:''}

    for i=0L, n_elements(sf_q)-1L do begin
        for ii=0L, n_elements(fields)-1L do begin
            targetfile = get_path(/envcount)+'envcount_cylr'+radius_string+'h'+height_string+'_thresh'+threshold_string+'_nbin'$
                +nbin_string+'_'+fields[ii]+'_'+sf_q[i]+litsuffix+envfname
            print, targetfile
            target  = mrdfits(targetfile, 1)
            data    = replicate(data_struct,n_elements(target))

            data.redshift   = target.redshift
            data.mass       = target.mass
            data.masslimit  = target.masslimit
            data.envcount   = target.envcount
            data.vmax       = target.vmax*fieldfactor[ii]
            data.vmaxavail  = target.vmaxavail*fieldfactor[ii]
            data.weight     = target.weight
            data.edgecut    = target.edgecut
            data.sfq        = sf_q[i]
            
            if (ii EQ 0L) then field_combine = data $
                else field_combine = [field_combine,data]
        endfor    
        if (i EQ 0L) then allgal_combine = field_combine $
            else allgal_combine = [allgal_combine,field_combine]
    endfor 
    print, 'All Galaxies and all fields combined=',n_elements(allgal_combine)

;    if keyword_set(prank) then begin 
;        lowthresh = 20
;        print, '     z_low      z_high      zero percent'
;        for k=0L,nzbins-1L do begin
;            edgecut = allgal_combine.edgecut eq 1
;            zindx   = allgal_combine.redshift ge zbins[k].zlo and allgal_combine.redshift lt zbins[k].zup
;            zdata   = allgal_combine[where(edgecut and zindx, zcount)]
;
;            zeroenv = zdata.envcount eq 0 
;            zeroenvdata = zdata[where(zeroenv, zerocount)]
;
;            zeropercent = 100.0*float(zerocount)/float(zcount)
;            print, zbins[k].zlo, zbins[k].zup, zeropercent
;            lowthresh = lowthresh>zeropercent
;        endfor
;        print, lowthresh
;        lowthresh = float(ceil(lowthresh))
;        print, 'Low environment percentage rank threshold=', lowthresh
;    endif 

    for iii=0L,nzbins-1L do begin
        edgecut = allgal_combine.edgecut eq 1
        zindx   = allgal_combine.redshift ge zbins[iii].zlo and allgal_combine.redshift lt zbins[iii].zup
        if keyword_set(prank) then begin
            zindx_data = allgal_combine[where(edgecut and zindx, zindx_count)]
            zindx_data.prank = get_percentage_rank(zindx_data.envcount)

            for j=0L,n_elements(sf_q)-1L do begin
                sfqindx = zindx_data.sfq eq sf_q[j]
                
                ;if keyword_set(literature) then minmass = get_min_masslimit(iii,sf_q[j],/literature) $
                ;    else minmass = get_min_masslimit(iii,sf_q[j])
                smcomp = zindx_data.mass gt zindx_data.masslimit

                for jj=0L,n_elements(hilow)-1L do begin 
                    lowthresh = 20.0>min(zindx_data.prank)
;                    print, lowthresh
                    if jj eq 0 then prankindx = zindx_data.prank gt 80.0
                    if jj eq 1 then prankindx = zindx_data.prank gt lowthresh and zindx_data.prank le 80.0
                    if jj eq 2 then prankindx = zindx_data.prank le lowthresh

                    prank_data = zindx_data[where(sfqindx and smcomp and prankindx, prankcount)]
                    print, sf_q[j],zbin[iii],hilow[jj],100.0*float(prankcount)/float(zindx_count),' percent, min=',min(prank_data.envcount),',max=',max(prank_data.envcount)
                    
                    hist = hist_smf(prank_data.mass,prank_data.vmaxavail,xmin=8.75,xmax=12.0,weight=prank_data.weight,binsize=bsize)
                    hist_dim = size(hist,/dimensions)
                    
                    if keyword_set(primus) then fname='smf_primus_cylr'+radius_string+'h'+height_string+$
                        '_thresh'+threshold_string+'_nbin'+nbin_string+litsuffix+strtrim(sf_q[j])+'_prank_'+$
                        hilow[jj]+'_'+zbin[iii]+'z.dat'
                    if keyword_set(sdss) then fname='smf_sdss_cylr'+radius_string+'h'+height_string+$
                        '_thresh'+threshold_string+'_nbin'+nbin_string+litsuffix+strtrim(sf_q[j])+'_prank_'+$
                        hilow[jj]+'_'+zbin[iii]+'z.dat'

;                    print, '****SMF:',fname
                    openw, lun, get_path(/smf)+'hist_'+fname, /get_lun
                        for l=0L, hist_dim[1]-1L do printf,lun, hist[0,l],hist[1,l],hist[2,l],format='(f,f,f)'
                    free_lun, lun
                    openw, lun, get_path(/smf)+fname, /get_lun
                        for l=0L, hist_dim[1]-1L do printf,lun, hist[0,l],alog10(hist[1,l]/bsize),$
                            0.434*(hist[2,l]/hist[1,l]),format='(f,f,f)'
                    free_lun, lun
                endfor
            endfor
        endif else begin       
            smcomp = allgal_combine.mass GT allgal_combine.masslimit
            zbin_data = allgal_combine[where(edgecut AND zindx AND smcomp,zbin_count)]
            for k=0L,n_elements(sf_q)-1L do begin
                sfqindx = zbin_data.sfq EQ sf_q[k]
                for kk=0L,n_elements(hilow)-1L do begin
                    highenvthresh = 1.5
                    if kk eq 0 then env_indx = zbin_data.envcount GT highenvthresh 
                    if kk eq 1 then env_indx = zbin_data.envcount GT 0.0 AND $ 
                        zbin_data.envcount LE highenvthresh 
                    if kk eq 2 then env_indx = zbin_data.envcount EQ 0.0
                    
                    finaldata = zbin_data[where(sfqindx AND env_indx,finaldata_count)]
                    print, zbin[iii], sf_q[k], hilow[kk],float(finaldata_count)/float(zbin_count),max(finaldata.envcount)
                    print, min(finaldata.mass), max(finaldata.mass), min(finaldata.vmaxavail), max(finaldata.vmaxavail), min(finaldata.weight), max(finaldata.weight)

                    hist = hist_smf(finaldata.mass,finaldata.vmaxavail, xmin=8.75,xmax=12.0,weight=finaldata.weight,binsize=bsize)
                    hist_dim = size(hist, /dimensions)
                    
                    if keyword_set(primus) then fname='smf_primus_cylr'+radius_string+'h'+height_string+$
                        '_thresh'+threshold_string+'_nbin'+nbin_string+litsuffix+strtrim(sf_q[k])+'_'+$
                        hilow[kk]+'_'+zbin[iii]+'z.dat'
                    if keyword_set(sdss) then fname='smf_sdss_cylr'+radius_string+'h'+height_string+$
                        '_thresh'+threshold_string+'_nbin'+nbin_string+litsuffix+strtrim(sf_q[k])+'_'+$
                        hilow[kk]+'_'+zbin[iii]+'z.dat'
                    
                    print, fname 
                    openw, lun, get_path(/smf)+'hist_'+fname, /get_lun
                        for l=0L, hist_dim[1]-1L do printf,lun, hist[0,l],hist[1,l],hist[2,l],format='(f,f,f)'
                    free_lun, lun
                    openw, lun, get_path(/smf)+fname, /get_lun
                        for l=0L, hist_dim[1]-1L do printf,lun, hist[0,l],alog10(hist[1,l]/bsize),$
                            0.434*(hist[2,l]/hist[1,l]),format='(f,f,f)'
                    free_lun, lun
                endfor 
            endfor
        endelse
    endfor 
end
