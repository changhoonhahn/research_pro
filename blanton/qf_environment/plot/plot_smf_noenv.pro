pro plot_smf_noenv, primus_run, sdss_run, noedge=noedge
; Directories 
    smfdir  = '/global/data/scr/chh327/primus/data/smf/'
    pardir  = get_path(/qfenvpara) 
    figdir  = '/home/users/hahn/research/figures/primus/qf/'
    
    fidmass = 10.5          ; Fiducial mass
    bin_string = '0.10'     ; Mass bin width

; Parameter file
    parameters = yanny_readone(pardir+'zerod_environment_parameters.par', hdr=hdr)
    para        = where(parameters.run eq primus_run)
    param       = parameters[para]
    cylrad      = strtrim(string(param.cylrad),2)
    cylheight   = strtrim(string(param.cylheight),2)
    thresh      = strtrim(string(param.thresh),2)

; SMF Plot
    if keyword_set(noedge) then noedge_flag = '_noedge'$ 
        else noedge_flag = ''
    psfile = figdir+'fig_smf_cylr'+cylrad+'h'+cylheight+'_thresh'+thresh+'_bin'+bin_string+'_noenv'+noedge_flag+'.ps'
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.8
        
    xrange = [8.25,12.00]
    yrange = [-6.0,-1.5]

    z_label = ['0.0375-0.145','0.2-0.4','0.4-0.6','0.6-0.8']
    nohat = 0
    allpsym = 16 & allsymsize=1.3

    surveys = ['sdss', 'primus']
    sfq = ['active','quiescent']

    for iiii=0L, n_elements(sfq)-1L do begin 
        if (sfq[iiii] eq 'active') then begin 
            loadct, 24, /silent, file='/home/users/hahn/research/pro/im_pro/etc/fsc_brewer.tbl'
            sdss_zcolor = [250] 
            primus_zcolor = [223,197,170]
        endif
        if (sfq[iiii] eq 'quiescent') then begin 
            loadct, 10, /silent, file='/home/users/hahn/research/pro/im_pro/etc/fsc_brewer.tbl'
            sdss_zcolor = [210] 
            primus_zcolor = [160,110,60] 
        endif

        if (iiii EQ 0L) then begin  
            plot, [0], [0], /nodata, position=pos[*,iiii], xsty=1, ysty=1, $
                yrange=yrange, xrange=xrange, xtickinterval=0.5, ytickinterval=1.0, $
                color=im_color('black',255)
            im_legend, z_label, /left, /bottom, box=0,$
                color=[sdss_zcolor,primus_zcolor], psym=replicate(16,n_elements(z_label)), $
                symsize=replicate(1.5,n_elements(z_label))*1.7, $
                symthick=8, spacing=2.7, charsize=2.1
        endif else begin
            plot, [0], [0], /nodata, /noerase, position=pos[*,iiii], xsty=1, ysty=1, $
                yrange=yrange, xrange=xrange, xtickinterval=0.5,ytickinterval=1.0, $
                ytickname=replicate(' ',10), color=im_color('black',255)
            if keyword_set(noedge) then im_legend, 'No Edge Galaxies', /left, /bottom, box=0, charsize=2.1 $
                else im_legend, 'Edge Galaxies', /left, /bottom, box=0, charsize=2.1 
        endelse  

        im_legend, sfq[iiii], /right, /top, box=0, charsize=2.1

        for i=0L,n_elements(surveys)-1L do begin  
            if (surveys[i] EQ 'sdss') then begin run = sdss_run & zcolor=sdss_zcolor & zbin =  ['nobinz'] & litsuffix='_' & endif 
            if (surveys[i] EQ 'primus') then begin run = primus_run & zcolor=primus_zcolor & zbin $
                    = ['0204z','0406z','0608z'] & litsuffix='_lit_' & endif

            para        = where(parameters.run eq run)
            param       = parameters[para]
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
          
            for ii=0L,n_elements(zbin)-1L do begin  
; Import SMF file
                smf_file = 'smf_'+surveys[i]+'_cylr'+rad_string+'h'+h_string+'_thresh'+thresh_string+$
                        litsuffix+sfq[iiii]+'_'+zbin[ii]+'_bin'+bin_string+'_noenv'+noedge_flag+'.fits'
                smf_data = mrdfits(smfdir+smf_file,1)
                print, smf_file
; Limits on the SMF
                limits = where(smf_data.limit EQ 1)


                smf_lower_err = alog10(smf_data.phi[limits])-(smf_data.phierr[limits]/smf_data.phi[limits]/alog(10))
                smf_upper_err = alog10(smf_data.phi[limits])+(smf_data.phierr[limits]/smf_data.phi[limits]/alog(10))
                polyfill, [smf_data.mass[limits], reverse(smf_data.mass[limits])],[smf_lower_err,reverse(smf_upper_err)], $
                    /data, color=im_color(zcolor[ii],255), noclip=0, /fill;, orientation=135, spacing=0.1
;                if (surveys[i] EQ 'sdss') then begin 
;                    sdss_data = smf_data
;                endif 

            endfor 
        endfor
        jm_mf = mrdfits('/global/data/scr/chh327/primus/mfdata/mf_'+sfq[iiii]+'_supergrid01_bin0.10_evol_sdss.fits.gz',1)
        djs_oplot, jm_mf.mass, alog10(jm_mf.phi), linestyle=0, linewidth=2, color=im_color('red')
        
        jm_mf = mrdfits('/global/data/scr/chh327/primus/mfdata/mf_'+sfq[iiii]+'_supergrid01_bin0.10_evol_lit.fits.gz',1)
        for i=0L,2L do begin
            djs_oplot, jm_mf[i].mass, alog10(jm_mf[i].phi), linestyle=1, linewidth=2, color=im_color('black')
        endfor
    endfor 
; Plot ytitle: 
    xyouts, pos[0,0]-0.06, pos[0,1],textoidl('log (\Phi/Mpc^{-3} dex^{-1})'), align=0.5,$
        /normal, orientation=90, color=im_color('black',255)
    xyouts, pos[0,1], pos[1,0]-0.125,textoidl('log (M_{*}/M_{\odot})'), align=0.5,$
        /normal, orientation=0, color=im_color('black',255)
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf

;############################################################
; Difference between John's SMF and the edge corrected SMF: 
;############################################################
; Parameter file
    parameters = yanny_readone(pardir+'zerod_environment_parameters.par', hdr=hdr)
    para        = where(parameters.run eq primus_run)
    param       = parameters[para]
    cylrad      = strtrim(string(param.cylrad),2)
    cylheight   = strtrim(string(param.cylheight),2)
    thresh      = strtrim(string(param.thresh),2)

    psfile = figdir+'fig_delta_smf_cylr'+cylrad+'h'+cylheight+'_thresh'+thresh+'_bin'+bin_string+'_noenv'+noedge_flag+'.ps'
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.8

    yrange = [-1.0, 1.0]
    xrange = [8.25,12.00]

    sdss_zcolor = [60]
    primus_zcolor = [110,160,210]
    z_label = ['0.0375-0.145','0.2-0.4','0.4-0.6','0.6-0.8']
    nohat = 0
    allpsym = 16 & allsymsize=1.3

    surveys = ['sdss', 'primus']
    sfq = ['active','quiescent']
        
    for iiii=0L, n_elements(sfq)-1L do begin 
        if (sfq[iiii] eq 'active') then begin 
            loadct, 24, /silent, file='/home/users/hahn/research/pro/im_pro/etc/fsc_brewer.tbl'
            sdss_zcolor = [250] 
            primus_zcolor = [223,197,170]
        endif
        if (sfq[iiii] eq 'quiescent') then begin 
            loadct, 10, /silent, file='/home/users/hahn/research/pro/im_pro/etc/fsc_brewer.tbl'
            sdss_zcolor = [210] 
            primus_zcolor = [160,110,60] 
        endif
; Configure plot legends and etc
        if (iiii EQ 0L) then begin  
            plot, [0], [0], /nodata, position=pos[*,iiii], xsty=1, ysty=1, $
                yrange=yrange, xrange=xrange, xtickinterval=0.5, ytickinterval=1.0, $
                color=im_color('black',255)
            im_legend, z_label, /left, /bottom, box=0,$
                color=[sdss_zcolor,primus_zcolor], psym=replicate(16,n_elements(z_label)), $
                symsize=replicate(1.5,n_elements(z_label))*1.7, $
                symthick=8, spacing=2.7, charsize=2.1
        endif else begin
            plot, [0], [0], /nodata, /noerase, position=pos[*,iiii], xsty=1, ysty=1, $
                yrange=yrange, xrange=xrange, xtickinterval=0.5,ytickinterval=1.0, $
                ytickname=replicate(' ',10), color=im_color('black',255)
            if keyword_set(noedge) then im_legend, 'No Edge Galaxies', /left, /bottom, box=0, charsize=2.1 $
                else im_legend, 'Edge Galaxies', /left, /bottom, box=0, charsize=2.1 
        endelse  

        im_legend, sfq[iiii], /right, /top, box=0, charsize=2.1

; Read in John's SMFs: 
        for i=0L,n_elements(surveys)-1L do begin  
            if (surveys[i] EQ 'sdss') then begin 
                run = sdss_run & zcolor=sdss_zcolor & zbin =  ['nobinz'] & litsuffix='_' 
                jm_mf = mrdfits('/global/data/scr/chh327/primus/mfdata/mf_'+sfq[iiii]+'_supergrid01_bin0.10_evol_sdss.fits.gz',1)
            endif 
            if (surveys[i] EQ 'primus') then begin 
                run = primus_run & zcolor=primus_zcolor & zbin $
                    = ['0204z','0406z','0608z'] & litsuffix='_lit_'
                jm_mf = mrdfits('/global/data/scr/chh327/primus/mfdata/mf_'+sfq[iiii]+'_supergrid01_bin0.10_evol_lit.fits.gz',1)
            endif

            para        = where(parameters.run eq run)
            param       = parameters[para]
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
          
            for ii=0L,n_elements(zbin)-1L do begin  
; Import SMF file
                smf_file = 'smf_'+surveys[i]+'_cylr'+rad_string+'h'+h_string+'_thresh'+thresh_string+$
                        litsuffix+sfq[iiii]+'_'+zbin[ii]+'_bin'+bin_string+'_noenv'+noedge_flag+'.fits'
                smf_data = mrdfits(smfdir+smf_file,1)
                print, smf_file
; Limits on the SMF
                limits = where(smf_data.limit EQ 1)

                jm_smf = jm_mf[ii]
                if array_equal(jm_smf.mass, smf_data.mass) then begin 
                    del_phi = alog10(smf_data.phi)-alog10(jm_smf.phi)
                endif else begin
                    print, "ERROR NOT HAPPY"
                    stop 
                endelse 
                djs_oplot, smf_data.mass[limits], del_phi[limits], $
                    psym=symcat(16,thick=5),color=im_color(zcolor[ii])
                djs_oplot, smf_data.mass[limits], del_phi[limits], $ 
                    linestyle=0,linewidth=3,color=im_color(zcolor[ii])
                print, min(del_phi[limits]), max(del_phi[limits])
            endfor 
        endfor
    endfor 
    
    xyouts, pos[0,0]-0.06, pos[0,1],textoidl('\Delta log (\Phi/Mpc^{-3} dex^{-1})'), align=0.5,$
        /normal, orientation=90, color=im_color('black',255)
    xyouts, pos[0,1], pos[1,0]-0.125,textoidl('log (M_{*}/M_{\odot})'), align=0.5,$
        /normal, orientation=0, color=im_color('black',255)
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
end 
