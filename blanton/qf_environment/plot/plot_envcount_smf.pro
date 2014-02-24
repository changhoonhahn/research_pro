pro plot_envcount_smf, primus_run, sdss_run, twoenv=twoenv
    fidmass = 10.5          ; Fiducial mass
    bin_string = '0.25'     ; Mass bin width

; Parameter file
    parameters = yanny_readone(get_path(/qfenvpara)+'zerod_environment_parameters.par', hdr=hdr)
    para        = where(parameters.run eq primus_run)
    param       = parameters[para]
    cylrad      = strtrim(string(param.cylrad),2)
    cylheight   = strtrim(string(param.cylheight),2)
    thresh      = strtrim(string(param.thresh),2)
    param_label  = 'R = '+cylrad+' H = '+cylheight

; SMF Plot
    figdir  = '/home/users/hahn/research/figures/primus/smf/'
    psfile = figdir+'fig_smf_cylr'+cylrad+'h'+cylheight+'_thresh'+thresh+'_bin'+bin_string+'.eps'
    if keyword_set(twoenv) then begin 
        im_plotconfig, 2, pos, psfile=psfile, charsize=1.8
    endif else begin 
        im_plotconfig, 14, pos, psfile=psfile, charsize=1.8, $
            xmargin=[1.1,0.3], width=3.6*[1,1,1], height=2.6*[1,1]
    endelse 
        
    xrange = [8.25,12.25]
    yrange = [-6.0,-1.75]
    
    z_label = ['0.0375-0.145','0.2-0.4','0.4-0.6','0.6-0.8']
    nohat = 0
    allpsym = 16 & allsymsize=1.3

    surveys = ['sdss', 'primus']
    if keyword_set(twoenv) then envbin_string = '2envbin' $
        else envbin_string = '3envbin'
    if keyword_set(twoenv) then envbin = ['lowenv', 'hienv'] $
        else envbin = ['lowenv', 'midenv', 'hienv'] 
    envbin_label = ['Low Env', 'High Env']
    sfq = ['active','quiescent']

    for iiii=0L, n_elements(sfq)-1L do begin 
        for iii=0L,n_elements(envbin)-1L do begin   
; Color scheme for the SMFs
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
; Initializing the plot
            if (iiii EQ 0L AND iii EQ 0L) then begin  
                plot, [0], [0], /nodata, position=pos[*,iiii*n_elements(envbin)+iii], xsty=1, ysty=1, $
                    yrange=yrange, xrange=xrange, xtickinterval=0.5, ytickinterval=1.0, xtickname=replicate(' ',10),$
                    color=im_color('black',255)
                im_legend, z_label, /left, /bottom, box=0,$
                    color=[sdss_zcolor,primus_zcolor], psym=replicate(16,n_elements(z_label)), $
                    symsize=replicate(1.5,n_elements(z_label))*1.7, $
                    symthick=8, spacing=2.7, charsize=2.1
            endif else if (iii NE 0L AND iiii NE 1L) then begin 
                plot, [0], [0], /nodata, /noerase, position=pos[*,iiii*n_elements(envbin)+iii], xsty=1, ysty=1, $
                    yrange=yrange, xrange=xrange, xtickinterval=0.5,ytickinterval=1.0, ytickname=replicate(' ',10), $
                    xtickname=replicate(' ',10), color=im_color('black',255)
            endif else if (iii NE 0L AND iiii EQ 1L) then begin 
                plot, [0], [0], /nodata, /noerase, position=pos[*,iiii*n_elements(envbin)+iii], xsty=1, ysty=1, $
                    yrange=yrange, xrange=xrange, xtickinterval=0.5,ytickinterval=1.0, ytickname=replicate(' ',10),$
                    color=im_color('black',255)
                im_legend, textoidl(param_label), /left, /bottom, box=0, charsize=2.1
            endif else begin
                plot, [0], [0], /nodata, /noerase, position=pos[*,iiii*n_elements(envbin)+iii], xsty=1, ysty=1, $
                    yrange=yrange, xrange=xrange, xtickinterval=0.5,ytickinterval=1.0,color=im_color('black',255)
            endelse  
            if (iiii EQ 0L ) then im_legend, envbin_label[iii], /right, /top, box=0, charsize=2.1
            

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
                    smfdir  = '/global/data/scr/chh327/primus/data/smf/'
                    smf_file = 'smf_'+surveys[i]+'_cylr'+rad_string+'h'+h_string+'_thresh'+thresh_string+$
                        litsuffix+sfq[iiii]+'_'+envbin[iii]+'_'+zbin[ii]+'_bin'+bin_string+'_'+$
                        envbin_string+'.fits'
                    smf_data = mrdfits(smfdir+smf_file,1)
                    print, smf_file
    ; Limits on the SMF
                    limits = where(smf_data.limit EQ 1)

    ;                oploterror, qf_data.mass[limits], qf_data.qf[limits], qf_data.err[limits], $
    ;                    psym=symcat(allpsym,thick=5), symsize=allsymsize, errthick=5, color=im_color(zcolor[ii]), $
    ;                    errcolor=im_color(zcolor[ii]), /hibar, nohat=nohat
    ;                oploterror, qf_data.mass[limits], qf_data.qf[limits], qf_data.err[limits], $
    ;                    psym=3, symsize=allsymsize, errthick=5, color=im_color(zcolor[ii]), $
    ;                    errcolor=im_color(zcolor[ii]), /lobar, nohat=nohat

                    smf_lower_err = alog10(smf_data.phi[limits])-(smf_data.phierr_cv[limits]/smf_data.phi[limits]/alog(10))
                    smf_upper_err = alog10(smf_data.phi[limits])+(smf_data.phierr_cv[limits]/smf_data.phi[limits]/alog(10))
;                    polyfill, [smf_data.mass[limits], reverse(smf_data.mass[limits])],[smf_lower_err,reverse(smf_upper_err)], $
;                        /data, color=im_color(zcolor[ii],255), noclip=0, /fill;, orientation=45, spacing=0.1
                    polyfill, [smf_data.mass[limits], reverse(smf_data.mass[limits])],[smf_lower_err,reverse(smf_upper_err)], $
                        /data, color=im_color(zcolor[ii],255), noclip=0, /fill;, orientation=135, spacing=0.1
                endfor
            endfor 
        endfor
    endfor 
; Plot ytitle: 
    xyouts, pos[0,0]-0.06, pos[1,0],textoidl('log (\Phi/Mpc^{-3} dex^{-1})'), align=0.5,$
        /normal, orientation=90, color=im_color('black',255)
    xyouts, pos[3,0]+0.04, pos[1,0]+0.2,textoidl('Starforming'), align=0.5,$
        /normal, orientation=90, color=im_color('black',255)
    xyouts, pos[3,0]+0.04, pos[1,2]+0.2,textoidl('Quiescent'), align=0.5,$
        /normal, orientation=90, color=im_color('black',255)
    xyouts, pos[0,1], pos[1,2]-0.1,textoidl('log (M_{*}/M_{\odot})'), align=0.5,$
        /normal, orientation=0, color=im_color('black',255)
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
end 
