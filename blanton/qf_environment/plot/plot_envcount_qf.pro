pro plot_envcount_qf, primus_run, sdss_run, twoenv=twoenv
; plots the quiescent fraction for galaxies in redshift z = 0-1 in different
; environments
    fidmass = 10.5
    bin_string = '0.25'
; Parameter file
    pardir  = get_path(/qfenvpara) 
    parameters = yanny_readone(pardir+'zerod_environment_parameters.par', hdr=hdr)
    para        = where(parameters.run eq primus_run)
    param       = parameters[para]
    cylrad      = strtrim(string(param.cylrad),2)
    cylheight   = strtrim(string(param.cylheight),2)
    thresh      = strtrim(string(param.thresh),2)
    param_label = 'R = '+cylrad+' H = '+cylheight

; QF Plot
    figdir  = '/home/users/hahn/research/figures/primus/qf/'
    psfile = figdir+'fig_qf_cylr'+cylrad+'h'+cylheight+'_thresh'+thresh+'_bin'+bin_string+'.eps'
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.8
        
    xrange = [8.25,12.25]
    yrange = [0.0,1.0]

    loadct, 10, /silent, file='/home/users/hahn/research/pro/im_plot/fsc_brewer.tbl'
    sdss_zcolor = [210]
    primus_zcolor = [160,110,60]
    z_label = ['0.0375-0.145','0.2-0.4','0.4-0.6','0.6-0.8']
    nohat = 0
    allpsym = 16 & allsymsize=1.3

    surveys = ['sdss', 'primus']
    if keyword_set(twoenv) then envbin_string = '2envbin' $
        else envbin_string = '3envbin'
    if keyword_set(twoenv) then envbin = ['lowenv', 'hienv'] $
        else envbin = ['lowenv', 'midenv', 'hienv'] 

    fit_A = []
    fit_B = [] 
    sig_A = []
    sig_B = []
    for iii=0L,n_elements(envbin)-1L do begin   
        if (iii EQ 0L) then begin  
            plot, [0], [0], /nodata, position=pos[*,iii], xsty=1, ysty=1, $
                yrange=yrange, xrange=xrange, ytitle='Quiescent Fraction', $
                xtickinterval=0.5, color=im_color('black',255)
            im_legend, z_label, /left, /top, box=0,$
                color=[sdss_zcolor,primus_zcolor], psym=replicate(16,n_elements(z_label)), $
                symsize=replicate(1.5,n_elements(z_label))*1.7, $
                symthick=8, spacing=2.7, charsize=2.1
        endif else begin
            plot, [0], [0], /nodata, /noerase, position=pos[*,iii], xsty=1, ysty=1, $
                yrange=yrange, xrange=xrange, ytickname=replicate(' ',10), $
                xtickinterval=0.5, color=im_color('black',255)
            im_legend, param_label, /left, /top, box=0, charsize=2.1
        endelse  

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
; import QF file
                qf_file = get_path(/qfenvcount)+'qf_'+surveys[i]+'_cylr'+rad_string+$
                    'h'+h_string+'_thresh'+thresh_string+litsuffix+envbin[iii]+'_'+$
                    zbin[ii]+'_bin'+bin_string+'_'+envbin_string+'.fits'
                print, qf_file
                qf_data = mrdfits(qf_file,1)
                ; impose stellar mass limits on the quiescent fraction
                limits = where(qf_data.limit EQ 1 AND qf_data.qf GT 0.0 and qf_data.qf LT 1.0)

                qf_lower_err = qf_data.qf[limits]-qf_data.err_cv[limits]
                qf_upper_err = qf_data.qf[limits]+qf_data.err_cv[limits]
                polyfill, [qf_data.mass[limits], reverse(qf_data.mass[limits])],[qf_lower_err,reverse(qf_upper_err)], $
                    /data, color=im_color(zcolor[ii],255), noclip=0, /fill

; linear fit for the QF
                fitlimits = where(qf_data.limit EQ 1 AND qf_data.qf GT 0.0 and qf_data.qf LT 1.0, fitlimcount)
                fitexy, qf_data.mass[fitlimits]-fidmass, qf_data.qf[fitlimits], fitA, fitB, $
                    x_sig=replicate(0.0,fitlimcount), y_sig=qf_data.err_cv[fitlimits], sigma_A_B

                print, 'y-int= ', fitA, 'slope= ', fitB, '  errors= ', sigma_A_B
                fit_A = [fit_A, fitA]
                fit_B = [fit_B, fitB]
                sig_A = [sig_A, sigma_A_B[0]]
                sig_B = [sig_B, sigma_A_B[1]]
            endfor
        endfor 
    endfor
    print, fit_A
    print, sig_A
    xyouts, pos[0,1], pos[1,0]-0.15,textoidl('log (M_{*}/M_{\odot})'), align=0.5,$
        /normal, orientation=0, color=im_color('black',255) 
    xyouts, pos[0,0]+0.2, pos[2,1],'Low Env', align=0.5,$
        /normal, orientation=0, color=im_color('black',255) 
    xyouts, pos[0,1]+0.2, pos[2,1],'High Env', align=0.5,$
        /normal, orientation=0, color=im_color('black',255) 
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
    
    para        = where(parameters.run eq primus_run)
    param       = parameters[para]
    cylrad      = strtrim(string(param.cylrad),2)
    cylheight   = strtrim(string(param.cylheight),2)
    thresh      = strtrim(string(param.thresh),2)
    ; QF Fit Plot
    psfile = figdir+'fig_qffit_cylr'+cylrad+'h'+cylheight+'_thresh'+thresh+'_bin'+bin_string+'.eps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8 
    
    xrange = [0.0,1.0]
    yrange = [0.0,1.0]

    highpsym = 16 & highcolor = 'black' & highsymsize = 1.5
    midpsym = 15 & midcolor = 'black' & midsymsize = 1.5
    lowpsym = 14 & lowcolor = 'black' & lowsymsize = 1.5

    plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle='Redshift '+textoidl('z'), $
        ytitle='Quiescent Fraction at '+textoidl('M_{fid}'), xtickinterval=0.1,$
        color=im_color('black',255)
    im_legend, ['High Density','Low Density'], /right, /top, box=0,$
        color=[highcolor,lowcolor], psym=[highpsym,lowpsym], $
        symsize=[highsymsize,lowsymsize]*1.7, $
        symthick=8, spacing=2.7, charsize=2.1
    im_legend, param_label, /top, /left, box=0, charsize=2.1
    
    zlo = [0.0375,0.2,0.4,0.6]
    zup = [0.145,0.4,0.6,0.8]
    zval = [0.09125,0.3,0.5,0.7]
    for iii=0L,n_elements(envbin)-1L do begin
        if (envbin[iii] EQ 'lowenv') then envpsym = lowpsym 
        if (envbin[iii] EQ 'midenv') then envpsym = midpsym
        if (envbin[iii] EQ 'hienv') then envpsym = highpsym
        
        ; QF fit data 
        qffit_dir = '/global/data/scr/chh327/primus/data/qf/qffit/'
        openw, lun, qffit_dir+'qffitdata_cylr'+cylrad+'h'+cylheight+'_thresh'+thresh+$
            '_bin'+bin_string+'_'+envbin[iii]+'.dat',/get_lun

        for ii=0L,n_elements(zval)-1L do begin
            qffit = fit_A[iii*n_elements(zval)+ii]
;            qffit_sig = sig_A[iii*n_elements(zval)+ii]+sig_B[iii*n_elements(zval)+ii]*fidmass
;            qffit_sig = sqrt((sig_A[iii*n_elements(zval)+ii])^2+(sig_B[iii*n_elements(zval)+ii])^2)
            qffit_sig = sig_A[iii*n_elements(zval)+ii]
            print, qffit_sig
            oploterror, zval[ii], qffit, qffit_sig, $
                psym=symcat(envpsym,thick=5), symsize=1.5, errthick=5, color=im_color('black',255), $
                errcolor=im_color('black',255), /hibar, nohat=nohat
            oploterror, zval[ii], qffit, qffit_sig, $
                psym=3, symsize=1.5, errthick=5, color=im_color('black',255), $
                errcolor=im_color('black',255), /lobar, nohat=nohat
            printf,lun,zval[ii],zlo[ii],zup[ii],qffit,qffit_sig
        endfor
        free_lun,lun
    endfor
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
end 
