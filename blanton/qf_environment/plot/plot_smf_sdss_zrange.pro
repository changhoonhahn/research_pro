pro plot_smf_sdss_zrange, vmaxcorr=vmaxcorr
; Directories 
    smfdir  = '/global/data/scr/chh327/primus/data/smf/'
    pardir  = get_path(/qfenvpara) 
    figdir  = '/home/users/hahn/research/figures/primus/'
    
    fidmass = 10.5          ; Fiducial mass
    bin_string = '0.25'     ; Mass bin width

; SMF Plot
    if keyword_set(vmaxcorr) then psfile = figdir+'fig_smf_sdss_mfdata_edp_jm_zrangecomp_vmaxcorrection.ps' $
        else psfile = figdir+'fig_smf_sdss_mfdata_edp_jm_zrangecomp_novmaxcorrection.ps'
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.8
        
    xrange = [8.25,11.75]
    yrange = [-6.0,-1.5]

    sdss_zcolor = [60,160]
    z_label = ['0.01-0.2','0.0375-0.145']
    nohat = 0
    allpsym = 16 & allsymsize=1.3

    sfq = ['active','quiescent']
    for iiii=0L, n_elements(sfq)-1L do begin 
        if keyword_set(vmaxcorr) then vmaxlegend = 'Vmax Corrected' $
            else vmaxlegend = 'Vmax Not Corrected'
        if (iiii EQ 0L) then begin  
            djs_plot, [0], [0], /nodata, position=pos[*,iiii], xsty=1, ysty=1, $
                yrange=yrange, xrange=xrange, xtickinterval=0.5, ytickinterval=1.0
            im_legend, z_label, /left, /bottom, box=0,$
                color=[sdss_zcolor], psym=replicate(16,n_elements(z_label)), $
                symsize=replicate(1.5,n_elements(z_label))*1.7, $
                symthick=8, spacing=2.7, charsize=2.1
        endif else begin
            djs_plot, [0], [0], /nodata, /noerase, position=pos[*,iiii], xsty=1, ysty=1, $
                yrange=yrange, xrange=xrange, xtickinterval=0.5,ytickinterval=1.0, $
                ytickname=replicate(' ',10) 
            im_legend, vmaxlegend, /right, /bottom, box=0,charsize=2.1
        endelse  

        im_legend, sfq[iiii], /right, /top, box=0, charsize=2.1

        zcolor=sdss_zcolor & zbin =  ['jm','edp'] 
        if keyword_set(vmaxcorr) then vmaxlabel = '_vmaxcorr' $
            else vmaxlabel = '_novmaxcorr'
        for ii=0L,n_elements(zbin)-1L do begin  
; Import SMF file
            smf_file = 'smf_mfdata_'+sfq[iiii]+'_supergrid01_sdss_nobinz_'+zbin[ii]+'zrange.fits'
            if (zbin[ii] EQ 'edp') then smf_file = 'smf_mfdata_'+sfq[iiii]+'_supergrid01_sdss_nobinz_'+zbin[ii]+'zrange'+vmaxlabel+'.fits'
            smf_data = mrdfits(smfdir+smf_file,1)
            print, smf_file
; Limits on the SMF
            limits = where(smf_data.limit EQ 1)

            smf_lower_err = alog10(smf_data.phi[limits])-(smf_data.phierr[limits]/smf_data.phi[limits]/alog(10))
            smf_upper_err = alog10(smf_data.phi[limits])+(smf_data.phierr[limits]/smf_data.phi[limits]/alog(10))
;                    polyfill, [smf_data.mass[limits], reverse(smf_data.mass[limits])],[smf_lower_err,reverse(smf_upper_err)], $
;                        /data, color=im_color(zcolor[ii],255), noclip=0, /fill;, orientation=45, spacing=0.1
            polyfill, [smf_data.mass[limits], reverse(smf_data.mass[limits])],[smf_lower_err,reverse(smf_upper_err)], $
                /data, color=im_color(zcolor[ii]), noclip=0, /fill;, orientation=135, spacing=0.1

            if (zbin[ii] EQ 'jm') then mysmf = smf_data
        endfor
        jm_mf = mrdfits('/global/data/scr/chh327/primus/mfdata/mf_'+sfq[iiii]+'_supergrid01_bin0.10_evol_sdss.fits.gz',1)
        djs_oplot, jm_mf.mass, alog10(jm_mf.phi),psym=symcat(16, thick=5), color=im_color('red')
    endfor 
; Plot ytitle: 
    xyouts, pos[0,0]-0.06, pos[0,1],textoidl('log (\Phi/Mpc^{-3} dex^{-1})'), align=0.5,$
        /normal, orientation=90
    xyouts, pos[0,1], pos[1,0]-0.125,textoidl('log (M_{*}/M_{\odot})'), align=0.5,$
        /normal, orientation=0
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
end 
