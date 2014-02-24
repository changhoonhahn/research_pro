pro plot_primus_sfr_mstar
; Plot the SFR versus Mstar relationship for PRIMUS galaxies
; Mass range and redshift range hardcoded
    mass = fltarr(15)  
    mass_min = 9.5
    mass_max = 11.5
    mass[0] = mass_min  
    for i=1L,n_elements(mass)-1L do mass[i] = mass_min+float(i)*(mass_max-mass_min)/float(n_elements(mass))
    
    SFR = fltarr(n_elements(mass)) 
    sig = fltarr(n_elements(mass)) 
    SSFR = fltarr(n_elements(mass)) 
    for i=0L,n_elements(mass)-1L do begin 
        SFR[i] = get_primus_sfr_mstar(mass[i],0.9,sigma=mass_sigma,ssfr=ssfr_tmp,/silent)
        sig[i] = mass_sigma
        SSFR[i] = ssfr_tmp
    endfor 
    
    for i=0L,n_elements(mass)-1L do print, mass[i],SFR[i],sig[i],SSFR[i] 
; Set up plot
    psfile = '/home/users/hahn/research/figures/tinker/primus_sfr_mstar.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8

    xrange = [min(mass), max(mass)] 
    yrange = [0.0, 2.0]
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle=textoidl('log (M_{*}/M_{\odot})'), $
        ytitle=textoidl('log (SFR [M_{\odot} yr^{-1}])') 
    im_legend, 'PRIMUS '+textoidl('z=0.9'), /right, /top, box=0, charsize=2
    djs_oplot, mass, SFR, psym=symcat(16,thick=2.0), $
        color=im_color('blue')
    djs_oplot, mass, SFR+sig, line=2, thick=4, color=im_color('blue')
    djs_oplot, mass, SFR-sig, line=2, thick=4, color=im_color('blue')

    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
    
    psfile = 'primus_ssfr_mstar.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8

    xrange = [min(mass), max(mass)] 
    yrange = [-11.0,-7.0]
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle=textoidl('log (M_{*}/M_{\odot})'), $
        ytitle=textoidl('log (SSFR [yr^{-1}])') 
    djs_oplot, mass, SSFR, psym=symcat(16,thick=2.0), $
        color=im_color('black')

    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
return
end
