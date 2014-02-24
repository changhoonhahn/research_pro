pro plot_edp_absmag_redshift
    figdir = '/home/users/hahn/research/figures/primus/'
    psfile = figdir+'fig_EDP_absmag_redshift.eps'
    im_plotconfig, 30, pos, psfile=psfile, charsize=1.8

    xrange = [0.0375,0.145]
    yrange = [-12.0,-25.7]

    envdir = '/global/data/scr/chh327/primus/data/EDP/'
    env_file = ['EDP-sdss-z00375_0145-numden.fits','EDP-primus-z0210-numden.fits']
    zbin_lo = [0.0375,0.2,0.4,0.6]
    zbin_hi = [0.145,0.4,0.6,0.8]

    targetdir = '/global/data/scr/chh327/primus/data/target/'
    target_file = ['target_cylr2h50_thresh75_nbin5_all_sdss_lit.fits', $
        'target_cylr2h50_thresh75_nbin20_all_lit.fits'] 

    envpsym = 16 & envcolor = 'red' & envsymsize = 1.4
    targpsym = 14 & targcolor = 'black' & targsymsize = 1.5

    nohat = 0

; SDSS
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, ytitle=textoidl('M_{r}'), xtickinterval=0.05 
    im_legend, ['EDP','Galaxy'], /right, /bottom, box=0,$
        color=[envcolor,targcolor], psym=[envpsym,targpsym], $
        symsize=[envsymsize,targsymsize]*1.7, $
        symthick=8, spacing=2.7, charsize=1.5

    sdss_targ = mrdfits(targetdir+target_file[0],1)
    targ_masslimit  = where(sdss_targ.mass GT sdss_targ.masslimit, targ_count)
    targ_redshift   = sdss_targ[targ_masslimit].redshift 
    targ_mg_01      = sdss_targ[targ_masslimit].mg_01
    targ_mr_01      = sdss_targ[targ_masslimit].mr_01
    djs_oplot, targ_redshift, targ_mr_01, psym=symcat(targpsym,thick=2), $ 
        color=im_color(targcolor)
   
    sdss_env = mrdfits(envdir+env_file[0],1) 
    env_mg_01 = sdss_env.absmag[3]
    env_mr_01 = sdss_env.absmag[4]
    env_redshift = sdss_env.redshift
    djs_oplot, env_redshift , env_mr_01, psym=symcat(envpsym,thick=3), $
        color=im_color(envcolor)

    primus_targ = mrdfits(targetdir+target_file[1],1)
    primus_env = mrdfits(envdir+env_file[1],1)

    for i=0L,n_elements(zbin_lo)-2L do begin 
        xrange = [zbin_lo[i+1L],zbin_hi[i+1L]]
        print, xrange
        djs_plot, [0], [0], /nodata, /noerase, position=pos[*,i+1], xsty=1, ysty=1, $
            yrange=yrange, xrange=xrange, ytickname=replicate(' ',10), $
            ytitle='', xtickinterval=0.1

        targ_cut  = where(primus_targ.redshift GT zbin_lo[i+1L] AND $
            primus_targ.redshift LT zbin_hi[i+1L], cutcount)
        targ_redshift   = primus_targ[targ_cut].redshift 
        targ_mg_01      = primus_targ[targ_cut].mg_01
        targ_mr_01      = primus_targ[targ_cut].mr_01
        djs_oplot, targ_redshift, targ_mr_01, psym=symcat(targpsym,thick=2), $ 
            color=im_color(targcolor)
        print, cutcount
       
        env_cut  = where(primus_env.redshift GT zbin_lo[i+1L] AND $
            primus_env.redshift LT zbin_hi[i+1L], cutcount)
        env_mg_01 = primus_env[env_cut].absmag[3]
        env_mr_01 = primus_env[env_cut].absmag[4]
        env_redshift = primus_env[env_cut].redshift
        djs_oplot, env_redshift , env_mr_01, psym=symcat(envpsym,thick=3), $
            color=im_color(envcolor)
        print, cutcount
    endfor
    print, pos
    xyouts, pos[0,2], pos[1,0]-0.15,'Redshift '+textoidl('(z)'), align=0.5,$
        /normal, orientation=0, color=im_color('black')
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
end 
