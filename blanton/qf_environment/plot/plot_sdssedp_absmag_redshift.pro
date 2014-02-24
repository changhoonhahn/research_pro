pro plot_sdssedp_absmag_redshift, plt=plt
    figdir = '/home/users/hahn/research/figures/primus/'
    psfile = figdir+'fig_SDSS_EDP_absmag_redshift.ps'

    xrange = [0.0375,0.145]
    yrange = [-12.0,-25.7]

    envdir = '/global/data/scr/chh327/primus/data/EDP/'
    env_file = ['EDP-sdss-z00375_0145-numden.fits','EDP-primus-z0210-numden.fits']
    zbin_lo = [0.0375,0.2,0.4,0.6]
    zbin_hi = [0.145,0.4,0.6,0.8]

    targetdir = '/global/data/scr/chh327/primus/data/target/'
    target_file = ['target_cylr2h50_thresh75_nbin5_all_sdss_lit.fits', $
        'target_cylr2h50_thresh75_nbin20_all_lit.fits'] 

    envpsym = 16 & envcolor = 'red' & envsymsize = 0.5
    targpsym = 14 & targcolor = 'black' & targsymsize = 1.5

    nohat = 0

    sdss_target = mrdfits(targetdir+target_file[0],1)
    sdss_env = mrdfits(envdir+env_file[0],1)
    kcorr = mrdfits('/global/data/scr/chh327/primus/science/mf/kcorrect/'+$
        'sdss_kcorr.fits.gz',1)
    kcorredp = mrdfits('/global/data/scr/chh327/primus/science/mf/kcorrect/edp_kcorrect/'+$
        'sdss_kcorr.fits.gz',1)

    targ_min_absmag = min(sdss_target.mg_01, imin) 
    print, sdss_target[imin].ra, sdss_target[imin].dec, sdss_target[imin].redshift, sdss_target[imin].mg_01

    envmatch = where(sdss_env.zprimus EQ sdss_target[imin].redshift, nenvmatch)
;    GT sdss_target[imin].redshift-0.01 AND $
;        sdss_env.zprimus LT sdss_target[imin].redshift+0.01, nenvmatch)
    print, nenvmatch
    print, sdss_env[envmatch].ra, sdss_env[envmatch].dec, sdss_env[envmatch].zprimus, sdss_env[envmatch].mg 

    kcorrmatch = where(kcorr.k_zobj EQ sdss_target[imin].redshift, nkcorrmatch)
    kcorredpmatch = where(kcorredp.k_zobj EQ sdss_target[imin].redshift, nkcorredpmatch)
    print, 'kcorr', nkcorrmatch, '  kcorr_edp', nkcorredpmatch
    print, kcorr[kcorrmatch].k_zobj, kcorredp[kcorredpmatch].k_zobj
    print, kcorr[kcorrmatch].k_maggies
    print, kcorredp[kcorredpmatch].k_maggies
    print, kcorr[kcorrmatch].k_ivarmaggies
    print, kcorredp[kcorredpmatch].k_ivarmaggies
    print, kcorr[kcorrmatch].k_coeffs
    print, kcorredp[kcorredpmatch].k_coeffs
    print, kcorr[kcorrmatch].k_fnuvugrizjhk_absmag_01
    print, kcorredp[kcorredpmatch].k_fnuvugrizjhk_absmag_01

; SDSS
    if keyword_set(plt) then begin 
        im_plotconfig, 1, pos, psfile=psfile, charsize=1.8

        djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
            yrange=yrange, xrange=xrange, xtitle='Redshift '+textoidl('(z)'), $
            ytitle=textoidl('M_{g}'), xtickinterval=0.05 ;xtickname=replicate(' ',10), 
        im_legend, ['Galaxy Sample','EDP'], /right, /bottom, box=0,$
            color=[envcolor,targcolor], psym=[envpsym,targpsym], $
            symsize=[envsymsize,targsymsize]*1.7, $
            symthick=8, spacing=2.7, charsize=1.5

        sdss_targ = mrdfits(targetdir+target_file[0],1)
        targ_masslimit  = where(sdss_targ.mass GT sdss_targ.masslimit, targ_count)
        targ_redshift   = sdss_targ[targ_masslimit].redshift 
        targ_mg_01      = sdss_targ[targ_masslimit].mg_01
        djs_oplot, targ_redshift, targ_mg_01, psym=symcat(targpsym,thick=2), $ 
            color=im_color(targcolor)
       
        sdss_env = mrdfits(envdir+env_file[0],1) 
        env_mg_01 = sdss_env.mg
        env_redshift = sdss_env.zprimus
        djs_oplot, env_redshift , env_mg_01, psym=symcat(envpsym,thick=3), $
            color=im_color(envcolor)

        im_plotconfig, /psclose, psfile=psfile, pdf=pdf
    endif 
return 
end 
