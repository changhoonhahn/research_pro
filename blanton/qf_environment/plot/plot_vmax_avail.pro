pro plot_vmax_avail, run_primus, run_sdss, n_ransack, n_random
; Import parameter file: 
    parpath = get_path(/repo)
    parameters = yanny_readone(parpath+'zerod_environment_parameters.par', hdr=hdr)
    run_index   = where(parameters.run eq run_primus)
    param       = parameters[run_index]
    cylrad      = param.cylrad
    cylheight   = param.cylheight
    primus_nbin = param.nbin
    threshold   = param.thresh

    run_index   = where(parameters.run eq run_sdss)
    param       = parameters[run_index]
    sdss_nbin = param.nbin

    radius_string       = strtrim(string(cylrad),2)
    height_string       = strtrim(string(cylheight),2)
    threshold_string    = strtrim(string(threshold),2)
    primus_nbin_string  = strtrim(string(primus_nbin),2)
    sdss_nbin_string    = strtrim(string(sdss_nbin),2)

    rsk_string = strtrim(string(n_ransack),2)
    rnd_string = strtrim(string(n_random),2)

    figdir = '/home/users/hahn/research/figures/primus/'
    psfile = figdir+'fig_vcomov_avail.ps'
    im_plotconfig, 6, pos, psfile=psfile, charsize=1.8
    
    nohat = 0

; Import Vmaxavail files
    vmax_dir = '/global/data/scr/chh327/primus/data/vmax_avail/'
    vmax_primus_file = vmax_dir+'vmax_avail_cylr'+radius_string+'h'+height_string+$
        '_thresh'+threshold_string+'_nbin'+primus_nbin_string+$
        '_es1_cosmos_cfhtls_xmm_cdfs_xmm_swire_ran'+rnd_string+'_rsk'+rsk_string+'.fits'
    vmax_sdss_file = vmax_dir+'vmax_avail_cylr'+radius_string+'h'+height_string+$
        '_thresh'+threshold_string+'_nbin'+sdss_nbin_string+$
        '_sdss_ran'+rnd_string+'_rsk'+rsk_string+'.fits'

    vmax_primus = mrdfits(vmax_primus_file, 1)
    print, vmax_primus_file
    vmax_sdss = mrdfits(vmax_sdss_file, 1)
    print, vmax_sdss_file

; Set up the plot: 
; SDSS Panel
    primuspsym = 14 & primuscolor = 'black' & primussymsize = 1.5
    sdsspsym = 16 & sdsscolor = 'red' & sdsssymsize = 1.4
    
    xrange = [0.0,0.2]
    yrange = [0.0, lf_comvol(0.2)]

    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle='Redshift '+textoidl('(z)'), $
        ytitle=textoidl('V_{comov,avail}'), xtickinterval=0.1 
    im_legend, ['SDSS','PRIMUS',textoidl('V_{comvo}')], linestyle=[0,0,1], $
        /left, /top, box=0, color=[sdsscolor,primuscolor,'black'], $ 
        psym=[sdsspsym,primuspsym,0], symsize=[sdsssymsize,primussymsize,0]*1.7, $
        symthick=8, spacing=2.7, charsize=1.5
    djs_oplot, vmax_sdss.z, vmax_sdss.v_max_avail, psym=symcat(sdsspsym,thick=2), $ 
        color=im_color(sdsscolor)
    djs_oplot, vmax_sdss.z, lf_comvol(vmax_sdss.z), line=2, thick=4, color=im_color('black')
    
    print, vmax_sdss.v_max_avail/lf_comvol(vmax_sdss.z)

; PRIMUS panel
    xrange = [0.2,1.0]
    yrange = [0.0, lf_comvol(1.0)]
    djs_plot, [0], [0], /nodata, position=pos[*,1], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle='Redshift '+textoidl('(z)'), $
        ytitle=textoidl('V_{comov,avail}'), xtickinterval=0.1 

    djs_oplot, vmax_primus.z, vmax_primus.v_max_avail, psym=symcat(primuspsym,thick=2), $ 
        color=im_color(primuscolor)
    djs_oplot, vmax_primus.z, lf_comvol(vmax_primus.z), line=2, thick=4, color=im_color('black')
    print, vmax_primus.v_max_avail/lf_comvol(vmax_primus.z)

    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
end 
