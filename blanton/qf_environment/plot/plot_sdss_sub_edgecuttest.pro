pro plot_sdss_sub_edgecuttest
    figdir = '/home/users/hahn/research/figures/primus/'
    psfile = figdir+'fig_sdss_sub_edgecuttest_radec.ps'

    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8
; Set up the plot: 
; SDSS Panel
    edgepsym = 14 & edgecolor = 'red' & edgesymsize = 1.5
    sdsspsym = 16 & sdsscolor = 'black' & sdsssymsize = 1.4
    
    xrange = [100.0,150.0] 
    yrange = [35.0,45.0] 

    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle='RA', ytitle='Dec'
    im_legend, ['SDSS','Edge'], /left, /top, box=0, color=[sdsscolor,edgecolor], $ 
        psym=[sdsspsym,edgepsym], symsize=[sdsssymsize,edgesymsize]*1.7, $
        symthick=8, spacing=2.7, charsize=1.5
    nohat = 0
    sfq = ['active','quiescent']

    ransack = mrdfits('/global/data/scr/chh327/primus/data/ransack/ransack_sdss_edpfield_1000000.fits',1)
    windowindx = where(ransack.ra GT xrange[0] AND ransack.ra LT xrange[1] AND $
        ransack.dec GT yrange[0] AND ransack.dec LT yrange[1])
    ransack_sub = ransack[windowindx]
    
    djs_oplot, ransack_sub.ra, ransack_sub.dec, psym=symcat(sdsspsym,thick=2), $ 
        color=im_color('blue')
; Import Environment Count files for both active and quiescent
    envcount_dir = '/global/data/scr/chh327/primus/data/envcount/' 
    envcount_file = 'EDGECUTTEST_envcount_cylr2h25_thresh75_nbin2_sdss.fits' 
    envcount_data = mrdfits(envcount_dir+envcount_file, 1)

    djs_oplot, envcount_data.ra, envcount_data.dec, psym=symcat(sdsspsym,thick=2), $ 
        color=im_color(sdsscolor)

; Only keeping galaxies with edgecut = 0 (galaxies on the edge)
    edgecut_indx = where(envcount_data.edgecut EQ 0, edgecut_count)
    print, edgecut_count, n_elements(envcount_data)
    edge_sdss_ra = envcount_data[edgecut_indx].ra
    edge_sdss_dec = envcount_data[edgecut_indx] .dec
    
    djs_oplot, edge_sdss_ra, edge_sdss_dec, psym=symcat(edgepsym,thick=2), $ 
        color=im_color(edgecolor)

    aroundgal = mrdfits('/home/users/hahn/research/pro/blanton/qf_environment/sdss_ransack_aroundgalaxy.fits',1)
    djs_oplot, aroundgal.ra, aroundgal.dec, psym=symcat(edgepsym,thick=2), $ 
        color=im_color('green')

    singgal = mrdfits('/home/users/hahn/research/pro/blanton/qf_environment/sdss_target_singlegalaxy.fits',1)
    djs_oplot, [singgal.ra], [singgal.dec], psym=symcat(edgepsym,thick=3), $ 
        color=im_color('coral')
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
end 
