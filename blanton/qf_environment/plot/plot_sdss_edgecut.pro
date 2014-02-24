pro plot_sdss_edgecut, run_sdss
; Checks that the edgecuts on the SDSS environment count data 
; chh - Feb 3, 2014

; Import parameter file: 
    parpath = get_path(/repo)
    parameters = yanny_readone(parpath+'zerod_environment_parameters.par', hdr=hdr)
    run_index   = where(parameters.run eq run_sdss)
    param       = parameters[run_index]
    cylrad      = param.cylrad
    cylheight   = param.cylheight
    threshold   = param.thresh

    radius_string       = strtrim(string(cylrad),2)
    height_string       = strtrim(string(cylheight),2)
    threshold_string    = strtrim(string(threshold),2)

    figdir = '/home/users/hahn/research/figures/primus/'
    psfile = figdir+'fig_sdss_edgecut_radec.ps'

    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8
; Set up the plot: 
; SDSS Panel
    edgepsym = 16 & edgecolor = 'red' & edgesymsize = 1.5
    sdsspsym = 16 & sdsscolor = 'black' & sdsssymsize = 1.4
    
    xrange = [0.0,360.0] 
    yrange = [-15.0,75.0] 

    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle='RA', ytitle='Dec'
    im_legend, ['SDSS','Edge'], /left, /top, box=0, color=[sdsscolor,edgecolor], $ 
        psym=[sdsspsym,edgepsym], symsize=[sdsssymsize,edgesymsize]*1.7, $
        symthick=8, spacing=2.7, charsize=1.5
    nohat = 0
    sfq = ['active','quiescent']
    sdss_ra = []
    sdss_dec = []
    sdss_edgecut = []
    for i=0L,n_elements(sfq)-1L do begin
; Import Environment Count files for both active and quiescent
        envcount_dir = '/global/data/scr/chh327/primus/data/envcount/' 
        envcount_file = 'EDGECUTTEST_envcount_cylr'+radius_string+'h'+height_string+$
            '_thresh'+threshold_string+'_sdss_'+sfq[i]+'_EDP-sdss-z00375_0145-numden.fits'

        envcount_data = mrdfits(envcount_dir+envcount_file, 1)
        print, envcount_dir+envcount_file
        sdss_ra = [sdss_ra, envcount_data.ra]
        sdss_dec = [sdss_dec, envcount_data.dec]
        sdss_edgecut = [sdss_edgecut, envcount_data.edgecut]
    endfor 

    djs_oplot, sdss_ra, sdss_dec, psym=symcat(sdsspsym,thick=2), $ 
        symsize=0.5, color=im_color(sdsscolor)

; Only keeping galaxies with edgecut = 0 (galaxies on the edge)
    edgecut_indx = where(sdss_edgecut EQ 0, edgecut_count)
    edge_sdss_ra = sdss_ra[edgecut_indx]
    edge_sdss_dec = sdss_dec[edgecut_indx] 
    
    djs_oplot, edge_sdss_ra, edge_sdss_dec, psym=symcat(edgepsym,thick=2.0), $ 
        symsize=0.25, color=im_color(edgecolor)

    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
end 
