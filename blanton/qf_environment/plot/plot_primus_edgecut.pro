pro plot_primus_edgecut, run_primus
; Checks that the edgecuts on the SDSS environment count data 
; chh - Feb 3, 2014

; Import parameter file: 
    parpath = get_path(/repo)
    parameters = yanny_readone(parpath+'zerod_environment_parameters.par', hdr=hdr)
    run_index   = where(parameters.run eq run_primus)
    param       = parameters[run_index]
    cylrad      = param.cylrad
    cylheight   = param.cylheight
    threshold   = param.thresh

    radius_string       = strtrim(string(cylrad),2)
    height_string       = strtrim(string(cylheight),2)
    threshold_string    = strtrim(string(threshold),2)

    figdir = '/home/users/hahn/research/figures/primus/'
    psfile = figdir+'fig_primus_edgecut_radec.ps'

    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8
; Set up the plot: 
; SDSS Panel
    edgepsym = 7 & edgecolor = 'red' & edgesymsize = 1.5
    noedgepsym = 14 & noedgecolor = 'blue' & noedgesymsize = 1.5
    primuspsym = 16 & primuscolor = 'black' & primussymsize = 1.4
    
    xrange = [32.0,38.0] 
    yrange = [-8.0,-3.0] 

    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle='RA', ytitle='Dec'
    im_legend, ['SDSS','Edge','Not Edge'], /left, /top, box=0, color=[primuscolor,edgecolor,noedgecolor], $ 
        psym=[primuspsym,edgepsym,noedgepsym], symsize=[primussymsize,edgesymsize,noedgesymsize]*1.7, $
        symthick=8, spacing=2.7, charsize=1.5
    nohat = 0
    sfq = ['active','quiescent']
    primus_ra = []
    primus_dec = []
    primus_edgecut = []

; Import Environment Count files for both active and quiescent and combine
    for i=0L,n_elements(sfq)-1L do begin
        envcount_dir = '/global/data/scr/chh327/primus/data/envcount/' 
        envcount_file = 'EDGECUTTEST_envcount_cylr'+radius_string+'h'+height_string+$
            '_thresh'+threshold_string+'_'+sfq[i]+'_lit_EDP-primus-z0210-numden.fits'
        envcount_data = mrdfits(envcount_dir+envcount_file, 1)
        print, envcount_dir+envcount_file
        primus_ra = [primus_ra, envcount_data.ra]
        primus_dec = [primus_dec, envcount_data.dec]
        primus_edgecut = [primus_edgecut, envcount_data.edgecut]
    endfor 

    djs_oplot, primus_ra, primus_dec, psym=symcat(primuspsym,thick=2), $ 
        color=im_color(primuscolor)

; Only keeping galaxies with edgecut = 1 (galaxies not on the edge)
    noedge_indx = where(primus_edgecut EQ 1, noedge_count)
    noedge_primus_ra = primus_ra[noedge_indx]
    noedge_primus_dec = primus_dec[noedge_indx] 
    
    djs_oplot, noedge_primus_ra, noedge_primus_dec, psym=symcat(noedgepsym,thick=1), $ 
        color=im_color(noedgecolor)

; Only keeping galaxies with edgecut = 0 (galaxies on the edge)
    edge_indx = where(primus_edgecut EQ 0, edge_count)
    edge_primus_ra = primus_ra[edge_indx]
    edge_primus_dec = primus_dec[edge_indx] 
    
    djs_oplot, edge_primus_ra, edge_primus_dec, psym=symcat(edgepsym,thick=1), $ 
        psymsize=0.1, color=im_color(edgecolor)

    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
end 
