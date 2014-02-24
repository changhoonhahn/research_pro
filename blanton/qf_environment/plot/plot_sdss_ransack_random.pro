pro plot_sdss_ransack_random
; Checks the RA and Dec coverage of the ransack-EDPfield and random-targetfield
; chh - Feb 3, 2014

    figdir = '/home/users/hahn/research/figures/primus/'
    psfile = figdir+'fig_sdss_ransack_random_radec.ps'

    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8
; Set up the plot: 
; SDSS Panel
    randompsym = 14 & randomcolor = 'red' & randomsymsize = 1.5
    ransackpsym = 16 & ransackcolor = 'black' & ransacksymsize = 1.4
    
    xrange = [0.0,360.0] 
    yrange = [-15.0,75.0] 

    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle='RA', $
        ytitle='Dec'
    im_legend, ['Random','Ransack'], /left, /top, box=0, color=[randomcolor,ransackcolor], $ 
        psym=[randompsym,ransackpsym], symsize=[randomsymsize,ransacksymsize]*1.7, $
        symthick=8, spacing=2.7, charsize=1.5
    nohat = 0
    
    ransack_dir = '/global/data/scr/chh327/primus/data/ransack/' 
    ransack_data = mrdfits(ransack_dir+'ransack_sdss_edpfield_1000000.fits',1)
    random_data = mrdfits(ransack_dir+'random_sdss_targetfield_1000000.fits',1)
    
    djs_oplot, ransack_data.ra, ransack_data.dec, psym=symcat(ransackpsym,thick=2), $ 
        color=im_color(ransackcolor)

    djs_oplot, random_data.ra, random_data.dec, psym=symcat(randompsym,thick=2), $ 
        color=im_color(randomcolor)

    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
end 
