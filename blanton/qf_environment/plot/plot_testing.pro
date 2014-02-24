pro plot_testing
    psfile = 'fig_testing.ps'
    im_plotconfig, 12, pos, psfile=psfile, charsize=1.8
;    , xmargin=[1.3,0.4], $
;        width=6.8, height=[5,2.5], charsize=1.8

    xrange = [9,12.0]
    yrange = [-7.0,-1.5]
    maxis1 = range(xrange[0]+0.15,12.5,100)
    mingal = 0

    sfpsym = 16 & sfcolor = 'dodger blue' & sffitcolor = 'black' & sfline = 5 & sfsymsize = 1.4
    qqpsym = 14 & qqcolor = 'firebrick' & qqfitcolor = 'black' & qqline = 3 & qqsymsize = 1.5
    allpsym = 15 & allcolor = 'black' & allfitcolor = 'black' & allline = 0 & allsymsize = 1.3

    nohat = 0

    print, pos

    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, xtitle='', xtickname=replicate(' ',10), $
        ytitle=mfplot_phititle(),xtickinterval=1

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
        yrange=yrange, xrange=xrange, ytickname=replicate(' ',10), $
        ytitle=mfplot_phititle(), xtitle=textoidl('log (\Phi \Phi \Phi)'), xtickinterval=1
    im_legend, ['All','Quiescent','Star-Forming'], /left, /bottom, box=0, $
        color=[allcolor,qqcolor,sfcolor], psym=[allpsym,qqpsym,sfpsym], $
        symsize=[allsymsize,qqsymsize,sfsymsize]*1.7, $
        symthick=8, spacing=2.7, charsize=2.1

end 
