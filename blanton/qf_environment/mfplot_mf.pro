pro outline_mf, mass, phimin, phimax, color=color, line=line
    color = 'black'
    nm = n_elements(zbins)
    djs_oplot, [mass[0],mass[0]], [phimin[0],phimax[0]], line=line, color=im_color(color)
    djs_oplot, [mass[nm-1],mass[nm-1]], [phimin[nm-1],phimax[nm-1]], line=line, thick=8, color=im_color(color)
    for ii = 0, nm-1 do begin
       djs_oplot, [mass[ii],mass[ii]], [phimax[ii],phimax[ii]], line=line, color=im_color(color)
       djs_oplot, [mass[ii],mass[ii]], [phimin[ii],phimin[ii]], line=line, color=im_color(color)
    endfor
;    for ii = 1, nm-1 do begin
;       djs_oplot, [zbins[ii].zlo,zbins[ii].zlo], $
;         [phi[ii-1]+phierr[ii-1],phi[ii]+phierr[ii]], line=0, color=im_color(color)
;       djs_oplot, [zbins[ii].zlo,zbins[ii].zlo], $
;         [phi[ii-1]-phierr[ii-1],phi[ii]-phierr[ii]], line=0, color=im_color(color)
;    endfor
return
end

pro mfplot_mf, noevol=noevol, pdf=pdf
; jm10aug25ucsd - plot the MF results

    mfpath = mf_path()
    mfspath = mf_path(/mfs)
    qapath = mf_path(/qaplots)
    if keyword_set(pdf) then begin
       paperpath = qapath
       suffix = '.ps'
    endif else begin
       paperpath = mf_path(/paper)
       suffix = '.eps'
    endelse

; the "no-evolution" plots are just for comparison; always write them
; out as QAplots
    if keyword_set(noevol) then begin
       esuffix = 'noevol' 
       paperpath = qapath
       suffix = '.ps'
    endif else begin
       esuffix = 'evol'
    endelse

; read the supergrid parameter file    
    h100 = mf_h100()

    massrange1 = [8.7,12.0]
    phirange1 = [-5.2,-1.7]
    allmassbelow = 1.0
    sfmassbelow = 1.0
    qqmassbelow = 0.5
    
; ---------------------------------------------------------------------------
; Figure 3: sdss - all, quiescent, and active samples
    psfile = paperpath+'mf_sdss'+suffix
    im_plotconfig, 6, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.8, height=[5,2.5], charsize=1.8

    xrange = [9,12.0]
    yrange = [-7.0,-1.5]
    maxis1 = range(xrange[0]+0.15,12.5,100)
    mingal = 0

    sfpsym = 16 & sfcolor = 'dodger blue' & sffitcolor = 'black' & sfline = 5 & sfsymsize = 1.4
    qqpsym = 14 & qqcolor = 'firebrick' & qqfitcolor = 'black' & qqline = 3 & qqsymsize = 1.5
    allpsym = 15 & allcolor = 'black' & allfitcolor = 'black' & allline = 0 & allsymsize = 1.3

    nohat = 0
    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle='', xtickname=replicate(' ',10), $
      ytitle=mfplot_phititle(), xtickinterval=1
    im_legend, ['All','Quiescent','Star-Forming'], /left, /bottom, box=0, $
      color=[allcolor,qqcolor,sfcolor], psym=[allpsym,qqpsym,sfpsym], $
      symsize=[allsymsize,qqsymsize,sfsymsize]*1.7, $
      symthick=8, spacing=2.7, charsize=2.1
    
; calculate the "quenching mass" for the paper
    qq = alog10(mf_schechter(maxis1,schechter=read_mf_vmax('sdss',/quie,/best)))
    aa = alog10(mf_schechter(maxis1,schechter=read_mf_vmax('sdss',/act,/best)) )
    splog, 'quenching mass ', interpol(maxis1,qq-aa,0.0)

; integrate the MF    
    totrho = read_mf_vmax('sdss',/rho)
    sfrho = read_mf_vmax('sdss',/active,/rho)
    qqrho = read_mf_vmax('sdss',/quiescent,/rho)
    splog, 'total, sf, qq, sf_fration, qq_fraction ', 10^totrho.rho, 10^sfrho.rho, $
      10^qqrho.rho, 10^sfrho.rho/10^totrho.rho, 10^qqrho.rho/10^totrho.rho
    
; now make the plot    
    all = read_mf_vmax('sdss',/log)
    allfit = read_mf_vmax('sdss',/bestfit)
    good = where(all.limit eq 1 and all.number gt mingal and $
      all.mass gt xrange[0] and all.mass lt xrange[1])
    oploterror, all.meanmass[good], all.phi[good], all.phierr_upper_stat[good], $
      psym=symcat(allpsym,thick=5), symsize=allsymsize, errthick=5, $
      color=im_color(allcolor), errcolor=im_color(allcolor), /hibar, nohat=nohat
    oploterror, all.meanmass[good], all.phi[good], all.phierr_lower_stat[good], psym=3, $
      symsize=allsymsize, errthick=5, color=im_color(allcolor), $
      errcolor=im_color(allcolor), /lobar, nohat=nohat

    qq = read_mf_vmax('sdss',/quiescent,/log)
    qqfit = read_mf_vmax('sdss',/quiescent,/bestfit)
    good = where(qq.limit eq 1 and qq.number gt mingal and $
      qq.mass gt xrange[0] and qq.mass lt xrange[1])
    oploterror, qq.meanmass[good], qq.phi[good], qq.phierr_upper_stat[good], $
      psym=symcat(qqpsym,thick=5), symsize=qqsymsize, errthick=5, $
      color=im_color(qqcolor), errcolor=im_color(qqcolor), /hibar, nohat=nohat
    oploterror, qq.meanmass[good], qq.phi[good], qq.phierr_lower_stat[good], psym=3, $
      symsize=qqsymsize, errthick=5, color=im_color(qqcolor), $
      errcolor=im_color(qqcolor), /lobar, nohat=nohat

    sf = read_mf_vmax('sdss',/active,/log)
    sffit = read_mf_vmax('sdss',/active,/bestfit)
    good = where(sf.limit eq 1 and sf.number gt mingal and $
      sf.mass gt xrange[0] and sf.mass lt xrange[1])
    oploterror, sf.meanmass[good], sf.phi[good], sf.phierr_upper_stat[good], $
      psym=symcat(sfpsym,thick=5), symsize=sfsymsize, errthick=5, $
      color=im_color(sfcolor), errcolor=im_color(sfcolor), /hibar, nohat=nohat
    oploterror, sf.meanmass[good], sf.phi[good], sf.phierr_lower_stat[good], psym=3, $
      symsize=sfsymsize, errthick=5, color=im_color(sfcolor), $
      errcolor=im_color(sfcolor), /lobar, nohat=nohat

; ratio of quiescent to 'all' galaxies    
    good = where(all.limit eq 1 and all.number gt mingal and $
      qq.limit eq 1 and qq.number gt mingal and $
      qq.mass gt xrange[0] and qq.mass lt xrange[1])
    meanmass = qq.meanmass[good]
    ratio = 10^(qq.phi[good]-all.phi[good])
    ratio_err = alog(10)*qq.phierr[good]*ratio

    yrange = [0,1.2]

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytitle='Quiescent Fraction', xtickinterval=1, ytickinterval=0.5
    oploterror, meanmass, ratio, meanmass*0+0.1/2, ratio_err, $
      psym=symcat(allpsym,thick=5), symsize=allsymsize, errthick=5
;   djs_oplot, !x.crange, 0.5*[1,1]
    
;   djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=5, ysty=5, $
;     yrange=yrange, xrange=xrange
;   polyfill, [meanmass,reverse(meanmass)], [ratio-ratio_err,reverse(ratio+ratio_err)], $
;     /data, color=im_color('grey70'), noclip=0, /fill
;   djs_oplot, meanmass, ratio, psym=symcat(allpsym,thick=5), symsize=1.0, $
;     color=im_color('grey20')
;   djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
;     yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
;     ytitle='Quiescent Fraction', xtickinterval=1, ytickinterval=0.5
    
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf

; ---------------------------------------------------------------------------
; Figure 4: sdss - compare with the literature - all galaxies only
    litpsym = 15 & litcolor = 'black' & litpsize = 1.1

    bell_color = 'firebrick' & bell_psym = 9
    baldry_color = 'dodger blue' & baldry_psym = 6
    li_color = 'forest green' & li_psym = 5
    cole_color = 'orange' & cole_psym = 4
    bern_color = 'tan' & bern_psym = 11

    cole = read_01cole()
    bell = read_03bell_smf()
    baldry = read_12baldry()
    li = read_09li()
    bern = read_10bernardi()
    
    xrange = [8.8,12.2]
;   xrange = [7.9,12.2]
    yrange = [-7.2,-1.3]
    residrange = 1.49*[-1,1]
    maxis1 = range(xrange[0]+0.15,13,100)
    maxis_bern = range(10.5,13,100)
    maxis_cole = range(9.3,12.0,100)

; single-panel plot showing just all the galaxies                                 
    psfile = paperpath+'mf_sdss_lit'+suffix
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.8, height=6.0, charsize=2.0
    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytitle=mfplot_phititle(), xtickinterval=1
    im_legend, 'All Galaxies', /right, /top, box=0, charsize=2

;   label = ['Cole+01','Bell+03','Li & White 09','Baldry+12','This paper']
    label = ['SDSS-GALEX (This paper)','Cole+01','Bell+03','Li & White 09','Baldry+12']
    im_legend, label, /left, /bottom, box=0, psym=[litpsym,cole_psym,bell_psym,li_psym,baldry_psym], $
      pspacing=1.9, charsize=2.0, color=[litcolor,cole_color,bell_color,li_color,baldry_color], $
      symthick=6, symsize=[1.7,1.9,1.7,1.7,1.7,1.9]
;   label = ['SDSS-GALEX (This paper)','Cole+01','Bell+03','Li & White 09','Bernardi+10','Baldry+12']
;   im_legend, label, /left, /bottom, box=0, psym=[litpsym,cole_psym,bell_psym,li_psym,bern_psym,baldry_psym], $
;     pspacing=1.9, charsize=2.0, color=[litcolor,cole_color,bell_color,li_color,bern_color,baldry_color], $
;     symthick=6, symsize=[1.7,1.9,1.7,1.7,1.7,1.9]

; overplot the data from the literature
    nohat = 0
    
    oploterror, cole.mass, cole.phi, cole.phierr, symsize=1.4, $
      psym=symcat(cole_psym,thick=3), color=im_color(cole_color), $
      errcolor=im_color(cole_color), nohat=nohat

    oploterror, bell.mass, bell.phi, bell.phierr, symsize=1.2, $
      psym=symcat(bell_psym,thick=3), color=im_color(bell_color), $
      errcolor=im_color(bell_color), nohat=nohat

    oploterror, baldry.mass, baldry.phi, baldry.phierr, symsize=1.2, $
      psym=symcat(baldry_psym,thick=3), color=im_color(baldry_color), $
      errcolor=im_color(baldry_color), nohat=nohat

    oploterror, li.mass, li.phi, li.phierr, symsize=1.2, $
      psym=symcat(li_psym,thick=3), color=im_color(li_color), $
      errcolor=im_color(li_color), nohat=nohat
    
;    oploterror, bern.mass, bern.phi, bern.phierr, symsize=1.2, $
;      psym=symcat(bern_psym,thick=3), color=im_color(bern_color), $
;      errcolor=im_color(bern_color), nohat=nohat
;;   mfoplot_lit, maxis_bern, /bernardi, /log
    
    mfdata = read_mf_vmax('sdss',/log)
    good = where(mfdata.limit eq 1,nthese)

    oploterror, mfdata.meanmass[good], mfdata.phi[good], mfdata.phierr_upper[good], $
      psym=symcat(litpsym,thick=8), symsize=litpsize, errthick=6, $
      color=im_color(litcolor,100+0), errcolor=im_color(litcolor,100), /hibar, nohat=nohat
    oploterror, mfdata.meanmass[good], mfdata.phi[good], mfdata.phierr_lower[good], $
      psym=3, symsize=litpsize, errthick=5, $
      color=im_color(litcolor,100+0), errcolor=im_color(litcolor,100), /lobar, nohat=nohat

    im_plotconfig, /psclose, psfile=psfile, pdf=pdf


; ---------------------------------------------------------------------------
; Figure 5: sdss - compare with the literature - qq/red, sf/blue
    xrange = [8.8,12.0]
    yrange = [-7.0,-1.5]
    residrange = 1.49*[-1,1]
    maxis1 = range(xrange[0]+0.15,13,100)
    maxis_ber = range(10.5+0.25,13,100)
    maxis_cole = range(9.3,12.0,100)
    maxis = range(xrange[0]+0.15,xrange[1],100)

    qqpsym = 14 & qqcolor = 'black' & qqpsize = 1.0
    sfpsym = 16 & sfcolor = 'black' & sfpsize = 1.0

    psfile = paperpath+'mf_sdss_qq_sf_lit'+suffix
    im_plotconfig, 6, pos, psfile=psfile, charsize=1.8, height=4.2*[1,1], $
      yspace=0.1, width=5

; star-forming/blue
    bell_color = 'dodger blue' & bell_psym = 6
    baldry_color = 'midnight blue' & baldry_psym = 9

    bell = read_03bell_smf(/blue)
    baldry = read_12baldry(/blue)

    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle='', xtickname=replicate(' ',10), $
      ytitle=mfplot_phititle(), xtickinterval=1
    im_legend, ['Star-Forming/Blue'], /right, /top, $
      box=0, charsize=1.7, margin=0, spacing=1.7
    label = ['SDSS-GALEX (This Paper)','Baldry+12','Bell+03']
    im_legend, label, /left, /bottom, box=0, psym=[sfpsym,baldry_psym,bell_psym], $
      pspacing=1.9, charsize=1.4, color=[sfcolor,baldry_color,bell_color], $
      symthick=6, symsize=[1.2,1.0,1.0];, margin=0

    mfdata = read_mf_vmax('sdss',/log,/active)
    good = where(mfdata.limit eq 1 and $
      mfdata.mass gt xrange[0] and mfdata.mass lt xrange[1])

    oploterror, mfdata.meanmass[good], mfdata.phi[good], mfdata.phierr_upper_stat[good], $
      psym=symcat(sfpsym,thick=8), symsize=sfpsize, errthick=6, $
      color=im_color(sfcolor,100+0), errcolor=im_color(sfcolor,100), /hibar, nohat=nohat
    oploterror, mfdata.meanmass[good], mfdata.phi[good], mfdata.phierr_lower_stat[good], $
      psym=3, symsize=sfpsize, errthick=5, $
      color=im_color(sfcolor,100+0), errcolor=im_color(sfcolor,100), /lobar, nohat=nohat

    oploterror, bell.mass, bell.phi, bell.phierr, $
      psym=symcat(bell_psym,thick=5), color=im_color(bell_color), $
      errcolor=im_color(bell_color), nohat=nohat
    oploterror, baldry.mass, baldry.phi, baldry.phierr, $
      psym=symcat(baldry_psym,thick=5), color=im_color(baldry_color), $
      errcolor=im_color(baldry_color), nohat=nohat

; quiescent/red
    bell_color = 'firebrick' & bell_psym = 6
    baldry_color = 'orange' & baldry_psym = 9

    bell = read_03bell_smf(/red)
    baldry = read_12baldry(/red)

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytitle=mfplot_phititle(), xtickinterval=1
    im_legend, ['Quiescent/Red'], /right, /top, $
      box=0, charsize=1.7, margin=0
    label = ['SDSS-GALEX (This Paper)','Baldry+12','Bell+03']
    im_legend, label, /left, /bottom, box=0, psym=[qqpsym,baldry_psym,bell_psym], $
      pspacing=1.9, charsize=1.4, color=[qqcolor,baldry_color,bell_color], $
      symthick=6, symsize=[1.3,1.0,1.0];, margin=0

    mfdata = read_mf_vmax('sdss',/log,/quiescent)
    good = where(mfdata.limit eq 1 and $
      mfdata.mass gt xrange[0] and mfdata.mass lt xrange[1])

    oploterror, mfdata.meanmass[good], mfdata.phi[good], mfdata.phierr_upper_stat[good], $
      psym=symcat(qqpsym,thick=8), symsize=qqpsize, errthick=6, $
      color=im_color(qqcolor,100+0), errcolor=im_color(qqcolor,100), /hibar, nohat=nohat
    oploterror, mfdata.meanmass[good], mfdata.phi[good], mfdata.phierr_lower_stat[good], $
      psym=3, symsize=qqpsize, errthick=5, $
      color=im_color(qqcolor,100+0), errcolor=im_color(qqcolor,100), /lobar, nohat=nohat

    oploterror, bell.mass, bell.phi, bell.phierr, $
      psym=symcat(bell_psym,thick=5), color=im_color(bell_color), $
      errcolor=im_color(bell_color), nohat=nohat
    oploterror, baldry.mass, baldry.phi, baldry.phierr, $
      psym=symcat(baldry_psym,thick=5), color=im_color(baldry_color), $
      errcolor=im_color(baldry_color), nohat=nohat

    im_plotconfig, /psclose, psfile=psfile, pdf=pdf

; --------------------------------------------------
; Figure 6: primus - fiducial supergrid, 'all' MF comparing the
; individual fields 
    fields = get_mf_fields()

    psfile = paperpath+'mf_all_byfield'+suffix
    plot_binned_mf, noevol=noevol, psfile=psfile, field=fields, supergrid=1 ; paper plot

    psfile = qapath+'mf_qq_byfield.ps'
    plot_binned_mf, noevol=noevol, psfile=psfile, field=fields, supergrid=1, /quiescent, /pdf ; qaplot

    psfile = qapath+'mf_sf_byfield.ps'
    plot_binned_mf, noevol=noevol, psfile=psfile, field=fields, supergrid=1, /active, /pdf ; qaplot

; ---------------------------------------------------------------------------
; Figure 8: sdss+primus - 'all' stellar mass function in seven
; redshift bins 
    xrange = [8.7,12.2]
    yrange = [-6.8,-1.01]

    zbins_sdss = mf_zbins(/sdss)
    zbins = mf_zbins(nzbins)
    maxis1 = range(xrange[0]+0.2,12.5,100)
    scolor = 'black' & sline = 5
    mingal = 0
    
    psfile = paperpath+'mf_all_sevenpanel'+suffix
    im_plotconfig, 18, pos, psfile=psfile, charsize=1.8, $
      xmargin=[1.1,0.3], width=2.4*[1,1,1,1], height=2.6*[1,1]

; SDSS
    mf_sdss = read_mf_vmax('sdss',noevol=noevol,/log)
    mffit_sdss = read_mf_vmax('sdss',/bestfit,noevol=noevol)

    plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle='', xtickname=replicate(' ',10), $
      ytitle='', xtickinterval=1, xminor=4
    legend, 'z='+strtrim(string(zbins_sdss.zlo,format='(F12.2)'),2)+$
      '-'+strtrim(string(zbins_sdss.zup,format='(F12.2)'),2), $
      /right, /top, box=0, margin=-0.1, charsize=1.5

    good_sdss = where(mf_sdss.limit eq 1 and mf_sdss.number ge mingal and mf_sdss.mass lt xrange[1])
    polyfill, [mf_sdss.meanmass[good_sdss],reverse(mf_sdss.meanmass[good_sdss])],$
      [mf_sdss.phi_lower_stat[good_sdss],reverse(mf_sdss.phi_upper_stat[good_sdss])], $
      /data, color=im_color('tan'), noclip=0, /fill
    djs_oplot, mf_sdss.meanmass[good_sdss], mf_sdss.phi[good_sdss], $
      psym=-symcat(15,thick=3), symsize=0.75, line=0, thick=4
;   djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=mffit_sdss)), $
;     line=sline, thick=4, color=im_color(scolor)

; PRIMUS    
    mf = read_mf_vmax(noevol=noevol,/avgfield,/log)
    for iz = 0, nzbins-1 do begin
       if (iz le 1) then begin
          xtickname = replicate(' ',10)
          xtitle = ''
       endif else begin
          delvarx, xtickname
          xtitle = mfplot_masstitle()
       endelse
       if (iz eq 3) then delvarx, ytickname else ytickname = replicate(' ',10)
;      if odd(iz) then delvarx, ytickname else ytickname = replicate(' ',10)
       
       plot, [0], [0], /nodata, /noerase, position=pos[*,iz+1], xsty=1, ysty=1, $
         yrange=yrange, xrange=xrange, ytickname=ytickname, $
         xtitle=xtitle, ytitle='', xtickinterval=1, xminor=4, xtickname=xtickname
       legend, 'z='+strtrim(string(zbins[iz].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[iz].zup,format='(F12.2)'),2), $
         /right, /top, box=0, margin=-0.1, charsize=1.5
;      djs_oplot, 10*[1,1], !y.crange, line=5
       
;      minmass = xrange[0]>10.0
       minmass = xrange[0]>(min(mf[iz].mass[where(mf[iz].limit)])-allmassbelow) ; 1 dex below
       good = where(mf[iz].number gt mingal and mf[iz].mass gt minmass)
;      good = where(mf[iz].limit eq 1 and mf[iz].number gt mingal)
       polyfill, [mf[iz].meanmass[good],reverse(mf[iz].meanmass[good])],$
         [mf[iz].phi_lower_stat[good],reverse(mf[iz].phi_upper_stat[good])], $
         /data, color=im_color('tan'), noclip=0, /fill
;      djs_oplot, mf[iz].meanmass[good], mf[iz].phi[good], $
;        psym=symcat(15,thick=3), symsize=0.75

       above = where(mf[iz].limit[good] eq 1 and mf[iz].mass[good] gt minmass)
       below = where(mf[iz].limit[good] eq 0 and mf[iz].mass[good] gt minmass,nbelow)
       djs_oplot, mf[iz].meanmass[good[above]], mf[iz].phi[good[above]], $
         psym=symcat(15,thick=3), symsize=0.75
       if nbelow gt 0 then begin
          djs_oplot, mf[iz].meanmass[good[below]], mf[iz].phi[good[below]], $
            psym=symcat(6,thick=3), symsize=0.75
;         plotsym, 2, 1.2, thick=4
;         djs_oplot, mf[iz].meanmass[good[below]], mf[iz].phi[good[below]], psym=8
       endif
       djs_oplot, mf_sdss.meanmass[good_sdss], mf_sdss.phi[good_sdss], line=0, thick=4
    endfor 

; ytitle    
    xyouts, pos[0,0]-0.06, pos[1,0], mfplot_phititle(), align=0.5, $
      /normal, orientation=90
;   xyouts, pos[0,0]-0.07, pos[1,4], mfplot_phititle(), align=0.5, $
;     /normal, orientation=90
    
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf

; ---------------------------------------------------------------------------
; Figure 10: sdss+primus - quiescent and SF stellar mass function in
; seven redshift bins 
    xrange = [8.7,12.2]
    yrange = [-5.4,-1.5]
;   yrange = [-7.0,-1.5]

    zbins_sdss = mf_zbins(/sdss)
    zbins = mf_zbins(nzbins)
    maxis1 = range(xrange[0]+0.2,12.5,100)
    scolor = 'black' & sline = 5

    qqcolor = 'firebrick'  & qqpolycolor = 'coral'
    qqfitcolor = 'dark red' & qqline = 0
    allfitcolor = 'black' & qqline = 0

    sfcolor = 'navy' & sfpolycolor = 'powder blue'
    sffitcolor = 'navy'     & sfline = 5

    mingal = 2 ; 3

; for the paper, calculate how the quenching mass varies with redshift
    zz = [zbins_sdss.zbin,zbins.zbin]
    mquench = alog10(10^10.4*(1+zz)^1.5)
    niceprint, zz, mquench
;    mquench = fltarr(7)
;    qq = read_mf_vmax('sdss',/quiescent,/silent,/log)
;    sf = read_mf_vmax('sdss',/active,/silent,/log)
;    mquench[0] = interpol(qq.mass,qq.phi-sf.phi,0.0)
;    ploterror, qq.mass, qq.phi-sf.phi, qq.phierr, xr=[10,11.3], $
;      psym=6, /trad, yrange=[-0.3,0.4], xsty=3, ysty=3
;    qq = read_mf_vmax(/avgfield,/quiescent,/silent,/log)
;    sf = read_mf_vmax(/avgfield,/active,/silent,/log)
;    for iz = 0, nzbins-1 do begin
;       oploterror, qq[iz].mass, qq[iz].phi-sf[iz].phi, qq[iz].phierr, line=iz
;       cc = get_kbrd(1)
;;      mquench[iz+1] = interpol(qq[iz].mass,qq[iz].phi-sf[iz].phi,0.0)
;    endfor
    
; now make the plot    
    psfile = paperpath+'mf_qq_sf_sevenpanel'+suffix
    im_plotconfig, 18, pos, psfile=psfile, charsize=1.7, $
      xmargin=[1.1,0.3], width=2.4*[1,1,1,1], height=2.6*[1,1]

; ----------
; SDSS
    plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle='', xtickname=replicate(' ',10), $
      ytitle='', xtickinterval=1, xminor=4, ytickinterval=1
    legend, 'z='+strtrim(string(zbins_sdss.zlo,format='(F12.2)'),2)+$
      '-'+strtrim(string(zbins_sdss.zup,format='(F12.2)'),2), $
      /right, /top, box=0, charsize=1.4, margin=0;-0.3

    qq = read_mf_vmax('sdss',noevol=noevol,/log,/quiescent,/silent)
    sf = read_mf_vmax('sdss',noevol=noevol,/log,/active,/silent)
    qqfit = read_mf_vmax('sdss',/bestfit,noevol=noevol,/quiescent,/silent)
    sffit = read_mf_vmax('sdss',/bestfit,noevol=noevol,/active,/silent)
    allfit = read_mf_vmax('sdss',/bestfit,noevol=noevol,/silent)

; star-forming
    gd = where(sf.limit eq 1 and sf.number ge mingal and sf.mass lt xrange[1],ngd)
;   gd = where(sf.limit eq 1 and sf.number gt mingal)
    mm = [sf.meanmass[gd[0]]-sf[0].binsize/2,sf.meanmass[gd],sf.meanmass[gd[ngd-1]]+sf[0].binsize/2]
    polyfill, [mm,reverse(mm)],$
      [[sf.phi_lower_stat[gd[0]],sf.phi_lower_stat[gd],sf.phi_lower_stat[gd[ngd-1]]],$
      reverse([sf.phi_upper_stat[gd[0]],sf.phi_upper_stat[gd],sf.phi_upper_stat[gd[ngd-1]]])], $
      /data, color=im_color(sfpolycolor), noclip=0, /fill
    djs_oplot, sf.meanmass[gd], sf.phi[gd], line=sfline, thick=4, $
      psym=-symcat(16,thick=4), symsize=0.75, color=im_color(sfcolor)

; quiescent    
    gd = where(qq.limit eq 1 and qq.number gt mingal and qq.mass lt xrange[1])
;   gd = where(qq.limit eq 1 and qq.number gt mingal); and qq.number gt 3)
    polyfill, [qq.meanmass[gd],reverse(qq.meanmass[gd])],$
      [qq.phi_lower_stat[gd],reverse(qq.phi_upper_stat[gd])], $
      /data, color=im_color(qqpolycolor), noclip=0, /line_fill, orientation=45, spacing=0.1
    polyfill, [qq.meanmass[gd],reverse(qq.meanmass[gd])],$
      [qq.phi_lower_stat[gd],reverse(qq.phi_upper_stat[gd])], $
      /data, color=im_color(qqpolycolor), noclip=0, /line_fill, orientation=135, spacing=0.1
;   polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], thick=4, $
;     /data, color=im_color(qqpolycolor), noclip=0, /line_fill, orientation=45, spacing=0.1
;   polyfill, [qq.meanmass[gd],reverse(qq.meanmass[gd])],$
;     [qq.phi_lower_stat[gd],reverse(qq.phi_upper_stat[gd])], thick=4, $
;     /data, color=im_color(qqpolycolor), noclip=0, /line_fill, orientation=135, spacing=0.1
    djs_oplot, qq.meanmass[gd], qq.phi[gd], line=qqline, thick=4, $
      psym=-symcat(14,thick=4), symsize=0.9, color=im_color(qqcolor)

;   djs_oplot, mquench[0]*[1,1], !y.crange, line=5
    
; schechter fits    
;   djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=allfit)), $
;     line=allline, thick=5, color=im_color(allfitcolor)
;   djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=qqfit)), $
;     line=qqline, thick=5, color=im_color(qqfitcolor)
;   djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=sffit)), $
;     line=sfline, thick=5, color=im_color(sffitcolor)

; ----------
; PRIMUS    
    qq = read_mf_vmax(noevol=noevol,/avgfield,/log,/quiescent,/silent)
    sf = read_mf_vmax(noevol=noevol,/avgfield,/log,/active,/silent)
    for iz = 0, nzbins-1 do begin
       if (iz le 1) then begin
;      if (iz le 3) then begin
          xtickname = replicate(' ',10)
          xtitle = ''
       endif else begin
          delvarx, xtickname
          xtitle = mfplot_masstitle()
       endelse
       if (iz eq 3) then delvarx, ytickname else ytickname = replicate(' ',10)
;      if odd(iz) then delvarx, ytickname else ytickname = replicate(' ',10)
       
       plot, [0], [0], /nodata, /noerase, position=pos[*,iz+1], xsty=1, ysty=1, $
         yrange=yrange, xrange=xrange, ytickname=ytickname, $
         xtitle=xtitle, ytitle='', xtickinterval=1, xminor=4, xtickname=xtickname, $
         ytickinterval=1
       legend, 'z='+strtrim(string(zbins[iz].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[iz].zup,format='(F12.2)'),2), $
         /right, /top, box=0, charsize=1.4, margin=0;-0.3

; star-forming 
       minmass = xrange[0]>(min(sf[iz].mass[where(sf[iz].limit)])-sfmassbelow) ; 1 dex below
;      gd = where(sf[iz].limit eq 1 and sf[iz].number gt mingal)
       gd = where(sf[iz].number gt mingal and sf[iz].mass gt minmass)
       polyfill, [sf[iz].meanmass[gd],reverse(sf[iz].meanmass[gd])],$
         [sf[iz].phi_lower_stat[gd],reverse(sf[iz].phi_upper_stat[gd])], $
         /data, color=im_color(sfpolycolor), noclip=0, /fill

       above = where(sf[iz].limit[gd] eq 1 and sf[iz].mass[gd] gt minmass)
       below = where(sf[iz].limit[gd] eq 0 and sf[iz].mass[gd] gt minmass,nbelow)
       djs_oplot, sf[iz].meanmass[gd[above]], sf[iz].phi[gd[above]], $
         psym=symcat(16,thick=3), symsize=0.75, color=im_color(sfcolor)
       if nbelow gt 0 then begin
          djs_oplot, sf[iz].meanmass[gd[below]], sf[iz].phi[gd[below]], $
            psym=symcat(9,thick=3), symsize=0.75, color=im_color(sfcolor)
       endif
;      djs_oplot, sf[iz].meanmass[gd], sf[iz].phi[gd], $
;        psym=symcat(16,thick=3), symsize=0.75, color=im_color(sfcolor)

; quiescent       
       minmass = xrange[0]>(min(qq[iz].mass[where(qq[iz].limit)])-qqmassbelow) ; 1 dex below
;      gd = where(qq[iz].limit eq 1 and qq[iz].number gt mingal)
       gd = where(qq[iz].number gt mingal and qq[iz].mass gt minmass)
       polyfill, [qq[iz].meanmass[gd],reverse(qq[iz].meanmass[gd])],$
         [qq[iz].phi_lower_stat[gd],reverse(qq[iz].phi_upper_stat[gd])], $
         /data, color=im_color(qqpolycolor), noclip=0, /line_fill, orientation=45, spacing=0.1
       polyfill, [qq[iz].meanmass[gd],reverse(qq[iz].meanmass[gd])],$
         [qq[iz].phi_lower_stat[gd],reverse(qq[iz].phi_upper_stat[gd])], $
         /data, color=im_color(qqpolycolor), noclip=0, /line_fill, orientation=135, spacing=0.1
;      polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], thick=4, $
;        /data, color=im_color(qqpolycolor), noclip=0, /line_fill, orientation=45, spacing=0.1
;      polyfill, [qq[iz].meanmass[gd],reverse(qq[iz].meanmass[gd])],$
;        [qq[iz].phi_lower_stat[gd],reverse(qq[iz].phi_upper_stat[gd])], thick=4, $
;        /data, color=im_color(qqpolycolor), noclip=0, /line_fill, orientation=135, spacing=0.1

       above = where(qq[iz].limit[gd] eq 1 and qq[iz].mass[gd] gt minmass)
       below = where(qq[iz].limit[gd] eq 0 and qq[iz].mass[gd] gt minmass,nbelow)
       djs_oplot, qq[iz].meanmass[gd[above]], qq[iz].phi[gd[above]], $
         psym=symcat(14,thick=3), symsize=0.9, color=im_color(qqcolor)
       if nbelow gt 0 then begin
          djs_oplot, qq[iz].meanmass[gd[below]], qq[iz].phi[gd[below]], $
            psym=symcat(4,thick=3), symsize=0.9, color=im_color(qqcolor)
       endif
;      djs_oplot, qq[iz].meanmass[gd], qq[iz].phi[gd], $
;        psym=symcat(14,thick=3), symsize=0.9, color=im_color(qqcolor)

; SDSS fits       
;      djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=allfit)), $
;        line=allline, thick=5, color=im_color(allfitcolor)
;      djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=qqfit)), $
;        line=qqline, thick=7, color=im_color(qqfitcolor)
;      djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=sffit)), $
;        line=sfline, thick=7, color=im_color(sffitcolor)

; SDSS MFs       
       sdsssf = read_mf_vmax('sdss',noevol=noevol,/log,/active,/silent)
       sdssgd = where(sdsssf.limit eq 1 and sdsssf.number gt mingal)
       djs_oplot, sdsssf.meanmass[sdssgd], sdsssf.phi[sdssgd], $
         color=im_color(sfcolor), line=sfline, thick=4

       sdssqq = read_mf_vmax('sdss',noevol=noevol,/log,/quiescent,/silent)
       sdssgd = where(sdssqq.limit eq 1 and sdssqq.number gt mingal)
       djs_oplot, sdssqq.meanmass[sdssgd], sdssqq.phi[sdssgd], $
         color=im_color(qqcolor), line=qqline, thick=4
;      djs_oplot, sdssqq.meanmass[sdssgd], sdssqq.phi[sdssgd], $
;        psym=symcat(4,thick=4), symsize=0.75, color=im_color(qqcolor)

;      djs_oplot, mquench[iz+1]*[1,1], !y.crange, line=5
;      djs_oplot, 10.5*[1,1], !y.crange, line=5
    endfor 

; ytitle    
    xyouts, pos[0,0]-0.07, pos[1,0], mfplot_phititle(), align=0.5, $
      /normal, orientation=90
;   xyouts, pos[0,0]-0.07, pos[1,4], mfplot_phititle(), align=0.5, $
;     /normal, orientation=90
    
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf

; ---------------------------------------------------------------------------
; Figure 11: sdss+primus - 'qq' and 'sf' stellar mass functions in a
; 2-panel plot 
    xrange = massrange1 
    yrange = phirange1
    maxis1 = range(xrange[0]+0.2,12.5,100)
    mingal = 3
    
    sdss_zbins = mf_zbins(/sdss)
    zbins = mf_zbins(nzbins)

;   zlabel = ['z=',replicate('   ',7)]+strtrim(string([sdss_zbins.zlo,zbins.zlo],format='(F12.2)'),2)+$
;     '-'+strtrim(string([sdss_zbins.zup,zbins.zup],format='(F12.2)'),2)
    zlabel = ['<z>=',replicate('      ',7)]+strtrim(string([sdss_zbins.zbin,zbins.zbin],format='(F12.2)'),2)
;   zlabel = ' '+strtrim(string([sdss_zbins.zbin,zbins.zbin],format='(F12.2)'),2)

    yoff = 0.15
    wid = 0.2 & hei = 0.1
    x1 = xrange[0]+0.2 & y1 = -4.0
    
    psfile = paperpath+'mf_qq_sf_onepanel'+suffix
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.8, $
      height=4.3, xspace=0.1

; sf; table 13, 22, and 24 also looks good
    sdss_col = 'black'
;   loadct, 1, /silent
    loadct, 24, /silent, file=getenv('IMPRO_DIR')+'/etc/fsc_brewer.tbl'
    col = reverse(findgen(nzbins)/(nzbins-1)*(250-170)+170)
;   col = reverse(findgen(nzbins)/(nzbins-1)*(230-130)+130)

    plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytitle=mfplot_phititle(), xtickinterval=1, ytickinterval=1, $
      color=im_color('black',255)
    im_legend, 'Star-Forming', /right, /top, box=0, charsize=1.6, margin=0

    polyfill, [x1,x1+wid,x1+wid,x1],[y1,y1,y1+hei,y1+hei], $
      /data, color=im_color(sdss_col,255), noclip=0, /fill
    xyouts, x1+wid*1.15, y1+hei*0.2, zlabel[0], align=0.0, /data, $
      charsize=1.3, color=im_color('black',255)
    for iz = 0, nzbins-1 do begin
       polyfill, [x1,x1+wid,x1+wid,x1],[y1,y1,y1+hei,y1+hei]-yoff*(iz+1), $
         /data, color=im_color(col[iz]), noclip=0, /fill
       xyouts, x1+wid*1.15, y1+hei*0.2-yoff*(iz+1), zlabel[iz+1], align=0.0, $
         /data, charsize=1.3, color=im_color('black',255)
    endfor
;   xyouts, x1+wid*1.5, y1+hei*2.2, '<z>', align=0.0, /data, charsize=1.3, color=im_color('black',255)
      
; PRIMUS    
    mfdata = read_mf_vmax(noevol=noevol,/avgfield,/log,/active)
    for iz = 0, nzbins-1 do begin
       good = where(mfdata[iz].limit eq 1 and mfdata[iz].number gt mingal)
       phimass = mfdata[iz].meanmass[good]
       phimin = mfdata[iz].phi_lower_stat[good]
       phimax = mfdata[iz].phi_upper_stat[good]
       polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
         /data, color=im_color(col[iz]), noclip=0, /fill
;      if odd(iz) then begin
;         polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
;           /data, color=im_color(col[iz]), noclip=0, /line_fill, $
;           spacing=0.1, orientation=45
;      endif else begin
;         polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
;           /data, color=im_color(col[iz]), noclip=0, /fill
;      endelse

; draw some outlines
;      outline_mf, phimass, phimin, phimax, line=iz
    endfor 

; SDSS
    mfdata = read_mf_vmax('sdss',noevol=noevol,/log,/active)
    mffit_sdss = read_mf_vmax('sdss',/bestfit,noevol=noevol,/active)

    good = where(mfdata.limit eq 1 and mfdata.number gt mingal)
    phimass = mfdata.meanmass[good]
    phimin = mfdata.phi_lower_stat[good]
    phimax = mfdata.phi_upper_stat[good]
    polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
      /data, color=im_color(sdss_col,255), noclip=0, /fill

; qq ; 10, 15
    sdss_col = 'black'
    loadct, 10, /silent, file=getenv('IMPRO_DIR')+'/etc/fsc_brewer.tbl'
    col = reverse(findgen(nzbins)/(nzbins-1)*(210-60)+60)
;   col = reverse(findgen(nzbins)/(nzbins-1)*(240-120)+120)
;   col = ['orange1','orange2','orange3','orange4','tomato1','tomato2']

    plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
      ytitle='', ytickname=replicate(' ',10), xtickinterval=1, ytickinterval=1, $
      color=im_color('black',255)
    im_legend, 'Quiescent', /right, /top, box=0, charsize=1.6, margin=0

    polyfill, [x1,x1+wid,x1+wid,x1],[y1,y1,y1+hei,y1+hei], $
      /data, color=im_color(sdss_col,255), noclip=0, /fill
    xyouts, x1+wid*1.15, y1+hei*0.2, zlabel[0], align=0.0, /data, $
      charsize=1.3, color=im_color('black',255)
    for iz = 0, nzbins-1 do begin
       polyfill, [x1,x1+wid,x1+wid,x1],[y1,y1,y1+hei,y1+hei]-yoff*(iz+1), $
         /data, color=im_color(col[iz],255), noclip=0, /fill
       xyouts, x1+wid*1.15, y1+hei*0.2-yoff*(iz+1), zlabel[iz+1], align=0.0, $
         /data, charsize=1.3, color=im_color('black',255)
    endfor
      
; PRIMUS    
    mfdata = read_mf_vmax(noevol=noevol,/avgfield,/log,/quiescent)
    for iz = 0, nzbins-1 do begin
       good = where(mfdata[iz].limit eq 1 and mfdata[iz].number gt mingal)
       phimass = mfdata[iz].meanmass[good]
       phimin = mfdata[iz].phi_lower_stat[good]
       phimax = mfdata[iz].phi_upper_stat[good]
       polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
         /data, color=im_color(col[iz],255), noclip=0, /fill
    endfor 

; SDSS
    mfdata = read_mf_vmax('sdss',noevol=noevol,/log,/quiescent)
    mffit_sdss = read_mf_vmax('sdss',/bestfit,noevol=noevol,/quiescent)

    good = where(mfdata.limit eq 1 and mfdata.number gt mingal)
    phimass = mfdata.meanmass[good]
    phimin = mfdata.phi_lower_stat[good]
    phimax = mfdata.phi_upper_stat[good]
    polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
      /data, color=im_color(sdss_col,255), noclip=0, /fill

    im_plotconfig, /psclose, psfile=psfile, pdf=pdf

; ---------------------------------------------------------------------------
; Figre 13: primus - compare the sf and qq MF with the literature 
    zbins = mf_zbins(nzbins,/lit)

    xrange = [8.7,12.1]
    yrange = [-5.5,-1.3]
    maxis = range(xrange[0],xrange[1],100)

    nohat = 1
    mingal = 3
    primuscolor = 'tan'
    
; assemble the literature preferences here
    drory_color = 'purple'         & drory_psym = 4
    pozzetti_color = 'dodger blue' & pozzetti_psym = 6
    borch_color = 'firebrick'      & borch_psym = 5
    ilbert_color = 'orange'  & ilbert_psym = 7
    
    psfile = paperpath+'mf_qq_sf_lit'+suffix
    im_plotconfig, 18, pos, psfile=psfile, charsize=1.7, height=2.8*[1,1], yspace=0.1

; #########################    
; top row = quiescent galaxies
    mfdata = read_mf_vmax(noevol=noevol,/avgfield,/log,/lit,/quiescent)
    for iz = 0, nzbins-1 do begin
       if (iz eq 0) then delvarx, ytickname else ytickname = replicate(' ',10)
       
       plot, [0], [0], /nodata, noerase=(iz gt 0), position=pos[*,iz], xsty=1, ysty=1, $
         yrange=yrange, xrange=xrange, xtickname=replicate(' ',10), ytickname=ytickname, $
         xtitle='', ytitle='', xtickinterval=1
       legend, 'z='+strtrim(string(zbins[iz].zlo,format='(F12.1)'),2)+$
         '-'+strtrim(string(zbins[iz].zup,format='(F12.1)'),2), $
         /right, /top, box=0, margin=0, charsize=1.3


; overplot the PRIMUS data
       good = where(mfdata[iz].limit eq 1 and mfdata[iz].number gt mingal)
       polyfill, [mfdata[iz].meanmass[good],reverse(mfdata[iz].meanmass[good])],$
         [mfdata[iz].phi_lower_stat[good],reverse(mfdata[iz].phi_upper_stat[good])], $
         /data, color=im_color(primuscolor), noclip=0, /fill

;      good = where(mfdata[iz].limit eq 1)
;      oploterror, mfdata[iz].meanmass[good], mfdata[iz].phi[good], $
;        mfdata[iz].phierr_upper_stat[good], $
;        psym=symcat(15,thick=3), symsize=0.7, nohat=nohat, /hibar
;      oploterror, mfdata[iz].meanmass[good], mfdata[iz].phi[good], $
;        mfdata[iz].phierr_lower_stat[good], $
;        psym=3, nohat=nohat, /lobar

       case iz of
          0: begin
             drory = read_09drory(1,/quiescent)
             oploterror, drory.mass, drory.phi, drory.phierr, $
               psym=symcat(drory_psym,thick=6), color=im_color(drory_color), $
               errcolor=im_color(drory_color), nohat=nohat

             pozzetti = read_10pozzetti(1,/type1)
             oploterror, pozzetti.mass, pozzetti.phi, pozzetti.phierr, $
               psym=symcat(pozzetti_psym,thick=4), color=im_color(pozzetti_color), $
               errcolor=im_color(pozzetti_color), nohat=nohat

             borch = read_06borch(1,/red)
             oploterror, borch.mass, borch.phi, borch.phierr, $
               psym=symcat(borch_psym,thick=4), color=im_color(borch_color), $
               errcolor=im_color(borch_color), nohat=nohat

             ilbert = read_10ilbert(1,/quiescent)
             oploterror, ilbert.mass, ilbert.phi, ilbert.phierr, $
               psym=symcat(ilbert_psym,thick=4), color=im_color(ilbert_color), $
               errcolor=im_color(ilbert_color), nohat=nohat
;            ilbert = get_mf_literature(/ilbert,/nolog,/quiescent)
;            djs_oplot, maxis, alog10(mf_schechter(maxis,schechter=ilbert[0])), $
;              line=ilbert[0].line, thick=6, color=im_color(ilbert[0].color1)

             legend, 'Quiescent', /left, /bottom, box=0, charsize=1.3, margin=0
             
          end
          1: begin
             drory = read_09drory(2,/quiescent)
             oploterror, drory.mass, drory.phi, drory.phierr, $
               psym=symcat(drory_psym,thick=5), color=im_color(drory_color), $
               errcolor=im_color(drory_color), nohat=nohat

             pozzetti = read_10pozzetti(2,/type1)
             oploterror, pozzetti.mass, pozzetti.phi, pozzetti.phierr, $
               psym=symcat(pozzetti_psym,thick=5), color=im_color(pozzetti_color), $
               errcolor=im_color(pozzetti_color), nohat=nohat

             borch = read_06borch(2,/red)
             oploterror, borch.mass, borch.phi, borch.phierr, $
               psym=symcat(borch_psym,thick=4), color=im_color(borch_color), $
               errcolor=im_color(borch_color), nohat=nohat

             ilbert = read_10ilbert(2,/quiescent)
             oploterror, ilbert.mass, ilbert.phi, ilbert.phierr, $
               psym=symcat(ilbert_psym,thick=4), color=im_color(ilbert_color), $
               errcolor=im_color(ilbert_color), nohat=nohat

          end
          2: begin
             drory = read_09drory(3,/quiescent)
             oploterror, drory.mass, drory.phi, drory.phierr, $
               psym=symcat(drory_psym,thick=5), color=im_color(drory_color), $
               errcolor=im_color(drory_color), nohat=nohat

             pozzetti = read_10pozzetti(3,/type1)
             oploterror, pozzetti.mass, pozzetti.phi, pozzetti.phierr, $
               psym=symcat(pozzetti_psym,thick=5), color=im_color(pozzetti_color), $
               errcolor=im_color(pozzetti_color), nohat=nohat

             borch = read_06borch(3,/red)
             oploterror, borch.mass, borch.phi, borch.phierr, $
               psym=symcat(borch_psym,thick=4), color=im_color(borch_color), $
               errcolor=im_color(borch_color), nohat=nohat

             ilbert = read_10ilbert(3,/quiescent)
             oploterror, ilbert.mass, ilbert.phi, ilbert.phierr, $
               psym=symcat(ilbert_psym,thick=4), color=im_color(ilbert_color), $
               errcolor=im_color(ilbert_color), nohat=nohat
             
          end
          3: begin
             drory = read_09drory(4,/quiescent)
             oploterror, drory.mass, drory.phi, drory.phierr, $
               psym=symcat(drory_psym,thick=5), color=im_color(drory_color), $
               errcolor=im_color(drory_color), nohat=nohat

             pozzetti = read_10pozzetti(4,/type1)
             oploterror, pozzetti.mass, pozzetti.phi, pozzetti.phierr, $
               psym=symcat(pozzetti_psym,thick=5), color=im_color(pozzetti_color), $
               errcolor=im_color(pozzetti_color), nohat=nohat

             borch = read_06borch(4,/red)
             oploterror, borch.mass, borch.phi, borch.phierr, $
               psym=symcat(borch_psym,thick=4), color=im_color(borch_color), $
               errcolor=im_color(borch_color), nohat=nohat

             ilbert = read_10ilbert(4,/quiescent)
             oploterror, ilbert.mass, ilbert.phi, ilbert.phierr, $
               psym=symcat(ilbert_psym,thick=4), color=im_color(ilbert_color), $
               errcolor=im_color(ilbert_color), nohat=nohat

          end
          else: 
       endcase
    endfor 

; #########################    
; bottom row = sf galaxies
    mfdata = read_mf_vmax(noevol=noevol,/avgfield,/log,/lit,/active)
    for iz = 0, nzbins-1 do begin
       if (iz eq 0) then delvarx, ytickname else ytickname = replicate(' ',10)

       plot, [0], [0], /nodata, /noerase, position=pos[*,iz+4], xsty=1, ysty=1, $
         yrange=yrange, xrange=xrange, ytickname=ytickname, $
         xtitle='', ytitle='', xtickinterval=1

; overplot the PRIMUS data
       good = where(mfdata[iz].limit eq 1 and mfdata[iz].number gt mingal)
       polyfill, [mfdata[iz].meanmass[good],reverse(mfdata[iz].meanmass[good])],$
         [mfdata[iz].phi_lower_stat[good],reverse(mfdata[iz].phi_upper_stat[good])], $
         /data, color=im_color(primuscolor), noclip=0, /fill

;      good = where(mfdata[iz].limit eq 1)
;      oploterror, mfdata[iz].meanmass[good], mfdata[iz].phi[good], $
;        mfdata[iz].phierr_upper_stat[good], $
;        psym=symcat(15,thick=3), symsize=0.7, nohat=nohat, /hibar
;      oploterror, mfdata[iz].meanmass[good], mfdata[iz].phi[good], $
;        mfdata[iz].phierr_lower_stat[good], $
;        psym=3, nohat=nohat, /lobar
       
       case iz of
          0: begin
             drory = read_09drory(1,/sf)
             oploterror, drory.mass, drory.phi, drory.phierr, $
               psym=symcat(drory_psym,thick=6), color=im_color(drory_color), $
               errcolor=im_color(drory_color), nohat=nohat

             pozzetti = read_10pozzetti(1,/type234)
             oploterror, pozzetti.mass, pozzetti.phi, pozzetti.phierr, $
               psym=symcat(pozzetti_psym,thick=4), color=im_color(pozzetti_color), $
               errcolor=im_color(pozzetti_color), nohat=nohat

             borch = read_06borch(1,/blue)
             oploterror, borch.mass, borch.phi, borch.phierr, $
               psym=symcat(borch_psym,thick=4), color=im_color(borch_color), $
               errcolor=im_color(borch_color), nohat=nohat

             ilbert = read_10ilbert(1,/sf)
             oploterror, ilbert.mass, ilbert.phi, ilbert.phierr, $
               psym=symcat(ilbert_psym,thick=4), color=im_color(ilbert_color), $
               errcolor=im_color(ilbert_color), nohat=nohat

             legend, 'Star-Forming', /left, /bottom, box=0, charsize=1.3, margin=0

          end
          1: begin
             drory = read_09drory(2,/sf)
             oploterror, drory.mass, drory.phi, drory.phierr, $
               psym=symcat(drory_psym,thick=5), color=im_color(drory_color), $
               errcolor=im_color(drory_color), nohat=nohat

             pozzetti = read_10pozzetti(2,/type234)
             oploterror, pozzetti.mass, pozzetti.phi, pozzetti.phierr, $
               psym=symcat(pozzetti_psym,thick=5), color=im_color(pozzetti_color), $
               errcolor=im_color(pozzetti_color), nohat=nohat

             borch = read_06borch(2,/blue)
             oploterror, borch.mass, borch.phi, borch.phierr, $
               psym=symcat(borch_psym,thick=4), color=im_color(borch_color), $
               errcolor=im_color(borch_color), nohat=nohat

             ilbert = read_10ilbert(2,/sf)
             oploterror, ilbert.mass, ilbert.phi, ilbert.phierr, $
               psym=symcat(ilbert_psym,thick=4), color=im_color(ilbert_color), $
               errcolor=im_color(ilbert_color), nohat=nohat

          end
          2: begin
             drory = read_09drory(3,/sf)
             oploterror, drory.mass, drory.phi, drory.phierr, $
               psym=symcat(drory_psym,thick=5), color=im_color(drory_color), $
               errcolor=im_color(drory_color), nohat=nohat

             pozzetti = read_10pozzetti(3,/type234)
             oploterror, pozzetti.mass, pozzetti.phi, pozzetti.phierr, $
               psym=symcat(pozzetti_psym,thick=5), color=im_color(pozzetti_color), $
               errcolor=im_color(pozzetti_color), nohat=nohat

             borch = read_06borch(3,/blue)
             oploterror, borch.mass, borch.phi, borch.phierr, $
               psym=symcat(borch_psym,thick=4), color=im_color(borch_color), $
               errcolor=im_color(borch_color), nohat=nohat

             ilbert = read_10ilbert(3,/sf)
             oploterror, ilbert.mass, ilbert.phi, ilbert.phierr, $
               psym=symcat(ilbert_psym,thick=4), color=im_color(ilbert_color), $
               errcolor=im_color(ilbert_color), nohat=nohat

          end
          3: begin
             drory = read_09drory(4,/sf)
             oploterror, drory.mass, drory.phi, drory.phierr, $
               psym=symcat(drory_psym,thick=5), color=im_color(drory_color), $
               errcolor=im_color(drory_color), nohat=nohat

             pozzetti = read_10pozzetti(4,/type234)
             oploterror, pozzetti.mass, pozzetti.phi, pozzetti.phierr, $
               psym=symcat(pozzetti_psym,thick=5), color=im_color(pozzetti_color), $
               errcolor=im_color(pozzetti_color), nohat=nohat

             borch = read_06borch(4,/blue)
             oploterror, borch.mass, borch.phi, borch.phierr, $
               psym=symcat(borch_psym,thick=4), color=im_color(borch_color), $
               errcolor=im_color(borch_color), nohat=nohat

             ilbert = read_10ilbert(4,/sf)
             oploterror, ilbert.mass, ilbert.phi, ilbert.phierr, $
               psym=symcat(ilbert_psym,thick=4), color=im_color(ilbert_color), $
               errcolor=im_color(ilbert_color), nohat=nohat

             label = ['PRIMUS','Borch+06','Drory+09','Pozzetti+10','Ilbert+10']
             color = ['tan',borch_color,drory_color,pozzetti_color,ilbert_color]
             psym = [15,borch_psym,drory_psym,pozzetti_psym,ilbert_psym]
             im_legend, label, /left, /bottom, box=0, color=color, $
               psym=psym, charsize=1.1, margin=0, symthick=5
          end
          else: 
       endcase
    endfor 

; render the titles    
    xyouts, pos[0,5], pos[1,5]-0.1, mfplot_masstitle(), align=0.5, charsize=1.7, /norm
    xyouts, pos[0,7], pos[1,7]-0.1, mfplot_masstitle(), align=0.5, charsize=1.7, /norm
    xyouts, pos[0,0]-0.07, pos[1,0], mfplot_phititle(), align=0.5, orientation=90, charsize=1.7, /norm

    im_plotconfig, psfile=psfile, /psclose, pdf=pdf

; ---------------------------------------------------------------------------
; QAplot: primus - compare the 'all' MF with the literature
    supergrid = 1
    zbins = mf_zbins(nzbins,/lit)
    mfdata = read_mf_vmax(noevol=noevol,/avgfield,/log,/lit)
    mffit_sdss = read_mf_vmax('sdss',supergrid=supergrid,/bestfit,noevol=noevol)

    xrange = [8.7,12.2]
    yrange = [-5.9,-1.3]
    maxis = range(xrange[0]+0.15,xrange[1],100)

    nohat = 1
    mingal = 0
    primuscolor = 'tan'
    
; assemble the literature preferences here
    drory_color = 'purple'         & drory_psym = 4
    pozzetti_color = 'dodger blue' & pozzetti_psym = 6
    borch_color = 'firebrick'      & borch_psym = 5
    perez_color = 'forest green'   & perez_psym = 9
    ilbert_color = 'orange'  & ilbert_psym = 7
    
    psfile = qapath+'mf_all_lit.ps'
;   psfile = qapath+'mf_all_lit.ps'
    im_plotconfig, 5, pos, psfile=psfile, charsize=1.7, $
      width=[3.4,3.4], height=3.0*[1,1], xmargin=[1.3,0.4]

    for iz = 0, nzbins-1 do begin
       if (iz ge 2) then delvarx, xtickname else $
         xtickname = replicate(' ',10)
       if odd(iz) then ytickname = replicate(' ',10) else $
         delvarx, ytickname
       
       plot, [0], [0], /nodata, noerase=(iz gt 0), position=pos[*,iz], xsty=1, ysty=1, $
         yrange=yrange, xrange=xrange, xtickname=xtickname, ytickname=ytickname, $
         xtitle='', ytitle=ytitle, xtickinterval=1
       legend, 'z='+strtrim(string(zbins[iz].zlo,format='(F12.1)'),2)+$
         '-'+strtrim(string(zbins[iz].zup,format='(F12.1)'),2), $
         /right, /top, box=0, margin=0, charsize=1.5

; PRIMUS       
       good = where(mfdata[iz].limit eq 1 and mfdata[iz].number gt mingal)
       polyfill, [mfdata[iz].meanmass[good],reverse(mfdata[iz].meanmass[good])],$
         [mfdata[iz].phi_lower_stat[good],reverse(mfdata[iz].phi_upper_stat[good])], $
         /data, color=im_color(primuscolor), noclip=0, /fill
       
;      djs_oplot, maxis, alog10(mf_schechter(maxis,schechter=mffit_sdss)), line=0, thick=4
       
; plot the literature data       
       case iz of
          0: begin
             drory = read_09drory(1)
             oploterror, drory.mass, drory.phi, drory.phierr, $
               psym=symcat(drory_psym,thick=6), color=im_color(drory_color), $
               errcolor=im_color(drory_color), nohat=nohat

             pozzetti = read_10pozzetti(1)
             oploterror, pozzetti.mass, pozzetti.phi, pozzetti.phierr, $
               psym=symcat(pozzetti_psym,thick=4), color=im_color(pozzetti_color), $
               errcolor=im_color(pozzetti_color), nohat=nohat

             borch = read_06borch(1)
             oploterror, borch.mass, borch.phi, borch.phierr, $
               psym=symcat(borch_psym,thick=4), color=im_color(borch_color), $
               errcolor=im_color(borch_color), nohat=nohat

             perez = read_08perez(1)
             oploterror, perez.mass, perez.phi, perez.phierr, $
               psym=symcat(perez_psym,thick=4), color=im_color(perez_color), $
               errcolor=im_color(perez_color), nohat=nohat

             ilbert = read_10ilbert(1)
             oploterror, ilbert.mass, ilbert.phi, ilbert.phierr, $
               psym=symcat(ilbert_psym,thick=4), color=im_color(ilbert_color), $
               errcolor=im_color(ilbert_color), nohat=nohat
          end
          1: begin
             drory = read_09drory(2)
             oploterror, drory.mass, drory.phi, drory.phierr, $
               psym=symcat(drory_psym,thick=5), color=im_color(drory_color), $
               errcolor=im_color(drory_color), nohat=nohat

             pozzetti = read_10pozzetti(2)
             oploterror, pozzetti.mass, pozzetti.phi, pozzetti.phierr, $
               psym=symcat(pozzetti_psym,thick=5), color=im_color(pozzetti_color), $
               errcolor=im_color(pozzetti_color), nohat=nohat

             borch = read_06borch(2)
             oploterror, borch.mass, borch.phi, borch.phierr, $
               psym=symcat(borch_psym,thick=4), color=im_color(borch_color), $
               errcolor=im_color(borch_color), nohat=nohat

             perez = read_08perez(2)
             oploterror, perez.mass, perez.phi, perez.phierr, $
               psym=symcat(perez_psym,thick=4), color=im_color(perez_color), $
               errcolor=im_color(perez_color), nohat=nohat

             ilbert = read_10ilbert(2)
             oploterror, ilbert.mass, ilbert.phi, ilbert.phierr, $
               psym=symcat(ilbert_psym,thick=4), color=im_color(ilbert_color), $
               errcolor=im_color(ilbert_color), nohat=nohat
          end
          2: begin
             drory = read_09drory(3)
             oploterror, drory.mass, drory.phi, drory.phierr, $
               psym=symcat(drory_psym,thick=5), color=im_color(drory_color), $
               errcolor=im_color(drory_color), nohat=nohat

             pozzetti = read_10pozzetti(3)
             oploterror, pozzetti.mass, pozzetti.phi, pozzetti.phierr, $
               psym=symcat(pozzetti_psym,thick=5), color=im_color(pozzetti_color), $
               errcolor=im_color(pozzetti_color), nohat=nohat

             borch = read_06borch(3)
             oploterror, borch.mass, borch.phi, borch.phierr, $
               psym=symcat(borch_psym,thick=4), color=im_color(borch_color), $
               errcolor=im_color(borch_color), nohat=nohat

             perez = read_08perez(3)
             oploterror, perez.mass, perez.phi, perez.phierr, $
               psym=symcat(perez_psym,thick=4), color=im_color(perez_color), $
               errcolor=im_color(perez_color), nohat=nohat

             ilbert = read_10ilbert(3)
             oploterror, ilbert.mass, ilbert.phi, ilbert.phierr, $
               psym=symcat(ilbert_psym,thick=4), color=im_color(ilbert_color), $
               errcolor=im_color(ilbert_color), nohat=nohat

             label = ['PRIMUS','Borch+06','Perez-Gonzalez+08','Drory+09','Pozzetti+10','Ilbert+10']
             color = [primuscolor,borch_color,perez_color,drory_color,pozzetti_color,ilbert_color]
             psym = [15,borch_psym,perez_psym,drory_psym,pozzetti_psym,ilbert_psym]
             im_legend, label, /left, /bottom, box=0, color=color, $
               psym=psym, charsize=1.3, margin=0, symthick=5
          end
          3: begin
             drory = read_09drory(4)
             oploterror, drory.mass, drory.phi, drory.phierr, $
               psym=symcat(drory_psym,thick=5), color=im_color(drory_color), $
               errcolor=im_color(drory_color), nohat=nohat

             pozzetti = read_10pozzetti(4)
             oploterror, pozzetti.mass, pozzetti.phi, pozzetti.phierr, $
               psym=symcat(pozzetti_psym,thick=5), color=im_color(pozzetti_color), $
               errcolor=im_color(pozzetti_color), nohat=nohat

             borch = read_06borch(4)
             oploterror, borch.mass, borch.phi, borch.phierr, $
               psym=symcat(borch_psym,thick=4), color=im_color(borch_color), $
               errcolor=im_color(borch_color), nohat=nohat

             perez = read_08perez(4)
             oploterror, perez.mass, perez.phi, perez.phierr, $
               psym=symcat(perez_psym,thick=4), color=im_color(perez_color), $
               errcolor=im_color(perez_color), nohat=nohat

             ilbert = read_10ilbert(4)
             oploterror, ilbert.mass, ilbert.phi, ilbert.phierr, $
               psym=symcat(ilbert_psym,thick=4), color=im_color(ilbert_color), $
               errcolor=im_color(ilbert_color), nohat=nohat
          end
          else: 
       endcase

; overplot the PRIMUS data
;       good = where(mfdata[iz].limit eq 1)
;       oploterror, mfdata[iz].meanmass[good], mfdata[iz].phi[good], $
;         mfdata[iz].phierr_upper_stat[good], $
;         psym=symcat(15,thick=3), symsize=0.7, nohat=nohat, /hibar
;       oploterror, mfdata[iz].meanmass[good], mfdata[iz].phi[good], $
;         mfdata[iz].phierr_lower_stat[good], $
;         psym=3, nohat=nohat, /lobar
    endfor 
    xyouts, pos[2,2], pos[1,2]-0.1, mfplot_masstitle(), align=0.5, $
      charsize=1.7, /norm
    xyouts, pos[0,0]-0.1, pos[1,0], mfplot_phititle(), align=0.5, $
      orientation=90, charsize=1.7, /norm

    im_plotconfig, psfile=psfile, /psclose, /pdf

return
end


;; ---------------------------------------------------------------------------
;; sdss+primus - qq fraction
;    psfile = paperpath+'qq_fraction'+suffix
;    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
;      width=6.8, height=5.0, charsize=2.0
;    
;    sdss_zbins = mf_zbins(/sdss)
;    zbins = mf_zbins(nzbins)
;    zlabel = ['z=',replicate('   ',7)]+strtrim(string([sdss_zbins.zlo,zbins.zlo],format='(F12.2)'),2)+$
;      '-'+strtrim(string([sdss_zbins.zup,zbins.zup],format='(F12.2)'),2)
;
;    mingal = 15
;    xrange = [8.9,12.1]
;    yrange = [0,2.1]
;    loadct, 3, /silent
;    col = reverse(findgen(nzbins)/(nzbins-1)*(230-50)+50)
;
;    sdss_order = 3
;    order = [2,2,2,2,2,2]
;    
;    sdss_psym = 16 & sdss_col = 'forest green'
;    psym = [15,6,14,13,5,2]
;    
;    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
;      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
;      ytitle='Quiescent Fraction', xtickinterval=1, ytickinterval=0.2
;
;; legend    
;    yoff = 0.06
;    wid = 0.2 & hei = 0.05
;    x1 = xrange[0]+0.15 & y1 = yrange[1]*0.9
;    polyfill, [x1,x1+wid,x1+wid,x1],[y1,y1,y1+hei,y1+hei], $
;      /data, color=im_color(sdss_col,255), noclip=0, /fill
;    xyouts, x1+wid*1.15, y1+hei*0.2, zlabel[0], align=0.0, /data, charsize=1.3
;    for iz = 0, nzbins-1 do begin
;       polyfill, [x1,x1+wid,x1+wid,x1],[y1,y1,y1+hei,y1+hei]-yoff*(iz+1), $
;         /data, color=im_color(col[iz]), noclip=0, /fill
;       xyouts, x1+wid*1.15, y1+hei*0.2-yoff*(iz+1), zlabel[iz+1], $
;         align=0.0, /data, charsize=1.3
;    endfor
;    
;; sdss    
;    sdssall = read_mf_vmax('sdss',noevol=noevol,/log)
;    sdssqq = read_mf_vmax('sdss',noevol=noevol,/log,/quiescent)
;    sdssgood = where(sdssqq.limit eq 1 and sdssqq.number gt mingal)
;          
;; primus    
;    qq = read_mf_vmax(noevol=noevol,/avgfield,/log,/quiescent)
;;   for iz = 1, nzbins-1 do begin ; exclude the lowest redshift bin
;    for iz = 0, nzbins-1 do begin
;       good = where(qq[iz].limit eq 1 and qq[iz].number gt mingal and $
;         qq[iz].mass ge min(sdssqq.mass[sdssgood]))
;       phimass = qq[iz].mass[good]
;       sdssphi = interpol(sdssqq.phi[sdssgood],sdssqq.mass[sdssgood],phimass)
;
;       phimin = 10^(qq[iz].phi_lower_stat[good]-sdssphi)
;       phimax = 10^(qq[iz].phi_upper_stat[good]-sdssphi)
;       polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
;         /data, color=im_color(col[iz]), noclip=0, /fill
;;      djs_oplot, all[iz].meanmass[good], 10^(qq[iz].phi[good]-all[iz].phi[good]), $
;;        psym=-symcat(psym[iz],thick=5), color=im_color(color[iz])
;    endfor
;
;;; now fit polynomials to the ratio in each redshift bin and overlay
;;; those
;;    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
;;      yrange=[0,1.1], xrange=[8.9,12.1], xtitle=mfplot_masstitle(), $
;;      ytitle='Quiescent Fraction', xtickinterval=1, ytickinterval=0.2
;;
;;; sdss    
;;;   all = read_mf_vmax('sdss',noevol=noevol)
;;    allfit = read_mf_vmax('sdss',noevol=noevol,/bestfit)
;;    qqfit = read_mf_vmax('sdss',noevol=noevol,/bestfit,/quiescent)
;;    qq = read_mf_vmax('sdss',noevol=noevol,/quiescent)
;;    good = where(qq.limit eq 1 and qq.number gt mingal)
;;;   good = where(all.limit eq 1 and all.number gt mingal and $
;;;     qq.limit eq 1 and qq.number gt mingal)
;;;   rcoeff = poly_fit(qq.meanmass[good],qq.phi[good]/all.phi[good],sdss_order)
;;;     measure_errors=qq.phierr[good]/all.phi[good])
;;;   rcoeff = poly_fit(qq.meanmass[good],qq.phi[good]/all.phi[good],$
;;;     measure_errors=qq.phierr[good]/all.phi[good],sdss_order)
;;;   mm = range(min(qq.meanmass[good]),max(qq.meanmass[good]),50)
;;;   djs_oploterr, qq.meanmass[good], qq.phi[good]/all.phi[good], $
;;    allmodel = mf_schechter(qq.meanmass[good],schechter=allfit)
;;    qqmodel = mf_schechter(qq.meanmass[good],schechter=qqfit)
;;    djs_oploterr, qq.meanmass[good], qq.phi[good]/allmodel, $
;;      yerr=qq.phierr[good]/allmodel, psym=6, color=im_color(sdss_col)
;;;   djs_oplot, qq.meanmass[good], qqmodel/allmodel, line=0, thick=4
;;;   djs_oplot, mm, poly(mm,rcoeff), line=0, color=im_color(sdss_col)
;;    
;;; primus    
;;;   all = read_mf_vmax(noevol=noevol,/avgfield)
;;    allfit = read_mf_vmax(noevol=noevol,/avgfield,/bestfit)
;;    qqfit = read_mf_vmax(noevol=noevol,/avgfield,/bestfit,/quiescent)
;;    qq = read_mf_vmax(noevol=noevol,/avgfield,/quiescent)
;;;   for iz = 5, 5 do begin      ; exclude the lowest redshift bin
;;    for iz = 1, nzbins-1 do begin ; exclude the lowest redshift bin
;;;   for iz = 0, nzbins-1 do begin
;;       good = where(qq[iz].limit eq 1 and qq[iz].number gt mingal)
;;;      good = where(qq[iz].limit eq 1 and all[iz].limit eq 1 and $
;;;        qq[iz].number gt mingal and all[iz].number gt mingal)
;;;      rcoeff = poly_fit(qq[iz].meanmass[good],qq[iz].phi[good]/all[iz].phi[good],order[iz]);,$
;;;        measure_errors=qq[iz].phierr[good]/all[iz].phi[good])
;;;      mm = range(min(qq[iz].meanmass[good]),max(qq[iz].meanmass[good]),50)
;;       allmodel = mf_schechter(qq[iz].meanmass[good],schechter=allfit[iz])
;;       qqmodel = mf_schechter(qq[iz].meanmass[good],schechter=qqfit[iz])
;;       djs_oploterr, qq[iz].meanmass[good], qq[iz].phi[good]/allmodel, $
;;         yerr=qq[iz].phierr[good]/allmodel, psym=6, color=im_color(col[iz])
;;;       djs_oplot, qq[iz].meanmass[good], qqmodel/allmodel, line=0, thick=4, $
;; ;        color=im_color(col[iz]);, psym=symcat(16)
;;;      djs_oplot, mm, poly(mm,rcoeff), line=iz-1, color=im_color(col[iz])
;;    endfor
;
;    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
;
;stop    
;    
;
;; ---------------------------------------------------------------------------
;; sdss+primus - 'all' stellar mass function on one single panel 
;    xrange = massrange1 
;    yrange = phirange1
;    maxis1 = range(xrange[0]+0.15,13,100)
;
;    sdss_zbins = mf_zbins(/sdss)
;    zbins = mf_zbins(nzbins)
;
;    maxis1 = range(xrange[0]+0.2,12.5,100)
;    sdss_col = 'orange'
;    loadct, 0, /silent
;    col = reverse(findgen(nzbins)/(nzbins-1)*(220-70)+70)
;
;    mingal = 3
;    
;    psfile = paperpath+'mf_all_onepanel'+suffix
;    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8, height=6.0
;
;    plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
;      ytitle=mfplot_phititle(), xtickinterval=1, ytickinterval=1
;    im_legend, 'All Galaxies', /right, /top, box=0, charsize=1.8
;    zlabel = ['z=',replicate('   ',7)]+string([sdss_zbins.zbin,zbins.zbin],format='(F4.2)')
;;   zlabel = ['z=',replicate('   ',7)]+strtrim(string([sdss_zbins.zlo,zbins.zlo],format='(F12.2)'),2)+$
;;     '-'+strtrim(string([sdss_zbins.zup,zbins.zup],format='(F12.2)'),2)
;    yoff = 0.15
;    wid = 0.2 & hei = 0.1
;    x1 = xrange[0]+0.15 & y1 = -4.1
;    polyfill, [x1,x1+wid,x1+wid,x1],[y1,y1,y1+hei,y1+hei], $
;      /data, color=im_color(sdss_col,255), noclip=0, /fill
;    xyouts, x1+wid*1.15, y1+hei*0.2, zlabel[0], align=0.0, /data, charsize=1.3
;    for iz = 0, nzbins-1 do begin
;       polyfill, [x1,x1+wid,x1+wid,x1],[y1,y1,y1+hei,y1+hei]-yoff*(iz+1), $
;         /data, color=im_color(col[iz]), noclip=0, /fill
;       xyouts, x1+wid*1.15, y1+hei*0.2-yoff*(iz+1), zlabel[iz+1], align=0.0, /data, charsize=1.3
;    endfor
;      
;; PRIMUS    
;    mfdata = read_mf_vmax(noevol=noevol,/avgfield,/log)
;    mffit = read_mf_vmax(noevol=noevol,/avgfield,/bestfit)
;    for iz = 0, nzbins-1 do begin
;       good = where(mfdata[iz].limit eq 1 and mfdata[iz].number gt mingal)
;       phimass = mfdata[iz].meanmass[good]
;       phimin = mfdata[iz].phi_lower_stat[good]
;       phimax = mfdata[iz].phi_upper_stat[good]
;       polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
;         /data, color=im_color(col[iz]), noclip=0, /fill
;;      djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=mffit[iz]))
;    endfor 
;
;; SDSS
;    mfdata_sdss = read_mf_vmax('sdss',noevol=noevol,/log)
;    mffit_sdss = read_mf_vmax('sdss',/bestfit,noevol=noevol)
;
;    good = where(mfdata_sdss.limit eq 1 and mfdata_sdss.number gt mingal)
;    phimass = mfdata_sdss.meanmass[good]
;    phimin = mfdata_sdss.phi_lower_stat[good]
;    phimax = mfdata_sdss.phi_upper_stat[good]
;    polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
;      /data, color=im_color(sdss_col,255), noclip=0, /fill
;
;;; add an inset showing the Schechter fits
;;    inset_pos = [pos[0]+0.25,pos[1]+0.23,pos[0]+0.57,pos[1]+0.5]
;;    djs_plot, [0], [0], /nodata, /noerase, position=inset_pos, $
;;      xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
;;      charsize=1.3, xtickinterval=1, ytickinterval=1, $
;;      xtitle=mfplot_masstitle(/simple), ytitle=mfplot_phititle(/simple)
;    
;    good = where(mfdata_sdss.limit eq 1)
;    maxis1 = range(min(mfdata_sdss.meanmass[good]),max(mfdata_sdss.meanmass[good]),100)
;    djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=mffit_sdss)), $
;      color=im_color(sdss_col,255), thick=6
;    for iz = 0, nzbins-1 do begin
;       good = where(mfdata[iz].limit eq 1 and mfdata[iz].number gt mingal)
;       maxis1 = range(min(mfdata_sdss.meanmass[good]),max(mfdata_sdss.meanmass[good]),100)
;       djs_oplot, maxis1, alog10(mf_schechter(maxis1,schechter=mffit[iz])), $
;         color=im_color(col[iz]), thick=6
;    endfor 
;    
;    im_plotconfig, /psclose, psfile=psfile, pdf=pdf

;; ---------------------------------------------------------------------------
;; sdss+primus - qq fraction
;    psfile = paperpath+'qq_fraction'+suffix
;    im_plotconfig, 6, pos, psfile=psfile, xmargin=[1.3,0.4], $
;      width=6.8, height=[4.5,2.5], charsize=2.0
;    device, decomposed=0
;    
;    sdss_zbins = mf_zbins(/sdss)
;    zbins = mf_zbins(nzbins)
;;   zlabel = 'z='+string([zbins_sdss.zbin,zbins.zbin],format='(F4.2)')
;;   zlabel = ['z=',replicate('   ',6)]+$
;;     string([zbins_sdss.zbin,zbins.zbin],format='(F4.2)')
;    zlabel = ['z=',replicate('   ',7)]+strtrim(string([sdss_zbins.zlo,zbins.zlo],format='(F12.2)'),2)+$
;      '-'+strtrim(string([sdss_zbins.zup,zbins.zup],format='(F12.2)'),2)
;
;    mingal = 15
;    xrange = [8.9,12.1]
;    yrange = [0,1.1]
;    loadct, 3, /silent
;    col = reverse(findgen(nzbins)/(nzbins-1)*(230-50)+50)
;;   col = reverse(findgen(nzbins)/(nzbins-1)*(230-50)+50)
;
;    sdss_order = 3
;    order = [2,2,2,2,2,2]
;;   order = [2,1,1,1,1,1]
;    
;    sdss_psym = 16 & sdss_col = 'forest green'
;    psym = [15,6,14,13,5,2]
;;   symsize = [1,1,1.2
;;   col = ['dodger blue','spring green4','orange','firebrick','purple','powder blue']
;    
;    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
;      yrange=yrange, xrange=xrange, xtitle='', $
;      ytitle='Quiescent Fraction', xtickname=replicate(' ',10), $
;      xtickinterval=1, ytickinterval=0.2
;
;    yoff = 0.06
;    wid = 0.2 & hei = 0.05
;    x1 = xrange[0]+0.15 & y1 = yrange[1]*0.9
;    polyfill, [x1,x1+wid,x1+wid,x1],[y1,y1,y1+hei,y1+hei], $
;      /data, color=im_color(sdss_col,255), noclip=0, /fill
;    xyouts, x1+wid*1.15, y1+hei*0.2, zlabel[0], align=0.0, /data, charsize=1.3
;    for iz = 0, nzbins-1 do begin
;       polyfill, [x1,x1+wid,x1+wid,x1],[y1,y1,y1+hei,y1+hei]-yoff*(iz+1), $
;         /data, color=im_color(col[iz]), noclip=0, /fill
;       xyouts, x1+wid*1.15, y1+hei*0.2-yoff*(iz+1), zlabel[iz+1], $
;         align=0.0, /data, charsize=1.3
;    endfor
;    
;; sdss    
;    all = read_mf_vmax('sdss',noevol=noevol,/log)
;    qq = read_mf_vmax('sdss',noevol=noevol,/log,/quiescent)
;    good = where(all.limit eq 1 and all.number gt mingal and $
;      qq.limit eq 1 and qq.number gt mingal)
;    phimass = all.meanmass[good]
;    phimin = 10^(qq.phi_lower_stat[good]-all.phi[good])
;    phimax = 10^(qq.phi_upper_stat[good]-all.phi[good])
;    polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
;      /data, color=im_color(sdss_col,255), noclip=0, /fill
;;   djs_oplot, qq.meanmass[good], 10^(qq.phi[good]-all.phi[good]), $
;;     psym=-symcat(sdss_psym,thick=5), color=im_color(sdss_col)
;
;; primus    
;    all = read_mf_vmax(noevol=noevol,/avgfield,/log)
;    qq = read_mf_vmax(noevol=noevol,/avgfield,/log,/quiescent)
;    for iz = 1, nzbins-1 do begin ; exclude the lowest redshift bin
;;   for iz = 0, nzbins-1 do begin
;       good = where(qq[iz].limit eq 1 and all[iz].limit eq 1 and $
;         qq[iz].number gt mingal and all[iz].number gt mingal)
;       phimass = all[iz].meanmass[good]
;       phimin = 10^(qq[iz].phi_lower_stat[good]-all[iz].phi[good])
;       phimax = 10^(qq[iz].phi_upper_stat[good]-all[iz].phi[good])
;       polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
;         /data, color=im_color(col[iz]), noclip=0, /fill
;;      djs_oplot, all[iz].meanmass[good], 10^(qq[iz].phi[good]-all[iz].phi[good]), $
;;        psym=-symcat(psym[iz],thick=5), color=im_color(color[iz])
;    endfor
;
;; now fit polynomials to the ratio in each redshift bin and overlay
;; those
;    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
;      yrange=[0,1.1], xrange=[8.9,12.1], xtitle=mfplot_masstitle(), $
;      ytitle='Quiescent Fraction', xtickinterval=1, ytickinterval=0.2
;
;; sdss    
;;   all = read_mf_vmax('sdss',noevol=noevol)
;    allfit = read_mf_vmax('sdss',noevol=noevol,/bestfit)
;    qqfit = read_mf_vmax('sdss',noevol=noevol,/bestfit,/quiescent)
;    qq = read_mf_vmax('sdss',noevol=noevol,/quiescent)
;    good = where(qq.limit eq 1 and qq.number gt mingal)
;;   good = where(all.limit eq 1 and all.number gt mingal and $
;;     qq.limit eq 1 and qq.number gt mingal)
;;   rcoeff = poly_fit(qq.meanmass[good],qq.phi[good]/all.phi[good],sdss_order)
;;     measure_errors=qq.phierr[good]/all.phi[good])
;;   rcoeff = poly_fit(qq.meanmass[good],qq.phi[good]/all.phi[good],$
;;     measure_errors=qq.phierr[good]/all.phi[good],sdss_order)
;;   mm = range(min(qq.meanmass[good]),max(qq.meanmass[good]),50)
;;   djs_oploterr, qq.meanmass[good], qq.phi[good]/all.phi[good], $
;    allmodel = mf_schechter(qq.meanmass[good],schechter=allfit)
;    qqmodel = mf_schechter(qq.meanmass[good],schechter=qqfit)
;    djs_oploterr, qq.meanmass[good], qq.phi[good]/allmodel, $
;      yerr=qq.phierr[good]/allmodel, psym=6, color=im_color(sdss_col)
;;   djs_oplot, qq.meanmass[good], qqmodel/allmodel, line=0, thick=4
;;   djs_oplot, mm, poly(mm,rcoeff), line=0, color=im_color(sdss_col)
;    
;; primus    
;;   all = read_mf_vmax(noevol=noevol,/avgfield)
;    allfit = read_mf_vmax(noevol=noevol,/avgfield,/bestfit)
;    qqfit = read_mf_vmax(noevol=noevol,/avgfield,/bestfit,/quiescent)
;    qq = read_mf_vmax(noevol=noevol,/avgfield,/quiescent)
;;   for iz = 5, 5 do begin      ; exclude the lowest redshift bin
;    for iz = 1, nzbins-1 do begin ; exclude the lowest redshift bin
;;   for iz = 0, nzbins-1 do begin
;       good = where(qq[iz].limit eq 1 and qq[iz].number gt mingal)
;;      good = where(qq[iz].limit eq 1 and all[iz].limit eq 1 and $
;;        qq[iz].number gt mingal and all[iz].number gt mingal)
;;      rcoeff = poly_fit(qq[iz].meanmass[good],qq[iz].phi[good]/all[iz].phi[good],order[iz]);,$
;;        measure_errors=qq[iz].phierr[good]/all[iz].phi[good])
;;      mm = range(min(qq[iz].meanmass[good]),max(qq[iz].meanmass[good]),50)
;       allmodel = mf_schechter(qq[iz].meanmass[good],schechter=allfit[iz])
;       qqmodel = mf_schechter(qq[iz].meanmass[good],schechter=qqfit[iz])
;       djs_oploterr, qq[iz].meanmass[good], qq[iz].phi[good]/allmodel, $
;         yerr=qq[iz].phierr[good]/allmodel, psym=6, color=im_color(col[iz])
;;       djs_oplot, qq[iz].meanmass[good], qqmodel/allmodel, line=0, thick=4, $
; ;        color=im_color(col[iz]);, psym=symcat(16)
;;      djs_oplot, mm, poly(mm,rcoeff), line=iz-1, color=im_color(col[iz])
;    endfor
;
;    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
;


;; add an inset showing the M*-alpha covariance
;    scale = 100
;    nmonte = 500
;    inset_pos = [pos[0]+0.15,pos[1]+0.15,pos[0]+0.45,pos[1]+0.4]
;    djs_plot, [0], [0], /nodata, /noerase, position=inset_pos, $
;      xsty=1, ysty=1, xrange=scale*0.025*[-1,1], yrange=scale*0.025*[-1,1], $
;      xtitle='10^{-2} \Delta'+'log (M_{*}/M'+sunsymbol()+')', $
;      ytitle='10^{-2} \Delta\alpha', charsize=1.4, $
;      xtickinterval=1, ytickinterval=1
;    for jj = 0, n_elements(subsample)-1 do begin
;       mffit = read_mf_vmax('sdss',quiescent=(subsample[jj] eq 'quiescent'),$
;         active=(subsample[jj] eq 'active'),/bestfit)
;       covar = mffit.covar[1:2,1:2]*sqrt(mffit.chi2_dof)
;
;       rand = mrandomn(seed,covar,nmonte)
;       rand[*,0] = rand[*,0];+mffit.logmstar
;       rand[*,1] = rand[*,1];+mffit.alpha
;       ell1 = covar2ellipse(covar,nsigma=1.0)
;       ell3 = covar2ellipse(covar,nsigma=3.0)
;       tvellipse, scale*ell1.major, scale*ell1.minor, 0.0, 0.0, ell1.angle, $
;         thick=7, color=im_color(fitcolor[jj],101), /data, line=mfline[jj]
;;      tvellipse, scale*ell3.major, scale*ell3.minor, 0.0, 0.0, ell1.angle, $
;;        line=5, thick=4, color=im_color(fitcolor[jj],101), /data
;    endfor

;; ---------------------------------------------------------------------------
;; sdss - all, quiescent, and active samples
;    psfile = paperpath+'mf_sdss'+suffix
;    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
;      width=6.8, height=6.0, charsize=2.0
;
;    xrange = [8.8,12.2]
;    yrange = [-7.0,-1.5]
;    maxis1 = range(xrange[0]+0.15,12.5,100)
;    mingal = 0
;
;    sfpsym = 16 & sfcolor = 'dodger blue' & sffitcolor = 'black' & sfline = 5 & sfsymsize = 1.2
;    qqpsym = 14 & qqcolor = 'firebrick' & qqfitcolor = 'black' & qqline = 3 & qqsymsize = 1.4
;    allpsym = 15 & allcolor = 'black' & allfitcolor = 'black' & allline = 0 & allsymsize = 1.2
;    
;    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
;      ytitle=mfplot_phititle(), xtickinterval=1
;    im_legend, ['All','Quiescent','Star-Forming'], /left, /bottom, box=0, $
;      color=[allcolor,qqcolor,sfcolor], psym=[allpsym,qqpsym,sfpsym], $
;      symsize=[allsymsize,qqsymsize,sfsymsize]*1.7, $
;      symthick=8, spacing=2.4, charsize=2.1
;    
;; calculate the "quenching mass" for the paper
;    qq = alog10(mf_schechter(maxis1,schechter=read_mf_vmax('sdss',/quie,/best)))
;    aa = alog10(mf_schechter(maxis1,schechter=read_mf_vmax('sdss',/act,/best)) )
;    splog, 'quenching mass ', interpol(maxis1,qq-aa,0.0)
;
;; integrate the MF    
;    totrho = read_mf_vmax('sdss',/rho)
;    sfrho = read_mf_vmax('sdss',/active,/rho)
;    qqrho = read_mf_vmax('sdss',/quiescent,/rho)
;    splog, 'total, sf, qq, sf_fration, qq_fraction ', 10^totrho.rho, 10^sfrho.rho, $
;      10^qqrho.rho, 10^sfrho.rho/10^totrho.rho, 10^qqrho.rho/10^totrho.rho
;    
;; now make the plot    
;    all = read_mf_vmax('sdss',/log)
;    allfit = read_mf_vmax('sdss',/bestfit)
;    good = where(all.limit eq 1 and all.number gt mingal)
;    oploterror, all.meanmass[good], all.phi[good], all.phierr_upper[good], $
;      psym=symcat(allpsym,thick=5), symsize=allsymsize, errthick=5, $
;      color=im_color(allcolor), errcolor=im_color(allcolor), /hibar, /nohat
;    oploterror, all.meanmass[good], all.phi[good], all.phierr_lower[good], psym=3, $
;      symsize=allsymsize, errthick=5, color=im_color(allcolor), $
;      errcolor=im_color(allcolor), /lobar, /nohat
;
;    qq = read_mf_vmax('sdss',/quiescent,/log)
;    qqfit = read_mf_vmax('sdss',/quiescent,/bestfit)
;    good = where(qq.limit eq 1 and qq.number gt mingal)
;    oploterror, qq.meanmass[good], qq.phi[good], qq.phierr_upper[good], $
;      psym=symcat(qqpsym,thick=5), symsize=qqsymsize, errthick=5, $
;      color=im_color(qqcolor), errcolor=im_color(qqcolor), /hibar, /nohat
;    oploterror, qq.meanmass[good], qq.phi[good], qq.phierr_lower[good], psym=3, $
;      symsize=qqsymsize, errthick=5, color=im_color(qqcolor), $
;      errcolor=im_color(qqcolor), /lobar, /nohat
;
;    sf = read_mf_vmax('sdss',/active,/log)
;    sffit = read_mf_vmax('sdss',/active,/bestfit)
;    good = where(sf.limit eq 1 and sf.number gt mingal)
;    oploterror, sf.meanmass[good], sf.phi[good], sf.phierr_upper[good], $
;      psym=symcat(sfpsym,thick=5), symsize=sfsymsize, errthick=5, $
;      color=im_color(sfcolor), errcolor=im_color(sfcolor), /hibar, /nohat
;    oploterror, sf.meanmass[good], sf.phi[good], sf.phierr_lower[good], psym=3, $
;      symsize=sfsymsize, errthick=5, color=im_color(sfcolor), $
;      errcolor=im_color(sfcolor), /lobar, /nohat
;
;    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
;


;; ---------------------------------------------------------------------------
;; sdss+primus - 'sf' stellar mass function on one single panel 
;    xrange = massrange1 
;    yrange = phirange1
;
;    sdss_zbins = mf_zbins(/sdss)
;    zbins = mf_zbins(nzbins)
;
;    maxis1 = range(xrange[0]+0.2,12.5,100)
;    sdss_col = 'yellow2'
;    loadct, 1, /silent
;    col = reverse(findgen(nzbins)/(nzbins-1)*(230-130)+130)
;
;    mingal = 3
;    
;    psfile = paperpath+'mf_sf_onepanel'+suffix
;    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8, height=6.0
;
;    plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
;      ytitle=mfplot_phititle(), xtickinterval=1, ytickinterval=1
;    im_legend, 'Star-Forming Galaxies', /right, /top, box=0, charsize=1.8
;    zlabel = ['z=',replicate('   ',7)]+strtrim(string([sdss_zbins.zlo,zbins.zlo],format='(F12.2)'),2)+$
;      '-'+strtrim(string([sdss_zbins.zup,zbins.zup],format='(F12.2)'),2)
;
;    yoff = 0.15
;    wid = 0.2 & hei = 0.1
;    x1 = xrange[0]+0.15 & y1 = -4.1
;    polyfill, [x1,x1+wid,x1+wid,x1],[y1,y1,y1+hei,y1+hei], $
;      /data, color=im_color(sdss_col,255), noclip=0, /fill
;    xyouts, x1+wid*1.15, y1+hei*0.2, zlabel[0], align=0.0, /data, charsize=1.3
;    for iz = 0, nzbins-1 do begin
;       polyfill, [x1,x1+wid,x1+wid,x1],[y1,y1,y1+hei,y1+hei]-yoff*(iz+1), $
;         /data, color=im_color(col[iz]), noclip=0, /fill
;       xyouts, x1+wid*1.15, y1+hei*0.2-yoff*(iz+1), zlabel[iz+1], align=0.0, /data, charsize=1.3
;    endfor
;      
;; PRIMUS    
;    mfdata = read_mf_vmax(noevol=noevol,/avgfield,/log,/active)
;    for iz = 0, nzbins-1 do begin
;       good = where(mfdata[iz].limit eq 1 and mfdata[iz].number gt mingal)
;       phimass = mfdata[iz].meanmass[good]
;       phimin = mfdata[iz].phi_lower_stat[good]
;       phimax = mfdata[iz].phi_upper_stat[good]
;       polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
;         /data, color=im_color(col[iz]), noclip=0, /fill
;    endfor 
;
;; SDSS
;    mfdata = read_mf_vmax('sdss',noevol=noevol,/log,/active)
;    mffit_sdss = read_mf_vmax('sdss',/bestfit,noevol=noevol,/active)
;
;    good = where(mfdata.limit eq 1 and mfdata.number gt mingal)
;    phimass = mfdata.meanmass[good]
;    phimin = mfdata.phi_lower_stat[good]
;    phimax = mfdata.phi_upper_stat[good]
;    polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
;      /data, color=im_color(sdss_col,255), noclip=0, /fill
;
;    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
;
;; ---------------------------------------------------------------------------
;; sdss+primus - 'qq' stellar mass function on one single panel 
;    xrange = massrange1 
;    yrange = phirange1
;
;    sdss_zbins = mf_zbins(/sdss)
;    zbins = mf_zbins(nzbins)
;
;    maxis1 = range(xrange[0]+0.2,12.5,100)
;    sdss_col = 'forest green'
;    loadct, 3, /silent
;    col = reverse(findgen(nzbins)/(nzbins-1)*(230-50)+50)
;
;    mingal = 3
;    
;    psfile = paperpath+'mf_qq_onepanel'+suffix
;    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8, height=6.0
;
;    plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
;      ytitle=mfplot_phititle(), xtickinterval=1, ytickinterval=1
;    im_legend, 'Quiescent Galaxies', /right, /top, box=0, charsize=1.8
;    zlabel = ['z=',replicate('   ',7)]+strtrim(string([sdss_zbins.zlo,zbins.zlo],format='(F12.2)'),2)+$
;      '-'+strtrim(string([sdss_zbins.zup,zbins.zup],format='(F12.2)'),2)
;    yoff = 0.15
;    wid = 0.2 & hei = 0.1
;    x1 = xrange[0]+0.15 & y1 = -4.1
;    polyfill, [x1,x1+wid,x1+wid,x1],[y1,y1,y1+hei,y1+hei], $
;      /data, color=im_color(sdss_col), noclip=0, /fill
;    xyouts, x1+wid*1.15, y1+hei*0.2, zlabel[0], align=0.0, /data, charsize=1.3
;    for iz = 0, nzbins-1 do begin
;       polyfill, [x1,x1+wid,x1+wid,x1],[y1,y1,y1+hei,y1+hei]-yoff*(iz+1), $
;         /data, color=im_color(col[iz]), noclip=0, /fill
;       xyouts, x1+wid*1.15, y1+hei*0.2-yoff*(iz+1), zlabel[iz+1], align=0.0, /data, charsize=1.3
;    endfor
;      
;; PRIMUS    
;    mfdata = read_mf_vmax(noevol=noevol,/avgfield,/log,/quiescent)
;    for iz = 0, nzbins-1 do begin
;       good = where(mfdata[iz].limit eq 1 and mfdata[iz].number gt mingal)
;       phimass = mfdata[iz].meanmass[good]
;       phimin = mfdata[iz].phi_lower_stat[good]
;       phimax = mfdata[iz].phi_upper_stat[good]
;       polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
;         /data, color=im_color(col[iz]), noclip=0, /fill
;    endfor 
;
;; SDSS
;    mfdata = read_mf_vmax('sdss',noevol=noevol,/log,/quiescent)
;    mffit_sdss = read_mf_vmax('sdss',/bestfit,noevol=noevol,/quiescent)
;
;    good = where(mfdata.limit eq 1 and mfdata.number gt mingal)
;    phimass = mfdata.meanmass[good]
;    phimin = mfdata.phi_lower_stat[good]
;    phimax = mfdata.phi_upper_stat[good]
;    polyfill, [phimass,reverse(phimass)],[phimin,reverse(phimax)], $
;      /data, color=im_color(sdss_col), noclip=0, /fill
;
;    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
;


;;; ---------------------------------------------------------------------------
;;; sdss - compare with the literature - qq/red, sf/blue
;;    xrange = [8.8,12.2]
;;    yrange = [-7.0,-1.0]
;;    residrange = 1.49*[-1,1]
;;    maxis1 = range(xrange[0]+0.15,13,100)
;;    maxis_ber = range(10.5+0.25,13,100)
;;    maxis_cole = range(9.3,12.0,100)
;;    maxis = range(xrange[0]+0.15,xrange[1],100)
;;
;;    redpsym = 14 & redcolor = 'black' & redpsize = 0.8
;;    qqpsym = 15 & qqcolor = 'black' & qqpsize = 0.8
;;    sfpsym = 16 & sfcolor = 'black' & sfpsize = 0.8
;;    bluepsym = 17 & bluecolor = 'black' & bluepsize = 0.8
;;
;;; single-panel plot showing just all the galaxies                                 
;;    psfile = paperpath+'mf_sdss_qq_sf_lit'+suffix
;;    im_plotconfig, 5, pos, psfile=psfile, charsize=1.4, yspace=0.9, $
;;      height=[3.1,3.1], xmargin=[1.0,0.3], width=[3.4,3.4]
;;
;;; ##########
;;; red/g-r
;;    bell_color = 'firebrick' & bell_psym = 6
;;    baldry_color = 'orange' & baldry_psym = 9
;;
;;    bell = read_03bell_smf(/red)
;;    baldry = read_12baldry(/red)
;;
;;    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
;;      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
;;      ytitle=mfplot_phititle(), xtickinterval=1
;;    im_legend, '^{0.1}(g-r) Red Sequence', /right, /top, box=0, $
;;      charsize=1.2, margin=0
;;    label = ['Bell+03 - Red','Baldry+12 - Red','SDSS-GALEX - Red']
;;    im_legend, label, /left, /bottom, box=0, psym=[bell_psym,baldry_psym,redpsym], $
;;      pspacing=1.9, charsize=1.1, color=[bell_color,baldry_color,redcolor], $
;;      symthick=6, symsize=[1.0,1.0,1.5], margin=-0.2
;;
;;    mfdata = read_mf_vmax('sdss',/log,/gmr_red)
;;    good = where(mfdata.limit eq 1,nthese)
;;
;;    oploterror, mfdata.meanmass[good], mfdata.phi[good], mfdata.phierr_upper[good], $
;;      psym=symcat(redpsym,thick=8), symsize=redpsize, errthick=6, $
;;      color=im_color(redcolor,100+0), errcolor=im_color(redcolor,100), /hibar, /nohat
;;    oploterror, mfdata.meanmass[good], mfdata.phi[good], mfdata.phierr_lower[good], $
;;      psym=3, symsize=redpsize, errthick=5, $
;;      color=im_color(redcolor,100+0), errcolor=im_color(redcolor,100), /lobar, /nohat
;;
;;    oploterror, bell.mass, bell.phi, bell.phierr, $
;;      psym=symcat(bell_psym,thick=5), color=im_color(bell_color), $
;;      errcolor=im_color(bell_color), /nohat
;;    oploterror, baldry.mass, baldry.phi, baldry.phierr, $
;;      psym=symcat(baldry_psym,thick=5), color=im_color(baldry_color), $
;;      errcolor=im_color(baldry_color), /nohat
;;
;;;   params = {phistar:    3.25D-3, phiplus:    0.08D-3, logmstar:   10.72D, alpha:     -0.45, alphaplus: -1.45}
;;;   djs_oplot, maxis, alog10(mf_schechter_plus(maxis,schechter_plus=params)), color='red', line=5
;;
;;; ##########
;;; quiescent
;;    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
;;      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
;;      ytitle='', xtickinterval=1, ytickname=replicate(' ',10)
;;    im_legend, ['Quiescent based on SFR/M'], /right, /top, $
;;      box=0, charsize=1.2, margin=0
;;    label = ['Bell+03 - Red','Baldry+12 - Red','SDSS-GALEX - Quiescent']
;;    im_legend, label, /left, /bottom, box=0, psym=[bell_psym,baldry_psym,qqpsym], $
;;      pspacing=1.9, charsize=1.1, color=[bell_color,baldry_color,qqcolor], $
;;      symthick=6, symsize=[1.0,1.0,1.1], margin=-0.2
;;
;;    mfdata = read_mf_vmax('sdss',/log,/quiescent)
;;    good = where(mfdata.limit eq 1,nthese)
;;
;;    oploterror, mfdata.meanmass[good], mfdata.phi[good], mfdata.phierr_upper[good], $
;;      psym=symcat(qqpsym,thick=8), symsize=qqpsize, errthick=6, $
;;      color=im_color(qqcolor,100+0), errcolor=im_color(qqcolor,100), /hibar, /nohat
;;    oploterror, mfdata.meanmass[good], mfdata.phi[good], mfdata.phierr_lower[good], $
;;      psym=3, symsize=qqpsize, errthick=5, $
;;      color=im_color(qqcolor,100+0), errcolor=im_color(qqcolor,100), /lobar, /nohat
;;
;;    oploterror, bell.mass, bell.phi, bell.phierr, $
;;      psym=symcat(bell_psym,thick=5), color=im_color(bell_color), $
;;      errcolor=im_color(bell_color), /nohat
;;    oploterror, baldry.mass, baldry.phi, baldry.phierr, $
;;      psym=symcat(baldry_psym,thick=5), color=im_color(baldry_color), $
;;      errcolor=im_color(baldry_color), /nohat
;;
;;;   params = {phistar:    3.25D-3, phiplus:    0.08D-3, logmstar:   10.72D, alpha:     -0.45, alphaplus: -1.45}
;;;   djs_oplot, maxis, alog10(mf_schechter_plus(maxis,schechter_plus=params)), color='red', line=5
;;
;;; ##########
;;; blue
;;    bell_color = 'dodger blue' & bell_psym = 6
;;    baldry_color = 'midnight blue' & baldry_psym = 9
;;
;;    bell = read_03bell_smf(/blue)
;;    baldry = read_12baldry(/blue)
;;
;;    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,2], xsty=1, ysty=1, $
;;      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
;;      ytitle=mfplot_phititle(), xtickinterval=1
;;    im_legend, '^{0.1}(g-r) Blue Cloud', /right, /top, box=0, $
;;      charsize=1.2, margin=0
;;    label = ['Bell+03 - Blue','Baldry+12 - Blue','SDSS-GALEX - Blue']
;;    im_legend, label, /left, /bottom, box=0, psym=[bell_psym,baldry_psym,bluepsym], $
;;      pspacing=1.9, charsize=1.1, color=[bell_color,baldry_color,bluecolor], $
;;      symthick=6, symsize=[1.0,1.0,1.2], margin=-0.2
;;
;;    mfdata = read_mf_vmax('sdss',/log,/gmr_blue)
;;    good = where(mfdata.limit eq 1,nthese)
;;
;;    oploterror, mfdata.meanmass[good], mfdata.phi[good], mfdata.phierr_upper[good], $
;;      psym=symcat(bluepsym,thick=8), symsize=bluepsize, errthick=6, $
;;      color=im_color(bluecolor,100+0), errcolor=im_color(bluecolor,100), /hibar, /nohat
;;    oploterror, mfdata.meanmass[good], mfdata.phi[good], mfdata.phierr_lower[good], $
;;      psym=3, symsize=bluepsize, errthick=5, $
;;      color=im_color(bluecolor,100+0), errcolor=im_color(bluecolor,100), /lobar, /nohat
;;
;;    oploterror, bell.mass, bell.phi, bell.phierr, $
;;      psym=symcat(bell_psym,thick=5), color=im_color(bell_color), $
;;      errcolor=im_color(bell_color), /nohat
;;    oploterror, baldry.mass, baldry.phi, baldry.phierr, $
;;      psym=symcat(baldry_psym,thick=5), color=im_color(baldry_color), $
;;      errcolor=im_color(baldry_color), /nohat
;;
;;;   params = {phistar:    0.71D-3, logmstar:   10.72D, alpha:     -1.45}
;;;   djs_oplot, maxis, alog10(mf_schechter(maxis,schechter=params)), color='red', line=5
;;    
;;; ##########
;;; star-forming
;;    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,3], xsty=1, ysty=1, $
;;      yrange=yrange, xrange=xrange, xtitle=mfplot_masstitle(), $
;;      ytitle='', xtickinterval=1, ytickname=replicate(' ',10)
;;    im_legend, ['Star-Forming based','on SFR/M'], /right, /top, $
;;      box=0, charsize=1.2, margin=-0.1, spacing=1.7
;;    label = ['Bell+03 - Blue','Baldry+12 - Blue','SDSS-GALEX - Star-Forming']
;;    im_legend, label, /left, /bottom, box=0, psym=[bell_psym,baldry_psym,sfpsym], $
;;      pspacing=1.9, charsize=1.1, color=[bell_color,baldry_color,sfcolor], $
;;      symthick=6, symsize=[1.0,1.0,1.2], margin=-0.2
;;
;;    mfdata = read_mf_vmax('sdss',/log,/active)
;;    good = where(mfdata.limit eq 1,nthese)
;;
;;    oploterror, mfdata.meanmass[good], mfdata.phi[good], mfdata.phierr_upper[good], $
;;      psym=symcat(sfpsym,thick=8), symsize=sfpsize, errthick=6, $
;;      color=im_color(sfcolor,100+0), errcolor=im_color(sfcolor,100), /hibar, /nohat
;;    oploterror, mfdata.meanmass[good], mfdata.phi[good], mfdata.phierr_lower[good], $
;;      psym=3, symsize=sfpsize, errthick=5, $
;;      color=im_color(sfcolor,100+0), errcolor=im_color(sfcolor,100), /lobar, /nohat
;;
;;    oploterror, bell.mass, bell.phi, bell.phierr, $
;;      psym=symcat(bell_psym,thick=5), color=im_color(bell_color), $
;;      errcolor=im_color(bell_color), /nohat
;;    oploterror, baldry.mass, baldry.phi, baldry.phierr, $
;;      psym=symcat(baldry_psym,thick=5), color=im_color(baldry_color), $
;;      errcolor=im_color(baldry_color), /nohat
;;
;;;   params = {phistar:    0.71D-3, logmstar:   10.72D, alpha:     -1.45}
;;;   djs_oplot, maxis, alog10(mf_schechter(maxis,schechter=params)), color='red', line=5
;;    
;;    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
;;
;;stop    
;;    
