pro plot_sfr_dist, mr_cut, mass_cut
    datadir = '/global/data/scr/chh327/tinker/group_catalog/' 
    galdata_name    = 'clf_groups_M'+strtrim(string(mr_cut),2)+'_'+strmid(strtrim(string(mass_cut),2),0,3)+'_D360.galdata_corr.fits'
    probdata_name   = 'clf_groups_M'+strtrim(string(mr_cut),2)+'_'+strmid(strtrim(string(mass_cut),2),0,3)+'_D360.prob.fits'
    galdata = mrdfits(datadir+galdata_name,1)
    probdata= mrdfits(datadir+probdata_name,1)

    cen_lims = where(alog10(galdata.stellmass) GT 10.0 and alog10(galdata.stellmass) LT 10.5 and probdata.p_sat LE 0.5, cen_limcount)
    print, n_elements(galdata), cen_limcount
   
    central_gal = galdata[cen_lims]
    SSFR = central_gal.SFR-alog10(central_gal.stellmass)
   
    psfile = 'SSFR_hist.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8
    cgHistoplot, SSFR, BINSIZE=0.2, /Fill
    im_plotconfig, /psclose, psfile=psfile, pdf=pdf
end
