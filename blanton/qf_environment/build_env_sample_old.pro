;+
; NAME:
;   build_env_sample
;
; PURPOSE:
;   Defines an environmental sample using the constraints specified in definesample_index.par
;   only includes objects in 2mask fields. Note that the constraints must be first specified 
;   in definesample_index.par before running this procedure. 
;
; CALLING SEQUENCE:
;   define_environment, name=
;
; INPUTS:
;   All inputs are specified in define_sample_index.par
;
; OUTPUTS:
;   .fits file with name = name
;
; FUNCTIONS USED:
;   yanny_readone( filename, [ selectname, hdr=, enums=, structs=, /anonymous, stnames=, 
;                  /quick, errcode= ])
;   is_in_window( xyz= , ra= , dec= , polygons)
;
; PROCEDURES USED: 
;   read_mangle_polygons, infile, polygons, id [, unit=]
;
; MODIFICATION HISTORY:
; Oct 8 2011: Chang Hoon Hahn
		; Jan 2 2012: Chang Hoon Hahn
; May 25 2012: Chang Hoon Hahn
;--------------------------------------------------------
pro build_env_sample, run=run, primus=primus, sdss=sdss, numden=numden, parent=parent
    parpath = get_path(/repo)
    parameters = yanny_readone(parpath+'zerod_environment_parameters.par', hdr=hdr)
    para = where(parameters.run eq run)
    param= parameters[para]
    name = param.name
    zmin = param.zmin
    zmax = param.zmax
    dim = dblarr(5)


    fiddir = '/global/data/scr/chh327/primus/science/mf/isedfit/'
    kcorrdir = '/global/data/scr/chh327/primus/science/mf/kcorrect/'
    zeroddir = '/global/data/scr/chh327/primus/science/mf/2165/ubersample/'

    fiducial1 = mrdfits(fiddir+'xmm_swire_fsps_chab_charlot_sfhgrid01.fits.gz', 1)
    kcorr1    = mrdfits(kcorrdir+'xmm_swire_kcorr_v20.fits.gz',1)
    zerod1    = mrdfits(zeroddir+'xmm_swire_zerod_v20.fits.gz',1)
    dim[0] 	  = N_ELEMENTS(fiducial1) 

    fiducial2 = mrdfits(fiddir+'cdfs_fsps_chab_charlot_sfhgrid01.fits.gz', 1)
    kcorr2    = mrdfits(kcorrdir+'cdfs_kcorr_v20.fits.gz',1)
    zerod2    = mrdfits(zeroddir+'cdfs_zerod_v20.fits.gz',1)
    dim[1] 	  = N_ELEMENTS(fiducial2)  

    fiducial3 = mrdfits(fiddir+'cfhtls_xmm_fsps_chab_charlot_sfhgrid01.fits.gz', 1)
    kcorr3    = mrdfits(kcorrdir+'cfhtls_xmm_kcorr_v20.fits.gz',1)
    zerod3    = mrdfits(zeroddir+'cfhtls_xmm_zerod_v20.fits.gz',1)
    dim[2] 	  = N_ELEMENTS(fiducial3) 

    fiducial4 = mrdfits(fiddir+'cosmos_fsps_chab_charlot_sfhgrid01.fits.gz', 1)
    kcorr4    = mrdfits(kcorrdir+'cosmos_kcorr_v20.fits.gz',1)
    zerod4    = mrdfits(zeroddir+'cosmos_zerod_v20.fits.gz',1)
    dim[3]	  = N_ELEMENTS(fiducial4) 

    fiducial5 = mrdfits(fiddir+'es1_fsps_chab_charlot_sfhgrid01.fits.gz', 1)
    kcorr5    = mrdfits(kcorrdir+'es1_kcorr_v20.fits.gz',1)
    zerod5    = mrdfits(zeroddir+'es1_zerod_v20.fits.gz',1)
    dim[4] 	  = N_ELEMENTS(fiducial5) 

    targ = {class:' ',ra:0.D,dec:0.D,zconf:0,redshift:0.,g:0.,i:0.,fuv:0.,nuv:0.,mass:0.,age:0.,SFR:0.,mass_err:0.,age_err:0.,SFR_err:0.,inwindow_2mask:0L,comp:0L}
    dimension = total(dim)
    targs = replicate(targ, dimension)

    binstart = 0L 
    binend = binstart + dim[0] - 1		
    targs[binstart:binend].class=zerod1.zprimus_class
    targs[binstart:binend].ra = zerod1.ra
    targs[binstart:binend].dec = zerod1.dec
    targs[binstart:binend].zconf = zerod1.zprimus_zconf
    targs[binstart:binend].redshift = zerod1.zprimus
    targs[binstart:binend].g = transpose(kcorr1.k_fnuvugrizjhk_absmag_01[3,*])
    targs[binstart:binend].i = transpose(kcorr1.k_fnuvugrizjhk_absmag_01[5,*])
    targs[binstart:binend].fuv = transpose(kcorr1.k_fnuvugrizjhk_absmag_01[0,*])
    targs[binstart:binend].nuv = transpose(kcorr1.k_fnuvugrizjhk_absmag_01[1,*])
    targs[binstart:binend].mass = fiducial1.mass
    targs[binstart:binend].age = fiducial1.age
    targs[binstart:binend].SFR = fiducial1.SFR
    targs[binstart:binend].mass_err = fiducial1.mass_err
    targs[binstart:binend].age_err = fiducial1.age_err
    targs[binstart:binend].SFR_err = fiducial1.SFR_err
    targs[binstart:binend].inwindow_2mask = zerod1.inwindow_2mask
    binstart = binstart + dim[0] 

    binend = binstart + dim[1] - 1		
    targs[binstart:binend].class=zerod2.zprimus_class
    targs[binstart:binend].ra = zerod2.ra
    targs[binstart:binend].dec = zerod2.dec
    targs[binstart:binend].zconf = zerod2.zprimus_zconf
    targs[binstart:binend].redshift = zerod2.zprimus
    targs[binstart:binend].g = transpose(kcorr2.k_fnuvugrizjhk_absmag_01[3,*])
    targs[binstart:binend].i = transpose(kcorr2.k_fnuvugrizjhk_absmag_01[5,*])
    targs[binstart:binend].fuv = transpose(kcorr2.k_fnuvugrizjhk_absmag_01[0,*])
    targs[binstart:binend].nuv = transpose(kcorr2.k_fnuvugrizjhk_absmag_01[1,*])
    targs[binstart:binend].mass = fiducial2.mass
    targs[binstart:binend].age = fiducial2.age
    targs[binstart:binend].SFR = fiducial2.SFR
    targs[binstart:binend].mass_err = fiducial2.mass_err
    targs[binstart:binend].age_err = fiducial2.age_err
    targs[binstart:binend].SFR_err = fiducial2.SFR_err
    targs[binstart:binend].inwindow_2mask = zerod2.inwindow_2mask
    binstart = binstart + dim[1] 

    binend = binstart + dim[2] - 1	
    targs[binstart:binend].class=zerod3.zprimus_class	
    targs[binstart:binend].ra = zerod3.ra
    targs[binstart:binend].dec = zerod3.dec
    targs[binstart:binend].zconf = zerod3.zprimus_zconf
    targs[binstart:binend].redshift = zerod3.zprimus
    targs[binstart:binend].g = transpose(kcorr3.k_fnuvugrizjhk_absmag_01[3,*])
    targs[binstart:binend].i = transpose(kcorr3.k_fnuvugrizjhk_absmag_01[5,*])
    targs[binstart:binend].fuv = transpose(kcorr3.k_fnuvugrizjhk_absmag_01[0,*])
    targs[binstart:binend].nuv = transpose(kcorr3.k_fnuvugrizjhk_absmag_01[1,*])
    targs[binstart:binend].mass = fiducial3.mass
    targs[binstart:binend].age = fiducial3.age
    targs[binstart:binend].SFR = fiducial3.SFR
    targs[binstart:binend].mass_err = fiducial3.mass_err
    targs[binstart:binend].age_err = fiducial3.age_err
    targs[binstart:binend].SFR_err = fiducial3.SFR_err
    targs[binstart:binend].inwindow_2mask = zerod3.inwindow_2mask
    binstart = binstart + dim[2] 

    binend = binstart + dim[3] - 1	
    targs[binstart:binend].class=zerod4.zprimus_class	
    targs[binstart:binend].ra = zerod4.ra
    targs[binstart:binend].dec = zerod4.dec
    targs[binstart:binend].zconf = zerod4.zprimus_zconf
    targs[binstart:binend].redshift = zerod4.zprimus
    targs[binstart:binend].g = transpose(kcorr4.k_fnuvugrizjhk_absmag_01[3,*])
    targs[binstart:binend].i = transpose(kcorr4.k_fnuvugrizjhk_absmag_01[5,*])
    targs[binstart:binend].fuv = transpose(kcorr4.k_fnuvugrizjhk_absmag_01[0,*])
    targs[binstart:binend].nuv = transpose(kcorr4.k_fnuvugrizjhk_absmag_01[1,*])
    targs[binstart:binend].mass = fiducial4.mass
    targs[binstart:binend].age = fiducial4.age
    targs[binstart:binend].SFR = fiducial4.SFR
    targs[binstart:binend].mass_err = fiducial4.mass_err
    targs[binstart:binend].age_err = fiducial4.age_err
    targs[binstart:binend].SFR_err = fiducial4.SFR_err
    targs[binstart:binend].inwindow_2mask = zerod4.inwindow_2mask
    binstart = binstart + dim[3] 

    binend = binstart + dim[4] - 1	
    targs[binstart:binend].class=zerod5.zprimus_class	
    targs[binstart:binend].ra = zerod5.ra
    targs[binstart:binend].dec = zerod5.dec
    targs[binstart:binend].zconf = zerod5.zprimus_zconf
    targs[binstart:binend].redshift = zerod5.zprimus
    targs[binstart:binend].g = transpose(kcorr5.k_fnuvugrizjhk_absmag_01[3,*])
    targs[binstart:binend].i = transpose(kcorr5.k_fnuvugrizjhk_absmag_01[5,*])
    targs[binstart:binend].fuv = transpose(kcorr5.k_fnuvugrizjhk_absmag_01[0,*])
    targs[binstart:binend].nuv = transpose(kcorr5.k_fnuvugrizjhk_absmag_01[1,*])
    targs[binstart:binend].mass = fiducial5.mass
    targs[binstart:binend].age = fiducial5.age
    targs[binstart:binend].SFR = fiducial5.SFR
    targs[binstart:binend].mass_err = fiducial5.mass_err
    targs[binstart:binend].age_err = fiducial5.age_err
    targs[binstart:binend].SFR_err = fiducial5.SFR_err
    targs[binstart:binend].inwindow_2mask = zerod5.inwindow_2mask
    binstart = binstart + dim[4] 
    indx0 = where(targs.zconf ge 3 AND targs.zconf le 4 AND targs.class eq 'GALAXY' and targs.inwindow_2mask eq 1$ 
        AND targs.redshift lt 0.8 AND targs.redshift ge 0.2)
    smfenv = targs[indx0]

    if keyword_set(parent) then begin
        if keyword_set(primus) then begin
            dir = get_path(/envt)
            mwrfits, smfenv, dir+'parent_'+name, /create
        endif
    endif

    primuspoly = get_poly_area(/smf, /sr)
    h100 = mf_h100()
    
    bin1 = where(smfenv.redshift ge 0.2 AND smfenv.redshift lt 0.4)
    bin2 = where(smfenv.redshift ge 0.4 AND smfenv.redshift lt 0.6) 
    bin3 = where(smfenv.redshift ge 0.6 AND smfenv.redshift lt 0.8)

    high= smfenv[bin3]  
    high_comvol = (primuspoly/3.0)*(lf_comvol(0.8)-lf_comvol(0.6))*(1.0/h100)^3.0
    high_hist = hist_smf(high.i, high_comvol[0], binsize=0.01)
    high_cumul= total(high_hist[1,*], /cumulative)
    high_x = high_hist[0,*]
    high_n = value_locate(high_x,-21.0)
    high_lim = high_x[high_n+1L]

    mid = smfenv[bin2]
    mid_comvol = (primuspoly/3.0)*(lf_comvol(0.6)-lf_comvol(0.4))*(1.0/h100)^3.0
    mid_hist = hist_smf(mid.i, mid_comvol[0], binsize=0.01)
    mid_cumul = total(mid_hist[1,*], /cumulative)
    mid_x = mid_hist[0,*]
    mid_n = value_locate(mid_cumul, high_cumul[high_n])
    mid_lim = mid_x[mid_n+1L]

    low = smfenv[bin1]
    low_comvol = (primuspoly/3.0)*(lf_comvol(0.4)-lf_comvol(0.2))*(1.0/h100)^3.0
    low_hist = hist_smf(low.i, low_comvol[0], binsize=0.01)
    low_cumul = total(low_hist[1,*], /cumulative)
    low_x = low_hist[0,*]
    low_n = value_locate(low_cumul, high_cumul[high_n])
    low_lim = low_x[low_n+1L]


    limit = dblarr(n_elements(high_x))
    limit[*] = high_cumul[high_n]

            plot, low_x, low_cumul, xrange=[-26.0, -18.0], linestyle=0, thick=3
	        oplot, high_x, limit
	        oplot, mid_x, mid_cumul, linestyle=2, thick=3
	        oplot, high_x, high_cumul, linestyle=3, thick=3

    print, n_elements(bin1), n_elements(bin2), n_elements(bin3)	
    print, 'primus_numden:',n_elements(low)/low_comvol,n_elements(mid)/mid_comvol,n_elements(high)/high_comvol
    print, 'primus_n:',low_n, mid_n, high_n
    print, 'primus_comvol:', low_comvol[0], mid_comvol[0], high_comvol[0]
    print, 'primus_lim:',low_lim, mid_lim, high_lim
    print, 'primus_cumul:', low_cumul[low_n], low_cumul[low_n+1L], mid_cumul[mid_n], high_cumul[high_n]

    low_bin = where(smfenv.redshift ge 0.2 AND smfenv.redshift lt 0.4 AND smfenv.i le low_lim) 
    mid_bin = where(smfenv.redshift ge 0.4 AND smfenv.redshift lt 0.6 AND smfenv.i le mid_lim)
    high_bin= where(smfenv.redshift ge 0.6 AND smfenv.redshift lt 0.8 AND smfenv.i le high_lim)

    targs.comp = 0L 
    smfenv[low_bin].comp = 1L 
    smfenv[mid_bin].comp = 1L
    smfenv[high_bin].comp= 1L 

;We impose the constraints for our sample above from above. 
    complete = where(smfenv.comp eq 1L,completecount)
    print, 'Completeness:',float(completecount)/float(n_elements(smfenv))
    smf_env = smfenv[complete]

    print, "length=", n_elements(smf_env)
    
    if keyword_set(numden) then begin
        if keyword_set(primus) then begin 
            output = {ra:0.D,dec:0.D,zprimus:0.,Mi:0.}
            outputs = replicate(output, n_elements(smf_env))
            outputs.ra = smf_env.ra
            outputs.dec = smf_env.dec
            outputs.zprimus = smf_env.redshift
            outputs.Mi = smf_env.i 
        
            dir = get_path(/envt)
            mwrfits, outputs, dir+name, /create
        endif
    endif

    if keyword_set(sdss) then begin
        dir = get_path(/mfdata)
        sdss_data = mrdfits( dir+'mfdata_all_supergrid01_sdss.fits.gz',1)
        sdss_fiducial = mrdfits('/mount/moon1/ioannis/archive/primus/mf/isedfit/sdss_fsps_chab_charlot_sfhgrid01.fits.gz', 1)
        sdss_kcorr    = mrdfits('/mount/moon1/ioannis/archive/primus/mf/kcorrect/sdss_kcorr.fits.gz',1)
        sdss_zerod    = mrdfits('/mount/moon1/ioannis/research/projects/primus/mf/2165/ubersample/sdss_ubersample.fits.gz',1)
        sdss_dim    = N_ELEMENTS(sdss_fiducial)
        sdss = replicate({ra:0.D,dec:0.D,redshift:0.,g:0.,r:0.,i:0.,fuv:0.,nuv:0.,mass:0.,age:0.,SFR:0.,$
            mass_err:0.,age_err:0.,SFR_err:0.,comp:0L,weight:0.,vmax:0L}, sdss_dim)

        sdss.ra = sdss_zerod.ra
        sdss.dec = sdss_zerod.dec
        sdss.redshift = sdss_zerod.z
        sdss.g = sdss_kcorr.k_fnuvugrizjhk_absmag_01[3]
        sdss.r = sdss_kcorr.k_fnuvugrizjhk_absmag_01[4]
        sdss.i = sdss_kcorr.k_fnuvugrizjhk_absmag_01[5]
        sdss.fuv = sdss_kcorr.k_fnuvugrizjhk_absmag_01[0]
        sdss.nuv = sdss_kcorr.k_fnuvugrizjhk_absmag_01[1]
        sdss.mass = sdss_fiducial.mass
        sdss.age = sdss_fiducial.age
        sdss.SFR = sdss_fiducial.SFR
        sdss.mass_err = sdss_fiducial.mass_err
        sdss.age_err = sdss_fiducial.age_err
        sdss.SFR_err = sdss_fiducial.SFR_err
        sdss.weight = sdss_data.weight
        sdss.vmax = sdss_data.vmax_evol

        if keyword_set(parent) then begin
            dir = get_path(/envt)
            mwrfits, sdss, dir+'parent_sdss_'+name, /create
        endif
        sdss_hist = hist_smf(sdss.i, sdss.vmax, weight=sdss.weight, binsize=0.01)
        sdss_cumul = float(total(sdss_hist[1,*], /cumulative))
        sdss_x = sdss_hist[0,*]
        sdsspoly = get_poly_area(/sdss, /sr)
        sdss_comvol = (sdsspoly/3.0)*(lf_comvol(0.2)-lf_comvol(0.01))*(1.0/h100)^3.0
               
        sdss_n = value_locate(sdss_cumul, high_cumul[high_n])
        sdss_lim = sdss_x[sdss_n]
        print, 'sdss_numden:', sdss_dim/sdss_comvol
        print, 'sdss_comvol', sdss_comvol
        print, 'sdss_lim:',sdss_lim
        print, 'sdss_cumul',sdss_cumul[sdss_n]

        oplot, sdss_x, sdss_cumul
        sdss_bin = where(sdss.i le sdss_lim)
        sdss.comp = 0L
        sdss[sdss_bin].comp = 1L
        sdss_complete = where(sdss.comp eq 1L,completecount)
        sdss_env = sdss[sdss_complete]
        
        print, "length=", n_elements(sdss_env), 'fraction=', float(n_elements(sdss_env))/float(sdss_dim)
        
        if keyword_set(numden) then begin
            sdss_output = replicate({ra:0.D,dec:0.D,zprimus:0.,Mi:0.}, n_elements(sdss_env))
            sdss_output.ra = sdss_env.ra
            sdss_output.dec = sdss_env.dec
            sdss_output.zprimus = sdss_env.redshift
            sdss_output.Mi = sdss_env.i
            
            sdss_dir = get_path(/envt)
            mwrfits, sdss_output, sdss_dir+name, /create
        endif
    endif
end
