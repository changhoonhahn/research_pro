function get_vovervmax, zmin, zmax, area=area, zbins=zbins, $
  zprimus=zprimus
; see POST_VMAX; area in sr
;   vmax = (area/3.0)*(lf_comvol(zmax)-lf_comvol(zmin))
    nzvals = 100
    zvals = range(zbins.zlo,zbins.zup,nzvals)
    dz = zvals[1]-zvals[0]
    dh = im_light(/km)/100D
    dvdz = dh^3*dcomvoldz(zvals,0.3,0.7) ; =dvcomoving(z)
    
    comvolint = area*total(dvdz*dz,/cumulative)
    local_vmin = interpolate(comvolint,findex(zvals,zmin))
    local_vmax = interpolate(comvolint,findex(zvals,zmax))-local_vmin
    comvol = interpolate(comvolint,findex(zvals,zprimus))-local_vmin
    vovervmax = comvol/local_vmax
return, vovervmax
end

pro build_mf_data_chh, sdss=sdss, clobber=clobber, literature=literature
; jm11aug12ucsd - build all the data structures we are going to need
;    mfpath = '/global/data/scr/chh327/primus/mfdata/mfdata_chh/new/env_zrange/' ; z range of EDP
    mfpath = '/global/data/scr/chh327/primus/mfdata/mfdata_chh/new/'    ; z range of John's data 
    h100 = mf_h100()
    q0 = mf_q0(qz0=qz0)
    
    if keyword_set(sdss) then prefix = '_sdss' else prefix = ''
    if keyword_set(literature) then litsuffix = '_lit' else litsuffix = ''

; read the VAGC window function area to speed up the computations below
    if keyword_set(sdss) and (n_elements(vagc_area) eq 0) then $
      windowfile = get_mf_window_chh('sdss',area=vagc_area,/sr)

; specify the fields and redshift bins
    field = get_mf_fields(sdss=sdss)
;   field = field[where(strtrim(field,2) ne 'cosmos')]
    nfield = n_elements(field)
;    zbins = mf_zbins_chh(sdss=sdss,nzbins,literature=literature)
    zbins = mf_zbins(sdss=sdss, nzbins, literature=literature)
    struct_print, zbins

; select quiescent and star-forming galaxies using various criteria
   subsample = ['all','quiescent','active']
;   subsample = ['all','quiescent','active','red','blue']
;   subsample = ['all','quiescent','active','quiescent_nuvmr','active_nuvmr','red','blue']
;    subsample = ['all','quiescent','active','red','blue']
;    if keyword_set(sdss) then subsample = [subsample,'red_gmr','blue_gmr']
;   subsample = ['blue_gmr']
;   subsample = ['all','quiescent','active','red','blue','red_gmr',$
;     'blue_gmr','uvj_quiescent','uvj_active']

; hard-coded for now!
    if keyword_set(sdss) then absmaglimit = -17.0 else begin
       if keyword_set(literature) then $
         absmaglimit = [-17.5,-19.0,-20.0,-21.0] else $
           absmaglimit = [-17.5,-18.5,-19.0,-19.5,-20.5,-21.0]
    endelse
    
; read the supergrid parameter file    
    super = get_mf_supergrid(nsuper=nsuper,superstring=superstring,/models)
;   super = get_mf_supergrid(nsuper=nsuper,superstring=superstring)

; just do core subsamples and supergrid01 if /literature    
    if keyword_set(literature) then begin
       super = get_mf_supergrid(1,nsuper=nsuper,superstring=superstring)
       subsample = ['all','quiescent','active']
    endif

; gather all the quantities we will need; looping over each field
; supergrid, and galaxy sample
    allarea = dblarr(nfield)
;   for ii = 1, nfield-1 do begin
    for ii = 0, nfield-1 do begin
; read the full parent sample
       allparent = read_mf_parent_chh(field[ii])
       windowfile = get_mf_window_chh(field[ii],area=area,/sr)
       allarea[ii] = area ; [sr]
       filters = get_mf_filters_chh(field[ii])
       galex = strmatch(filters,'*galex*',/fold)
       nfilt = n_elements(filters)
       for jj = 0, nsuper-1 do begin
          allisedparent = read_mf_isedfit_chh(field[ii],supergrid=super[jj].supergrid,/parent)
; three selections: (1) ^{0.1}(NUV-r) vs mass; (2) ^{0.1}(u-r) vs
; Mr; and (3) U-V vs V-J
          sfr = allisedparent.sfr100_50+0.30
;          sfr = allisedparent.sfr100_50
          mass = allisedparent.mass_50
          sfrm = sfr-mass+9.0 ; [Gyr^-1]
          nuv = allparent.k_fnuvugrizjhk_absmag_01[1]
          mr = allparent.k_fnuvugrizjhk_absmag_01[4]
          mu = allparent.k_fnuvugrizjhk_absmag_01[2]
          mg = allparent.k_fnuvugrizjhk_absmag_01[3]
          nuvmr = nuv-mr
          umr = mu-mr
          gmr = mg-mr

          bigmu = allparent.k_fnuvubvrijhk_absmag_00[2]
          bigmv = allparent.k_fnuvubvrijhk_absmag_00[4]
          bigmj = allparent.k_fnuvubvrijhk_absmag_00[7]
          bigumv = bigmu-bigmv
          bigvmj = bigmv-bigmj

;         for ss = 6, 6 do begin
          delvarx, indx
          for ss = 0, n_elements(subsample)-1 do begin
; pick the right sample
             if subsample[ss] eq 'all' then indx = lindgen(n_elements(allparent))
             if subsample[ss] eq 'quiescent' then indx = mf_select_sfr_quiescent_chh(sfr,mass,z=allparent.zprimus)                   ; quiescent
             if subsample[ss] eq 'active' then indx = mf_select_sfr_quiescent_chh(sfr,mass,z=allparent.zprimus,/active)              ; star-forming

             parent = allparent[indx]
             isedparent = allisedparent[indx]
             ngal = n_elements(parent)

             print, subsample[ss], 'N_GAL = ',ngal 
; read the limits file
             if keyword_set(sdss) then begin
                masslim = 8.95
;               masslim = 8.775 ; 9.0
;               masslim = 8.975 ; 9.0
                masslim_50 = masslim
             endif else begin
                limitsfile = '/global/data/scr/chh327/primus/mfdata/'+$
                    'limits'+prefix+'_'+subsample[ss]+$
                  '_supergrid'+superstring[jj]+litsuffix+'.fits.gz'
                print, limitsfile 
                limits = mrdfits(limitsfile,1)
                limits1 = limits[where(strtrim(limits.field,2) eq field[ii])]
                masslim = limits1.masslim_95
                masslim_50 = limits1.masslim_50
             endelse
             
; initialize the output data structure
             mfdata = replicate({area: area, ra: 0D, dec: 0D, parent_id: 0L, isedfit_id: 0L, $
               masslimit: 0.0, masslimit_50: 0.0, minfitmass: 0.0, absmaglimit: 0.0, $
               z: 0.0, zprimus: 0.0, zspec: -1.0, ised_chi2: 0.0, k_maggies: fltarr(nfilt), k_ivarmaggies: fltarr(nfilt), $
               k_chi2: 0.0, k_coeffs: fltarr(5), galex: 0, mass: 0.0, k_mass: 0.0, $
               targ_mag: 0.0, r50: -999.0, sb50: -999.0, mb_00: 0.0, mr_01: 0.0, mg_01:0.0, $
               umb_00: 0.0, nuvmr_01: 0.0, nuvmr_00: 0.0, umr_01: 0.0, gmr_01: 0.0, rmj_01: 0.0, rmj_00: 0.0, $
               bigumv: 0.0, bigvmj: 0.0, mb_00_evol: 0.0, mr_01_evol: 0.0, $
               weight: 0.0, zmin_evol:0.0, zmin_noevol:0.0, zmax_evol:0.0, zmax_noevol:0.0,$
               vmax_evol: 0D, vmax_noevol: 0D, vovervmax_evol: 0D, vovervmax_noevol: 0D, $
               age: 0.0, tau: 0.0, av: 0.0, metallicity: 0.0, sfr: 0.0, sfr100: 0.0, sfrm100: 0.0},ngal)
             mfdata.parent_id = indx
             mfdata.ra = parent.ra
             mfdata.dec = parent.dec
             mfdata.z = parent.zprimus
             mfdata.k_chi2 = parent.k_chi2
             mfdata.k_coeffs = parent.k_coeffs
             mfdata.k_mass = parent.k_mass
             mfdata.k_maggies = parent.k_maggies
             mfdata.k_ivarmaggies = parent.k_ivarmaggies
             mfdata.weight = parent.final_weight
             mfdata.r50 = parent.r50
             mfdata.sb50 = parent.sb50
             mfdata.targ_mag = parent.targ_mag

; get the original spectroscopic redshifts so I can calculate the
; outlier rate
             if tag_exist(parent,'zspec_flag') then begin
                iszspec = where(parent.zspec_flag,nzspec,comp=iszprimus)
                if (nzspec ne 0L) then begin
                   mfdata.zprimus = parent.zprimus
                   mfdata[iszspec].zspec = parent[iszspec].zprimus
                   mfdata[iszspec].zprimus = parent[iszspec].zprimus_original
                endif
             endif
             
             mfdata.galex = total(parent.k_ivarmaggies[galex] gt 0,1) gt 0
             
             mfdata.mb_00 = parent.k_fnuvubvrijhk_absmag_00[3]
             mfdata.mr_01 = parent.k_fnuvugrizjhk_absmag_01[4]
             mfdata.mg_01 = parent.k_fnuvugrizjhk_absmag_01[4]
;            mfdata.mz_06 = parent.k_ugrizjhk_absmag_06[4]

             mfdata.umb_00 = parent.k_fnuvubvrijhk_absmag_00[2]-parent.k_fnuvubvrijhk_absmag_00[3]
             mfdata.nuvmr_01 = parent.k_fnuvugrizjhk_absmag_01[1]-parent.k_fnuvugrizjhk_absmag_01[4]
             mfdata.nuvmr_00 = parent.k_fnuvubvrijhk_absmag_00[1]-parent.k_fnuvubvrijhk_absmag_00[5]
             mfdata.gmr_01 = parent.k_fnuvugrizjhk_absmag_01[3]-parent.k_fnuvugrizjhk_absmag_01[4]
             mfdata.umr_01 = parent.k_fnuvugrizjhk_absmag_01[2]-parent.k_fnuvugrizjhk_absmag_01[4]
             mfdata.rmj_01 = parent.k_fnuvugrizjhk_absmag_01[4]-parent.k_fnuvugrizjhk_absmag_01[7]
             mfdata.rmj_00 = parent.k_fnuvubvrijhk_absmag_00[5]-parent.k_fnuvubvrijhk_absmag_00[7]
;            mfdata.umz_06 = parent.k_ugrizjhk_absmag_06[0]-parent.k_ugrizjhk_absmag_06[4]
             
             mfdata.bigumv = parent.k_fnuvubvrijhk_absmag_00[2]-parent.k_fnuvubvrijhk_absmag_00[4]
             mfdata.bigvmj = parent.k_fnuvubvrijhk_absmag_00[4]-parent.k_fnuvubvrijhk_absmag_00[7]
             
             mfdata.mb_00_evol = k_evolve(mfdata.mb_00,mfdata.z,q0,0.0,qz0)
             mfdata.mr_01_evol = k_evolve(mfdata.mr_01,mfdata.z,q0,0.0,qz0)
;            mfdata.mz_06_evol = k_evolve(mfdata.mz_06,mfdata.z,q0,0.0,qz0)


             mfdata.isedfit_id = isedparent.isedfit_id
             mfdata.ised_chi2 = isedparent.chi2
             mfdata.mass = isedparent.mass_50
             mfdata.age = isedparent.age_50
             mfdata.tau = isedparent.tau_50
             mfdata.av = isedparent.av_50
             mfdata.sfr = isedparent.sfr_50
             mfdata.sfr100 = isedparent.sfr100_50+0.30  ;Adjustment made to compensate for SFR discrepancies in mfdata
;             mfdata.sfr100 = isedparent.sfr100_50
             mfdata.sfrm100 = isedparent.sfr100_50-isedparent.mass_50+9.0 ; [Gyr^-1]
             mfdata.metallicity = isedparent.Z_50
            
;             if keyword_set(sdss) then mfdata = mfdata[where(mfdata.z GE zbins.zlo AND mfdata.z LT zbins.zup)]

             mfdatafile = mfpath+'mfdata_'+subsample[ss]+'_supergrid'+$
               superstring[jj]+'_'+field[ii]+litsuffix+'.fits'
             if im_file_test(mfdatafile+'.gz',clobber=clobber) then continue

; loop on redshift to derive Vmax for each object
             for iz = 0, nzbins-1 do begin
                print, zbins[iz].zlo, zbins[iz].zup
                these = where((parent.zprimus ge zbins[iz].zlo) and $
                  (parent.zprimus lt zbins[iz].zup))
                sample = parent[these]
                ised = isedparent[these]
                mfdata[these].masslimit = masslim[iz]
                mfdata[these].masslimit_50 = masslim_50[iz]
                mfdata[these].absmaglimit = absmaglimit[iz]

                zmin_evol = sample.zmin_evol>zbins[iz].zlo
                zmax_evol = sample.zmax_evol<zbins[iz].zup
                zmin_noevol = sample.zmin_noevol>zbins[iz].zlo
                zmax_noevol = sample.zmax_noevol<zbins[iz].zup

                mfdata[these].zmin_evol = sample.zmin_evol
                mfdata[these].zmax_evol = sample.zmax_evol 
                mfdata[these].zmin_noevol = sample.zmin_noevol
                mfdata[these].zmax_noevol = sample.zmax_noevol 
                mfdata[these].vmax_evol = (area/3.0)*(lf_comvol(zmax_evol)-$
                  lf_comvol(zmin_evol))*(1.0/h100)^3.0 ; h=1-->h=0.7
                mfdata[these].vmax_noevol = (area/3.0)*(lf_comvol(zmax_noevol)-$
                  lf_comvol(zmin_noevol))*(1.0/h100)^3.0

                mfdata[these].vovervmax_noevol = get_vovervmax(zmin_noevol,$
                  zmax_noevol,area=area,zbins=zbins[iz],zprimus=sample.zprimus)
                mfdata[these].vovervmax_evol = get_vovervmax(zmin_evol,$
                  zmax_evol,area=area,zbins=zbins[iz],zprimus=sample.zprimus)
             endfor 
             im_mwrfits, mfdata, mfdatafile, clobber=clobber
          endfor
       endfor
    endfor 

; now loop back through and merge all the fields together
    if (keyword_set(sdss) eq 0) then begin
       totarea = total(allarea,/double)
       for jj = 0, nsuper-1 do begin
          for ss = 0, n_elements(subsample)-1 do begin
             for ii = 0, nfield-1 do begin
                mfdatafile = mfpath+'mfdata_'+subsample[ss]+'_supergrid'+$
                  superstring[jj]+'_'+field[ii]+litsuffix+'.fits.gz'
                mfdata = mrdfits(mfdatafile,1)
; weight each field by the fractional area
                mfdata.vmax_evol = (totarea/allarea[ii])*mfdata.vmax_evol
                mfdata.vmax_noevol = (totarea/allarea[ii])*mfdata.vmax_noevol

                mfdata = struct_addtags(replicate({field: field[ii]},$
                  n_elements(mfdata)),mfdata)
                mfdata = struct_trimtags(temporary(mfdata),except=['k_maggies','k_ivarmaggies'])
                if (ii eq 0) then bigmfdata = mfdata else $
                  bigmfdata = [bigmfdata,mfdata]
             endfor
             bigmfdatafile = mfpath+'mfdata_'+subsample[ss]+'_supergrid'+$
               superstring[jj]+litsuffix+'.fits'
             if im_file_test(bigmfdatafile+'.gz',clobber=clobber) then continue
             im_mwrfits, bigmfdata, bigmfdatafile, clobber=clobber
          endfor
       endfor
    endif

return
end

