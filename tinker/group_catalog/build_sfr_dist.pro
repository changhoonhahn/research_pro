pro build_sfr_dist, silent=silent
; Builds a sample of galaxies with SFR distribution using SFR main sequence
; from PRIMUS SFR-M* relation 
    massbin = mass_bins(nmassbin,minmass=minmass, maxmass=maxmass,/simple)
;    massbin = massbin[where(massbin.massbin EQ 10.25,nmassbin)] ; limit to one mass bin 
    struct_print, massbin

    z = 0.9
    mbin = 0.01
; import SMF file used to normalize the number of galaxies in mass bins
    smf_dir = '/global/data/scr/chh327/tinker/group_catalog/smf/'
    smf_file = 'smf_mfdata_all_supergrid01_lit_0810z_mbin'+strmid(strtrim(string(mbin),2),0,4)+'_john.fits' ; currently hardcoded
    smf_data = mrdfits(smf_dir+smf_file,1)
    print, smf_dir+smf_file

    total_mass = []
    total_SSFR = []
; randomly chose number of active star-forming field galaxies for the smallest mass bin
; the rest of the number of galaxies is scaled based on this number
    ngal_active = 100000L   
    ngal_total = 0L 
    ngal_bin = 1L   ; just some arbitrary non-zero number
    for i=0L,n_elements(smf_data.mass)-1L do begin 
        if ((i EQ 0L) OR (ngal_bin NE 0)) then begin 
            bin_qf = get_qf(smf_data.mass[i], z)  ; Quiescent fraction of bin as a function of Mstar and z
            print, smf_data.mass[i]
            if (i EQ 0L) then begin
                ngal_quiescent = long((ngal_active*bin_qf)/(1.0-bin_qf))
                print, i, '     qf=', bin_qf, '     ngalQ=', ngal_quiescent, '      ngalSF=', ngal_active
            endif else begin 
                ngal_quiescent = long(ngal_bin*bin_qf)
                ngal_active = long(ngal_bin - ngal_quiescent)
                print, i, '     qf=', bin_qf, '     ngalQ=', ngal_quiescent, '      ngalSF=', ngal_active
            endelse 
; Active galaxies
            bin_mean_SFR = get_primus_sfr_mstar(smf_data.mass[i], z, sigma=bin_sigma_SFR, ssfr=bin_ssfr, silent=silent)
            bin_SFR = bin_sigma_SFR*randomn(rseed,ngal_active,/normal)+bin_mean_SFR   ; log-normal sampling of SFR based on observed PRIMUS mean SFR and sigma SFR
            bin_active_mass = replicate(smf_data.mass[i],ngal_active)  ;massbin[i].masssigma*randomn(rseed,ngal_active,/normal)+massbin[i].massbin
            bin_active_SSFR = bin_SFR-bin_active_mass
; Quiescent galaxies   
            if (ngal_quiescent NE 0L) then bin_quiescent_mass = replicate(smf_data.mass[i], ngal_quiescent) $  ;0.15*randomn(rseed,ngal_quiescent,/normal)+massbin[i].massbin
                else bin_quiescent_mass = [] 
            if (ngal_quiescent NE 0L) then bin_quiescent_SSFR = 0.2*randomn(rseed,ngal_quiescent,/normal)-12.0 $    ; log-normal distribution of quiescent SSFR with mean ~10^-12 and 0.2 dex scatter
                else bin_quiescent_SSFR = []     ; if statement included because quiescent fraction is effectively 0 for low masses
            
            total_mass = [total_mass, bin_active_mass, bin_quiescent_mass]
            total_SSFR = [total_SSFR, bin_active_SSFR, bin_quiescent_SSFR]
            ngal_total = ngal_total + ngal_active + ngal_quiescent
        endif 
        
        if (i LT n_elements(smf_data.mass)-1L) then begin 
            if (i EQ 0) OR (ngal_bin NE 0) then bin_smf_phi = smf_data.phi[i]
            next_bin_smf_phi = smf_data.phi[i+1L]
                ; bin_smf_phi = smf_data[where((smf_data.mass GE massbin[i].masslo) and (smf_data.mass LT massbin[i].massup),smf_bin_match)].phi 
                ; next_bin_smf_phi = smf_data[where((smf_data.mass GE massbin[i+1L].masslo) and (smf_data.mass LT massbin[i+1L].massup),smf_next_bin_match)].phi 
                ; if (smf_bin_match NE 1) then STOP 
            ngal_bin = ngal_active+ngal_quiescent
            if (next_bin_smf_phi NE 0) then ngal_bin = ngal_bin*(next_bin_smf_phi/bin_smf_phi)  ; scale the next bin ngal by SMF phi ratio
        endif 
        ;print, 'i =', i, '      N_quiescent = ', ngal_quiescent,'     N_active = ',ngal_active, '   ngal_total =', ngal_total, n_elements(total_mass)
    endfor

    output = replicate({mass:0.0, ssfr:0.0},ngal_total)
    output.mass = total_mass
    output.ssfr = total_SSFR
    output_file = 'ssfr_dist_z'+strmid(strtrim(string(z),2),0,4)+'_mbin'+strmid(strtrim(string(mbin),2),0,4)+'.fits' 
    output_dir  = '/global/data/scr/chh327/tinker/group_catalog/ssfr_dist/'
    mwrfits, output, output_dir+output_file, /create
end
