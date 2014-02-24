function get_primus_sfr_mstar, mass, z, sigma=sigma, ssfr=ssfr, $
    silent=silent
; Obtains the mean and sigma of SFR for specific z and mass ranges.
; Also computes the SSFR for the mid-mass of the bin
; Determine the mass range and z range: 
    massbin = get_massbin(mass, masslo=masslo, massup=massup, delta_mass=0.1, silent=silent)
    zbin = get_zbin(z, zlo=zlo, zup=zup, /primus)

; Import PRIMUS active galaxies 
    primus_envcount_dir = '/global/data/scr/chh327/primus/data/envcount/'
    primus_active_file = 'envcount_cylr2h25_thresh75_active_lit_EDP-primus-z0210-numden.fits'
    primus_active = mrdfits(primus_envcount_dir+primus_active_file,1,/silent)

; Impose mass and z ranges and select only field galaxies:
    filter = where((primus_active.mass GT masslo) and (primus_active.mass LT massup) and $
        (primus_active.redshift GT zlo) and (primus_active.redshift LT zup), n_filter)
    ; and (primus_active.envcount EQ 0), n_filter) 
    ;print, n_filter
    primus_active_bin = primus_active[filter]

    SFR = mean(primus_active_bin.SFR)
    sigma = stddev(primus_active_bin.SFR)
    ssfr = SFR-massbin
;    print, SFR, sigma, ssfr
return, SFR
end
