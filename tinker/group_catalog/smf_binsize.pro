function smf_binsize, bin=bin, minmass=minmass, maxmass=maxmass
; Obtains the binsize for build_smf_primus SMFs
    if n_elements(bin) eq 0 then binsize = 0.1 else binsize = bin
    minmass = 9.0 ; shift by a half-bin
    maxmass = 12.0
return, binsize
end
