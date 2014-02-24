function mass_bins, nmassbin, minmass=minmass, maxmass=maxmass, $
    simple=simple, delta_mass=delta_mass
; Sets the mass bins to be used
    if keyword_set(simple) then begin 
        mmin = [ 9.5, 10.0, 10.5, 11.0] 
        mmax = [ 10.0, 10.5, 11.0, 11.5] 
    endif 
    if keyword_set(delta_mass) then begin 
        min = 9.5
        max = 11.5
        nbin = ceil((max-min)/delta_mass)
        mmin = min+delta_mass*float(bindgen(nbin)) 
        mmax = mmin+delta_mass 
    endif
    nmassbin = n_elements(mmin)
    massbins = replicate({massbin:0.0, masslo:0.0 , massup: 0.0, masssigma:0.0}, nmassbin)

    minmass = min(mmin)
    maxmass = max(mmax)
    
    massbins.masslo = mmin 
    massbins.massup = mmax 
    massbins.massbin = (mmax-mmin)/2.0+mmin
    massbins.masssigma = (mmax-mmin)/2.0
return, massbins
end
