function get_massbin, mass, masslo=masslo, massup=massup,$
    simple=simple, delta_mass=delta_mass, silent=silent
; Finds the mass bins that the mass belongs to
    massbins = mass_bins(nmassbins,minmass=minmass,maxmass=maxmass,simple=simple,delta_mass=delta_mass)
    if (keyword_set(silent) EQ 0) then struct_print, massbins
    
    i_match = where((round(10000.0*massbins.masslo) LE round(10000.0*mass)) AND $
        (round(10000.0*massbins.massup) GT round(10000.0*mass)),n_match)    ; to deal with pesky precision issues
    if (n_match NE 1) then begin
        print, mass, "OUT OF RANGE"  
        STOP     ; Sanity check to ensure that there's not more than one match
    endif

    masslo = massbins[i_match].masslo
    massup = massbins[i_match].massup
    return, massbins[i_match].massbin
end 
