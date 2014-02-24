function get_edgecut,ran_ra,ran_dec, target,zmin=zmin,zmax=zmax,rad=rad,thresh=thresh,$
    primus=primus,sdss=sdss
    if (n_elements(ran_ra) eq 0 OR n_elements(ran_dec) EQ 0 OR n_elements(target) EQ 0) then print, 'INPUT ERROR'

    n_ransack = n_elements(ran_ra)
    print, 'Nransack=', n_ransack, ' R_aperture=',rad

    outputs=replicate({edgecut:0L}, n_elements(target))

; Angular separation of the matching radius defined at z_min in degrees
    r_angsep_max = (float(rad)/(3000.0*angdidis(zmin,0.3,0.7)))*(180.0/!PI)

;Total area of all the polygons for the survey
    if keyword_set(primus) then totalarea = get_poly_area(/primus)
    if keyword_set(sdss) then totalarea = get_poly_area(/sdss,/edp)

; Spherematch cannot handle running on the entire PRIMUS sample so it is divded into bins 
    if keyword_set(primus) then nbin = ceil(float(n_elements(target))/25000.0)
    if keyword_set(sdss) then nbin = 1L 
    print, 'Number of bins in edgecut', nbin
    bin_z = fltarr(nbin+1L) 
    bin_z[0] = zmin
    for j=1L,nbin do begin 
        bin_z[j]=bin_z[j-1]+float(zmax-zmin)/float(nbin)
    endfor 

    for i=0L, nbin-1L do begin
        bin_index = where(target.redshift GE bin_z[i] AND target.redshift LT bin_z[i+1], bin_count)
        print, 'Number of galaxies in bin #'+strtrim(string(i+1L),2)+' = ', bin_count
        bin = target[bin_index]

        r_angsep = (float(rad)/(3000.0*angdidis(bin.redshift, 0.3, 0.7)))*(180.0/!PI)
     	
        spherematch, ran_ra, ran_dec, bin.ra, bin.dec, r_angsep_max, m_ran, m_targ,$
            bindist12, maxmatch=0
    
        r_angsep_m_targ = r_angsep[m_targ]
        print, 'n[r_angsepz] = ', n_elements(r_angsep_m_targ), '  n[m_targ] = ', n_elements(m_targ), n_elements(bindist12)
        print, 'min[r_angsep_m_targ] = ', min(r_angsep_m_targ), '   max[r_angsep_m_targ] = ', max(r_angsep_m_targ)
        ikeep = where(bindist12 lt r_angsep_m_targ, nkeep)
        print, 'number of matches to keep=', nkeep
    
        mkeep = m_targ[ikeep]
        isort   = sort(mkeep)
        sorted  = mkeep[sort(mkeep)]
        iuniq   = uniq(mkeep[isort])

        nmatch  = lonarr(n_elements(bin))
        threshold = fltarr(n_elements(bin))

        istart  = 0L
        for m=0L, n_elements(iuniq)-1L do begin
            iend    = iuniq[m]
            icurr   = isort[istart:iend]
            nmatch[mkeep[icurr[0]]] = n_elements(icurr)
            if (m EQ 0L) then begin
                mwrfits, bin[mkeep[icurr[0]]], 'sdss_bin_singlegalaxy.fits', /create
                around_gal = replicate({ra:0.,dec:0.},n_elements(icurr)) 
                around_gal.ra = ran_ra[m_ran[ikeep[icurr]]]
                around_gal.dec = ran_dec[m_ran[ikeep[icurr]]]
                mwrfits, around_gal, 'sdss_ransack_aroundgalaxy.fits', /create
            endif
            istart  = iend+1L
        endfor 
        
        threshold = (float(n_ransack)/totalarea)*(r_angsep^2*!PI)*thresh
        print, 'n_threshold=',n_elements(threshold),'   mean threshold=', mean(threshold)

; If nmatch < threshold then point is on the edge: edgecut = 0 
; If nmatch > threshold then point is NOT on the edge: edgecut = 1 
        bin_output = lonarr(n_elements(bin))
        bin_output[where(nmatch GE threshold)] = 1L
        outputs[bin_index].edgecut = bin_output 

    endfor
    print, 'mean edgecut=',mean(outputs.edgecut)
    return, outputs
end 
; We determine the (number of random points within r_angsep)/(expected number of random points): 
;        randense = (float(nmatch)/(r_angsep[i]^2*!PI))/(float(ransack_num)/totalarea)
