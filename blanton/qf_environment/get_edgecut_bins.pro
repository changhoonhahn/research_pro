function get_edgecut,ran_ra,ran_dec, target,zmin=zmin,zmax=zmax,rad=rad,nbin=nbin,thresh=thresh,$
    primus=primus,sdss=sdss,crude=crude
    if (n_elements(ran_ra) eq 0 OR n_elements(ran_dec) EQ 0 OR n_elements(target) EQ 0) then print, 'INPUT ERROR'

    n_ransack = n_elements(ran_ra)
    print, 'Nransack=', n_ransack, ' R_aperture=',rad

; Divide redshift range into bins: 
    b=dblarr(nbin+1L)
    z0=zmin
    for j=0L,nbin do begin
        b[j]=z0
        z0=z0+(zmax-zmin)/float(nbin)
    endfor

    outputs=replicate({edgecut:0L}, n_elements(target))

; Angular separation of the matching radius defined at the lower redshift of the bin in degrees
; Lower redshift limit of the bin is selected so that the angular radius is larger 
    r_angsep=dblarr(nbin)
    for k=0L,nbin-1L do begin
        r_angsep[k]=(float(rad)/(3000.0*angdidis(b[k],0.3,0.7)))*(180.0/!PI)
    endfor

;Total area of all the polygons for the survey
    if keyword_set(primus) then totalarea = get_poly_area(/primus)
    if keyword_set(sdss) then totalarea = get_poly_area(/sdss,/edp)

    for i=0L,nbin-1 do begin
        
    	zbinindx = where(target.redshift ge b[i] and target.redshift lt b[i+1], zbincount)
     	bin = target[zbinindx]
        print, 'z_mid = ',b[i],'  bin dimensions = ',zbincount,'    r_angsep = ',r_angsep[i] 

        if keyword_set(crude) then bin_r_angsep = r_angsep[i] $
            else bin_r_angsep = (float(rad)/(3000.0*angdidis(bin.redshift, 0.3, 0.7)))*(180.0/!PI)
     	
        spherematch, ran_ra, ran_dec, bin.ra, bin.dec, r_angsep[i], m_ran, m_targ,$
            bindist12, maxmatch=0
    
        if keyword_set(crude) then begin 
            mkeep = m_targ
        endif else begin 
            r_angsep_m_targ = bin_r_angsep[m_targ]
            print, 'n[r_angsepz] = ', n_elements(r_angsep_m_targ), '  n[m_targ] = ', n_elements(m_targ), n_elements(bindist12)
            print, 'min[r_angsep_m_targ] = ', min(r_angsep_m_targ), '   max[r_angsep_m_targ] = ', max(r_angsep_m_targ)
            ikeep = where(bindist12 lt r_angsep_m_targ, nkeep)
            print, 'number of matches to keep=', nkeep
            mkeep = m_targ[ikeep]
        endelse 

        isort   = sort(mkeep)
        sorted  = mkeep[sort(mkeep)]
        iuniq   = uniq(mkeep[isort])

        nmatch  = lonarr(zbincount)
        threshold = fltarr(zbincount)
 
        istart  = 0L
        for m=0L, n_elements(iuniq)-1L do begin
            iend    = iuniq[m]
            icurr   = isort[istart:iend]
            nmatch[mkeep[icurr[0]]] = n_elements(icurr)
            if (m EQ 0L AND i EQ 0L) then begin
                mwrfits, bin[mkeep[icurr[0]]], 'sdss_target_singlegalaxy.fits', /create
                around_gal = replicate({ra:0.,dec:0.},n_elements(icurr)) 
                if keyword_set(crude) then begin 
                    around_gal.ra = ran_ra[m_ran[icurr]]
                    around_gal.dec = ran_dec[m_ran[icurr]]
                endif else begin
                    around_gal.ra = ran_ra[m_ran[ikeep[icurr]]]
                    around_gal.dec = ran_dec[m_ran[ikeep[icurr]]]
                endelse 
                mwrfits, around_gal, 'sdss_ransack_aroundgalaxy.fits', /create
            endif
            istart  = iend+1L
        endfor 
        
        threshold = (float(n_ransack)/totalarea)*(bin_r_angsep^2*!PI)*thresh

; If nmatch < threshold then point is on the edge: edgecut = 0 
; If nmatch > threshold then point is NOT on the edge: edgecut = 1 
        edgecut=lonarr(zbincount)
        edgecut[where(nmatch GE threshold)] = 1

        print, 'n_threshold=',n_elements(threshold),'   mean threshold=', mean(threshold)
        print, 'mean edgecut=',mean(edgecut)

        outputs[zbinindx].edgecut = edgecut
    endfor
    return, outputs
end 
; We determine the (number of random points within r_angsep)/(expected number of random points): 
;        randense = (float(nmatch)/(r_angsep[i]^2*!PI))/(float(ransack_num)/totalarea)
