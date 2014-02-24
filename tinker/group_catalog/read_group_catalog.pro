pro read_group_catalog, mr_cut, mass_cut,groups=groups, prob=prob, galdata=galdata
    file_dir = '/global/data/scr/chh327/tinker/group_catalog/'
    filename = 'clf_groups_M'+strtrim(string(mr_cut),1)+'_'+strtrim(mass_cut,1)+'_D360'
    
    if keyword_set(prob) then begin 
        file_ext = '.prob'
        readcol, file_dir+filename+file_ext, v1, id_gal, id_group, id_cent, M_r, P_sat, $
            stellmass, v8, v9, v10, proj_r_rhalo, proj_r_rad, angrad_halo,$
            format='A,UL,UL,UL,F,F,F,F,F,F,F,F,F'

        output = replicate({id_gal:0, id_group:0, id_cent:0, M_r:0.0, P_sat:0.0, $
            stellmass:0.0, proj_r_rhalo:0.0, proj_r_rad:0.0, angrad_halo:0.0},n_elements(id_gal))
        output.id_gal       = id_gal 
        output.id_group     = id_group
        output.id_cent      = id_cent
        output.M_r          = M_r
        output.P_sat        = P_sat
        output.stellmass    = stellmass
        output.proj_r_rhalo = proj_r_rhalo
        output.proj_r_rad   = proj_r_rad
        output.angrad_halo  = angrad_halo

        output_name = file_dir+filename+file_ext+'.fits'
        print, output_name
        mwrfits, output, output_name, /create  
    endif

    if keyword_set(groups) then begin 
        file_ext = '.groups'
        readcol, file_dir+filename+file_ext, v1, id_group, id_cent, mass_group, v5, n_sat, $
            stellmass_group, stellmass_cent, v9, ra_cent, dec_cent, cz_cent, v13, v14, $
            format='A,UL,UL,F,F,F,F,F,F,F,F,F,F,F'

        output = replicate({id_group:0, id_cent:0, mass_group:0.0, n_sat:0, stellmass_group:0.0, $
            stellmass_cent:0.0, ra_cent:0.0, dec_cent:0.0, cz_cent:0.0},n_elements(id_group))
        output.id_group         = id_group
        output.id_cent          = id_cent
        output.mass_group       = mass_group
        output.n_sat            = n_sat
        output.stellmass_group  = stellmass_group 
        output.stellmass_cent   = stellmass_cent
        output.ra_cent          = ra_cent
        output.dec_cent         = dec_cent
        output.cz_cent          = cz_cent

        output_name = file_dir+filename+file_ext+'.fits'
        print, output_name
        mwrfits, output, output_name, /create  
    endif
    if keyword_set(galdata) then begin  
        file_ext = '.galdata_corr'
        readcol, file_dir+filename+file_ext, v1, id_gal, M_r, M_g, cz, Dn4000, SFR, H_delta_EW, $
            stellmass, ra, dec, v_disp, sn, n_sersic, $
            format='A,UL,F,F,F,F,F,F,F,F,F,F,F,F'

        output = replicate({id_gal:0, M_r:0.0, M_g:0.0, cz:0.0, Dn4000:0.0, SFR:0.0, H_delta_EW:0.0, $
            stellmass:0.0, ra:0.0, dec:0.0, v_disp:0.0, sn:0.0, n_sersic:0.0},n_elements(id_gal))
        output.id_gal       = id_gal
        output.M_r          = M_r
        output.M_g          = M_g
        output.cz           = cz
        output.Dn4000       = Dn4000
        output.SFR          = SFR
        output.H_delta_EW   = H_delta_EW
        output.stellmass    = stellmass
        output.ra           = ra
        output.dec          = dec
        output.v_disp       = v_disp
        output.sn           = sn
        output.n_sersic     = n_sersic
        
        output_name = file_dir+filename+file_ext+'.fits'
        print, output_name
        mwrfits, output, output_name, /create  
    endif 
end
