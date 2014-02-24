pro build_env_sample
    dir = get_path(/envt)
    envdir      = '/global/data/scr/chh327/primus/env/'
    parentdir   = '/global/data/scr/chh327/primus/science/mf/2165/parent_v20/'

    fields = ['xmm_swire','cdfs','cfhtls_xmm','cosmos','es1']
    primus_ngal = 0L 
    for i=0L,n_elements(fields)-1L do begin 
        primus_field    = mrdfits(parentdir+'mf_parent_'+fields[i]+'.fits.gz',1)
        field_ngal      = n_elements(primus_field)
        if (i EQ 0L) then primus_tags = struct_selecttags(primus_field[0],$
            select_tags=['ra','dec','zprimus','zmin_evol','zmax_evol',$
            'K_FNUVUGRIZJHK_ABSMAG_01','targ_weight'])
        primus_parent   = replicate(primus_tags,primus_ngal+field_ngal) 

        if (i GT 0L) then begin 
            primus_parent[0L:primus_ngal-1L] = tmp_primus_parent
        endif 
        primus_parent[primus_ngal:primus_ngal+field_ngal-1L] = struct_selecttags(primus_field,$
            select_tags=['ra','dec','zprimus','zmin_evol','zmax_evol','K_FNUVUGRIZJHK_ABSMAG_01'$
            ,'targ_weight'])

        tmp_primus_parent = primus_parent
        primus_ngal = primus_ngal+field_ngal 
    endfor
    primus_parent = struct_addtags(temporary(primus_parent),replicate({vmax_evol:0.,comp:0L},primus_ngal)) 
    ; SDSS
    sdss_parent = mrdfits('/global/data/scr/chh327/primus/science/mf/2165/parent_v20/edp_parent/mf_parent_sdss.fits.gz',1)
    sdss_ngal   = n_elements(sdss_parent) 
    sdss_parent = struct_addtags(temporary(sdss_parent),replicate({vmax_evol:0.,comp:0L},sdss_ngal)) 

    print, ' PRIMUS Parent Sample = ', n_elements(primus_parent)
    print, ' SDSS Parent Sample = ', n_elements(sdss_parent)

    h100 = mf_h100()

    sdss_zbins = mf_zbins_chh(sdss_nzbins, zmin=sdss_zmin,zmax=sdss_zmax,/sdss)
    primus_zbins = mf_zbins_chh(primus_nzbins, zmin=primus_zmin,zmax=primus_zmax,/literature)
    
    redbin = [0.9,0.7,0.5,0.3,0.1] 
    for j=0L,n_elements(redbin)-1L do begin 
        if (redbin[j] EQ 0.1) then begin 
            bin_zmin    = sdss_zmin
            bin_zmax    = sdss_zmax
            zbin_indx   = sdss_parent.zprimus GE sdss_zmin AND sdss_parent.zprimus LT sdss_zmax
            bin         = sdss_parent[where(zbin_indx,zbin_count)]
            polyarea    = get_poly_area(/sdss, /edp, /sr) 
        endif else begin 
            bin_zmin = primus_zbins[primus_nzbins-j-1].zlo
            bin_zmax = primus_zbins[primus_nzbins-j-1].zup
            zbin_indx   = primus_parent.zprimus GE bin_zmin AND primus_parent.zprimus LT bin_zmax
            bin         = primus_parent[where(zbin_indx,zbin_count)]
            polyarea    = get_poly_area(/primus, /sr) 
        endelse 
        print, j, bin_zmin, bin_zmax, zbin_count
        zmin_evol = bin.zmin_evol>bin_zmin
        zmax_evol = bin.zmax_evol<bin_zmax
        bin.vmax_evol = (polyarea/3.0)*(lf_comvol(zmax_evol)-lf_comvol(zmin_evol))*(1.0/h100)^3.0 

        bin_r           = bin.K_FNUVUGRIZJHK_ABSMAG_01[4]
        bin_vmax_evol   = bin.vmax_evol
        if (redbin[j] EQ 0.1) then bin_weights = bin.final_weight else bin_weights = bin.targ_weight

        openw,lun,dir+'env_zbin'+strtrim(string(redbin[j]),1)+'.dat',/get_lun 
        for ii=0L,n_elements(bin)-1L do printf,lun,bin_r[ii],bin[ii].zprimus,bin_vmax_evol[ii]
        free_lun,lun

        bin_hist    = hist_smf(bin_r,bin_vmax_evol,weight=bin_weights, binsize=0.01)
        bin_cumul   = total(bin_hist[1,*], /cumulative)
        bin_x       = bin_hist[0,*]

        print, 'Bin Redshift range = [',bin_zmin,',',bin_zmax,']'
        print, 'Bin Number Density = ', max(bin_cumul)

        if (j eq 0) then begin
            bin_n = value_locate(bin_x, -20.75)
            high_x = bin_x
            high_n = bin_n 
            high_cumul = bin_cumul
            print, 'High Redshift Number Density limit=',high_cumul[high_n]
        endif else begin 
            bin_n = value_locate(bin_cumul, high_cumul[high_n])
        endelse 
        bin_absmag_lim = bin_x[bin_n+1L]

        print, ' Bin Absolute Magnitude Limit = ', bin_absmag_lim
        print, ' Bin Number Density Limit = ', bin_cumul[bin_n]

        if (redbin[j] EQ 0.1) then begin 
            absmag_lim_indx = sdss_parent.K_FNUVUGRIZJHK_ABSMAG_01[4] le bin_absmag_lim
            sdss_parent[where(zbin_indx AND absmag_lim_indx)].comp = 1L 
        endif else begin 
            absmag_lim_indx = primus_parent.K_FNUVUGRIZJHK_ABSMAG_01[4] le bin_absmag_lim 
            primus_parent[where(zbin_indx AND absmag_lim_indx)].comp = 1L
        endelse 
        bin_absmag_lim_indx = bin_r le bin_absmag_lim 
        bin = bin[where(bin_absmag_lim_indx)]
        openw,lun,dir+'env_zbin'+strtrim(string(redbin[j]),1)+'_lim.dat',/get_lun 
        for ii=0L,n_elements(bin)-1L do  begin 
            printf,lun,bin[ii].K_FNUVUGRIZJHK_ABSMAG_01[4],bin[ii].zprimus,bin[ii].vmax_evol
        endfor
        free_lun,lun
    endfor
;We impose the constraints for our sample above from above. 
    primus_complete_indx    = where(primus_parent.comp eq 1L,primus_completecount)
    sdss_complete_indx      = where(sdss_parent.comp eq 1L,sdss_completecount)
    print, 'PRIMUS Completeness:',float(primus_completecount)/float(n_elements(primus_parent))
    print, 'SDSS Completeness:',float(sdss_completecount)/float(n_elements(sdss_parent))
    primus_parent = primus_parent[primus_complete_indx]
    sdss_parent = sdss_parent[sdss_complete_indx] 

    outputs = replicate({ra:0.D,dec:0.D,redshift:0.,weight:0.,absmag:fltarr(10)}, $
        n_elements(primus_parent))
    outputs.ra      = primus_parent.ra
    outputs.dec     = primus_parent.dec
    outputs.redshift= primus_parent.zprimus
    outputs.weight  = primus_parent.targ_weight
    outputs.absmag  = primus_parent.K_FNUVUGRIZJHK_ABSMAG_01 
    mwrfits, outputs, dir+'EDP-primus-z0210-numden.fits', /create

    sdss_output = replicate({ra:0.D,dec:0.D,redshift:0.,weight:0.,absmag:fltarr(10)},$
        n_elements(sdss_parent))
    sdss_output.ra      = sdss_parent.ra
    sdss_output.dec     = sdss_parent.dec
    sdss_output.redshift= sdss_parent.zprimus
    sdss_output.weight  = sdss_parent.final_weight
    sdss_output.absmag  = sdss_parent.K_FNUVUGRIZJHK_ABSMAG_01
    mwrfits, sdss_output, dir+'EDP-sdss-z00375_0145-numden.fits', /create
end
