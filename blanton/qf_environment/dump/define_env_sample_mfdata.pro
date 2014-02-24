;+
; NAME:
;   define_env_sample_smf_new
;
; PURPOSE:
;   Defines an environmental sample using the constraints specified in definesample_index.par
;   only includes objects in 2mask fields. Note that the constraints must be first specified 
;   in definesample_index.par before running this procedure. 
;
; CALLING SEQUENCE:
;   define_environment, name=
;
; INPUTS:
;   All inputs are specified in define_sample_index.par
;
; OUTPUTS:
;   .fits file with name = name
;
; FUNCTIONS USED:
;   yanny_readone( filename, [ selectname, hdr=, enums=, structs=, /anonymous, stnames=, 
;                  /quick, errcode= ])
;   is_in_window( xyz= , ra= , dec= , polygons)
;
; PROCEDURES USED: 
;   read_mangle_polygons, infile, polygons, id [, unit=]
;
; MODIFICATION HISTORY:
; Oct 8 2011: Chang Hoon Hahn
		; Jan 2 2012: Chang Hoon Hahn
; May 25 2012: Chang Hoon Hahn
;--------------------------------------------------------
pro define_env_sample_smf_new, name=name 
        dim = dblarr(5)

        data0 = mrdfits('/mount/moon1/ioannis/research/projects/primus/mf/2165/mfs_v20/mfdata_all_supergrid01_xmm_swire.fits.gz',1)
        dim[0]    = N_ELEMENTS(data0)

        data1 = mrdfits('/mount/moon1/ioannis/research/projects/primus/mf/2165/mfs_v20/mfdata_all_supergrid01_cdfs.fits.gz',1)
        dim[1]    = N_ELEMENTS(data1)

        data2 = mrdfits('/mount/moon1/ioannis/research/projects/primus/mf/2165/mfs_v20/mfdata_all_supergrid01_cfhtls_xmm.fits.gz',1)
        dim[2]    = N_ELEMENTS(data2)

        data3 = mrdfits('/mount/moon1/ioannis/research/projects/primus/mf/2165/mfs_v20/mfdata_all_supergrid01_cosmos.fits.gz',1)
        dim[3]    = N_ELEMENTS(data3)

        data4 = mrdfits('/mount/moon1/ioannis/research/projects/primus/mf/2165/mfs_v20/mfdata_all_supergrid01_es1.fits.gz',1)
        dim[4]    = N_ELEMENTS(data4)

        targ = {class:' ',ra:0.D,dec:0.D,redshift:0.,mr_01:0.,mb_00:0.,mass:0.,age:0.,SFR:0., comp:0L}
        dimension = total(dim)
        targs = replicate(targ, dimension)

        binstart = 0L
                binend = binstart + dim[0] - 1L
                targs[binstart:binend].ra = data0.ra
                targs[binstart:binend].dec = data0.dec
                targs[binstart:binend].redshift = data0.zprimus
                targs[binstart:binend].mr_01 = data0.mr_01
                targs[binstart:binend].mb_00 = data0.mb_00
                targs[binstart:binend].mass = data0.mass
                targs[binstart:binend].age = data0.age
                targs[binstart:binend].SFR = data0.SFR
                binstart = binstart + dim[0]

                binend = binstart + dim[1] - 1L
                targs[binstart:binend].ra = data1.ra
                targs[binstart:binend].dec = data1.dec
                targs[binstart:binend].redshift = data1.zprimus
                targs[binstart:binend].mr_01 = data1.mr_01
                targs[binstart:binend].mb_00 = data1.mb_00
                targs[binstart:binend].mass = data1.mass
                targs[binstart:binend].age = data1.age
                targs[binstart:binend].SFR = data1.SFR
                binstart = binstart + dim[1]

                binend = binstart + dim[2] - 1L
                targs[binstart:binend].ra = data2.ra
                targs[binstart:binend].dec = data2.dec
                targs[binstart:binend].redshift = data2.zprimus
                targs[binstart:binend].mr_01 = data2.mr_01
                targs[binstart:binend].mb_00 = data2.mb_00
                targs[binstart:binend].mass = data2.mass
                targs[binstart:binend].age = data2.age
                targs[binstart:binend].SFR = data2.SFR
                binstart = binstart + dim[2]

                binend = binstart + dim[3] - 1L
                targs[binstart:binend].ra = data3.ra
                targs[binstart:binend].dec = data3.dec
                targs[binstart:binend].redshift = data3.zprimus
                targs[binstart:binend].mr_01 = data3.mr_01
                targs[binstart:binend].mb_00 = data3.mb_00
                targs[binstart:binend].mass = data3.mass
                targs[binstart:binend].age = data3.age
                targs[binstart:binend].SFR = data3.SFR
                binstart = binstart + dim[3]

                binend = binstart + dim[4] - 1L
                targs[binstart:binend].ra = data4.ra
                targs[binstart:binend].dec = data4.dec
                targs[binstart:binend].redshift = data4.zprimus
                targs[binstart:binend].mr_01 = data4.mr_01
                targs[binstart:binend].mb_00 = data4.mb_00
                targs[binstart:binend].mass = data4.mass
                targs[binstart:binend].age = data4.age
                targs[binstart:binend].SFR = data4.SFR
                binstart = binstart + dim[4]


print, 'dimension =', dimension

;We are only interested in data with zconf [3,4] and class='GALAXY', so we impose "indx0"
	indx0 = where(targs.redshift lt 0.8 AND targs.redshift ge 0.2) 
	smfenv = targs[indx0]
	save, smfenv, filename='smf_env_galaxies.sav'

	bin1 = where(smfenv.redshift ge 0.2 AND smfenv.redshift lt 0.4)
	bin2 = where(smfenv.redshift ge 0.4 AND smfenv.redshift lt 0.6) 
	bin3 = where(smfenv.redshift ge 0.6 AND smfenv.redshift lt 0.8)

        high= smfenv[bin3]  
        high_hist = histogram(high.mr_01, binsize=0.01, min=-26.0)
        high_comvol = lf_comvol(0.8)-lf_comvol(0.6)
        high_cumul= float(total(high_hist, /cumulative))/float(high_comvol[0])
        high_x = dblarr(n_elements(high_hist))
        for i=0L,n_elements(high_hist)-1L do begin
                high_x[i] = -26.0+0.01*float(i)
        endfor
        high_n = value_locate(high_x,-20.75)
        high_lim = high_x[high_n+1L]
        print, '-20.75 array value=', value_locate(high_x, -20.75), high_cumul[high_n]

        mid = smfenv[bin2]
        mid_hist = histogram(mid.mr_01, binsize=0.01, min=-26.0)
        mid_comvol = lf_comvol(0.6)-lf_comvol(0.4)
	mid_cumul= float(total(mid_hist, /cumulative))/float(mid_comvol[0])
        mid_x = dblarr(n_elements(mid_hist))
	for i=0L,n_elements(mid_hist)-1L do mid_x[i] = -26.0+0.01*float(i)
	mid_n = value_locate(mid_cumul, high_cumul[high_n])
	mid_lim = mid_x[mid_n+1L]
;        mid_lim = -26.0+0.01*(value_locate(mid_cumul, complim)+1L)


	low = smfenv[bin1]
	low_hist = histogram(low.mr_01, binsize=0.01, min=-26.0) 
	low_comvol = lf_comvol(0.4)-lf_comvol(0.2)
        low_cumul= float(total(low_hist, /cumulative))/float(low_comvol[0])
	low_x = dblarr(n_elements(low_hist))
	for i=0L,n_elements(low_hist)-1L do low_x[i] = -26.0+0.01*float(i)
	low_n = value_locate(low_cumul, high_cumul[high_n])
	low_lim = low_x[low_n+1L]
;	low_lim = -26.0+0.01*(value_locate(low_cumul, complim)+1L)

	limit = dblarr(n_elements(high_x))
	limit[*] = high_cumul[high_n]

	set_plot, 'z'
	device, set_resolution=[640,680]
        plot, low_x, low_cumul, xrange=[-26.0, -18.0], linestyle=0, thick=3
	oplot, high_x, limit
	oplot, mid_x, mid_cumul, linestyle=2, thick=3
	oplot, high_x, high_cumul, linestyle=3, thick=3
	legend, linestyle=[0,2,3], ['Low Redshift', 'Mid Redshift', 'High Redshfit'], font=50
	write_png,'smf_numden_completeness.png', tvrd()

	print, n_elements(bin1), n_elements(bin2), n_elements(bin3)	
	print, low_n, mid_n, high_n
	print, low_comvol[0], mid_comvol[0], high_comvol[0]
	print, low_lim, mid_lim, high_lim
	print, low_cumul[low_n], mid_cumul[mid_n], high_cumul[high_n]
	
	low_bin = where(smfenv.redshift ge 0.2 AND smfenv.redshift lt 0.4 AND smfenv.mr_01 le low_lim) 
	mid_bin = where(smfenv.redshift ge 0.4 AND smfenv.redshift lt 0.6 AND smfenv.mr_01 le mid_lim)
	high_bin= where(smfenv.redshift ge 0.6 AND smfenv.redshift lt 0.8 AND smfenv.mr_01 le high_lim)
	
	targs.comp = 0L 
	smfenv[low_bin].comp = 1L 
	smfenv[mid_bin].comp = 1L
	smfenv[high_bin].comp= 1L 

;We impose the constraints for our sample above from above. 
 	complete = where(smfenv.comp eq 1L,completecount)
	print, 'Completeness:',float(completecount)/float(n_elements(smfenv))
	smf_env = smfenv[complete]

 	print, "length=", n_elements(smf_env)
	
	output = {ra:0.D,dec:0.D,zprimus:0.,mr_01:0.}
 	outputs = replicate(output, n_elements(smf_env))
 	outputs.ra = smf_env.ra
 	outputs.dec = smf_env.dec
 	outputs.zprimus = smf_env.redshift
	outputs.mr_01 = smf_env.mr_01 
	
	mwrfits, outputs, name, /create
	return
end
