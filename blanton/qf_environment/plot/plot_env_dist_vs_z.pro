pro plot_env_dist_vs_z, target=target, env=env
;Reading in target file:
	datum = mrdfits(target, 1)
	indx = where(datum.edgecut eq 1)
	data = datum[indx]
;Reading in environment file:
	envn = mrdfits(env, 1)

;Defining the bin indexes: 
	indx0 = where(data.redshift ge 0.2 AND data.redshift lt 0.4)
	indx1 = where(data.redshift ge 0.4 AND data.redshift lt 0.6)
	indx2 = where(data.redshift ge 0.6 AND data.redshift lt 0.8)
	indx3 = where(envn.zprimus ge 0.2 AND envn.zprimus lt 0.4)
	indx4 = where(envn.zprimus ge 0.4 AND envn.zprimus lt 0.6)
	indx5 = where(envn.zprimus ge 0.6 AND envn.zprimus lt 0.8)	
;Computing the fraction of the sky the polygons cover. 
	totarea = get_total_poly_area()/(129600.0/!PI)
	print, totarea

;Defining the target and environment bins	
	bin1 = data[indx0]
	bin2 = data[indx1]
	bin3 = data[indx2]
	
	env1 = envn[indx3]
	env2 = envn[indx4]
	env3 = envn[indx5]

	print, 'size of environment =', size(envn)
	print, 'size of target =', size(data)


;BIN 1: 
	print, '***************************BIN 1*********************************'
	maxcount1 = max(bin1.count)
	print, maxcount1
	hist1 = HISTOGRAM(bin1.count, BINSIZE = 1, min = 0, max=maxcount1)
	histsize1 = n_elements(hist1)	

;Fraction histogram:
	tot1 = float(total(hist1))
	hist1 = hist1/tot1

;Computing histogram total:
	histtot1 = 0.0
	for i = 0L, histsize1-1 do begin
		histtot1 = histtot1+hist1[i]*tot1*float(i)
	endfor 
	print, 'histtot=', histtot1

	encvol1 = (4.0*!PI/3.0)*(lf_comvol(0.4) - lf_comvol(0.2))*totarea

;Evaluating expected # in cylinder: (total number of environment galaxies in bin1)/(total enclosed volume in bin1)*(400*!PI)
	total1 = n_elements(env1)
	avg_n1= ((float(total1))/encvol1)*(!PI*50.0)

;X value of histogram:
	xval1 = DBLARR(histsize1)
	for i= 0L, histsize1-1L do begin
		xval1[i] = float(i)/avg_n1
	endfor 	

	print, 'TOTAL ENVIRONMENT GALAXIES IN BIN 1 =', total1
	print, 'TOTAL TARGET GALAXIES IN BIN 1 =', tot1
	print, 'TOTAL ENCLOSED VOLUME IN BIN 1 =', encvol1
	print, 'EXPECTED NUMBER OF GALAXIES IN CYLINDER =', avg_n1
; (HISTOGRAM TOTAL)/(TOTAL NUMBER OF TARGET GALAXIES IN BIN 1)
	print, 'AVERAGE ENVCOUNT FOR TARGET GALAXIES =', (float(histtot1))/(float(tot1))
	print, '(TOTAL # IN CYLINDERS)/(TOTAL # EXPECTED)', (float(histtot1))/(avg_n1*float(tot1))

;BIN 2: 
	print, '***************************BIN 2*********************************'
	maxcount2 = max(bin2.count)
	print, maxcount2
	hist2 = HISTOGRAM(bin2.count, BINSIZE = 1, min = 0, max=maxcount2)
	histsize2 = n_elements(hist2)	

;Fraction histogram:
	tot2 = total(hist2)
	hist2 = hist2/tot2

;Computing histogram total:
	histtot2 = 0.0
	for i = 0L, histsize2-1L do begin
		histtot2 = histtot2+hist2[i]*tot2*float(i)
	endfor 
	print, 'histtot=', histtot2

	encvol2 = (4.0*!PI/3.0)*(lf_comvol(0.6) - lf_comvol(0.4))*totarea

;Evaluating expected # in cylinder: (total number of environment galaxies in bin1)/(total enclosed volume in bin1)*(400*!PI)
	total2 = n_elements(env2)
	avg_n2= ((float(total2))/encvol2)*(!PI*50.0)

;X value of histogram:
	xval2 = DBLARR(histsize2)
	for i= 0L, histsize2-1L do begin
		xval2[i] = float(i)/avg_n2
	endfor 	

	print, 'TOTAL ENVIRONMENT GALAXIES IN BIN 1 =', total2
	print, 'TOTAL TARGET GALAXIES IN BIN 1 =', tot2
	print, 'TOTAL ENCLOSED VOLUME IN BIN 1 =', encvol2
	print, 'EXPECTED NUMBER OF GALAXIES IN CYLINDER =', avg_n2
; (HISTOGRAM TOTAL)/(TOTAL NUMBER OF TARGET GALAXIES IN BIN 1)
	print, 'AVERAGE ENVCOUNT FOR TARGET GALAXIES =', (float(histtot2))/(float(tot2))
	print, '(TOTAL # IN CYLINDERS)/(TOTAL # EXPECTED)', (float(histtot2))/(avg_n2*float(tot2))

;BIN 1: 
	print, '***************************BIN 3*********************************'
	maxcount3 = max(bin3.count)
	print, maxcount3
	hist3 = HISTOGRAM(bin3.count, BINSIZE = 1, min=0, max=maxcount3)
	histsize3 = n_elements(hist3)	

;Fraction histogram:
	tot3 = total(hist3)
	hist3 = hist3/tot3
	
;Computing histogram total:
	histtot3 = 0.0
	for i = 0L, histsize3-1L do begin
		histtot3 = histtot3+hist3[i]*tot3*float(i)
	endfor 
	print, 'histtot=', histtot3

	encvol3 = (4.0*!PI/3.0)*(lf_comvol(0.8) - lf_comvol(0.6))*totarea

;Evaluating expected # in cylinder: (total number of environment galaxies in bin1)/(total enclosed volume in bin1)*(400*!PI)
	total3 = n_elements(env3)
	avg_n3= ((float(total3))/encvol3)*(!PI*50.0)

;X value of histogram:
	xval3 = DBLARR(histsize3)
	for i= 0L, histsize3-1L do begin
		xval3[i] = float(i)/avg_n3
	endfor 	

	print, 'TOTAL ENVIRONMENT GALAXIES IN BIN 1 =', total3
	print, 'TOTAL TARGET GALAXIES IN BIN 1 =', tot3
	print, 'TOTAL ENCLOSED VOLUME IN BIN 1 =', encvol3
	print, 'EXPECTED NUMBER OF GALAXIES IN CYLINDER =', avg_n3
; (HISTOGRAM TOTAL)/(TOTAL NUMBER OF TARGET GALAXIES IN BIN 1)
	print, 'AVERAGE ENVCOUNT FOR TARGET GALAXIES =', (float(histtot3))/(float(tot3))
	print, '(TOTAL # IN CYLINDERS)/(TOTAL # EXPECTED)', (float(histtot3))/(avg_n3*float(tot3))


	print, 'TOTAL TARGET GALAXIES =', tot1+tot2+tot3
	print, histsize1, histsize2, histsize3
	print, 'HISTSIZE SUM =', histsize1+histsize2+histsize3
	print, 'TOTAL ENVIRONMENT GALAXIES =', total1+total2+total3
	print, 'Integral of distribution check:', total(hist1), total(hist2), total(hist3)
;	print, hist1
;	print, hist2
;	print, hist3
;	save, hist1, hist2, hist3, filename='hist.sav'

	dir = get_path(/envdist)
	namestring = strtrim(target,2)
	fname1 = dir+'plot_envdist_vs_z_'+namestring+'_lowz_hist.dat'
        fname2 = dir+'plot_envdist_vs_z_'+namestring+'_midz_hist.dat'
	fname3 = dir+'plot_envdist_vs_z_'+namestring+'_highz_hist.dat'
	fname4 = dir+'plot_envdist_vs_z_'+namestring+'_overplot_hist.dat'
	
	openw,lun,fname1,/get_lun
	for i=0L,n_elements(xval1)-1L do printf,lun,xval1[i],hist1[i],format='(f,f)'
	free_lun,lun
        openw,lun,fname2,/get_lun
        for i=0L,n_elements(xval2)-1L do printf,lun,xval2[i],hist2[i],format='(f,f)'
        free_lun,lun
        openw,lun,fname3,/get_lun
        for i=0L,n_elements(xval3)-1L do printf,lun,xval3[i],hist3[i],format='(f,f)'
        free_lun,lun
        openw,lun,fname4,/get_lun
        for i=0L,n_elements(hist3)-1L do printf,lun,hist1[i],hist2[i],hist3[i],format='(f,f,f)'
        free_lun,lun

	return
end
