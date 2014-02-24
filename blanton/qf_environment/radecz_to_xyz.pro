;+
; NAME:
;	radecz_to_xyz
;
; PURPOSE:
;	Simple program to convert coordinates in RA, DEC, and redshift into x, y, and z coordinates
;
; CALLING SEQUENCE:
;	radecz_to_xyz(ra=ra, dec=dec, redshift=redshift)
;
; INPUTS:
;	ra - Array with RA values
;	dec - Array with DEC values
; 	redshift - Array with redshift values
;
; OUTPUTS:
;	output - An array with three columns which correspond to x, y, and z respectively.
; 
;---------------------------------------------------------------------
function radecz_to_xyz, ra=ra, dec=dec, redshift=redshift
	n = n_elements(ra) 

	r = 3000.0*comdis(redshift, 0.3, 0.7)
	ra = ra*!PI/180.0
	dec = dec*!PI/180.0

	theta = ra
	phi = 0.5*!PI - dec

	x = r*cos(theta)*sin(phi)
	y = r*sin(theta)*sin(phi)
	z = r*cos(phi)

	output = dblarr(3, n)
	output(0,*) = x
	output(1,*) = y
	output(2,*) = z 

	return, output
end
