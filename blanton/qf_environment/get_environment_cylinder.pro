;+
; NAME:
;   get_environment_cylinder
;
; PURPOSE:
;   Finds the number of objects in the environmental sample that is within a cylinder of specified 
;   radius and height centered around each of the objects in the target sample using matchnd.
;
; CALLING SEQUENCE:
;   get_environment_cylinder(envfile=envfile, target=target)0	
;
; INPUTS:
;   envfile - .fits file name with environmental objects
;   target - structure with (x, y, z) coordinates of target objects
;
; OUTPUTS:
;   nmatch - number of environment objects inside the cylinder centered around each of the target objects
;
; FUNCTIONS USED: 
;   comdis(z, OmegaM, OmegaL, weq=weq)
;
; MODIFICATION HISTORY:
; Oct 8 2011: Chang Hoon Hahn
; Jan 2 2012: Chang Hoon Hahn
; May 25 2012: Chang Hoon Hahn
;--------------------------------------------------------

function get_environment_cylinder, envfile=envfile, targ=targ, rad=rad, h=h
;read in the environment sample file 
	envdata=mrdfits(envfile, 1) 
	
        r_targ = 3000.0*comdis(targ.redshift, 0.3, 0.7)
	phi_targ = targ.ra
	theta_targ = 90.0 - targ.dec
	
	r_env = 3000.0*comdis(envdata.redshift, 0.3, 0.7)
	phi_env = envdata.ra 
	theta_env = 90.0 - envdata.dec

	angles_to_xyz, r_targ, phi_targ, theta_targ, x_targ, y_targ, z_targ
	angles_to_xyz, r_env, phi_env, theta_env, x_env, y_env, z_env

	n1 = n_elements(r_targ)
	n2 = n_elements(r_env)

	print, 'target size =', n1,' environment size =', n2

	target = dblarr(3, n1)	
	target(0,*) = x_targ
	target(1,*) = y_targ
	target(2,*) = z_targ

	envpt = dblarr(3, n2)
	envpt(0,*) = x_env
	envpt(1,*) = y_env
	envpt(2,*) = z_env

;	target = radecz_to_xyz(ra=targ.ra, dec=targ.dec, redshift=targ.redshift)
;	envpt = radecz_to_xyz(ra=data.ra, dec=data.dec, redshift=data.zprimus)
;	print, 'target=',target
;	print, 'envpt=',envpt	
	
	cyl_radius = float(rad)
 	cyl_height = 0.5*float(h)
 	radius = sqrt(cyl_radius^2+cyl_height^2)
 	innerr = 0.0001
 	innerrsqr = innerr^2

;We use matchnd
;        if (n_elements(envpt) GT n_elements(target)) then begin
        matchnd, envpt, target, radius, m1=m1, m2=m2, maxmatch=0
;        endif else begin
;            matchnd, target, envpt, radius, m2=m2, m1=m1, maxmatch=0
;        endelse 
	dimtarg = size(target,/dimensions)
 	ntargets = dimtarg(1)
	print, 'ntargets', ntargets

 	isort   = sort(m2)
 	iuniq   = uniq(m2[isort])
	istart  =0L
	nmatch  = dblarr(ntargets)
 	for i=0L, n_elements(iuniq)-1L do begin
            iend    = iuniq[i]
            icurr   = isort[istart:iend]
;     	    nmatch[m2[icurr[0]]] = n_elements(icurr)
            nmatch[m2[icurr[0]]] = total(envdata[m1[icurr]].weight)
;            print, n_elements(icurr), total(envdata[m1[icurr]].weight)
            istart  = iend+1L

            targ_x  = target[0,m2[icurr[0]]]
            targ_y  = target[1,m2[icurr[0]]]
            targ_z  = target[2,m2[icurr[0]]]
            
            n_outside = 0L 
            for m=0L, n_elements(icurr)-1L do begin
;Project the environment point "vector" onto the target point "vector"
                env_x = envpt[0,m1[icurr[m]]]
                env_y = envpt[1,m1[icurr[m]]]
                env_z = envpt[2,m1[icurr[m]]]
                env_w = envdata[m1[icurr[m]]].weight
;                print, 'Weight of environment galaxy:',env_w

                distsqr=(targ_x-env_x)^2+(targ_y-env_y)^2+(targ_z-env_z)^2

                targ_mag= targ_x^2+targ_y^2+targ_z^2
                targ_dot_env= targ_x*env_x+targ_y*env_y+targ_z*env_z
                proj = targ_dot_env/targ_mag

                proj_x = proj*targ_x
                proj_y = proj*targ_y
                proj_z = proj*targ_z

; The distance between the target and the projection along the height of the cylinder
                proj_h = sqrt((proj_x-targ_x)^2+(proj_y-targ_y)^2+(proj_z-targ_z)^2)

; The distance between the target and the projection along the radius of the cyclinder
                proj_r = sqrt((env_x-proj_x)^2+(env_y-proj_y)^2+(env_z-proj_z)^2)

; If the point is within the inner sphere then subtract from nmatch since it's probably the same galaxy
; If outside the cylinder then subtract from nmatch 
                if (distsqr LT innerrsqr OR proj_h GT cyl_height OR proj_r GT cyl_radius) then begin
                    nmatch[m2[icurr[0]]]=nmatch[m2[icurr[0]]]-env_w
                    n_outside = n_outside+1L 
                endif
            endfor
;            print, 'galaxies outside cylinder', n_outside 
            if (n_outside EQ n_elements(icurr)) then nmatch[m2[icurr[0]]] = 0.0
    endfor
    return, nmatch
end
