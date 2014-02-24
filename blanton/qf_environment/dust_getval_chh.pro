function dust_getval_chh, gall, galb, infile=infile, skipline=skipline, $
 outfile=outfile, map=map_in, interp=interp_in, noloop=noloop, $
 verbose=verbose, ipath=ipath, bhpath=bhpath

   if (NOT keyword_set(infile) $
    AND (N_elements(gall) EQ 0 AND N_elements(galb) EQ 0) $
    ) then begin
      print, 'Must set either coordinates or INFILE'
      return, -1
   endif

   bitnames = [ $
    [ '' , '' ], $
    [ '' , '' ], $
    [ 'OK     ' , 'asteroi' ], $
    [ 'OK     ' , 'glitch ' ], $
    [ 'OK     ' , 'source ' ], $
    [ 'OK     ' , 'no_list' ], $
    [ 'OK     ' , 'big_obj' ], $
    [ 'OK     ' , 'no_IRAS' ] ]

   ; Convert map name to upper case; default to Ebv
   if keyword_set(map_in) then map = strupcase(map_in) else map = 'EBV'

   if keyword_set(interp_in) then interp = interp_in
   if (map EQ 'MASK') then interp=0B

   ; If INFILE is defined, then read galactic coordinates from that file
   if (keyword_set(infile)) then $
    readcol, infile, gall, galb, format='F,F', skipline=skipline

   if keyword_set(ipath) then begin
       dum = findfile(ipath+'SFD*.fits', count=ct)
       if (ct EQ 0) then begin
           message, 'Bad file path!'
           return, -1
       endif
   endif
        
   if (NOT keyword_set(ipath)) then begin
      dust_dir = getenv('DUST_DIR')
      if (dust_dir NE '') then ipath = dust_dir+'/maps/' $
       else ipath = '/home/users/hahn/research/data/sfd-dust-maps/'
   endif
   dum = findfile(ipath+'SFD*.fits', count=ct)
   if (ct EQ 0) then begin
       message, 'No data files found in path'
       return, -1
   endif

   if (map EQ 'BH' AND NOT keyword_set(bhpath)) then begin
      bhpath = byte(1, 0) ? dust_dir+'/BHdat.i686/' : dust_dir+'/BHdat.sun4/'
   endif

   case strupcase(map) of
   'BH': begin
      value = bh_rdfort( gall, galb, bhpath=bhpath )
      end
   'I100': begin
      value = wcs_getval( $
       ['SFD_i100_4096_ngp.fits', 'SFD_i100_4096_sgp.fits'], $
       gall, galb, path=ipath, interp=interp, noloop=noloop, verbose=verbose)
      end
   'I60': begin
      value = wcs_getval( $
       ['SFD_i60_4096_ngp.fits', 'SFD_i60_4096_sgp.fits'], $
       gall, galb, path=ipath, interp=interp, noloop=noloop, verbose=verbose)
      end
   'I25': begin
      value = wcs_getval( $
       ['SFD_i25_4096_ngp.fits', 'SFD_i25_4096_sgp.fits'], $
       gall, galb, path=ipath, interp=interp, noloop=noloop, verbose=verbose)
      end
   'I12': begin
      value = wcs_getval( $
       ['SFD_i12_4096_ngp.fits', 'SFD_i12_4096_sgp.fits'], $
       gall, galb, path=ipath, interp=interp, noloop=noloop, verbose=verbose)
      end
   'X': begin
      value = wcs_getval( $
       ['SFD_xmap_ngp.fits', 'SFD_xmap_sgp.fits'], $
       gall, galb, path=ipath, interp=interp, noloop=noloop, verbose=verbose)
      end
   'T': begin
      value = wcs_getval( $
       ['SFD_temp_ngp.fits', 'SFD_temp_sgp.fits'], $
       gall, galb, path=ipath, interp=interp, noloop=noloop, verbose=verbose)
      end
   'IX': begin
      value = wcs_getval( $
       ['SFD_dust_4096_ngp.fits', 'SFD_dust_4096_sgp.fits'], $
       gall, galb, path=ipath, interp=interp, noloop=noloop, verbose=verbose)
      value = value / 0.0184
      end
   'EBV': begin
      value = wcs_getval( $
       ['SFD_dust_4096_ngp.fits', 'SFD_dust_4096_sgp.fits'], $
       gall, galb, path=ipath, interp=interp, noloop=noloop, verbose=verbose)
      end
   'MASK': begin
      value = wcs_getval( $
       ['SFD_mask_4096_ngp.fits', 'SFD_mask_4096_sgp.fits'], $
       gall, galb, path=ipath, interp=interp, noloop=noloop, verbose=verbose)
      end
   else: begin
      message, 'Valid map names: Ebv, BH, I100, X, T, IX, mask'
      end
   endcase

   ; If OUTFILE is defined, then write to output file
   if (keyword_set(outfile)) then begin
      get_lun, olun
      openw, olun, outfile
      if (map EQ 'MASK') then begin
         for igal=0, N_elements(gall)-1 do begin
            bits = djs_int2bin(value[igal],ndigit=8)
            printf, olun, format='(f8.3,f8.3,7a8)', $
             gall[igal], galb[igal], $
             strcompress(string(bits[0]+2*bits[1])+'hcons'), $
             bitnames[indgen(6)*2+4+bits[2:7]]
         endfor
      endif else begin
         for igal=0, N_elements(gall)-1 do begin
            printf, olun, format='(f8.3,f8.3,f12.5)', $
             gall[igal], galb[igal], value[igal]
         endfor
      endelse
      close, olun
      free_lun, olun
   endif

   return, value
end
;------------------------------------------------------------------------------
