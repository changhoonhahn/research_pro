function mf_select_uvj_quiescent, umv, vmj, z=z, box=box, active=active
; jm12jan04ucsd - select quiescent and actively star-forming galaxies 
;   using the U-V vs V-J color-color plot

;   v2ab = k_vega2ab(filterlist=['bessell_U','bessell_V','twomass_J']+'.par',/silent,/kurucz)
;   box = {umv_const: 1.3+(v2ab[0]-v2ab[1]), vmj_turn: 1.6+(v2ab[1]-v2ab[2]), vmj_slope: 0.88}

;   box = {umv_const: 0.5, vmj_turn: 1.6, vmj_slope: 0.88}
    box = {umv_const: 1.3, vmj_turn: 0.69, vmj_red: 1.6, vmj_slope: 0.88}
    
    ngal = n_elements(vmj)
    if (ngal eq 0L) then return, -1 ; just return BOX
    if (n_elements(umv) ne ngal) then begin
       splog, 'Dimensions of VMJ and UMV do not agree!'
       return, -1
    endif

    all = lindgen(ngal)
    iquiescent = where((umv gt poly(vmj,[box.vmj_turn,box.vmj_slope])>box.umv_const) and $
      (vmj lt box.vmj_red))

    iquiescent = where((umv gt poly(vmj-box.vmj_turn,$
      [box.umv_const,box.vmj_slope])>box.umv_const) and $
      (vmj lt box.vmj_red))
    iactive = cmset_op(all,'and',/not2,iquiescent)
    
;   iactive = where((umv gt poly(vmj-box.vmj_turn,$
;     [box.umv_const,box.vmj_slope])>box.umv_const) and $
;     (vmj lt box.vmj_red))
    
    if keyword_set(active) then indx = iactive else indx = iquiescent

;   djs_plot, vmj, umv, psym=3, xsty=3, ysty=3, xr=[0,3], yr=[0,3]
;   djs_oplot, vmj[indx], umv[indx], psym=3, color='red'

;   vmjaxis = range(0.0,box.vmj_red,50)
;   djs_oplot, vmjaxis, poly(vmjaxis-box.vmj_turn,[box.umv_const,box.vmj_slope])>box.umv_const, color='red'
;   djs_oplot, box.vmj_red*[1,1], [poly(box.vmj_red-box.vmj_turn,[box.umv_const,box.vmj_slope]),!y.crange[1]], color='red'
    
return, indx
end


;function mf_select_uvj_quiescent, rmj, nuvmr, box=box, active=active
;; jm11sep01ucsd - select quiescent and actively star-forming galaxies 
;;   using the NUV-r vs r-J color-color plot
;
;; define the boundaries of the selection box (see
;; mfplot_select_quiescent)     
;    box = {nuvmr_const: 4.0, rmj_turn: 1.2, rmj_slope: 7.5}
;    
;    ngal = n_elements(rmj)
;    if (ngal eq 0L) then return, -1 ; just return BOX
;    if (n_elements(nuvmr) ne ngal) then begin
;       splog, 'Dimensions of RMJ and NUVMR do not agree!'
;       return, -1
;    endif
;
;    iquiescent = where(nuvmr gt poly(rmj-box.rmj_turn,$
;      [box.nuvmr_const,box.rmj_slope])>box.nuvmr_const,comp=iactive)
;    if keyword_set(active) then indx = iactive else $
;      indx = iquiescent
;    
;return, indx
;end
