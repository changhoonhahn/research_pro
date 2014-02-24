function get_mf_vagc_chh, sample=sample, letter=letter, $
  poststr=poststr, zminmax=zminmax
; jm10sep10ucsd - define the VAGC sample; ZMINMAX is needed by both
; MF_ISEDFIT and BUILD_MF_UBERSAMPLE

; POSTSTR 0 -30.   -10. 1.e-3 0.5   0.1 0.3 0.7 0.0  0.1  0.00 -10.0  10.0
;    zminmax = [0.01,0.2]
;   zminmax = [0.01,0.25]
;    zminmax = [0.06,0.145]
    zminmax = [0.0375,0.145]
    sample = 'dr72'
    letter = 'bsafe'
    poststr = '0'
    vagc_sample = sample+letter+poststr

;; POSTSTR/25: -17>Mr>-24; 0.01<z<0.25
;    zminmax = [0.01,0.25]
;    sample = 'dr72'
;    letter = 'bsafe'
;    poststr = '25'
;    vagc_sample = sample+letter+poststr
    
;; POSTSTR/35: -17>Mr>-24; 0.033<z<0.25
;    zminmax = [0.033,0.25]
;    sample = 'dr72'
;    letter = 'bsafe'
;    poststr = '35'
;    vagc_sample = sample+letter+poststr
    
return, vagc_sample
end
    
