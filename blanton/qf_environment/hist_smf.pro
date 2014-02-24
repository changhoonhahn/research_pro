function hist_smf, array, vmax, xmin=xmin, xmax=xmax, weight=weight, binsize=binsize
    arrsize = n_elements(array)
    if (arrsize EQ 1) then begin
        print, "NOT AN ARRAY"
        stop
    endif 
    if n_elements(xmin) eq 0 then x_min=min(array) 
    if n_elements(xmin) ne 0 then x_min=xmin 
    if n_elements(xmax) eq 0 then x_max=max(array) 
    if n_elements(xmax) ne 0 then x_max=xmax
    xsize = CEIL((x_max-x_min)/float(binsize))+1L

    if xsize eq 0L then xsize = 1L 
    if xsize eq 1L then xsize = 2L

    xarr    = dblarr(xsize)
    xarr[0] = x_min
    for i=0L,xsize-1L do xarr[i] = xarr[0]+float(i)*binsize

    arr     = {mass:0., vmax:0., weight:0.}
    arrays  = replicate(arr, arrsize)
    arrays.mass = array
    arrays.vmax = vmax
    if n_elements(weight) eq 0 then begin
        arrays.weight = 1.0
    endif else begin
        arrays.weight = weight
    endelse 

    hist    = dblarr(xsize-1L)
    hist_err= dblarr(xsize-1L)
	
    for i=0L, xsize-2L do begin
        if i eq xsize-2L then begin
            indx = where(arrays.mass ge xarr[i] and arrays.mass le xarr[i+1L],indxcount)
        endif else begin 
            indx = where(arrays.mass ge xarr[i] and arrays.mass lt xarr[i+1L],indxcount)
        endelse 
        bin = arrays[indx] 
        if (indxcount EQ 0L) then begin         
            w = 0.0
        endif else begin 
            w = float(bin.weight)/float(bin.vmax)
        endelse 
        hist[i] = total(w)
        hist_err[i] = sqrt(total(w^2))
    endfor

    x = xarr[0:xsize-2L]+0.5*float(binsize)

    output = dblarr(3, xsize-1L)
    output[0,0:(xsize-2L)] = x
    output[1,0:(xsize-2L)] = hist
    output[2,0:(xsize-2L)] = hist_err
    return, output
end
