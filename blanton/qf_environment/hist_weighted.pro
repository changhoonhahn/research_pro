function hist_weighted, array, xmin=xmin, xmax=xmax, weight=weight, binsize=binsize
	arrsize = n_elements(array)
    if n_elements(xmin) eq 0 then x_min=min(array) else x_min=xmin
    if n_elements(xmax) eq 0 then x_max=max(array) else x_max=xmax 
    xsize = CEIL((x_max-x_min)/float(binsize))+1L

    if xsize eq 0L then xsize = 1L 
    if xsize eq 1L then xsize = 2L

    xarr = dblarr(xsize)
	xarr[0] = x_min
	for i=0L,xsize-1L do begin
		xarr[i]=xarr[0]+float(i)*binsize
	endfor

	arrays = replicate({arr:0., weight:0.}, arrsize)
	arrays.arr = array
    if n_elements(weight) eq 0 then begin
        arrays.weight = 1.0
    endif else begin
        arrays.weight = weight
    endelse 

    hist = dblarr(xsize-1L)
	hist_err = dblarr(xsize-1L)
	
    for i=0L, xsize-2L do begin
        if i eq xsize-2L then begin
            indx = where(arrays.arr ge xarr[i] and arrays.arr le xarr[i+1L])
        endif else begin 
    		indx = where(arrays.arr ge xarr[i] and arrays.arr lt xarr[i+1L])
        endelse

		bin = arrays[indx] 	
        if (n_elements(bin) eq 1L) then begin
            if (bin.arr lt xarr[i] or bin.arr gt xarr[i+1L]) then w = 0.0
            if (bin.arr ge xarr[i] and bin.arr le xarr[i+1L]) then w = bin.weight
        endif else begin 
	    	w = bin.weight
        endelse 

        hist[i] = total(w)
        hist_err[i] = sqrt(total(w^2))
	endfor

	if arrsize eq 1L then begin
        x = xarr[0]
    endif else begin
        x = xarr[0:xsize-2L]+0.5*float(binsize)
    endelse 

	output = dblarr( 3, xsize-1L)
	output[0,0:(xsize-2L)] = x
	output[1,0:(xsize-2L)] = hist
	output[2,0:(xsize-2L)] = hist_err
    return, output
end
