function get_percentage_rank, array
	arrsize = n_elements(array)
	p_rank = dblarr(arrsize)

	for i = 0L, arrsize-1L do begin 
		tmp = where(array le array[i], count)
		percent = float(count)/float(arrsize-1)*100.0
		p_rank[i] = percent
	endfor 
	return, p_rank
end
