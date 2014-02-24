pro ransack_dat2fits, fname
    readcol,fname,ra,dec,rand
    
    ransack=replicate({ra:0.,dec:0.},n_elements(ra))
    ransack.ra=ra
    ransack.dec=dec

    mwrfits,ransack,strmid(fname,0,strpos(fname,'.dat'))+'.fits',/create
end
