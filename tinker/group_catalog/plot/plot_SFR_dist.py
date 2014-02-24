import numpy as np
import pylab as py
import pyfits as fits
from matplotlib import rc
import yanny as yan
import sys

rc('text', usetex=True)
rc('font', family='serif')

mr_cut = str(int(sys.argv[1]))
mass_cut = str(sys.argv[2])

file_dir = '/global/data/scr/chh327/tinker/group_catalog/'
file_name = 'clf_groups_M'+mr_cut+'_'+mass_cut+'_D360.galdata_corr.fits'

file_data = fits.open(file_dir+file_name)
SFR = file_data[1].data.field('SFR')

fig1 = py.figure(1)
sfrhist = fig1.add_subplot(111)
sfr_hist = sfrhist.hist(SFR, normed=True)

py.show()
