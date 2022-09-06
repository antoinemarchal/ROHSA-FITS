import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

#Open data
path = "../data/"
fitsname = "GHIGLS_DFN_Tb.fits"
hdu = fits.open(path+fitsname)
hdr = hdu[0].header
cube = hdu[0].data[0][100:400,:32,:32]

hdu0 = fits.PrimaryHDU(cube, header=hdr)
hdulist = fits.HDUList([hdu0])
hdulist.writeto(path + "GHIGLS_DFN_Tb_32.fits", clobber=True)
