import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

#Open data
path = "../data/"
fitsname = "GHIGLS_DFN_Tb_32.fits"
hdu = fits.open(path+fitsname)
hdr = hdu[0].header
cube = hdu[0].data

fitsname = "GHIGLS_DFN_Tb_32_gauss_run_0_gfort.fits"
hdu = fits.open(path+fitsname)
gfort = hdu[0].data

fitsname = "GHIGLS_DFN_Tb_32_gauss_run_0_ifort.fits"
hdu = fits.open(path+fitsname)
ifort = hdu[0].data

