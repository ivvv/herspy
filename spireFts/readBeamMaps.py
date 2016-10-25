# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 14:43:03 2016

@author: ivaltchanov
"""

from astropy.io import fits
import astropy.units as u
from spectral_cube import SpectralCube


file = "/Users/ivaltchanov/slw-beam-profile.fits"
#xx = fits.open(file)
xx = SpectralCube.read(file,hdu=1)
#file = "/Users/ivaltchanov/ssw-beam-profile.fits.gz"
#yy = fits.open(file)
