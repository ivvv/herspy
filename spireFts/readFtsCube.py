# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 15:00:02 2016

@author: ivaltchanov
"""

from astropy.io import fits
from astropy.wcs import WCS

import astropy.units as u
from spectral_cube import SpectralCube


#file = "/Users/ivaltchanov/slw-beam-profile.fits"
#file = "/Users/ivaltchanov/Tmp/HerschelData/ivaltcha255775486/1342204895/level2/HR_SLW_cube_convol/hspirespectrometer1342204895_spg_SLW_convol_HR_20ssc_1461253623981.fits.gz"
file = "/Users/ivaltchanov/Tmp/HerschelData/pacsRangeSpec/hpacs1342186983_20hps3deqibs_00_1452690668125.fits.gz"
xx = fits.open(file)
xx['image'].header
cube = xx['image'].data
header = xx['image'].header
if ('PACS' in xx[0].header["INSTRUME"]):
    header["CTYPE3"] = "WAVE"
elif ('SPIRE' in xx[0].header["INSTRUME"]):
    header["CTYPE3"] = "FREQ"
#

yy = SpectralCube(data=cube,wcs=WCS(header))
#
yy[10,:,:].quicklook()

#file = "/Users/ivaltchanov/ssw-beam-profile.fits.gz"
#yy = fits.open(file)
