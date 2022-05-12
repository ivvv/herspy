#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 14:14:25 2017

@author: ivaltchanov
"""

from astropy.io import fits
from astropy.visualization import *

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm

beam_url = "http://archives.esac.esa.int/hsa/legacy/ADP/PSF/SPIRE/SPIRE-P"
#
pixes = {"PSW": "6", "PMW": "10", "PLW": "14"}
waves = {"PSW": "250", "PMW": "350", "PLW": "500"}

ax = plt.figure(figsize=(10,10))

url2 = "%s/0x5000241aL_PSW_bgmod10_6arcsec.fits.gz"%(beam_url)
hdu2 = fits.open(url2)
norm2 = ImageNormalize(hdu2['image'].data,interval=PercentileInterval(99.5),stretch=AsinhStretch())
#norm2 = ImageNormalize(hdu2['image'].data,interval=PercentileInterval(99),stretch=ContrastBiasStretch(0.5,0.5))
plt.imshow(hdu2['image'].data, cmap=cm.gray_r, norm=norm2, origin="lower")
#plt.get_xaxis().set_ticks([])
#plt.get_yaxis().set_ticks([])    
#plt.savefig('test_beams.png',dpi=300)
plt.show()
