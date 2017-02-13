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

f, ax = plt.subplots(2, 3)

for i,band in enumerate(["PSW","PMW","PLW"]):
    url1 = "%s/0x5000241aL_%s_bgmod10_1arcsec.fits.gz"%(beam_url,band)
    url2 = "%s/0x5000241aL_%s_bgmod10_%sarcsec.fits.gz"%(beam_url,band,pixes[band])
    hdu1 = fits.open(url1)
    hdu2 = fits.open(url2)
    norm1 = ImageNormalize(hdu1['image'].data,interval=PercentileInterval(99),stretch=ContrastBiasStretch(0.5,0.5))
    norm2 = ImageNormalize(hdu2['image'].data,interval=PercentileInterval(99),stretch=ContrastBiasStretch(0.5,0.5))
    #norm1 = ImageNormalize(hdu1['image'].data,interval=AsymmetricPercentileInterval())
    #norm2 = ImageNormalize(hdu2['image'].data,interval=AsymmetricPercentileInterval())
    #norm1 = ImageNormalize(hdu1['image'].data,interval=PercentileInterval(99.9),stretch=LogStretch())
    #norm2 = ImageNormalize(hdu2['image'].data,interval=PercentileInterval(99.9),stretch=LogStretch())
    ax[0,i].imshow(hdu1['image'].data, cmap=cm.gray_r, norm=norm1, origin="lower")
    ax[0,i].get_xaxis().set_ticks([])
    ax[0,i].get_yaxis().set_ticks([])    
    ax[0,i].annotate('%s µm (1\" pixels)'%waves[band], textcoords='axes fraction',fontsize=8, xy=(0.1,0.9),xytext=(0.1,0.9))
    ax[1,i].imshow(hdu2['image'].data, cmap=cm.gray_r, norm=norm2, origin="lower")
    ax[1,i].get_xaxis().set_ticks([])
    ax[1,i].get_yaxis().set_ticks([])
    ax[1,i].annotate('%s µm (%s\" pixels)'%(waves[band],pixes[band]), textcoords='axes fraction', fontsize=8, xy=(0.1,0.9),xytext=(0.1,0.9))
    break
#plt.savefig('test_beams.png',dpi=300)
plt.show()
