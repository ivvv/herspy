# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 14:10:53 2017

@author: ivaltchanov
"""
import os
import zipfile
from astropy.io import fits

def readPhotBeam(jarFile,array='PSW',variant='fine'):
    #
    #
    wdir = os.path.dirname(jarFile)
    beamFile = 'SCalPhotBeamProf_%s_%s_v6.fits'%(array,variant)
    #
    with zipfile.ZipFile(jarFile, 'r') as zf:
        lst = zf.infolist()
        for zi in lst:
            fn = zi.filename
            if (beamFile in fn):
                hdu= fits.open(wdir + '/' + fn)
    return hdu
#
