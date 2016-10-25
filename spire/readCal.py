# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 15:59:54 2016

@author: ivaltchanov
"""

import zipfile
import numpy as np
from astropy.io import fits
import os
from astropy.wcs import WCS
import math

def bin_ndarray(ndarray, new_shape, operation='sum'):
    """
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.
    Number of output dimensions must match number of input dimensions.
    Example
    -------
    >>> m = np.arange(0,100,1).reshape((10,10))
    >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    >>> print(n)
    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]
    """
    if not operation.lower() in ['sum', 'mean', 'average', 'avg']:
        raise ValueError("Operation {} not supported.".format(operation))
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d, c in zip(new_shape,
                                                   ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        if operation.lower() == "sum":
            ndarray = ndarray.sum(-1*(i+1))
        elif operation.lower() in ["mean", "average", "avg"]:
            ndarray = ndarray.mean(-1*(i+1))
    return ndarray
    
#
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

jarFile = '/Users/ivaltchanov/SPCAL/14_3/spire_cal_14_3.jar'
#
xfwhm = {'PSW': 18.9, 'PMW': 25.8, 'PLW': 38.3}
yfwhm = {'PSW': 18.0, 'PMW': 24.6, 'PLW': 35.2}
coef = math.sqrt(8.0*math.log(2.0))
#
for arr in ["PSW","PMW","PLW"]:
    xsig = xfwhm[arr]/coef
    ysig = yfwhm[arr]/coef
    # 2D GAussian integral to infinity
    garray = 2.0*math.pi*xsig*ysig
    xBeam = readPhotBeam(jarFile,array=arr,variant='fine')
    xdelt = (xBeam[0].header["CDELT2"]*3600.0)**2 # arcsec^2
    yBeam = readPhotBeam(jarFile,array=arr,variant='nominal')
    ydelt = (yBeam[0].header["CDELT2"]*3600.0)**2 # arcsec^2
    area1 = xdelt*np.sum(xBeam[1].data)
    area6 = ydelt*np.sum(yBeam[1].data)
    print (arr,area1,area6 ,area6/area1,garray,garray/area6 )
