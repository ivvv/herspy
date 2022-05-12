#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 11:13:27 2017

@author: ivaltchanov
"""
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from spectral_cube import SpectralCube
import tarfile
import gzip
import requests
import os
import io

import numpy as np

import matplotlib.pyplot as plt

tmpDir = os.path.expanduser('~') + '/Tmp/'


def getSpireFtsCube(obsid, apod = False, cache=True):
    """
    Using the HTTP access to HAIO, retrieve the level-2 products in a tar file
    and extract only the requested fits files

    cubeType can be 'convol' or 'naive'
    apod controls wheter to extract the apodized version (if True)
    
    simple caching is impemented, checking of a tar file with the same name (OBSID based) already exists
    at the default folder.
    
    """
    vers = '_spg_'
    if (apod): 
        vers = '_spgApod_'
    #
    tarFile = tmpDir + "%i_level2.tar"%obsid
    if (os.path.isfile(tarFile) and cache):
        print ("Found an already existing tar file for OBSID %i. Will Use it"%obsid)
    else:
        haioRequest = "http://archives.esac.esa.int/hsa/aio/jsp/product.jsp?PROTOCOL=HTTP&OBSERVATION_ID=%i&PRODUCT_LEVEL=Level2"%obsid
        print ("Downloading level-2 data from the Herschel Science Archive. May take a while... be patient")
        r = requests.get(haioRequest)
        with open(tarFile, "wb") as tmp:
            tmp.write(r.content)
    # now read the downloaded tar file.
    isCube = False
    cube = {}
    with tarfile.open(tarFile,'r') as tar:
        for member in tar.getmembers():
            if (('_20ssc_' in member.name) and (vers in member.name)):
                ctype = 'naive'
                if ('convol' in member.name):
                    ctype = 'convol'
                isCube = True
                f=tar.extractfile(member)
                tmp1 = io.BytesIO(f.read())
                tmp2 = gzip.open(tmp1)
                xx = fits.open(tmp2)
                detx = xx[0].header['DETECTOR']
                res = xx[0].header['PROC_RES']
                # comined key for the output dictionary
                xkey = "%s_%s_%s"%(detx,res,ctype)
                cubeData = xx['image'].data
                header = xx['image'].header
                header["CTYPE3"] = "FREQ"
                cube[xkey]= SpectralCube(data=cubeData,wcs=WCS(header))
    tar.close()
    if (not isCube):
        print ("Level-2 file from the Herschel Science Archive has no cube FITS file.")
        print ("Are you sure %i is a FTS spectral mapping observation?"%obsid)
        return None
    return cube
#
#
#
obsid = 1342204895
#result = getSpireFtsCube(obsid,cubeType='naive')
result = getSpireFtsCube(obsid, apod=True)
# visualize the cube as an image at 1200 GHz
#
#%%
ix = result["SSW_HR_convol"].closest_spectral_channel(1200 * u.GHz)
result["SSW_HR_convol"][ix,:,:].quicklook()
#
# now extract 1-d spectrum from the central pixels
#
fig = plt.figure(figsize=(13,5))
cspec = {}
for array in ["SSW","SLW"]:
    nx = result["%s_HR_convol"%array].header["NAXIS1"]
    ny = result["%s_HR_convol"%array].header["NAXIS2"]
    cspec[array] = result["%s_HR_convol"%array][:,int(nx/2),int(ny/2)]
    plt.plot(cspec[array].spectral_axis/1.0e9*u.GHz,cspec[array].array,label="%s pixel(%i,%i)"%(array,nx/2,ny/2))
    pass
#
#%%
# let's plot all spectral form SLW
fig = plt.figure(figsize=(13,5))
xshape = result["SLW_HR_convol"].shape
for i in np.arange(xshape[1]):
    for j in np.arange(xshape[2]):
        plt.plot(result["SLW_HR_convol"].spectral_axis/1.0e9*u.GHz,\
                 result["SLW_HR_convol"][:,i,j].array)
               
plt.xlabel('Frequency (GHz)')
plt.ylabel('Flux (W/m2/Hz/sr)')
#plt.legend(loc='upper left')





