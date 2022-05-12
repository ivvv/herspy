# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 16:39:16 2016

@author: ivaltchanov
"""

#
import os
import io
import math
import zipfile
import gzip
import tarfile
import requests
import numpy as np

import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.modeling import models, fitting

from spectral_cube import SpectralCube

#%%
# global
home = os.path.expanduser('~')
#
# set the tmpDir to a sutable temprary place
#
for idir in [home+ "/tmp",home+"/Tmp",home+"/temp",home+"/Temp","/tmp"]:
    if (os.path.isdir(idir)):
        tmpDir = idir
        break
#
# In[7]:
def getSpireFtsLevel2(obsid, what='spss', cache=True):
    """
    Using the HTTP access to HAIO, retrieve all level-2 products in tar.gz file
    and extract only the requested fits file
    then it puts the spectral structure in a dictionary
    """
    tarFile = tmpDir + "/%i_level2.tar"%obsid
    if (os.path.isfile(tarFile) and cache):
        print ("Found an already existing tar file for OBSID %i. Will Use it"%obsid)
    else:
        haioRequest = "http://archives.esac.esa.int/hsa/aio/jsp/product.jsp?PROTOCOL=HTTP&OBSERVATION_ID=%i&PRODUCT_LEVEL=Level2"%obsid
        print ("Downloading level-2 data from the Herschel Science Archive. May take a while... be patient")
        r = requests.get(haioRequest)
        with open(tarFile, "wb") as tmp:
            tmp.write(r.content)
    # now read the downloaded tar file.
    with tarfile.open(tarFile,'r') as tar:
        for member in tar.getmembers():
            if (what in member.name and '_spg_' in member.name):
                f=tar.extractfile(member)
                xx = fits.open(gzip.open(io.BytesIO(f.read())))
    spec = {}
    with xx as hdu:
        #
        for k in hdu:
            extname = k.name
            if ('S' in extname):
                spec[k.name] = {}
                spec[k.name]['wave'] = k.data["wave"]
                spec[k.name]['flux'] = k.data["flux"]
                spec[k.name]['fluxErr'] = k.data["error"]
    return spec
    
# In[7]:
def getSpireFtsPrecube(obsid):
    """
    Using the HTTP access to HAIO, retrieve the level-2 products in tar.gz file
    and extract only the requested fits files
    """
    tarFile = tmpDir + "/%i_level2.tar"%obsid
    if (os.path.isfile(tarFile) and cache):
        print ("Found an already existing tar file for OBSID %i. Will Use it"%obsid)
    else:
        haioRequest = "http://archives.esac.esa.int/hsa/aio/jsp/product.jsp?PROTOCOL=HTTP&OBSERVATION_ID=%i&PRODUCT_LEVEL=Level2"%obsid
        print ("Downloading level-2 data from the Herschel Science Archive. May take a while... be patient")
        r = requests.get(haioRequest)
        with open(tarFile, "wb") as tmp:
            tmp.write(r.content)
    # now read the downloaded tar file.
    isCube = 0
    precube = {}
    with tarfile.open(tarFile,'r') as tar:
        for member in tar.getmembers():
            if ('spectrum2d' in member.name and '_spg_' in member.name):
                isCube = 1
                f=tar.extractfile(member)
                xx = fits.open(gzip.open(io.BytesIO(f.read())))
                detx = xx[0].header['DETECTOR']
                precube[detx] = xx[1]
    if (not isCube):
        print ("Level-2 file from the Herschel Science Archive has no precube (sprectum2d) FITS file.")
        print ("Are you sure %i is a FTS spectral mapping observation?"%obsid)
        
    return precube

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
    tarFile = tmpDir + "/%i_level2.tar"%obsid
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
                xx = fits.open(gzip.open(io.BytesIO(f.read())))
                detx = xx[0].header['DETECTOR']
                res = xx[0].header['PROC_RES']
                # comined key for the output dictionary
                xkey = "%s_%s_%s"%(detx,res,ctype)
                cubeData = xx['image'].data
                header = xx['image'].header
                header["CTYPE3"] = "FREQ"
                cube[xkey]= SpectralCube(data=cubeData,wcs=WCS(header))
    if (not isCube):
        print ("Level-2 file from the Herschel Science Archive has no cube FITS file.")
        print ("Are you sure %i is a FTS spectral mapping observation?"%obsid)
        return None
    return cube


def extractSpectrum(precube, detector, jiggle=0):
    """
    Extract spectrum from a precube
    """
    array = detector[0:3]
    chk = (precube[array].data["detector"] == detector) & \
        (precube[array].data["jiggId"] == jiggle)
    nspec = len(chk.nonzero())
    if (nspec != 1):
        print ("More than one spectrum for detector %s and jiggle %i"%(detector,jiggle))
    ix = chk.nonzero()[0][0]
    spec = {}
    spec["ra"] = precube[array].data["longitude"][ix]
    spec["dec"] = precube[array].data["latitude"][ix]
    spec[detector] = {}
    spec[detector]["wave"] = precube[array].data["wave"][ix]
    spec[detector]["flux"] = precube[array].data["flux"][ix]
    spec[detector]["fluxErr"] = precube[array].data["error"][ix]
    return spec
    
            
    
    
# In[2]:
def readSpireSparseSpec(spssFile):
    #
    spec = {}
    with fits.open(spssFile) as hdu:
        #
        for k in hdu:
            extname = k.name
            if ('S' in extname):
                spec[k.name] = {}
                spec[k.name]['wave'] = k.data["wave"]
                spec[k.name]['flux'] = k.data["flux"]
                spec[k.name]['fluxErr'] = k.data["error"]
    return spec
#
def plotSpireSparseSpec(specIn, onlyCentral=True):
    """
    """
    central = ['SSWD4','SLWC3']
    plt.figure(figsize=(8,5))
    for det in specIn.keys():
        if (det in central):
            plt.plot(specIn[det]['wave'],specIn[det]['flux'],'k-')
        else:
            plt.plot(specIn[det]['wave'],specIn[det]['flux'],'c-')
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Flux density (Jy)')

# In[2]:

def getFtsBeam(jarFile, band='SSW'):
    # the FTS beam is a 3-D array, one frequency dimension and two spatial dimensions
    zf = zipfile.ZipFile(jarFile, 'r')
    try:
        lst = zf.infolist()
        for zi in lst:
            fn = zi.filename
            if fn.count('SCalSpecBeamProf_%s'%(band)):
                # QUESTION: can this be extracted on the fly? currently the files are extracted in the 
                # current folder
                bb = zf.extract(fn)
                beam = fits.open(bb)
                break
    finally:
        zf.close()
    #
    hh = beam['image'].header
    #
    nfq = hh['NAXIS3']
    fqStart = hh['CRVAL3']
    fqStep = hh['CDELT3']
    # the frequency axis
    fq = fqStart + np.arange(nfq)*fqStep
    #
    return fq,beam
# In[3]:

def fitBeamAndPlot(jarFile):
    """
    """
    fq = {}
    coef = math.sqrt(8.0*math.log(2.0))
    arc2sr = math.pow((3600.0*180.0/math.pi),2) # 1 sq.arcsec to steradian
    #
    beamFwhm = {}
    beamArea = {}
    plt.figure(figsize=(8,5))
    for arr in ["SSW","SLW"]:
        fq[arr], res = getFtsBeam(jarFile, band=arr)
        beamMap = res['image'].data
        beamArea[arr] = np.sum(beamMap,axis=(1,2))/arc2sr
        bshape = beamMap.shape
        nfq = bshape[0]
        #
        # now fit a 2-D Gaussian to the beam
        #
        x,y = np.mgrid[:bshape[1],:bshape[2]]
        # 
        # set up the initial values of the model
        sx = sy = 17.0/coef
        xc = yc = bshape[1]/2.0
        g_init = models.Gaussian2D(amplitude=1, x_mean=xc, y_mean=yc, x_stddev=sx, y_stddev=sy,theta=0)
        #
        fit_p = fitting.LevMarLSQFitter()
        #
        beamFwhm[arr] = np.zeros(nfq)
        for i in np.arange(nfq):
            p = fit_p(g_init, x, y, beamMap[i,:,:])
            # take the geometric mean of the fitted FWHM
            beamFwhm[arr][i] = coef*math.sqrt(p.x_stddev*p.y_stddev)
    # Plot the data with the best-fit FWHM
        plt.subplot(211)
        plt.plot(fq[arr], beamFwhm[arr], 'k-')
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('FWHM (arcsec)')
        plt.subplot(212)
        plt.plot(fq[arr], beamArea[arr], 'k-')    
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Solid angle (sr)')
        pass
    #
    return fq,beamArea,beamFwhm

#
def example_cube():
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
