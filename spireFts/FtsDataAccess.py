# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 16:39:16 2016

@author: ivaltchanov
"""

#
import matplotlib.pyplot as plt
import numpy as np
import math
from astropy.io import fits
import zipfile
import tarfile
from astropy.modeling import models, fitting
import requests

# In[7]:
def getSpireFtsLevel2(obsid, what='spss'):
    """
    Using the HTTP access to HAIO, retrieve all level-2 products in tar.gz file
    and extract only the requested fits file
    then it puts the spectral structure in a dictionary
    """
    tarFile = "%i_level2.tar"%obsid
    haioRequest = "http://archives.esac.esa.int/hsa/aio/jsp/product.jsp?PROTOCOL=HTTP&OBSERVATION_ID=%i&PRODUCT_LEVEL=Level2"%obsid
    print ("Downloading level-2 data from the Herschel Science Archive. May take a while... be patient")
    r = requests.get(haioRequest)
    with open(tarFile, "wb") as tmp:
        tmp.write(r.content)
    # now read the downloaded tar file.
    with tarfile.open(tarFile,'r') as tar:
        for member in tar.getmembers():
            if (what in member.name and '_spg_' in member.name):
                f=tar.extract(member)
                xx = fits.open(member.name)
    tar.close()
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
    tarFile = "%i_level2.tar"%obsid
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
                f=tar.extract(member)
                xx = fits.open(member.name)
                detx = xx[0].header['DETECTOR']
                precube[detx] = xx[1]
    tar.close()
    if (not isCube):
        print ("Level-2 file from the Herschel Science Archive has no precube (sprectum2d) FITS file.")
        print ("Are you sure %i is a FTS spectral mapping observation?"%obsid)
        
    return precube

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