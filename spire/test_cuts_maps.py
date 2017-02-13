# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 09:32:36 2016

@author: ivaltchanov
"""
import requests
import os, io
import gzip
import tarfile
import warnings
import numpy as np

from astropy.io.votable import parse_single_table
from astropy.io import fits
from astropy.visualization import (LogStretch, ImageNormalize, \
    PercentileInterval, AsymmetricPercentileInterval, SqrtStretch,\
    AsinhStretch)

#from sklearn.preprocessing import scale
from scipy import stats
#from sklearn.neighbors.kde import KernelDensity

import matplotlib.pyplot as plt
#import matplotlib.colors as colors
from matplotlib import cm

from PIL import Image, PngImagePlugin

# global
home = os.path.expanduser('~')
#
# set the tmpDir to a sutable temprary place
#
for idir in [home+ "/tmp",home+"/Tmp",home+"/temp",home+"/Temp","/tmp"]:
    if (os.path.isdir(idir)):
        tmpDir = idir
        break

def getPhotMaps(obsid, level='Level2', what="extd", cache=True):
    """
    Download VOTABLE with the metadata search
    Then save the products in a tar file and return the maps in a dictionary
    level - what level to extract, can be Level2, Level2_5 or Level3
    what - 'extd' for extended-source calibrated maps, 'psrc' for point-source
    cache - will check if votable or a tar files are already available in tmpDir
    
    """
    archiveMetaUrl = "http://archives.esac.esa.int/hsa/aio/jsp/metadata.jsp"
    archiveProdUrl = "http://archives.esac.esa.int/hsa/aio/jsp/product.jsp"
    #
    haioMetaRequest = "%s?RESOURCE_CLASS=OBSERVATION&INSTRUMENT=%%%%27SPIRE%%%%27&OBSERVATION_ID=%i"%(archiveMetaUrl,obsid)
    #print (haioMetaRequest)
    #
    voFile =  tmpDir + "/%i_metadata.vot"%obsid
    if (os.path.isfile(voFile)):
        print ("Found an already existing metadata VOTABLE file for OBSID %i. Will reuse it"%obsid)
    else:
        r = requests.get(haioMetaRequest)
        with open(voFile, "wb") as tmp:
            tmp.write(r.content)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        votable = parse_single_table(voFile, pedantic=False)
    #
    data = votable.array
    #
    procLevel = data['Obs_State'].data[0].decode()
    instMode = data["Instrument_mode"].data[0].decode()
    #
    if ('LEVEL3' in procLevel):
        print ("Observation %i has level-3 map"%obsid)
        levelx = 'Level3'
    elif ('LEVEL2_5' in procLevel):
        print ("Observation %i has level-2_5 map"%obsid)
        levelx = 'Level2_5'
    else:
        levelx = 'Level2'
    #
    print ("OBSID %i (%s) has %s"%(obsid,instMode,levelx) )
    print (" ===> you asked for %s"%level)
    #
    #
    tarFile = tmpDir + "/%i_%s_%s.tar"%(obsid,instMode,level)
    #
    #levelx = "Level2"
    if os.path.isfile(tarFile):
        print ("Found an already existing tar file for OBSID %i. Will Use it"%obsid)
    else:
        haioProdRequest = "%s?PROTOCOL=HTTP&OBSERVATION_ID=%i&PRODUCT_LEVEL=%s"%(archiveProdUrl,obsid, level)
        r = requests.get(haioProdRequest)
        if (b'No results' in r.content):
            print ("*** No results")
            return None
        else:
            with open(tarFile, "wb") as tmp:
                tmp.write(r.content)
    #
    maps = {}
    with tarfile.open(tarFile,'r') as tar:
        for member in tar.getmembers():
            if ((what in member.name) and (not 'diag' in member.name)):
                f=tar.extractfile(member)
                xx = fits.open(gzip.open(io.BytesIO(f.read())))
                band = xx[0].header["DETECTOR"]
                maps[band] = xx
    #
    return maps
#
#

#%%
#NGC2841, 1342193797, PACS 1342207168
#M101, 1342188750, PACS 1342198471
#NGC1291, 1342189508, PACS 1342189569
#NGC2915, 1342189510 => circus and sources, PACS 1342204428

#myObsid = 1342193797 # NGC2841
#myObsid = 1342188750 # M101
#myObsid = 1342189508 # NGC1291
#myObsid = 1342189510 # NGC2915
#myObsid = 1342182459 # M83
#myObsid = 1342182474 # B0420-014
#myObsid = 1342183050 # NGC6946
#myObsid = 1342183071 # GP-field1
#myObsid = 1342183496 # dark sky
myObsid = 1342185644 # Lockman Hole
#myObsid = 1342186840 # DR21

mapx = getPhotMaps(myObsid, level="Level2_5")
if (mapx == None):
    print ('You asked for level 2.5 but there isn\'t one. Will get level 2 instead')
    mapx = getPhotMaps(myObsid, level="Level2")
#
#
#%%
band = 'PSW'
wcsheader = mapx[band]['image'].header
if ('OBJECT' in mapx[band][0].header.keys()):
    target = mapx[band][0].header['OBJECT']
else:
    target = 'Unknown'
image = mapx[band]['image'].data
imx = image[~np.isnan(image)]
#mean_im = np.mean(imx)
#stdev_im = np.std(imx)
#print (mean_im,stdev_im)
min_image = np.min(imx)
max_image = np.max(imx)
#
#image_scaled = (image - mean_im)/stdev_im
#image_scaled = (image - min_image)/(max_image - min_image)
#
# now find the mode
#
mode = stats.mode(image,axis=None,nan_policy='omit')
# normalized mode
mode_x = (mode.mode[0] - min_image)/(max_image - min_image)
print ("Mode is at %f percentile point, will use it as lower cut)"%(mode_x))
#
#image_scaled = scale(image[~np.isnan(image)],with_mean=True, with_std=True, copy=True)
#n1, bins1, patches1 = plt.hist(imx.ravel(),256,[0,256],normed=1)
#n1, bins1, patches1 = plt.hist(image_scaled[~np.isnan(image_scaled)], 50, normed=1, facecolor='green', histtype='step',alpha=0.75)
#n1, bins1, patches1 = plt.hist(image_scaled[~np.isnan(image_scaled)], bins='auto', normed=1, facecolor='green', histtype='step',alpha=0.75)
#
#
#grad = np.gradient(n1)
#plt.plot(bins1[0:-1],grad,'r-')
# gradient peak
#g_peak = bins1[np.argmax(grad)]
#g_dip = bins1[np.argmin(grad)]
# histogram peak
#h_peak = bins1[np.argmax(n1)]
#h_peak = max(h_peak,0.1) # use 10% as the lower minimum, when the histogram peak is below 0.1

#low_perc = bins1[max(0,np.argmax(n1)-1)]
#
#imx = image_scaled[~np.isnan(image_scaled)]
#kde = stats.gaussian_kde(imx)
#density = kde(imx)
#peak = max(peak,0.01)
#low_perc = max(low_perc,0.01)
# set 1% as the minimum cut
#print ("Histo Peak at %f, will use it as lower percentile (i.e. remove everything below)"%(h_peak))
#
#peak_norm = (peak - np.min(bins1))/(np.max(bins1) - np.min(bins1))
#print ("Normalized peak is at %f"%(peak_norm))
#n1, bins1, patches1 = plt.hist(image[~np.isnan(image)], 50, normed=1, facecolor='green', histtype='step',alpha=0.75)
#
#%%
fig  = plt.figure(figsize=(10,10))

norm = ImageNormalize(image,interval=AsymmetricPercentileInterval(mode_x*100,99.5),stretch=AsinhStretch())
#norm = ImageNormalize(image,interval=PercentileInterval(98),stretch=SqrtStretch())
#norm = ImageNormalize(image,interval=PercentileInterval(99.5),stretch=LogStretch())
ax = fig.add_subplot(111)
ax.imshow(image, cmap=cm.gray, norm=norm, origin="lower", interpolation="nearest")
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])    
ax.annotate("%s, %i %s"%(target,myObsid,band), textcoords='axes fraction',fontsize=10, xy=(0.05,0.95),xytext=(0.05,0.95))
#plt.show()
pngFile = tmpDir + "/%i_%s_plot.png"%(myObsid,band)
plt.savefig(pngFile, dpi=300)
plt.close()

im = Image.open(pngFile)
meta = PngImagePlugin.PngInfo()
useMeta = ["NAXIS","NAXIS1","NAXIS2","CRPIX1","CRPIX2","CRVAL1","CRVAL2","CDELT1","CDELT2",
           "CTYPE1","CTYPE2","EQUINOX","CROTA2"]
for key in wcsheader.keys():
    if (key not in useMeta):
        continue
    print (key,wcsheader[key])
    meta.add_text(key, "%s"%wcsheader[key])
im.convert(mode='L').save(pngFile, "png", pnginfo=meta)

#f, ax = plt.subplots(1, 3)

#for i,band in enumerate(["PSW","PMW","PLW"]):
#    image = mapx[band]['image'].data
#    norm = ImageNormalize(image,interval=AsymmetricPercentileInterval(30,99.9),stretch=AsinhStretch())
#    #norm = ImageNormalize(image,interval=PercentileInterval(98),stretch=SqrtStretch())
#    #norm = ImageNormalize(image,interval=PercentileInterval(99.5),stretch=LogStretch())
#    ax[i].imshow(image, cmap=cm.gray_r, norm=norm, origin="lower")
#    ax[i].get_xaxis().set_ticks([])
#    ax[i].get_yaxis().set_ticks([])    
#    ax[i].annotate('%s'%band, textcoords='axes fraction',fontsize=8, xy=(0.1,0.9),xytext=(0.1,0.9))
#plt.show()
#plt.savefig(tmpDir + "/%i_plot.jpeg"%myObsid, dpi=300)
#plt.close()