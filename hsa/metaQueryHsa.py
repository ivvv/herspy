# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 09:32:36 2016

@author: ivaltchanov
"""
import requests
import os, io
from astropy.io.votable import parse_single_table
from astropy.io import fits
import gzip
import tarfile
import warnings

# global
home = os.path.expanduser('~')
#
# set the tmpDir to a sutable temprary place
#
for idir in [home+ "/tmp",home+"/Tmp",home+"/temp",home+"/Temp","/tmp"]:
    if (os.path.isdir(idir)):
        tmpDir = idir
        break

#%%
#obsid = 1342239959 # level-3 is big!!!
obsid = 1342193797 # level-3 is big!!!

archive = 'hsa'
archiveMetaUrl = "http://archives.esac.esa.int/%s/aio/jsp/metadata.jsp"%archive
archiveProdUrl = "http://archives.esac.esa.int/%s/aio/jsp/product.jsp"%archive

#haioMetaRequest = "%s?RESOURCE_CLASS=OBSERVATION&QUERY=%%28INSTRUMENT==%%27SPIRE%%27%%20AND%%20"%(archiveMetaUrl)
#haioMetaRequest += "OBSERVATION_ID==%i%%29"%(obsid)
haioMetaRequest = "%s?RESOURCE_CLASS=OBSERVATION&INSTRUMENT=%%%%27SPIRE%%%%27&OBSERVATION_ID=%i"%(archiveMetaUrl,obsid)
print (haioMetaRequest)
#
r = requests.get(haioMetaRequest)
voFile =  "/Users/ivaltchanov/Tmp/test.vot"
with open(voFile, "wb") as tmp:
    tmp.write(r.content)
#
#%%
#
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
#%%
#
tarFile = tmpDir + "/%i_%s_%s.tar"%(obsid,instMode,levelx)
#
levelx = "Level2"
if os.path.isfile(tarFile):
    print ("Found an already existing tar file for OBSID %i. Will Use it"%obsid)
else:
    haioProdRequest = "%s?PROTOCOL=HTTP&OBSERVATION_ID=%i&PRODUCT_LEVEL=%s"%(archiveProdUrl,obsid, levelx)
    r = requests.get(haioProdRequest)
    if (b'No results' in r.content):
        print ("*** No results")
    else:
        with open(tarFile, "wb") as tmp:
            tmp.write(r.content)
#
#%%
what = "extd"
maps = {}
with tarfile.open(tarFile,'r') as tar:
    for member in tar.getmembers():
        if ((what in member.name) and (not 'diag' in member.name)):
            f=tar.extractfile(member)
            xx = fits.open(gzip.open(io.BytesIO(f.read())))
            band = xx[0].header["DETECTOR"]
            maps[band] = xx
#
#

