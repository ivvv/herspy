# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 09:32:36 2016

@author: ivaltchanov
"""
import requests
from astropy.io.votable import parse_single_table
import warnings

#obsid = 1342239959 # level-3 is big!!!
obsid = 1342239959 # level-3 is big!!!

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
#-----------
#
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    votable = parse_single_table(voFile, pedantic=False)
data = votable.array
procLevel = data['Obs_State'].data[0]
if (b'LEVEL3' in procLevel):
    print ("Observation %i has level-3 map"%obsid)
    levelx = 'Level3'
elif (b'LEVEL2_5' in procLevel):
    print ("Observation %i has level-2_5 map"%obsid)
    levelx = 'Level2_5'
else:
    levelx = 'Level2'
#
#----------------
#
haioProdRequest = "%s?PROTOCOL=HTTP&OBSERVATION_ID=%i&PRODUCT_LEVEL=%s&COMPRESSION=TARGZ"%(archiveProdUrl,obsid, levelx)
r = requests.get(haioProdRequest)
tgzFile =  "/Users/ivaltchanov/Tmp/test.tgz"
with open(tgzFile, "wb") as tmp:
    tmp.write(r.content)
#
