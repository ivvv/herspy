#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 14:15:45 2017

Metadata query for HIFI observations

@author: ivaltchanov
"""
import os
import requests
import warnings
from astropy.io.votable import parse_single_table

def getHifiMeta(obsid, cache = False, useTemp=True):
    #
    # Only metadata search is performed and raNominal,decNominal are xtracted from the VOT
    #
    archive = 'hsa'
    archiveMetaUrl = "http://archives.esac.esa.int/%s/aio/jsp/metadata.jsp"%archive
    haioMetaRequest = "%s?RESOURCE_CLASS=OBSERVATION&INSTRUMENT=%%%%27HIFI%%%%27&OBSERVATION_ID=%i"%(archiveMetaUrl,obsid)
    #
    # VOTable to keep the metadata queryresults
    #
    if (useTemp):
        voFile = "/Users/ivaltchanov/Tmp/hifi_meta_temp.vot"
    else:
        voFile =  "/Users/ivaltchanov/Tmp/%i_hifi_meta.vot"%obsid
    if (cache and (not useTemp)):
        if (os.path.isfile(voFile)):
            print("Found an already existing VOTable with metadata: %s"%voFile)
    else:
        r = requests.get(haioMetaRequest)
        with open(voFile, "wb") as tmp:
            tmp.write(r.content)
    #
    # now red the table
    #
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        votable = parse_single_table(voFile, pedantic=False)
    #
    # extract the data
    data = votable.array
    #
    return data
