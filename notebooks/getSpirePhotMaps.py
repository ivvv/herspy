import io
import os
import requests
import tarfile
import gzip

from astropy.io import fits
#
# a function to access the Herschel Science Archive and download maps
#

def getSpirePhotMaps(obsid, level='Level2', what='extd', cache=True, tmpDir=None):
    """
    Direct download of the HSA tar file using the tap interface

    level - what level to extract, can be Level2, Level2_5 or Level3
    what - 'extd' for extended-source calibrated maps, 'psrc' for point-source
    cache - will check if votable or a tar files are already available in tmpDir

    """
    if (level == 'Level3'):
        lx = '3'
    elif (level == 'Level2_5'):
        lx = '2_5'
    else:
        lx = '2'
    #
    URLbase = 'http://archives.esac.esa.int/hsa/whsa-tap-server/data?'
    query = f'RETRIEVAL_TYPE=OBSERVATION&OBSERVATION_ID={obsid}&INSTRUMENT_NAME=SPIRE&PRODUCT_LEVEL={lx}'
    if (not os.path.isdir(tmpDir)):
        raise FileNotFound ("*** tmpDir is not set, please set it first. It will be used for cache and to save the tar files")
        return None
    tarFile = f'{tmpDir}/{obsid}_{level}.tar'
    #
    if os.path.isfile(tarFile) and cache:
        print (f'Found an already existing tar file for OBSID {obsid}. Will use it')
    else:
        print ("Downloading products from the Herschel Science Archive. Please wait")
        tapRequest = URLbase + query
        r = requests.get(tapRequest)
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
