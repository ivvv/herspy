import numpy as np
from astropy.wcs import WCS

def check_coverage(ra,dec,hdu):
    #
    # check and return the average coverage in a 3x3 box around (ra,dec) coordinates
    #
    # Note: it does not check if the point is outside the HDU['coverage'] image
    # in this case it will return 0 anyway.
    #
    try:
        cvrg = hdu['coverage']
    except:
        print ("Cannot read the coverage map")
        return 0.0
    #
    wcs = WCS(cvrg.header)
    (ny,nx) = cvrg.data.shape
    (xp,yp) = wcs.wcs_world2pix(ra,dec, 1)
    xp = xp.astype(int)
    yp = yp.astype(int)
    if (yp < 0 or xp < 0 or yp > ny or xp > nx):
        print ("Point outside the input image")
        return 0.0
    #
    # get a 3x3 pixels around each source and take the average coverage
    #
    y0 = max(0,yp - 1)
    y1 = min(ny,yp + 2)
    x0 = max(0,xp - 1)
    x1 = min(nx,xp + 2)
    wcov = cvrg.data[y0:y1,x0:x1]
    return np.mean(wcov)
#
