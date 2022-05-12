import numpy as np
from spherical_geometry import polygon, vector
#
# function to search for a point in footprints polygons
#
def whereis(raIn,decIn,wcs_corners):
    #
    # will identify all obsids where the input raIn,DecIn point is within the footprint polygon
    #
    # raI, decIn are the coordinates in degrees
    # wcs_corners is a astropy.table.Table objects with the obsid and the 4 footprint corners.
    #
    # returns the list of obsids
    #
    nm = len(wcs_corners)
    #
    obsids = wcs_corners["obsid"].data
    ra1 = wcs_corners["ra1"].data
    ra2 = wcs_corners["ra2"].data
    ra3 = wcs_corners["ra3"].data
    ra4 = wcs_corners["ra4"].data
    dec1 = wcs_corners["dec1"].data
    dec2 = wcs_corners["dec2"].data
    dec3 = wcs_corners["dec3"].data
    dec4 = wcs_corners["dec4"].data
    ra0 = ra1
    dec0 = dec1
    #
    # storing in a dictionary with SphericalPolygons
    spx = {}
    for i in np.arange(nm):
        iobs = obsids[i]
        ra_corn = [ra1[i],ra2[i],ra3[i],ra4[i],ra0[i]]
        dec_corn = [dec1[i],dec2[i],dec3[i],dec4[i],dec0[i]]
        spoly = polygon.SphericalPolygon.from_radec(ra_corn,dec_corn, degrees=True)
        spx[iobs] = spoly
    #
    # now checking if the input falls in a polygon
    #
    vec = vector.radec_to_vector(raIn,decIn)
    obsFound = []
    for j in spx.keys():
        if (spx[j].contains_point(vec)):
            obsFound.append(j)
        nfound = len(obsFound)
    return obsFound
#
