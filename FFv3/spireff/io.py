#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import io
import requests
import urllib3
import tarfile
import gzip

import numpy as np
from pandas import read_csv
from matplotlib.pyplot import imread

from astropy.io import fits
from astropy.wcs import WCS
from astropy.units import deg
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import pixel_to_pixel

verbose = False
def set_verbose(val):
    global _verboseprint, verbose
    verbose = val
    _verboseprint = print if verbose else lambda *a, **k: None
set_verbose(verbose)

http = urllib3.PoolManager()
spg_url = "http://archives.esac.esa.int/hsa/aio/jsp/"
hpdp_url = ("http://archives.esac.esa.int/hsa/legacy/HPDP/SPIRE/SPIRE-S/"
            "cal_targets/SPIREspec_calibrators/")
bgs_url = ("http://archives.esac.esa.int/hsa/legacy/HPDP/SPIRE/SPIRE-S/"
           "BKGS/spectra/")
archive = ('http://archives.esac.esa.int/hsa/legacy/HPDP/SPIRE/SPIRE-S/'
           'spectral_feature_catalogue/')

band_SLW, band_SSW = [447, 1018], [944, 1568] # [GHz]
beam = {"SLW": 0.00972222, "SSW": 0.00527778} # Nominal beam FWHM [deg]
resolution = {'LR': 25.0, 'MR': 7.2, 'HR': 1.18448225} # [GHz]

def frequency_beam(save_dir=None, save=False):
    output = {'SLW': None, 'SSW': None}
    files = {'SLW': 'fwhm_omega_SLW_HR.csv', 'SSW': 'fwhm_omega_SSW_HR.csv'}
    root = 'http://archives.esac.esa.int/hsa/legacy/ADP/PSF/SPIRE/SPIRE-S/'
    
    if save_dir is None:
        save_dir = os.getcwd()
    for ary, ext in files.items():
        save_file = os.path.join(save_dir, ext)
        if os.path.isfile(save_file):
            _verboseprint(f"file : {ext} already exists. Will use it.")
            output[ary] = read_csv(save_file).astype(float)
        else:
            f_beam = read_csv(root + ext, header=3)[3:].astype(float)
            output[ary] = f_beam
            if save: 
                f_beam.to_csv(save_file)
    return output


def _map_meta(vals, header):
    """
    Maps standard HSA keywords to the keywords in a fits file header.

    Parameters
    ----------
    vals : string or list of strings
        standard HSA keywords to be mapped.
    header : fits.header.Header
        Fits file header that `vals` are mapped to.

    Returns
    -------
    dict
        A dictionary with the form {vals: header_keyword}.

    Note : Many standared HSA keywords are too long for fits file header
        keywords and are thus replaced. In most cases the replacement keywords
        are standardized eg, "mapSampling" -> "MAPSAMPL", but this is not
        always the case. HSA keyword mapping is presented in the header
        hierarchy data. This information is used to reverse the mapping in
        this function.
    """
    if isinstance(vals, str):
        vals = [vals]
    return {v: k.replace('key.', '') for k, v, d in header.cards if v in vals}

def _keys_summary(obj, indent='--', _n=0):
    """
    A convenience function that prints the key structucture of a dictionary.

    Parameters
    ----------
    obj : dictionary or other
        The dictionary or subdictionary who's key sturcute is being
        investigated. If the object has no `.keys()` method, the recursive
        function exits.
    indent : string, optional
        Indentation annotation for the key hierarchy. The default is '--'.
    _n : int, optional
        Counter to keep track of recursion. Should not be manipulated by user.
        The default is 0.

    Returns
    -------
    None.
    """
    print(f"\n{' Summary ':_^15}") if _n == 0 else None
    for key in obj.keys():
        print(indent*_n + str(key) + (':' if _n == 0 else ''))
        try:
            obj_new = obj[key]
            _keys_summary(obj_new, _n=_n+1)
        except AttributeError:
            continue
    if _n == 0:
        print(f"{' End ':_^15}\n")

def _init_var(*args, **kwargs):
    """
    Converts all args and kwargs to iterables. Returns variables in the same
    order they were input.
    """
    args = list(args)
    for i, arg in enumerate(args):
        if arg is None:
            args[i] = []
        elif isinstance(arg, str):
            args[i] = [arg]
        else:
            try:
                _ = iter(arg)
            except TypeError:
                args[i] = [arg]

    for key, val in kwargs.items():
        if val is None:
            kwargs[key] = []
        elif isinstance(val, str):
            kwargs[key] = [val]
        else:
            try:
                _ = iter(val)
            except TypeError:
                kwargs[key] = [val]

    if len(kwargs) == 0 and len(args) > 0:
        return [*args]
    elif len(kwargs) > 0 and len(args) == 0:
        return kwargs
    else:
        return (*args, kwargs)

def _make_wcs(cont_table):
    """
    Generates a WCS object form a mapping continuum parameters table.

    Parameters
    ----------
    cont_table : HDUList
        Mapping continuum parameters HDUList FF product.

    Returns
    -------
    output : dict of WCS
        A WCS object derived from the row/column and dec/ra values of
        `cont_table`. Nominally, a WCS object is created for each SPIRE band.
        Format is {'SLW': wcs_slw, 'SSW': wcs_ssw}
    """

    output = dict()
    data = cont_table[1].data
    for ary in ['SLW', 'SSW']:
        ind_ary = np.where(data['array'] == ary)[0]
        data_ = data[ind_ary]
        ra, dec = data_['ra'][0], data_['dec'][0]
        row, col = data_['row'][0], data_['column'][0]
        #ra_lim = [np.max(data_['ra']), np.min(data_['ra'])]
        dec_lim = [np.min(data_['dec']), np.max(data_['dec'])]
        row_lim = [np.min(data_['row']), np.max(data_['row'])]
        col_lim = [np.min(data_['column']), np.max(data_['column'])]
        row_diff, col_diff = np.diff(row_lim), np.diff(col_lim)
        cdelt = np.diff(dec_lim)[0]/row_diff[0]

        w = WCS(naxis=2)
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.crpix = [col+1, row+1]
        w.wcs.crval = [ra, dec]
        w.wcs.cdelt = [-cdelt, cdelt]
        header = w.to_header()
        header['NAXIS1'] = int(col_diff+1)
        header['NAXIS2'] = int(row_diff+1)

        output[ary] = WCS(header)

    return output

def _pixel_to_pixel(pixels, wcs):
    """
    Convenience function that wraps astropy.wcs.utils.pixel_to_pixel and maps
    SLW cube pixels to nearest SSW cube pixels and vice versa.

    Parameters
    ----------
    pixels : dictionary
        Dictionary containing the SLW and SSW pixels to be mapped. Format
        should be {'SLW': [[cols], [rows]], 'SSW': [[cols], [rows]]} for
        multiple pixels, or {'SLW': [col, row], 'SSW': [col, row]} for
        a single pixel. Both SPIRE bands are not required.
    wcs : dictionary of WCS
        Dictionary of WCS for SLW and SSW bands. see `_make_wcs`.

    Returns
    -------
    output : dictionary
        Dictionary of mapped pixels. Format matches `pixels` with reciprocal
        SPIRE band keys.
    """
    swap = {'SLW': 'SSW', 'SSW': 'SLW'}
    output = dict()
    for ary, pix in pixels.items():
        pix = np.array(pix)
        output[swap[ary]] = np.round(pixel_to_pixel(wcs[ary], wcs[swap[ary]],
                                                *pix), 0).astype(int).tolist()
    return output

def _pixel_in_region(ra, dec, wcs, radius=None, flatten=True):
    """
    Identifies the SLW and SSW cube pixels, using the provided wcs variable,
    who's centers reside withing a radius about the ra, dec coordinates.

    Parameters
    ----------
    ra : float or list of floats
        Right ascension coordinates. Assumes units are degrees.
    dec : float or list of floats
        Declination coordinates. Assumes units are degrees.
    wcs : dictionary of WCS
        Dictionary of WCS for SLW and SSW bands. see `_make_wcs`.
    radius : astropy.units Quantity or None, optional
        Search radius about coordinates. If None, the radius is set to half
        the SPIRE SLW beam. The default is None.
    flatten : bool, optional
        If True, returns concatenated list of pixels for all search regions.
        If False, keeps the pixel lists seperate for each search region.
        The default is True.

    Returns
    -------
    output : dictionary
        Dictionary of pixels matching the search criteria. Format is
        {'SLW': [[cols], [rows]], 'SSW': [[[cols], [rows]]]} if `flatten`
        or
        {'SLW': [[[cols_1], [rows_1]], ..., [[cols_n], [rows_n]]],
         'SSW': [[[cols_1], [rows_1]], ..., [[cols_n], [rows_n]]]} if
        not `flatten`.

    Note : Length of ra and dec should match.

    """
    ra, dec = _init_var(ra, dec)
    if radius is None: radius = beam['SLW']/2.0 * deg

    output = dict()
    coord = SkyCoord(ra, dec, frame='icrs', unit='deg')
    for ary, w in wcs.items():
        output[ary] = []
        x, y = np.meshgrid(range(w.array_shape[0]), range(w.array_shape[1]))
        x, y = x.flatten(), y.flatten()
        search_coords = w.pixel_to_world(x, y)
        i_search, i_coord, *_ = coord.search_around_sky(search_coords, radius)
        if flatten:
            output[ary] += [list(x[i_search]), list(y[i_search])]
        else:
            for i in range(len(coord)):
                ind = i_search[i_coord == i]
                output[ary].append([list(x[ind]), list(y[ind])])

    return output

def _pixels_in_mapping_products(products):
    """
    Generates a dictionary of cube pixels in mapping products obtained from 
    the `parse_mapping` function.

    Parameters
    ----------
    products : dictionary
        Dictioanry returned from `parse_mapping` function.

    Returns
    -------
    output : dictionary
        Dictionary with the appropirate format for use as `pixels` argument in
        other functions in this library.

    """
    
    output = dict()
    for ary in ['SLW', 'SSW']:
        if ary in products:
            prod, pix = products[ary], [[],[]]
            for col in prod.keys():
                for row in prod[col].keys():
                    pix[0].append(col)
                    pix[1].append(row)
            output[ary] = pix
        else:
            output.update({ary: np.transpose([eval(key.replace(ary, '')) 
                        for key in products.keys() if ary in key]).tolist()})
        
    return output

def gen_frequency(cube):
    """
    Generates frequency axis for mapping spectral cube.

    Parameters
    ----------
    cube : fits HDUList
        HDUList for a spectral cube.

    Raises
    ------
    ValueError
        Cube must have an image hdu. This hdu contains the requed WCS
        information.

    Returns
    -------
    numpy array
        Spectral cube frequency axis.

    """
    if not isinstance(cube, fits.hdu.image.ImageHDU):
        try:
            cube = cube['image']
        except:
            raise ValueError('func : ``gen_frequency`` var : ``cube`` must '
                             'have ImageHDU.')

    header = cube.header
    return WCS(header).all_pix2world(0, 0, range(header['NAXIS3']), 0)[-1]

def parse_mapping(cat_table, cont_table, cube_hdu=None, pixels=None,
                  nearest=False, coords=None, radius=None, flatten=True):
    """
    Returns FF mapping products, and corresponding HSA spectra, meeting search
    criteria determined by input parameters.

    Parameters
    ----------
    cat_table : fits HDUList
        FF features catalogue.
    cont_table : fits HDUList
        FF continuum parameters.
    cube_hdu : dictionary of fits HDULists, optional
        Dictionary of HSA spectral cubes.
        Format is {'SLW': cube_sle, 'SSW': cube_ssw}.
        Does not need both SPIRE bands. If cubes are provided, output pixels
        will have spectra ("flux") and frequency axis ("wave") for all pixels
        in the associated cube. The default is None.
    pixels : dictionary of lists, optional
        Dictioinary of spectral cube pixels.
        Format must match `pixels` input for `_pixel_to_pixel`.
        The default is None.
    nearest : bool, optional
        If True, output will include nearest pixels from conjugate SPIRE band.
        In other words, the output of `_pixel_to_pixel` with input `pixels`
        will be included in the output. The default is False.
    coords : list or list of lists, optional
        Pixels from both SPIRE bands within a `radius` about `coord` will be
        included in the output. Format should be [[ras], [decs]] for multiple
        coordinates, or [ra, dec] for a single coordinate. Values should be
        floats with units of degrees. The default is None.
    radius : astropy.units Quantity or None, optional
        Search radius about coordinates. If None, the radius is set to half
        the SPIRE SLW beam. The default is None.
    flatten : boo, optional
        If True, output dictionary will have keys in the form "ary[col,row]".
        If False, output dictionary will have a cascading key structure in the
        form ["ary"][col][row]. "ary" implies a SPIRE band.
        The default is True.

    Returns
    -------
    output : dictionary
        Ouput dictionary for all cube pixels meeting the search criteria
        determined by input parameters.

    Note : Due to the large number of pixels in SPIRE cubes, output structure
        can be difficult to interpret. Calling `_key_summary(output)` can
        make the output more digestible.

    Caution : Output structure is ary -> column -> row, though the standard
        Python structure would be ary -> row -> column. We chose our structure
        to be more consistent with WCS objects, so that the indexing implies
        ary -> ra -> dec.
    """
    if radius is None: radius = beam['SLW']/2.0 * deg
    output = {'header': (cat_table['data'].header,)}
    data, cont_data = cat_table['data'].data, cont_table[1].data
    cat_slw = data[np.where(data['array'] == 'SLW')[0]]
    cat_ssw = data[np.where(data['array'] == 'SSW')[0]]
    con_slw = cont_data[np.where(cont_data['array'] == 'SLW')[0]]
    con_ssw = cont_data[np.where(cont_data['array'] == 'SSW')[0]]
    wcs = _make_wcs(cont_table)

    iter_pixels = {'SLW': [], 'SSW': []}
    if pixels is not None:
        pixels = {ary: [[item] if type(item) == int else list(item)
                        for item in vals] for ary, vals in pixels.items()}
        if nearest:
            nearest_pixels = _pixel_to_pixel(pixels, wcs)
            nearest_pixels.update({ary: [[], []] for ary in ['SLW', 'SSW'] \
                           if ary not in nearest_pixels})
            pixels.update({ary: [[], []] for ary in ['SLW', 'SSW'] \
                           if ary not in pixels})
            for ary in pixels:
                pixels[ary][0] += nearest_pixels[ary][0]
                pixels[ary][1] += nearest_pixels[ary][1]
        iter_pixels = {ary: iter_pixels[ary] + list(zip(*pixels[ary]))
                       for ary in pixels}

    if coords is not None:
        region_pixels = _pixel_in_region(*coords, wcs=wcs, radius=radius, 
                                         flatten=True)
        iter_pixels = {ary: iter_pixels[ary] + list(zip(*region_pixels[ary]))
                       for ary in region_pixels}

    if (pixels is None) and (coords is None):
        for ary, w in wcs.items():
            x, y = np.meshgrid(range(w.array_shape[0]), range(w.array_shape[1]))
            iter_pixels[ary] = list(zip(x.flatten(), y.flatten()))

    temp = {'SLW': [cat_slw, con_slw], 'SSW': [cat_ssw, con_ssw]}
    for ary, pix_values in iter_pixels.items():
        cat, con = temp[ary]
        if not flatten: output[ary] = dict()
        if (cube_hdu is not None) and (ary in cube_hdu):
            freq_axis = gen_frequency(cube_hdu[ary])
            flux_data = cube_hdu[ary]['image'].data
        for (col, row) in sorted(set(pix_values)):
            con_ind = np.where((row==con['row']) &
                               (col==con['column']))[0]
            cat_ind = np.where((row==cat['row']) &
                               (col==cat['column']))[0]
            if len(con_ind) == 0:
                print(f'WARNING : {ary}[{col},{row}] has no data...')
                continue
            if (not flatten) and (col not in output[ary]):
                output[ary][col] = dict()
            features = cat['frequency'][cat_ind]
            snrs = cat['SNR'][cat_ind]
            # vel = np.unique(cat['velocity'][cat_ind])[0]
            # vel_err = np.unique(cat['velocityError'][cat_ind])[0]
            cont_par = np.array(con[con_ind].item()[6::2], dtype=float)
            if flatten:
                output[f'{ary}[{col},{row}]'] = {'lines': {'freq': features,
                                            'snr': snrs}, 'con_par': cont_par}
            else:
                output[ary][col][row] = {'lines': {'freq': features,
                                            'snr': snrs}, 'con_par': cont_par}

            if (cube_hdu is not None) and (ary in cube_hdu):
                if flatten:
                    output[f'{ary}[{col},{row}]'].update({'wave': freq_axis,
                                         'flux': flux_data[:, row, col]})
                else:
                    output[ary][col][row].update({'wave': freq_axis,
                                         'flux': flux_data[:, row, col]})
    return output

def parse_sparse(cat_table, cont_table, spec_hdu=None, detectors=None):
    """
    Returns FF sparse products, and corresponding HSA spectra, meeting search
    criteria determined by input parameters.

    Parameters
    ----------
    cat_table : fits HDUList
        FF features catalogue.
    cont_table : fits HDUList
        FF continuum parameters.
    spec_hdu : fits HDUList, optional
        HSA spectral HDULis. If provided, output will contain spectra ("flux")
        and frequency axis ("wave") for the associated `detectors`.
        The default is None.
    detectors : TYPE, optional
        Array detectors to be included in the output. If None, output will
        contain all detectors in the `cat_table`. The default is None.

    Returns
    -------
    output : dictionary
        Dictionary for all detectors meeting the search criteria.
    """
    output = {'header': (cat_table['data'].header,)}
    cat_data = cat_table['data'].data
    cont_data = cont_table['data'].data
    ary_freq = {'SLW': None, 'SSW': None}

    if detectors is None:
        detectors = sorted(set(cat_data['detector']))
    elif isinstance(detectors, str):
        detectors = [detectors]

    for det in detectors:
        ary = 'SLW' if 'SLW' in det else 'SSW'
        prods = np.array([(cat_data['frequency'][i], cat_data['SNR'][i]) for \
                         i, d in enumerate(cat_data['detector']) if d == det])
        features = prods[:, 0]
        snrs = prods[:, 1]
        ind = np.where(cont_data['detector'] == det)[0]
        if len(ind) == 1:
            ind = ind[0]
            continuum_params = cont_data[ind][1::2]
        else:
            print(f"Warning : FF products not found for detector {det}")
            continue

        output[det] = {'lines': {'freq': features, 'snr': snrs},
                       'con_par': continuum_params}
        if spec_hdu is not None:
            spectrum = spec_hdu[det]
            freq, flux = spectrum.data['wave'], spectrum.data['flux']
            if ary_freq[ary] is None:
                ary_freq[ary] = freq
            else:
                freq = ary_freq[ary]
            output[det].update({'wave': freq, 'flux': flux})

    return output

def _get_ffpath(obsid, cal=None, obs_mod=None, off_axis=False, **kwargs):
    """
    Returns the URLs for the postcard, features catalogue and continuum files
    for a particular obsid.

    Parameters
    ----------
    obsid : int
         Herschel Science Archeive observation identification number.
    cal : string or list of strings, optional
        Data calibration type. options: 'point', 'ext'. The default is None.
    obs_mod : string, optional
        Observation sampling mode. options : 'intermediate', 'sparse', 'full'.
        The default is None.
    off_axis : bool, optional
        If True, output will contain off axis FF products paths.
        The default is False.
    kwargs : None
        Not implemented but catches extra kwargs.

    Returns
    -------
    obs_types : string
        Same as input.
    post_files : list of string
        Postcard file paths.
    cat_files : list of string
        Features catalogue file paths.
    cont_files : list of string
        Continuum paramters file paths.
    """
    cal, obs_mod = _init_var(cal, obs_mod)

    sparse_dir = {'pnt': f'{archive}HRpointProducts/',
                  'ext': f'{archive}HRextProducts/',
                  'offAx': f'{archive}HRextProducts-offAxis/'}
    map_dir = f'{archive}HRmapping/'

    obs_types, post_files, cat_files, cont_files = [], [], [], []
    if ('full' in obs_mod) or ('intermediate' in obs_mod):
        obs_types.append('mapping')
        post_files.append(f'{map_dir}postcards/{obsid}_postcard_comb_2x3.png')
        cat_files.append(f'{map_dir}featureCatalogues/{obsid}_featuresFound.fits')
        cont_files.append(f'{map_dir}continuumParameters/{obsid}_fittedContinuumParameters.fits')
    if 'sparse' in obs_mod:
        if off_axis:
            c = 'offAx'
            obs_types.append(f'sparse_{c}')
            post_files.append(f'{sparse_dir[c]}postcards/{obsid}_postcard_{c}.png')
            cat_files.append(f'{sparse_dir[c]}featureCatalogues/{obsid}_featuresFound_{c}.fits')
            cont_files.append(f'{sparse_dir[c]}continuumParameters/{obsid}_fittedContinuumParameters_{c}.fits')
        else:
            for c in cal:
                obs_types.append(f'sparse_{c}')
                if c == 'point': c = 'pnt';
                # PNG postcards are in a hidden folder
                post_files.append(f'{archive}.images/postcards_png/{obsid}_postcard_{c}.png')
                cat_files.append(f'{sparse_dir[c]}featureCatalogues/{obsid}_featuresFound_{c}.fits')
                cont_files.append(f'{sparse_dir[c]}continuumParameters/{obsid}_fittedContinuumParameters_{c}.fits')

    # now check if those files are available
    for i, post_file in enumerate(post_files):
        r = http.request('HEAD', post_file)
        if (r.status != 200):
            print(f"- Warning : postcard file <{post_file.replace(archive, '')}> not found")
            post_files[i] = None
    for i, cat_file in enumerate(cat_files):
        r = http.request('HEAD', cat_file)
        if (r.status != 200):
            print(f"- Warning : featureFinder table <{cat_file.replace(archive, '')}> not found")
            cat_files[i] = None
    for i, cont_file in enumerate(cont_files):
        r = http.request('HEAD', cont_file)
        if (r.status != 200):
            print(f"- Warning : continuum table <{cont_file.replace(archive, '')}> not found")
            cont_files[i] = None

    return obs_types, post_files, cat_files, cont_files

def _extract_spg(tar, cal=None, obs_mod=None, apod=False,
                cube_type=None, photo_band=None, diag=False, 
                **kwargs):
    """
    Extracts fits files meeting the criteria defined by the input parameters
    from the standard product generation (spg) requests response. This funciton
    is used by `get_observation` and should not be called by user directly.

    Parameters
    ----------
    tar : tar file
        Tar file obtained from HSA request response.
    cal : string or list of strings, optional
         Data calibration type. options: 'point', 'ext'. The default is None.
    obs_mod : string, optional
        Observation sampling mode. options : 'intermediate', 'sparse', 'full'.
        Not implemented. The default is None.
    apod : bool, optional
        If True, uses apodized spectra. The default is False.
    cube_type : string or list of strings, optional
        Spectral mapping cube types. options: 'naive', 'convol'
        The default is None.
    photo_band : sting or list of strings, optional
        SPIRE photometer bands. options: 'PLW', 'PMW', 'PSW'
        The default is None.
    diag : bool, optional
        If True, uses diag photometer maps. The default is False.
    kwargs : None
        Not implemented but catches extra kwargs.

    Returns
    -------
    name : list of strings
        Product identifier names.
    gz_files : list of gzip files
        Products for the desired spg fits files.
    """
    apod = '_apod' if apod else ''
    diag = 'diag' if diag else ''
    cal, obs_mod, cube_type, photo_band \
        = _init_var(cal, obs_mod, cube_type, photo_band)

    obs_keys = {'Photometer': {'point': (f'psrc{band}{diag}/' for band in photo_band),
                               'ext': (f'extd{band}{diag}/' for band in photo_band)},

                'sparse': {'point': (f'_point{apod}/',), 'ext': (f'_ext{apod}/',)},

                'full': {'naive': (f'HR_SLW_cube{apod}/', f'HR_SSW_cube{apod}/'),
                        'convol': (f'HR_SLW_cube_convol{apod}/', f'HR_SSW_cube_convol{apod}/')}}

    meta_keys = ['instrument', 'obsid', 'object', 'aot', 'creator',
                 'obsMode', 'mapSampling', 'photObsid000']

    name, gz_files = [], []
    include, exclude = [], []
    for member in tar.getmembers():
        if 'obs_' in  member.name:
            _verboseprint(member.name)
            meta = fits.getheader(gzip.open(tar.extractfile(member)))
            meta_map = _map_meta(meta_keys, meta)
            for key, value in meta_map.items():
                try: _verboseprint(f"--- {key:.<15}{meta[value]}");
                except KeyError: print(f'Warning : Missing {key}!');
            if meta[meta_map['aot']] == 'Spectrometer':
                if meta[meta_map['mapSampling']] == 'sparse':
                    include.extend([element for c in cal for element in
                                    obs_keys['sparse'][c]])
                elif meta[meta_map['mapSampling']] in ['full', 'intermediate']:
                    include.extend([element for c in cube_type for element in
                                    obs_keys['full'][c]])
            elif meta[meta_map['aot']] in ['Photometer', 'Parallel Mode']:
                include.extend([element for c in cal for element in
                                obs_keys['Photometer'][c]])
        if sum([key in member.name for key in include]) > 0 and \
            sum([key in member.name for key in exclude]) == 0:
            _verboseprint(member.name)
            gz_file = tar.extractfile(member)
            gz_files.append(io.BytesIO(gz_file.read()))
            name.append(member.name.split('/')[-2])

    return name, gz_files

def _is_hpdp(obsid):
    """
    Convenience function that returns True if FF used Highly Processed
    Data Products (HPDP) specra for this obsid.
    """
    r = requests.get(hpdp_url)
    i_start = int(r.text.find(str(obsid)))
    if i_start == -1:
        return False
    else:
        return True

def _is_bgs(obsid):
    """
    Convenience function that returns True if FF used BackGround Subctracted
    (BGS) specra for this obsid.
    """
    r = requests.get(bgs_url)
    i_start = int(r.text.find(str(obsid)))
    if i_start == -1:
        return False
    else:
        return True
    
def safecat(save_dir=None, save=False):
    """
    Returns the SpireAutomatedFeatureExtractionCATalogue (SAFECAT).

    Parameters
    ----------
    save_dir : string, optional
        Directory to save the product to. If None, uses working directory.
        The default is None.
    save : bool, optional
        If True, will save product into `save_dir`. The default is False.

    Returns
    -------
    fits HDUList
        SAFECAT product.
        
    Note : Features in the rest frequency catalogue with no associated radial 
        velocity estimate are given as their detected frequency. In most cases
        the radial velocity estimate is derived from 12CO features, and the 
        velocity estimate may not be appropriate for other lines in the 
        same spectrum.
    """
    
    if save_dir is None:
        save_dir = os.getcwd()
    
    save_file = os.path.join(save_dir, 'SAFECAT_V3.fits')
    url = archive + 'SAFECAT_V3.fits.gz'
    _verboseprint("Loading SAFECAT_V3.")

    if os.path.isfile(save_file):
        _verboseprint(f"-- File <{save_file.split('/')[-1]}> already exists... "
               "Will use it.")
        return fits.open(save_file)
    else:
        prod = fits.open(url)
        if save:
            prod.writeto(save_file)
        return prod

def get_observation(obsid, useFF=False, useHsa=False, level=2, spec_type=None,
                    obs_mod='sparse', cal='point', off_axis=False,
                    cube_type='convol', photo_band=('PLW','PMW','PSW'),
                    save=False, save_dir=None):
    """
    Returns FF and/or HSA products for a given obsid based on input parameters.

    Parameters
    ----------
    obsid : int
        SPIRE HSA observation identification number. Can represent both
        spectrometer or photometer observations.
    useFF : bool, optional
        If True, retreives FF products. Other input parameters are overwritten
        by their nominal values based on FF results. The default is False.
    useHsa : bool, optional
        If True, retreives HSA products. The default is False.
    level : int, optional
        Processing level of HSA products. Not tested for non-default values.
        The default is 2.
    spec_type : string, optional
        HSA product spectral type. options: 'spg', 'hpdp', 'bgs'
        The default is None.
    obs_mod : TYPE, optional
        HSA observation sampling mode. options: 'sparse', 'full'.
        The default is 'sparse'.
    cal : string or list of strings, optional
        HSA spectral calibration type. options: 'point', 'ext'
        The default is 'point'.
    off_axis : bool, optional
        If True, output will contain off axis FF products.
        The default is False.
    cube_type : string or list of strings, optional
        HSA spectral cube types. options: 'naive', 'convol'.
        The default is 'convol'.
    photo_band : string or list of strings, optional
        SPIRE photometer bands to include for SPIRE HSA photometer
        observations. The default is ('PLW','PMW','PSW').
    save : bool, optional
        If True, products will be saved to `save_dir`. The default is False.
    save_dir : TYPE, optional
        The directory products will be saved in. If None, will save in working
        directory. The default is None.

    Raises
    ------
    ValueError
        Must use valid `spec_type`.

    Returns
    -------
    output : dictionary
        A dictionary containing the FF and/or HSA products.

    Note : May need to implement 'intermediate' in `obs_mode`.

    """

    obsid = int(obsid)
    output = dict()

    if save_dir is None:
        save_dir = os.getcwd()

    if (spec_type is None) and useFF:
        # use `cwd` to avoid repeaded downloads since user defined `save_dir`
        # will often be an obsid specific directory
        save_files = (os.path.join(save_dir, 'ff-sparse-obs.csv'),
                      os.path.join(save_dir, 'ff-mapping-obs.csv'))
        #save_files = (os.path.join(os.getcwd(), 'ff-sparse-obs.csv'),
        #              os.path.join(os.getcwd(), 'ff-mapping-obs.csv'))
        mod_list_url = (archive + 'hrSparseObservations.csv',
                    archive + 'hrMappingObservations.csv')
        mod_names = ('sparse', 'full')

        for i, file in enumerate(save_files):
            if os.path.isfile(file):
                _verboseprint(f'File <{os.path.basename(file)}> already exists. '
                        'Will use it...')
                df = read_csv(file)
            else:
                df = read_csv(mod_list_url[i])
                if save:
                    df.to_csv(file)
            if obsid in df['obsid'].to_list():
                obs_mod = mod_names[i]
                if obs_mod == 'sparse':
                    spec_type = df['dataUsed'][df['obsid'] == obsid].values[0]
                    ff_cal = df['sourceExt'][df['obsid'] == obsid].values[0]
                    cal = 'point' if ff_cal == 'pointLike' else \
                    'ext' if ff_cal == 'extended' else ['point', 'ext']
                elif obs_mod == 'full':
                    spec_type, cal, off_axis = 'spg', 'ext', False
    elif (spec_type is None) and useHsa:
        spec_type = 'spg'
        
    _verboseprint(f'\n--- {obsid} ---')
    _verboseprint(f'spec_type = {spec_type}\nobs_mod = {obs_mod}\ncal = {cal}')
    cal, obs_mod, spec_types, cube_type, photo_band \
        = _init_var(cal, obs_mod, spec_type, cube_type, photo_band)
    obs_kwargs = {'cal': cal, 'obs_mod': obs_mod, 'cube_type': cube_type,
                  'off_axis': off_axis, 'photo_band': photo_band}
    if useFF:
        _verboseprint('\nRetrieving FF Files.')
        output['ff'] = dict()
        action = {'postcard': imread, 'catalogue': fits.open, 'continuum': fits.open}
        ff_files = list(_get_ffpath(obsid, **obs_kwargs))
        obs_types = ff_files.pop(0)
        for obs_type, *files in zip(obs_types, *ff_files):
            output['ff'][obs_type] = {'postcard': None, 'catalogue': None, 'continuum': None}
            for name, file in zip(action, files):
                if file is None:
                    continue
                func = action[name]
                save_file = os.path.join(save_dir, file.split('/')[-1])
                if os.path.isfile(save_file):
                    _verboseprint(f"-- File <{file.split('/')[-1]}> already exists... Will use it.")
                    output['ff'][obs_type][name] = func(save_file)
                else:
                    product = io.BytesIO(requests.get(file).content)
                    output['ff'][obs_type][name] = func(product)
                    if save:
                        with open(save_file, 'wb') as temp:
                            temp.write(product.getvalue())

    if useHsa:
        _verboseprint('\nRetrieving HSA Files.')
        for spec_type in spec_types:
            if spec_type.lower() in ('spg'):
                tar_file = os.path.join(save_dir, f"{obsid}_level{level}_{spec_type}.tar")
                tar_kwargs = {}
                if os.path.isfile(tar_file):
                    _verboseprint(f"-- Tar file <{tar_file}> already exits... Will use it.")
                    tar_kwargs['name'] = tar_file
                else:
                    request = (f"{spg_url}product.jsp?PROTOCOL=HTTP&"
                               f"OBSERVATION_ID={obsid}&PRODUCT_LEVEL=Level{level}"
                               "&INSTRUMENT=SPIRE")
                    _verboseprint('\nDownloading products from HSA. Please wait...')
                    product = io.BytesIO(requests.get(request).content)
                    tar_kwargs['fileobj'] = product
                    if save:
                        with open(tar_file, 'wb') as temp:
                            temp.write(product.getvalue())
                with tarfile.open(**tar_kwargs) as tar:
                    for name, prod in zip(*_extract_spg(tar, **obs_kwargs)):
                        output[name] = fits.open(gzip.open(prod))
            else:
                if spec_type.lower() in ('hpdp', 'calhpdp'):
                    r = requests.get(hpdp_url)
                    url = hpdp_url
                elif spec_type.lower() in ('bgs'):
                    r = requests.get(bgs_url)
                    url = bgs_url
                else:
                    raise ValueError(f"Error: Invalid spec_type {spec_type}.")
                i_start = int(r.text.find(str(obsid)))
                if i_start == -1:
                    _verboseprint(f"Error: Observation {obsid} not found in {spec_type}")
                    continue
                i_end = i_start + int(r.text[i_start:].find('.gz')) + 3
                save_file = os.path.join(save_dir, r.text[i_start:i_end])
                if os.path.isfile(save_file):
                    _verboseprint(f"-- Files <{save_file.replace(save_dir, '')}> "
                                       f"already exits... Will use it.")
                    output['spectrum'] = fits.open(gzip.open(save_file))
                else:
                    request = url + r.text[i_start:i_end]
                    print("\nDownloading... be patient")
                    gz_file = io.BytesIO(requests.get(request).content)
                    if save:
                        with open(save_file, 'wb') as temp:
                            temp.write(gz_file.getvalue())
                    output['spectrum'] = fits.open(gzip.open(gz_file))

    if verbose:
        _keys_summary(output)
    return output