#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from . import io as spio

from os import path
from pickle import load
from lmfit import Model
import warnings

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Polygon

from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

hipe_cm = load(open(path.join(path.dirname(__file__), "hipe_cm.pickle"), "rb"))

def _sinc(x, amp, cent, width):
    """
    Sinc or sine over argument function.

    Parameters
    ----------
    x : numpy array
        Independent variable.
    amp : float
        Amplitude of the sinc profile.
    cent : float
        Center of the sinc profile.
    width : float
        Width of the sinc profile.

    Returns
    -------
    model : numpy array
        Evaluated since function.
    """
    p2 = width/np.pi
    arg = (x - cent)/p2
    with warnings.catch_warnings():
        # ignore division by zero warning
        warnings.simplefilter("ignore")
        model = amp*np.sin(arg)/arg
    model[model != model] = amp
    return model

def SPIRE_ILS(bandwidth=None, res=spio.resolution['HR']):
    """
    Model of the SPIRE Instrument Line Shape (ILS), which is nominally a
    sinc profile.

    Parameters
    ----------
    bandwidth : iterable of length 2, optional
        Limits of the independent varialbe. If None, the limits of the SPIRE
        interferometer are used. The default is None.
    res : float, optional
        Resolution of the sinc profile. SPIRE operated in High Resolution (HR)
        and Low Resolution (LR) modes for science observations.
        The default is spio.resolution['HR'].

    Returns
    -------
    freq : numpy array
        Frequency axis in increments of `res`/pi.
    ils : numpy array
        Model ILS.

    """
    step = res/np.pi
    if bandwidth is None:
        SPIREBandwidth = [spio.band_SLW[0], spio.band_SSW[-1]]
        freq = np.arange(SPIREBandwidth[0], SPIREBandwidth[1], step)
    else:
        freq = np.arange(bandwidth[0],bandwidth[1], step)
    ils = _sinc(freq, 1.0, freq[int(len(freq)/2.0)], res)
    return freq, ils

def model_spectrum(products, res=spio.resolution['HR'], vary_width=False):
    """
    Model SPIRE spectra processed by the FF using the FF products. Essentially
    replicates the fitted models obtained by the FF routine.

    Parameters
    ----------
    products : dictionary of dictionaries
        A dictionary with FF and HSA products. Key structure should be:
            -header
            -ID_1 (sup product identifieer eg, "SLWC3", "SSW[4,3]")
            --lines
            ---freq (FF line frequencies)
            ---snr (FF snrs for lines)
            --con_par (FF fitted continuum parameters)
            --flux (HSA spectral flux of ID)
            --wave (HSA frequency axis of ID)
        Multiple 'ID' dictionaries can be provided in `products`. This
        structure is provided automatically by the `parse_sparse` and
        `parse_mapping(..., flatten=True)` functions in the io library.
    res : float, optional
        Resolution of sinc line shapes. The default is spio.resolution['HR'].
    vary_width : bool, optional
        If True, will allow the widths of the since models to vary. This may
        produce better results for partially resolved lines. Note the correct
        approach is to model such featues as Gaussian convolved sincs.
        The default is False.

    Returns
    -------
    products : dictionary of dictionaries
        Modified version of input `products`. Eeach 'ID' dictionary will
        receive a numpy array for the modeled continuum 'cont' and a
        lmfit.model.ModelResult 'model'.

    Note : This function fits the continuum subtracted flux. The continuum
        is determined by the FF continuum parameters in `products`.

    Caution : Accurate flux determination requires additional processing of
        the input spectrum before fitting occures.
    """
    for key, value in products.items():
        if key == 'header':
            continue
        features = value['lines']['freq']
        freq, flux = value['wave'], value['flux']
        continuum = np.polyval(value['con_par'][::-1], freq)
        flux_subcont = flux - continuum
        if len(features) == 0:
            print(f'func : `model_spectrum` {key} has no features. No model '
                  'will be added.')
            products[key].update({'cont': continuum})
            continue

        prefs = []
        totalModel = None
        for i, line in enumerate(features):
            ind = np.where(min(abs(line - freq)) == abs(line - freq))[0]
            pref = f"_{i}"
            prefs.append(pref)
            amp = flux_subcont[ind][0]
            newMod = Model(_sinc, prefix=pref)
            newMod.set_param_hint(pref+'amp', value=amp, vary=True)
            newMod.set_param_hint(pref+'cent', value=line, vary=False)
            # set vary to True for partially resolved features
            newMod.set_param_hint(pref+'width', value=res, vary=vary_width)
            if totalModel == None:
                totalModel = newMod
            else:
                totalModel += newMod

        products[key].update({'cont': continuum,
                            'model': totalModel.fit(flux_subcont, x=freq)})

    return products

def sky_image(hdu, foot_data=None, foot_color=None, detectors=None,
              annotate=False, pixels=None, fig=None, index=111, foot_kwargs=None, 
              fig_kwargs=None, ax_kwargs=None, img_kwargs=None, txt_kwargs=None):
    """
    Produces a sky image with SPIRE photometer data and adds footprint
    overlays from SPIRE spectrometer data.

    Parameters
    ----------
    hdu : fits HDUList
        SPIRE photometer HDUList.
    foot_data : fits HDUList or list of HDULists, optional
        HDULists for SPIRE HSA spectral data. The default is None.
    foot_color : dictionary, optional
        Dictionary of SPIRE band specific colors. Format is
        {'SLW': color, 'SSW': color}. The default is None.
    detectors : list, optional
        Detectors to use when `foot_data` has single pointing products.
        If None, all detectors in the products will be used.
        The default is None.
    annotate : bool, optional
        If True, detectors in array will be annotated. The default is False.
    pixels : dictionary, optional
        Dictionary of spectral cube pixels to include in footprint. 
        Format is {'SLW': [[cols], [rows]], 'SSW': [[cols], [rows]]} for 
        multiple pixels or {'SLW': [col, row], 'SSW': [col, row]} for a 
        single pixel. The default is None.
    fig : matplotlib figure, optional
        Plot axis will be added to this figure if provided. If None, a new
        figure will be created. The default is None.
    index : valid subplot index, optional
        Used to place plot at a particular location in the figure.
        The default is 111.
    foot_kwargs : dictionary, optional
        Additional footprint overlay kwargs passed to matplotlib patches
        constructor. The default is None.
    fig_kwargs : dictionary, optional
        Additional matplotlib figure kwargs passed to `figure` constructor.
        If None, the image colormap will be set to the HIPE default.
        The default is None.
    ax_kwargs : dictionary, optional
        Additional matplotlib axis kwargs passed to `add_subplot` constructor.
        The default is None.
    img_kwargs : dictionary, optional
        Additional matplotlib image kwargs passed to `imshow` constructor.
        The default is None.
    img_kwargs : dictionary, optional
        Additional matplotlib annotate kwargs passed to `annotate` constructor.
        The default is None.

    Returns
    -------
    fig : matplotlib Figure
        Resulting figure object.
    ax : matplotlib WCSAxesSubplot
        Resulting axes object.
    img : matplotlib AxesImage
        Resulting image objec.
    wcs : astropy WCS
        WCS object from input `hdu` photometer image.

    Note : All aspect of the output figure can be adjusted form the returned
        variables.
    """

    if ax_kwargs is None:
        ax_kwargs = dict()
    if fig_kwargs is None:
        fig_kwargs = dict()
    if foot_kwargs is None:
        foot_kwargs = dict()
    if txt_kwargs is None:
        txt_kwargs = dict()
    if img_kwargs is None:
        img_kwargs = {'cmap': hipe_cm}
    else:
        img_kwargs.setdefault('cmap', hipe_cm)

    header_img = hdu['Primary'].header
    units = header_img['SIGUNIT']
    image = hdu['image'].data
    wcs = WCS(header_img)
    if fig is None:
        fig = plt.figure(**fig_kwargs)
    ax = fig.add_subplot(index, projection=wcs, **ax_kwargs)
    img = ax.imshow(image, origin='lower', **img_kwargs)
    fig.colorbar(img, ax=ax, pad=0.0, label=f'[{units}]')

    if foot_data is not None:
        if type(foot_data) not in [list, tuple]:
            foot_data = [foot_data]
        if foot_color is None:
            foot_color = {'SLW': 'red', 'SSW': 'blue'}
        for i, foot_hdu in enumerate(foot_data):
            top_header = foot_hdu[0].header
            if top_header['MAPSAMPL'] in ['full', 'intermediate']:
                ra, dec = (top_header['RA_NOM'], top_header['DEC_NOM'])
                ary = top_header['DETECTOR']
                wcs_cube = WCS(foot_hdu['image']).dropaxis(2)
                # account for RA distortion from Dec
                dec_correction = 1.0/np.cos(np.deg2rad(dec))
                cube_limits = wcs_cube.calc_footprint() \
                    + np.array([[1,-1],[1,1],[-1,1],[-1,-1]]) \
                        * np.array([dec_correction, 1]) \
                            * wcs_cube.wcs.cdelt[1]/2.0
                # main cube footprint
                polygon = Polygon(cube_limits, fc='none',
                                  ec=foot_color[ary],
                                  transform=ax.get_transform('icrs'),
                                  label=ary, **foot_kwargs)
                ax.add_patch(polygon)
                
                # add pixels
                if (pixels is not None) and (ary in pixels):
                    pix_ = np.transpose(pixels[ary]).reshape(-1, 2)
                    d1 = np.array([0.5, -0.5])
                    d2 = np.array([0.5, 0.5])
                    lowerleft = wcs_cube.wcs_pix2world(pix_-d2, 0)
                    lowerright = wcs_cube.wcs_pix2world(pix_+d1, 0)
                    upperleft = wcs_cube.wcs_pix2world(pix_-d1, 0)
                    upperright = wcs_cube.wcs_pix2world(pix_+d2, 0)
                    for ll, ul, ur, lr in zip(lowerleft, upperleft, upperright, lowerright):
                        polygon = Polygon([ll, ul, ur, lr], fc='none',
                                        ec=foot_color[ary],
                                        transform=ax.get_transform('icrs'),
                                        **foot_kwargs)
                        ax.add_patch(polygon)                      

                # determin zoom in limits of image
                if i == 0:
                    offset = 1.5*spio.beam['SLW']
                    cent_coord = SkyCoord(ra, dec, frame='icrs', unit='deg')
                    cent_pix = np.round([*wcs.world_to_pixel(cent_coord)],
                                        0).astype(int)
                    img.set_clim(vmax=np.nanmax(
                                image[cent_pix[1]-5:cent_pix[1]+5,
                                      cent_pix[0]-5:cent_pix[0]+5])*1.2)
                    delta_ra, delta_dec = 0.0, 0.0
                    
                # account for RA distortion from Dec
                ra_offset = offset*dec_correction
                ra_max_, dec_min_ = (cube_limits[0] \
                                   + np.array([ra_offset, -offset]))
                ra_min_, dec_max_ = (cube_limits[2] \
                                   + np.array([-ra_offset, offset]))
                if dec_max_-dec_min_ > delta_dec:
                    dec_max, dec_min = dec_max_, dec_min_
                    delta_dec = dec_max_-dec_min_
                if ra_max_-ra_min_ > delta_ra:
                    ra_max, ra_min = ra_max_, ra_min_
                    delta_ra - ra_max_-ra_min_

            elif top_header['MAPSAMPL'] == 'sparse':
                arrays = {'SLW': dict(), 'SSW': dict(),
                          'beam': spio.beam}
                for table in foot_hdu:
                    try:
                        header = table.header
                        det = header['CHNLNAME']
                        ra = header['RA']
                        dec = header['DEC']
                        if 'SLW' in det:
                            arrays['SLW'][det] = [ra, dec]
                        elif 'SSW' in det:
                            arrays['SSW'][det] = [ra, dec]
                    except:
                        continue
                if detectors is None:
                    detectors = [det for ary in ('SLW', 'SSW') \
                                 for det in arrays[ary].keys()]
                for ary in ['SLW', 'SSW']:
                    set_label = True
                    for det, value in arrays[ary].items():
                        if det in detectors:
                            height = arrays['beam'][ary]
                            ellipse = Ellipse(value, height/np.cos(np.deg2rad(value[1])),
                                           height, ec=foot_color[ary], fc='none',
                                           transform=ax.get_transform('icrs'),
                                           **foot_kwargs)
                            if set_label: 
                                ellipse.set_label(ary)
                                set_label = False
                            ax.add_patch(ellipse)
                            if annotate:
                                ax.annotate(det[-2:], value+np.array([0, spio.beam[ary]/2]), 
                                            ha='center', va='bottom', 
                                            xycoords=ax.get_transform('icrs'), 
                                            c=foot_color[ary], **txt_kwargs)

                offset = 4.0*spio.beam['SLW']
                ra, dec = arrays['SLW']['SLWC3']
                # account for RA distortion from Dec
                ra_offset = offset/np.cos(np.deg2rad(dec))
                ra_min, ra_max = (ra - ra_offset), (ra + ra_offset)
                dec_min, dec_max = (dec - offset), (dec + offset)

        bottom_left_sky = SkyCoord(ra_max, dec_min, frame='icrs', unit='deg')
        top_right_sky = SkyCoord(ra_min, dec_max, frame='icrs', unit='deg')

        bottom_left_pix = wcs.world_to_pixel(bottom_left_sky)
        top_right_pix = wcs.world_to_pixel(top_right_sky)

        ax.set_xlim(bottom_left_pix[0], top_right_pix[0])
        ax.set_ylim(bottom_left_pix[1], top_right_pix[1])
        ax.set_xlabel('RA (J2000)', labelpad=0.4)
        ax.set_ylabel('Dec (J2000)', labelpad=-0.5)

    return fig, ax, img, wcs

def sparse_postcard(products, fig=None, index=111, fig_kwargs=None):
    """
    Generates a plot with a similar appearance as the FF postcards.

    Parameters
    ----------
    products : dictionary of dictionaries
        A dictionary with FF and HSA products. Key structure should be:
            -header
            -ID_1 (sup product identifieer eg, "SLWC3", "SSW[4,3]")
            --lines
            ---freq (FF line frequencies)
            ---snr (FF snrs for lines)
            --con_par (FF fitted continuum parameters)
            --flux (HSA spectral flux of ID)
            --wave (HSA frequency axis of ID)
        Multiple 'ID' dictionaries can be provided in `products`. This
        structure is provided automatically by the `parse_sparse` and
        `parse_mapping(..., flatten=True)` functions in the io library.
    fig : matplotlib figure, optional
        Plot axis will be added to this figure if provided. If None, a new
        figure will be created. The default is None.
    index : valid subplot index, optional
        Used to place plot at a particular location in the figure.
        The default is 111.
    fig_kwargs : dictionary, optional
        Additional matplotlib figure kwargs passed to `figure` constructor.
        The default is None.

    Returns
    -------
    fig : matplotlib Figure
        Resulting figure object.
    ax : matplotlib WCSAxesSubplot
        Resulting axes object.
    """

    if fig_kwargs is None:
        fig_kwargs = dict()

    header = products['header'][0]
    edge = header['EDGEMASK']
    source = header['OBJECT']
    try: flux_unit = header['FLXUNIT']
    except:
        fl = header.comments['MAX_CONT']
        flux_unit = fl[fl.find('[')+1:fl.find(']')]
    obsid = header['OBS_ID']

    colors = {'SLW': {'spec': (0.55, 0.0, 0.0), 'bar': (1.0, 0.0, 0.0),
                      'stripe': (0.5, 0.0, 0.0, 0.2)},
              'SSW': {'spec': (0.0, 0.0, 0.55), 'bar': (0.0, 0.7, 1.0),
                      'stripe': (0.0, 0.0, 0.5, 0.2)},
              'continuum': 'lime'}
    overlap = (spio.band_SSW[0], spio.band_SLW[-1])

    keys = list(products.keys())[1:]
    flux_lim = np.array([min([min(products[key]['flux']) for key in keys]),
            max([max(products[key]['flux']) for key in keys])])
    offset = (flux_lim[1] - flux_lim[0]) * 0.1
    ymin = flux_lim[0]

    if fig is None:
        fig = plt.figure(**fig_kwargs)
    ax = fig.add_subplot(index)

    for key, value in products.items():
        if key == 'header':
            continue
        ary = 'SLW' if 'SLW' in key else 'SSW'
        features = value['lines']['freq']
        snrs = products[key]['lines']['snr']
        freq, flux = value['wave'], value['flux']
        continuum = np.polyval(value['con_par'][::-1], freq)   
        
        ax.axvspan(freq[0], freq[0]+edge, color=colors[ary]['stripe'], ec=None, zorder=0)
        ax.axvspan(freq[-1]-edge, freq[-1], color=colors[ary]['stripe'], ec=None, zorder=0)

        ax.plot(freq, flux, color=colors[ary]['spec'], label=f'{key} ({len(features)})')
        if len(features) > 0:
            for i, line in enumerate(features):
                ind = np.where(min(abs(line - freq)) == abs(line - freq))[0]
                amp = flux[ind][0]
                sg = np.sign(snrs[i])
                ht = sg*np.log(abs(snrs[i]))*offset
                bars = (amp + sg*offset, amp+ht)
                if sg < 0.0:
                    ymin = flux_lim[0] - offset
                if ary == 'SLW' and line > overlap[0]:
                    ax.plot(np.ones(2)*line, bars, color=colors[ary]['bar'], linewidth=1.0)
                else:
                    ax.plot(np.ones(2)*line, bars, color=colors[ary]['bar'])

        ax.plot(freq, continuum, color=colors['continuum'])

    ylim = (ymin, flux_lim[1]*1.3)

    title = f'{source} Observation ID {obsid}'
    xlabel = 'Frequency [GHz]'
    ylabel = f'Flux Density [{flux_unit}]' if flux_unit == 'Jy' \
                else f'Brightness [{flux_unit}]'
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(430.0, 1580.0)
    ax.set_ylim(ylim)
    ax.set_xticks(np.arange(500, 1600, 100))
    ax.legend(loc=2)

    return fig, ax

