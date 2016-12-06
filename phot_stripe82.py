# -*- coding: utf-8 -*-
"""

Created on 31/10/2016

@Author: Carlos Eduardo Barbosa

"""
import os
from subprocess import call

import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io.ascii import SExtractor, BaseReader

from config import *


def visual_comparison(im):
    """ Images of the data products compared visually. """
    outfile = os.path.join(plots_dir, "visual", im.replace(".fits", ".png"))
    if os.path.exists(outfile):
        # print "\includegraphics[width=0.98\linewidth]" \
        # "{{figs/visual/{0}}}".format(im.replace(".fits", ".png"))
        return
    fig = plt.figure(1, figsize=(15,5))
    plt.clf()
    imgs = os.path.join(south_data_dir, im)
    imgn = os.path.join(cefca_data_dir, im.replace(".fits", "_prepro.fits"))
    datas = pf.getdata(imgs)
    datan = pf.getdata(imgn)
    h = pf.getheader(imgs)
    title = "Field: {0}, Band: {1}".format(h["object"], h["filter"])
    ps = np.percentile(datas, (1, 99))
    pn = np.percentile(datan, (1, 99))
    ax1 = plt.subplot(1,3,1)
    ax1.get_xaxis().set_ticks([])
    im1 = ax1.imshow(datas, vmin=pn[0], vmax=pn[1], origin="bottom",
                     cmap="cubehelix")
    plt.colorbar(im1, fraction=0.046, pad=0.04)
    ax2 = plt.subplot(1,3,2)
    im2 = ax2.imshow(datan, vmin=pn[0], vmax=pn[1], origin="bottom",
                     cmap="cubehelix")
    plt.colorbar(im2, fraction=0.046, pad=0.04)
    ax2.set_title(title)
    ax3 =  plt.subplot(1,3,3)
    res = 2. * (datas - datan) / (np.sqrt(datas + datan))
    pres = np.percentile(res[np.isfinite(res)], (10, 90))
    presmax = np.nanmax(np.sqrt(pres**2))
    im3 = ax3.imshow(res, vmin=-presmax, vmax=presmax,
                     origin="bottom", cmap="Spectral")
    plt.colorbar(im3, fraction=0.046, pad=0.04)
    xlabels = ["South pipeline", "CEFCA pipeline", "Residuals"]
    for i, ax in enumerate([ax1, ax2, ax3]):
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        ax.set_xlabel(xlabels[i])
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, im.replace(".fits", ".png")))
    plt.close(1)
    return

def run_sextractor(im, redo=False):
    """ Run sextractor for photometry. """
    catalog = im.replace(".fits", ".cat")
    if os.path.exists(catalog) and not redo:
        return
    call(["sextractor", im, "-c", config_files[0], "-PARAMETERS_NAME",
          config_files[1], "-FILTER_NAME", config_files[2], "-CATALOG_NAME",
          catalog])
    return

def load_sextractor_config_files():
    path = os.path.join(os.getcwd(), "tables")
    sdefault = os.path.join(path, "default.sex")
    sparams = os.path.join(path, "default.param")
    sconv = os.path.join(path, "default.conv")
    return (sdefault, sparams, sconv)

def match_catalog(im):
    sex = SExtractor()
    sexcat = sex.read(im.replace(".fits", ".cat"))
    band = pf.getval(im, "filter")
    exptime = pf.getval(im, "exptime")
    airmass = pf.getval(im, "airmass")
    ra = pf.getval(im, "ra")
    dec = pf.getval(im, "dec")
    center = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
    c = SkyCoord(ra=sexcat["ALPHA_J2000"], dec=sexcat["DELTA_J2000"])
    ##########################################################################
    # Reading catalog from file
    stcat = os.path.join(home, "stellar_catalog/Stripe82Stars_Ivezick.fit")
    hdulist = pf.open(stcat)
    catalog = hdulist[1].data
    cradec = SkyCoord(ra=catalog["RAJ2000"] * u.degree,
                      dec=catalog["DEJ2000"] * u.degree)
    # Trimming catalog considering only regions around the observed field
    idx = np.where(cradec.separation(center) < 1.6 * u.deg)[0]
    trimcat = catalog[idx]
    cradec = SkyCoord(ra=trimcat["RAJ2000"] * u.degree,
                      dec=trimcat["DEJ2000"] * u.degree)
    ##########################################################################
    idx, d2d, d3d = c.match_to_catalog_sky(cradec)
    data = np.column_stack((trimcat[idx]["RAJ2000"], trimcat[idx]["DEJ2000"],
                            trimcat[idx]["umag"], trimcat[idx]["e_umag"],
                            trimcat[idx]["gmag"], np.zeros(len(idx)),
                            trimcat[idx]["rmag"], trimcat[idx]["e_rmag"],
                            trimcat[idx]["imag"], trimcat[idx]["e_imag"],
                            trimcat[idx]["zmag"], trimcat[idx]["e_zmag"],
                            sexcat["MAG_AUTO"] + 2.5 * np.log10(exptime),
                           sexcat["MAGERR_AUTO"], np.ones_like(idx) * airmass))
    i = np.where(d2d < 1 * u.arcsec)[0]
    data = data[i]
    ###########################################################################
    # Cleaning data from photometry with large uncertainties
    data = data[data[:,13] < .1]
    ##########################################################################
    h = ["RAJ2000", "DEJ2000", "umag", "e_umag", "gmag", "e_gmag", "rmag",
         "e_rmag", "imag", "e_imag", "zmag", "e_zmag", "mag_inst",
         "e_mag_inst", "airmass"]
    newcat = "{}_{}".format(band, im.split("-")[2]).replace(".fits", ".cat")
    with open(newcat, "w") as f:
        f.write("#" + " ".join(h) + "\n")
        np.savetxt(f, data)
    return

def run_splus():
    os.chdir(south_data_dir)
    filenames = np.loadtxt(os.path.join(south_data_dir,
                            "fields_with_calibration_stars.txt"), dtype="S")
    for im in filenames:
        visual_comparison(im)
        run_sextractor(im)
        match_catalog(im)

def run_cefca():
    ###########################################################################
    def replace_header(im, redo=False):
        """ Change the header from original produced by CEFCA by the one
        produced by the SPLUS collaboration. This is necessary to run
        SExtractor and Swarp. """
        imgin = os.path.join(cefca_data_dir + "_bk", im)
        imgout = os.path.join(cefca_data_dir, im)
        imgh = os.path.join(south_data_dir, im.replace("_prepro", ""))
        if not os.path.exists(imgout) or redo:
            h = pf.getheader(imgh)
            data = pf.getdata(imgin)
            hdu = pf.PrimaryHDU(data, header=h)
            hdu.writeto(imgout, clobber=True)
   ###########################################################################
    filenames = np.loadtxt(os.path.join(south_data_dir,
                                        "fields_with_calibration_stars.txt"),
                           dtype="S")
    filenames = [x.replace(".fits", "_prepro.fits") for x in filenames]
    os.chdir(cefca_data_dir)
    for im in filenames:
        print im
        replace_header(im)
        run_sextractor(im)
        match_catalog(im)


if __name__ == "__main__":
    config_files = load_sextractor_config_files()
    run_splus()
    run_cefca()

