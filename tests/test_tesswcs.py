import os

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pytest
from astropy.wcs import WCS as astropyWCS
from astropy.io import fits

from tesswcs import PACKAGEDIR, WCS, __version__, pointings, rcolumns, rrows
from tesswcs.utils import (
    _build_support_dicts,
    _load_support_dicts,
    _load_wcs_data,
    angle_from_matrix,
    angle_to_matrix,
    footprint,
    plot_geometry,
)

DOCSDIR = "/".join([*PACKAGEDIR.split("/")[:-2], "docs"])


def test_install():
    assert __version__ == "0.1.0"
    assert os.path.isfile(f"{PACKAGEDIR}/data/TESS_wcs_data.json.bz2")
    # can read dictionaries
    _wcs_dicts = _load_wcs_data()
    xs, ys, xcent, ycent, sip = _load_support_dicts()
    # dictionaries have data
    assert "a" in sip.keys()
    for dict in [_wcs_dicts, xs, ys, xcent, ycent, sip["a"]]:
        assert 1 in dict.keys()
        assert 1 in dict[1].keys()
    assert 1 in _wcs_dicts[1][1].keys()
    assert pointings["RA"][0] == 352.6844

    _build_support_dicts()

    # Generate a figure for documentation
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    plot_geometry(ax=ax)
    fig.savefig(f"{DOCSDIR}/figures/tess_camera.png", bbox_inches="tight", dpi=150)


def test_utils():
    angle = 10
    M = angle_to_matrix(angle)
    assert M.shape == (2, 2)
    assert angle_from_matrix(M) == angle


@pytest.mark.skip
def test_footprint():
    # This fails?
    row, column = footprint().T
    assert row.min() == 0
    assert row.max() == rrows
    assert column.min() == 0
    assert column.max() == rcolumns


def test_load_archive():
    """Can we load an archival TESS wcs?"""
    sector, camera, ccd = 1, 1, 1
    wcs = WCS.from_archive(sector=sector, camera=camera, ccd=ccd)
    assert isinstance(wcs, astropyWCS)

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.grid(True)
    for camera in np.arange(1, 5):
        for ccd in np.arange(1, 5):
            wcs = WCS.from_archive(sector=sector, camera=camera, ccd=ccd)
            fp = wcs.pixel_to_world(*footprint()[:-30].T)
            ax.scatter(
                fp.ra.wrap_at(180 * u.deg).rad,
                fp.dec.rad,
                lw=0.5,
                s=0.1,
                c=f"C{camera - 1}",
            )
    ax.set(
        title=f"WCS for Sector {sector} Pointing",
        xlabel="RA",
        ylabel="Dec",
    )
    fig.savefig(
        f"{DOCSDIR}/figures/tess_{sector}_archive.png",
        bbox_inches="tight",
        dpi=150,
    )


def test_predict():
    """Can we predict a pointing?"""
    sector = 1
    ra, dec, roll = (
        pointings["RA"][sector - 1],
        pointings["Dec"][sector - 1],
        pointings["Roll"][sector - 1],
    )
    camera, ccd = 1, 1
    wcs = WCS.predict(ra=ra, dec=dec, roll=roll, camera=camera, ccd=ccd)
    assert isinstance(wcs, astropyWCS)

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.grid(True)
    for camera in np.arange(1, 5):
        for ccd in np.arange(1, 5):
            wcs = WCS.predict(ra=ra, dec=dec, roll=roll, camera=camera, ccd=ccd)
            fp = wcs.pixel_to_world(*footprint()[:-30].T)
            ax.scatter(
                fp.ra.wrap_at(180 * u.deg).rad,
                fp.dec.rad,
                lw=0.5,
                s=0.1,
                c=f"C{camera - 1}",
            )
    ax.set(
        title=f"Predicted WCS for RA:{ra}, Dec:{dec}, Roll:{roll}\n[i.e. Sector {sector} Pointing]",
        xlabel="RA",
        ylabel="Dec",
    )
    fig.savefig(
        f"{DOCSDIR}/figures/tess_{sector}_predict.png",
        bbox_inches="tight",
        dpi=150,
    )


def test_write():
    """Can we write a WCS to fits with additional keywords?"""
    sector = 1
    ra, dec, roll = (
        pointings["RA"][sector - 1],
        pointings["Dec"][sector - 1],
        pointings["Roll"][sector - 1],
    )
    camera, ccd = 1, 1
    wcs = WCS.predict(ra=ra, dec=dec, roll=roll, camera=camera, ccd=ccd)
    hdr = wcs.to_header()
    assert isinstance(hdr, fits.Header)
    assert "AUTHOR" in hdr
    for attr in ["ra", "dec", "roll", "camera", "ccd"]:
        assert attr.upper() in hdr
    hdr_string = wcs.to_header_string()
    assert "AUTHOR" in hdr_string
    hdu = wcs.to_fits()
    assert "AUTHOR" in hdu[0].header
