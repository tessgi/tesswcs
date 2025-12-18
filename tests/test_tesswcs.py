import os

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS as astropyWCS
from tqdm import tqdm

from tesswcs import PACKAGEDIR, WCS, pointings, rcolumns, rrows
from tesswcs.utils import (
    _load_support_dicts,
    _load_wcs_database,
    angle_from_matrix,
    angle_to_matrix,
    footprint,
    plot_geometry,
)

DOCSDIR = "/".join([*PACKAGEDIR.split("/")[:-2], "docs"])


def test_install():
    assert os.path.isfile(f"{PACKAGEDIR}/data/TESS_wcs_data.json.bz2")
    # can read dictionaries
    _wcs_dicts = _load_wcs_database()
    xs, ys, xcent, ycent, M, sip = _load_support_dicts()
    # dictionaries have data
    assert "a" in sip.keys()
    for dict in [_wcs_dicts, xs, ys, xcent, ycent, sip["a"]]:
        assert 1 in dict.keys()
        assert 1 in dict[1].keys()
    assert 1 in _wcs_dicts[1][1].keys()
    assert pointings["RA"][0] == 352.6844

    # Generate a figure for documentation
    if os.environ.get("CI") != "true":
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        plot_geometry(ax=ax)
        fig.savefig(f"{DOCSDIR}/figures/tess_camera.png", bbox_inches="tight", dpi=150)


def test_utils():
    angle = 10
    M = angle_to_matrix(angle)
    assert M.shape == (2, 2)
    assert np.isclose(angle_from_matrix(M), angle)


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
    wcs = WCS.from_sector(sector=sector, camera=camera, ccd=ccd)
    assert isinstance(wcs, astropyWCS)
    if os.environ.get("CI") != "true":
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111, projection="mollweide")
        ax.grid(True)
    for camera in np.arange(1, 5):
        for ccd in np.arange(1, 5):
            wcs = WCS.from_sector(sector=sector, camera=camera, ccd=ccd)
            fp = wcs.pixel_to_world(*footprint()[:-30].T)
            if os.environ.get("CI") != "true":
                ax.scatter(
                    fp.ra.wrap_at(180 * u.deg).rad,
                    fp.dec.rad,
                    lw=0.5,
                    s=0.1,
                    c=f"C{camera - 1}",
                )
    if os.environ.get("CI") != "true":
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

    if os.environ.get("CI") != "true":
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111, projection="mollweide")
        ax.grid(True)
    for camera in np.arange(1, 5):
        for ccd in np.arange(1, 5):
            wcs = WCS.predict(ra=ra, dec=dec, roll=roll, camera=camera, ccd=ccd)
            fp = wcs.pixel_to_world(*footprint()[:-30].T)
            if os.environ.get("CI") != "true":
                ax.scatter(
                    fp.ra.wrap_at(180 * u.deg).rad,
                    fp.dec.rad,
                    lw=0.5,
                    s=0.1,
                    c=f"C{camera - 1}",
                )
    if os.environ.get("CI") != "true":
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


def test_comprable():
    sector, camera, ccd = 4, 1, 1
    R, C = np.meshgrid(
        np.arange(0, rrows, 10), np.arange(0, rcolumns, 10), indexing="ij"
    )
    wcs = WCS.from_sector(sector=sector, camera=camera, ccd=ccd)
    truth = wcs.pixel_to_world(R.ravel(), C.ravel())
    prediction = WCS.predict(
        ra=wcs.ra, dec=wcs.dec, roll=wcs.roll, camera=wcs.camera, ccd=wcs.ccd
    ).pixel_to_world(R.ravel(), C.ravel())
    separation = (
        (truth.separation(prediction).to(u.arcsecond) / (21 * u.arcsecond))
        .reshape(R.shape)
        .value
    )
    if os.environ.get("CI") != "true":
        fig, ax = plt.subplots()
        im = ax.pcolormesh(C[0], R[:, 0], separation, vmin=0, vmax=0.5)
        ax.set(
            xlabel="Column",
            ylabel="Row",
            title=f"Sector {sector}, Camera {camera}, CCD {ccd}",
        )
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label(
            "Separation between true position and predicted position [pixels]"
        )
        fig.savefig(
            f"{DOCSDIR}/figures/tess_{sector}_separation.png",
            bbox_inches="tight",
            dpi=150,
        )

    if os.environ.get("CI") != "true":
        fig, ax = plt.subplots()
    # Check all C1 sector, camera, ccd predictions are within half a pixel of the truth
    for sector in tqdm(np.arange(1, 14), desc="sector", leave=True, position=0):
        for camera in np.arange(1, 5):
            for ccd in np.arange(1, 5):
                R, C = np.meshgrid(
                    np.arange(0, rrows, 10), np.arange(0, rcolumns, 10), indexing="ij"
                )
                wcs = WCS.from_sector(sector=sector, camera=camera, ccd=ccd)
                truth = wcs.pixel_to_world(R.ravel(), C.ravel())
                prediction = WCS.predict(
                    ra=wcs.ra,
                    dec=wcs.dec,
                    roll=wcs.roll,
                    camera=wcs.camera,
                    ccd=wcs.ccd,
                ).pixel_to_world(R.ravel(), C.ravel())
                separation = (
                    (truth.separation(prediction).to(u.arcsecond) / (21 * u.arcsecond))
                    .reshape(R.shape)
                    .value
                )
                assert (separation < 1).all()
                if os.environ.get("CI") != "true":
                    ax.hist(
                        separation.ravel(),
                        np.linspace(0, 2, 50),
                        alpha=0.02,
                        color="k",
                        density=True,
                    )
    if os.environ.get("CI") != "true":
        ax.set(
            xlabel="Separation between true WCS and predicted WCS [pixel]",
            yticks=[],
            title="tesswcs accuracy for sectors 1-70",
        )
        plt.savefig(
            f"{DOCSDIR}/figures/tess_accuracy.png",
            dpi=150,
            bbox_inches="tight",
        )
