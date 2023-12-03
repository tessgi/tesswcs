"""Convenience functions to help locate sources"""
import warnings

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time

from tesswcs import WCS, pointings

from . import wcs_dicts


def _process_input_parameters(coords, sector=None, cycle=None, time=None):
    if not isinstance(coords, (SkyCoord, list)):
        raise ValueError("Must pass a SkyCoord object")
    if isinstance(coords, list):
        if not isinstance(coords[0], SkyCoord):
            raise ValueError("Must pass a SkyCoord object")
        coords = SkyCoord(coords)
    if coords.isscalar:
        coords = SkyCoord([coords])

    num_args_set = sum(arg is not None for arg in [sector, cycle, time])
    if num_args_set > 1:
        raise ValueError("Only one of 'sector', 'cycle', or 'time' can be set")

    sector_times = [
        Time(pointings["Start"], format="jd"),
        Time(pointings["End"], format="jd"),
    ]
    sector_mask = np.ones(len(pointings), bool)
    if time is not None:
        sector_mask = (time > sector_times[0]) & (time < sector_times[1])
    elif sector is not None:
        sector_mask = np.in1d(pointings["Sector"], sector)
    elif cycle is not None:
        sector_mask = np.in1d(pointings["Cycle"], cycle)

    if sector_mask.sum() == 0:
        raise ValueError(
            "No past or planned TESS observations found in that 'sector', 'cycle', or 'time'"
        )

    return coords, sector_mask


def check_observability(
    coords: SkyCoord, sector: int = None, cycle: int = None, time: Time = None
) -> Table:
    """Checks whether an input set of SkyCoord objects are observable by TESS

    Parameters:
    -----------
    coords: astropy.coordinates.SkyCoord or List
        Set of sky coordinates to test.
    sector: int, list of ints, optional
        Optional list of sectors to narrow down search to
    cycle: int, list of ints, optional
        Optional list of cycles to narrow down search to
    time: astropy.time.Time
        Optional time to narrow down search. If a time is passed, will only search sector that encompasses that time.
    """
    coords, sector_mask = _process_input_parameters(coords, sector, cycle, time)

    observable = {}
    for sector, ra, dec, roll in pointings[sector_mask][
        ["Sector", "RA", "Dec", "Roll"]
    ]:
        observable[sector] = {}
        for camera in np.arange(1, 5):
            observable[sector][camera] = {}
            for ccd in np.arange(1, 5):
                observable[sector][camera][ccd] = []
                # wcs = WCS.predict(ra, dec, roll, camera, ccd)
                if sector in wcs_dicts.keys():
                    wcs = WCS.from_archive(sector, camera, ccd)
                else:
                    wcs = WCS.predict(ra, dec, roll, camera, ccd)
                for coord in coords:
                    onsilicon = wcs.footprint_contains(coord)
                    observable[sector][camera][ccd].append(onsilicon)

    sector = np.hstack(np.asarray(list(observable.keys()))[:, None] * np.ones(16, int))
    camera = np.hstack(
        [
            np.hstack(np.asarray(list(item.keys()))[:, None] * np.ones(4, int))
            for key, item in observable.items()
        ]
    )
    ccd = np.hstack(
        (
            [
                np.hstack(
                    [
                        [ccdkey for ccdkey, ccditem in camitem.items()]
                        for camkey, camitem in sectoritem.items()
                    ]
                )
                for sectorkey, sectoritem in observable.items()
            ]
        )
    )
    obs = np.vstack(
        [
            np.vstack(
                [
                    np.vstack([ccditem for ccdkey, ccditem in camitem.items()])
                    for camkey, camitem in sectoritem.items()
                ]
            )
            for sectorkey, sectoritem in observable.items()
        ]
    )
    observable = Table(
        [sector, camera, ccd, *obs.T],
        names=[
            "Sector",
            "Camera",
            "CCD",
            *[f"targ_{idx + 1:04}" for idx in range(len(coords))],
        ],
    )
    return observable


def get_pixel_locations(
    coords: SkyCoord, sector: int = None, cycle: int = None, time: Time = None
) -> Table:
    """Obtains the pixel locations of any sources in `coords` that fall on a TESS pixel.

    Parameters:
    -----------
    coords: astropy.coordinates.SkyCoord or List
        Set of sky coordinates to test.
    sector: int, list of ints, optional
        Optional list of sectors to narrow down search to
    cycle: int, list of ints, optional
        Optional list of cycles to narrow down search to
    time: astropy.time.Time
        Optional time to narrow down search. If a time is passed, will only search sector that encompasses that time.
    """
    coords, sector_mask = _process_input_parameters(coords, sector, cycle, time)
    target_ids, sectors, cameras, ccds, rows, columns = [], [], [], [], [], []
    for sector, ra, dec, roll in pointings[sector_mask][
        ["Sector", "RA", "Dec", "Roll"]
    ]:
        for camera in np.arange(1, 5):
            for ccd in np.arange(1, 5):
                # wcs = WCS.predict(ra, dec, roll, camera, ccd)
                if sector in wcs_dicts.keys():
                    wcs = WCS.from_archive(sector, camera, ccd)
                else:
                    wcs = WCS.predict(ra, dec, roll, camera, ccd)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    pix = wcs.world_to_pixel(coords)
                    row, col = np.atleast_1d(pix[0]), np.atleast_1d(pix[1])
                k = (
                    (row > 0)
                    & (col > 0)
                    & (row < (wcs.pixel_shape[0] + 1))
                    & (col < (wcs.pixel_shape[1] + 1))
                )
                if k.any():
                    target_ids.append(np.where(k)[0])
                    sectors.append(np.ones(k.sum(), int) * sector)
                    cameras.append(np.ones(k.sum(), int) * camera)
                    ccds.append(np.ones(k.sum(), int) * ccd)
                    rows.append(np.round(row[k], 1))
                    columns.append(np.round(col[k], 1))
    if len(target_ids) == 0:
        return Table(
            None, names=["Target Index", "Sector", "Camera", "CCD", "Row", "Column"]
        )
    target_ids, sectors, cameras, ccds, rows, columns = (
        np.hstack(target_ids),
        np.hstack(sectors),
        np.hstack(cameras),
        np.hstack(ccds),
        np.hstack(rows),
        np.hstack(columns),
    )
    pixel_locations = Table(
        [target_ids, sectors, cameras, ccds, rows, columns],
        names=["Target Index", "Sector", "Camera", "CCD", "Row", "Column"],
    )
    return pixel_locations
