"""Convenience functions to help locate sources"""

import warnings

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from tqdm import tqdm

from tesswcs import WCS, log, pointings

from .utils import _load_wcs_database


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


def get_observability_mask(wcs: WCS, coords: SkyCoord):
    """For given RAs and Decs, return an array of the same shape with True or False, where True indicates the point falls inside the WCS

    Parameters:
    -----------
    wcs: astropy.wcs.WCS
        Input WCS to check
    RA: np.ndarray
        RAs to check in degrees
    Dec: np.ndarray
        Decs to check in degrees

    Returns:
    --------
    mask: np.ndarray
        Array of booleans with same shape as `RA` and `Dec`. True where `RA` and `Dec` will fall on observable pixels.
    """
    # Firstly, use WCS to calculate a cone around the region of pixels
    # Calculate a "scale" matrix (we don't care about rotation"
    scale = np.asarray(
        [
            [np.hypot(*(wcs.wcs.pc * wcs.wcs.cdelt)[:, 0]), 0],
            [0, np.hypot(*(wcs.wcs.pc * wcs.wcs.cdelt)[:, 1])],
        ]
    )
    # Find the separation from crpix to the corners
    corners = (
        np.asarray(
            [
                [0, 0],
                [wcs._naxis[0], 0],
                [0, wcs._naxis[1]],
                [wcs._naxis[0], wcs._naxis[1]],
            ]
        )
        - wcs.wcs.crpix
    )
    # Calculate the radius in degrees to each corner, take the maximum
    rad = np.hypot(*corners.dot(scale).T).max()
    # This mask is true only where RAs and Decs fall closer than `rad` to the wcs crval.
    k = (coords.separation(SkyCoord(*wcs.wcs.crval, unit="deg")).degree < rad).reshape(
        coords.ra.shape
    )
    if k.any():
        # For those that do fall close, calculate their pixel value
        C, R = wcs._all_world2pix(
            np.asarray([coords.ra.deg[k], coords.dec.deg[k]]).T,
            0,
            tolerance=0.0001,
            maxiter=20,
            adaptive=False,
            detect_divergence=False,
            quiet=True,
        ).T
        # Update the mask with those pixels
        j = (
            (C >= 0)
            & (R >= 44)
            & (C <= (wcs._naxis[0] - 30))
            & (R <= (wcs._naxis[1] - 44))
        )
        k[k] = j
    return k


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
                wcs = WCS.from_sector(sector, camera, ccd)
                observable[sector][camera][ccd] = get_observability_mask(wcs, coords)

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
    wcs_dicts = _load_wcs_database()
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
                    # Pixels are indexed from 1
                    col, row = np.atleast_1d(pix[0]) + 1, np.atleast_1d(pix[1]) + 1

                k = (
                    (row > 0)
                    & (col > 0)
                    & (row < (wcs.pixel_shape[0] + 1))
                    & (col < (wcs.pixel_shape[1] + 1))
                )

                if k.any():
                    k[k] &= np.asarray([wcs.footprint_contains(c) for c in coords[k]])
                    target_ids.append(np.where(k)[0])
                    sectors.append(np.ones(k.sum(), int) * sector)
                    cameras.append(np.ones(k.sum(), int) * camera)
                    ccds.append(np.ones(k.sum(), int) * ccd)
                    rows.append(row[k])
                    columns.append(col[k])

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


def get_interesting_targets(sector, reflimit=50):
    """For a given sector, find any astrophysical targets that have a lot of references.

    This function queries Simbad for targets that fall within the given sector, and
    returns any targets that have up to `reflimit` references.

    Parameters:
    -----------
    sector: int
        Sector number to search
    reflimit: int
        Number of references a target must have to be considered interesting.

    Returns:
    --------
    result: pd.DataFrame
        DataFrame containing all interesting targets.
    """

    from astroquery.simbad import Simbad

    simbad = Simbad()
    simbad.add_votable_fields("nbref", "otypedef", "V", "alltypes")

    dfs = []
    for camera in tqdm(
        np.arange(1, 5), desc="Querying Cameras", position=0, leave=True
    ):
        for ccd in np.arange(1, 5):
            previous_level = log.level
            log.setLevel("CRITICAL")
            wcs = WCS.from_sector(sector, camera, ccd)
            log.setLevel(previous_level)
            coordinate = SkyCoord(*wcs.wcs.crval, unit=("deg", "deg"))
            radius = wcs.pixel_to_world(0, 0).separation(coordinate).deg
            result_table = simbad.query_region(
                coordinate, radius=f"{radius}d", criteria=f"nbref>{reflimit}"
            )
            k = wcs.footprint_contains(
                SkyCoord(*result_table.to_pandas()[["ra", "dec"]].values.T, unit="deg")
            )
            result_table = result_table[k]
            df = result_table.to_pandas()[
                ["main_id", "ra", "dec", "V", "otypedef.otype_longname", "nbref"]
            ].rename(
                {"main_id": "id", "V": "Vmag", "otypedef.otype_longname": "otype"},
                axis="columns",
            )
            df["camera"] = camera
            df["ccd"] = ccd
            df["planet"] = np.asarray(
                [
                    np.any(["Pl" in i0 for i0 in i.split("|")])
                    for i in result_table.to_pandas()["alltypes.otypes"]
                ]
            ).astype(bool)
            dfs.append(df)
    df = pd.concat(dfs).sort_values("nbref", ascending=False).reset_index(drop=True)
    return df
