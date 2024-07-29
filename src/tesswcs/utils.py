"""Utilities for getting WCS"""

import bz2
import json
import os

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import Angle, SkyCoord
from tqdm import tqdm

from . import PACKAGEDIR, log, pixel_corners, pixel_scale, pointings, rcolumns, rrows


def angle_to_matrix(theta):
    """Convert an angle theta in degrees to a rotation matrix"""
    cos_angle = np.cos(np.deg2rad(theta))
    sin_angle = np.sin(np.deg2rad(theta))
    return np.array([[cos_angle, -sin_angle], [sin_angle, cos_angle]])


def angle_from_matrix(matrix):
    """Convert a rotation matrix to an angle in degrees"""
    theta = np.rad2deg(np.arctan2(matrix[1, 0], matrix[0, 0]))
    return theta


def footprint(npoints=50, rows=2078, columns=2136):
    """Gets the column and row points for CCD edges"""
    column = np.hstack(
        [
            np.zeros(npoints),
            np.linspace(0, rows, npoints),
            np.linspace(0, rows, npoints),
            np.ones(npoints) * rows,
        ]
    )
    row = np.hstack(
        [
            np.linspace(0, columns, npoints),
            np.zeros(npoints),
            np.ones(npoints) * columns,
            np.linspace(0, columns, npoints),
        ]
    )
    return np.vstack([row, column]).T


def get_M(truth, approx):
    """Finds the transformation matrix that makes `approx` best fit `truth`"""
    truth = np.hstack([truth, np.ones((truth.shape[0], 1))])
    approx = np.hstack([approx, np.ones((approx.shape[0], 1))])
    a, b, e = np.linalg.solve(truth.T.dot(truth), truth.T.dot(approx[:, 0]))
    c, d, f = np.linalg.solve(truth.T.dot(truth), truth.T.dot(approx[:, 1]))

    M = np.linalg.inv(np.asarray([[a, b, e], [c, d, f], [0, 0, 1]]))
    return M


def _build_wcs_data():
    # This package isn't publicly available yet...
    try:
        from tessrip import Rip
    except ImportError:
        raise ImportError("Must have `tessrip` installed to build the wcs dictionary.")

    sep_dicts = {}
    for cycle, sector, ra, dec, roll, start, end in pointings.iterrows():
        sector = int(sector)
        sep_dicts[sector] = {}
        sep_dicts[sector] = {"ra": ra, "dec": dec, "roll": roll}
        c0 = SkyCoord(ra, dec, unit="deg", frame="icrs")
        for cam in tqdm(
            np.arange(1, 5), leave=True, position=0, desc=f"sector {sector}"
        ):
            cam = int(cam)
            sep_dicts[sector][cam] = {}
            wcss = []
            for ccd in np.arange(1, 5):
                wcss.append(Rip(sector, cam, ccd).wcs)
            for idx in range(4):
                crpix = [1045, 1001]
                crval = SkyCoord(
                    *wcss[idx].wcs_pix2world([crpix], 0)[0],
                    unit="deg",
                )
                sep = c0.separation(crval)
                pa = c0.position_angle(crval)
                # r = sep.deg / pixel_scale
                # phi = pa.rad
                # x, y = r * np.cos(phi), r * np.sin(phi)

                crval_sip = SkyCoord(
                    *wcss[idx].all_pix2world([crpix], 0)[0],
                    unit="deg",
                )
                sep_sip = c0.separation(crval_sip)
                pa_sip = c0.position_angle(crval_sip)

                sep_dicts[sector][cam][idx + 1] = {
                    "pa": pa.rad,
                    "sep": sep.deg,
                    "pa_sip": pa_sip.rad,
                    "sep_sip": sep_sip.deg,
                    "ccd_center_crval": [crval.ra.deg, crval.dec.deg],
                    "ccd_center_crval_sip": [crval_sip.ra.deg, crval_sip.dec.deg],
                    "ccd_center_crpix": crpix,
                    "crval0": list(wcss[idx].wcs.crval),
                    "crpix0": list(wcss[idx].wcs.crpix),
                    "cd": [list(wcss[idx].wcs.cd[jdx]) for jdx in range(2)],
                    "sip_a_order": wcss[idx].sip.a_order,
                    "sip_a": [
                        list(wcss[idx].sip.a[jdx])
                        for jdx in range(wcss[idx].sip.a.shape[0])
                    ],
                    "sip_b": [
                        list(wcss[idx].sip.b[jdx])
                        for jdx in range(wcss[idx].sip.b.shape[0])
                    ],
                    "sip_ap": [
                        list(wcss[idx].sip.ap[jdx])
                        for jdx in range(wcss[idx].sip.ap.shape[0])
                    ],
                    "sip_bp": [
                        list(wcss[idx].sip.bp[jdx])
                        for jdx in range(wcss[idx].sip.bp.shape[0])
                    ],
                }
                corners = SkyCoord(
                    [
                        SkyCoord(*wcss[idx].wcs_pix2world([[0, 0]], 0)[0], unit="deg"),
                        SkyCoord(
                            *wcss[idx].wcs_pix2world([[0, rcolumns]], 0)[0], unit="deg"
                        ),
                        SkyCoord(
                            *wcss[idx].wcs_pix2world([[rrows, 0]], 0)[0], unit="deg"
                        ),
                        SkyCoord(
                            *wcss[idx].wcs_pix2world([[rrows, rcolumns]], 0)[0],
                            unit="deg",
                        ),
                    ]
                )
                sep_dicts[sector][cam][idx + 1]["corner_pa"] = list(
                    c0.position_angle(corners).rad
                )
                sep_dicts[sector][cam][idx + 1]["corner_sep"] = list(
                    c0.separation(corners).deg
                )
                sep_dicts[sector][cam][idx + 1]["corner_ra"] = list(corners.ra.deg)
                sep_dicts[sector][cam][idx + 1]["corner_dec"] = list(corners.dec.deg)

                corners = SkyCoord(
                    [
                        SkyCoord(*wcss[idx].all_pix2world([[0, 0]], 0)[0], unit="deg"),
                        SkyCoord(
                            *wcss[idx].all_pix2world([[0, rcolumns]], 0)[0], unit="deg"
                        ),
                        SkyCoord(
                            *wcss[idx].all_pix2world([[rrows, 0]], 0)[0], unit="deg"
                        ),
                        SkyCoord(
                            *wcss[idx].all_pix2world([[rrows, rcolumns]], 0)[0],
                            unit="deg",
                        ),
                    ]
                )
                sep_dicts[sector][cam][idx + 1]["corner_sip_pa"] = list(
                    c0.position_angle(corners).rad
                )
                sep_dicts[sector][cam][idx + 1]["corner_sip_sep"] = list(
                    c0.separation(corners).deg
                )
                sep_dicts[sector][cam][idx + 1]["corner_sip_ra"] = list(corners.ra.deg)
                sep_dicts[sector][cam][idx + 1]["corner_sip_dec"] = list(
                    corners.dec.deg
                )

    filename = f"{PACKAGEDIR}/data/TESS_wcs_data.json.bz2"
    json_data = json.dumps(sep_dicts)
    with bz2.open(filename, "wt", encoding="utf-8") as f:
        f.write(json_data)
    _build_support_dicts()


def _load_wcs_data():
    log.debug("Loading WCS data from file")
    filename = f"{PACKAGEDIR}/data/TESS_wcs_data.json.bz2"
    with bz2.open(filename, "rt", encoding="utf-8") as f:
        wcs_dict = json.load(f)
    log.debug("Adding correct units to WCS data")
    _wcs_dicts = {}
    for sector in wcs_dict.keys():
        _wcs_dicts[int(sector)] = {
            "ra": wcs_dict[sector]["ra"],
            "dec": wcs_dict[sector]["dec"],
            "roll": wcs_dict[sector]["roll"],
        }
        for camera in "1234"[::-1]:
            _wcs_dicts[int(sector)][int(camera)] = {}
            for ccd in "1234"[::-1]:
                _wcs_dicts[int(sector)][int(camera)][int(ccd)] = {
                    "pa": Angle(wcs_dict[sector][camera][ccd]["pa"] * u.rad),
                    "sep": Angle(wcs_dict[sector][camera][ccd]["sep"] * u.deg),
                    "pa_sip": Angle(wcs_dict[sector][camera][ccd]["pa_sip"] * u.rad),
                    "sep_sip": Angle(wcs_dict[sector][camera][ccd]["sep_sip"] * u.deg),
                    "ccd_center_crval": wcs_dict[sector][camera][ccd][
                        "ccd_center_crval"
                    ],
                    "ccd_center_crval_sip": wcs_dict[sector][camera][ccd][
                        "ccd_center_crval_sip"
                    ],
                    "ccd_center_crpix": wcs_dict[sector][camera][ccd][
                        "ccd_center_crpix"
                    ],
                    "crval0": wcs_dict[sector][camera][ccd]["crval0"],
                    "crpix0": wcs_dict[sector][camera][ccd]["crpix0"],
                    "cd": np.asarray(wcs_dict[sector][camera][ccd]["cd"]),
                    "sip_a_order": wcs_dict[sector][camera][ccd]["sip_a_order"],
                    "sip_a": np.asarray(wcs_dict[sector][camera][ccd]["sip_a"]),
                    "sip_b": np.asarray(wcs_dict[sector][camera][ccd]["sip_b"]),
                    "sip_ap": np.asarray(wcs_dict[sector][camera][ccd]["sip_ap"]),
                    "sip_bp": np.asarray(wcs_dict[sector][camera][ccd]["sip_bp"]),
                    "corner_pa": Angle(
                        wcs_dict[sector][camera][ccd]["corner_pa"] * u.rad
                    ),
                    "corner_sep": Angle(
                        wcs_dict[sector][camera][ccd]["corner_sep"] * u.deg
                    ),
                    "corner": SkyCoord(
                        wcs_dict[sector][camera][ccd]["corner_ra"],
                        wcs_dict[sector][camera][ccd]["corner_dec"],
                        unit="deg",
                    ),
                    "corner_sip_pa": Angle(
                        wcs_dict[sector][camera][ccd]["corner_pa"] * u.rad
                    ),
                    "corner_sip_sep": Angle(
                        wcs_dict[sector][camera][ccd]["corner_sep"] * u.deg
                    ),
                    "corner_sip": SkyCoord(
                        wcs_dict[sector][camera][ccd]["corner_sip_ra"],
                        wcs_dict[sector][camera][ccd]["corner_sip_dec"],
                        unit="deg",
                    ),
                }
    return _wcs_dicts


def _build_support_dicts():
    """This function creates dictionaries of the best fitting CCD coordinates and SIP polynomials"""
    _wcs_dicts = _load_wcs_data()
    # Contains the corners of each CCD w.r.t the boresight
    xs, ys = {
        1: {1: [], 2: [], 3: [], 4: []},
        2: {1: [], 2: [], 3: [], 4: []},
        3: {1: [], 2: [], 3: [], 4: []},
        4: {1: [], 2: [], 3: [], 4: []},
    }, {
        1: {1: [], 2: [], 3: [], 4: []},
        2: {1: [], 2: [], 3: [], 4: []},
        3: {1: [], 2: [], 3: [], 4: []},
        4: {1: [], 2: [], 3: [], 4: []},
    }

    # Contains the center of each CCD w.r.t the boresight
    xcent, ycent = {
        1: {1: [], 2: [], 3: [], 4: []},
        2: {1: [], 2: [], 3: [], 4: []},
        3: {1: [], 2: [], 3: [], 4: []},
        4: {1: [], 2: [], 3: [], 4: []},
    }, {
        1: {1: [], 2: [], 3: [], 4: []},
        2: {1: [], 2: [], 3: [], 4: []},
        3: {1: [], 2: [], 3: [], 4: []},
        4: {1: [], 2: [], 3: [], 4: []},
    }
    log.debug("Building camera/CCD data.")
    for camera in np.arange(1, 5):
        for ccd in np.arange(1, 5):
            for sector in _wcs_dicts.keys():
                phi, r = (
                    _wcs_dicts[sector][camera][ccd]["corner_pa"]
                    + _wcs_dicts[sector]["roll"] * u.deg,
                    _wcs_dicts[sector][camera][ccd]["corner_sep"].deg,
                )
                x, y = r * np.cos(phi.rad), r * np.sin(phi.rad)

                xs[camera][ccd].append(x)
                ys[camera][ccd].append(y)

                phi, r = (
                    _wcs_dicts[sector][camera][ccd]["pa"]
                    + _wcs_dicts[sector]["roll"] * u.deg,
                    _wcs_dicts[sector][camera][ccd]["sep"].deg,
                )
                x, y = r * np.cos(phi.rad), r * np.sin(phi.rad)
                xcent[camera][ccd].append(x)
                ycent[camera][ccd].append(y)

            xs[int(camera)][int(ccd)] = list(
                np.mean(xs[camera][ccd], axis=0)
            )  # / pixel_scale
            ys[int(camera)][int(ccd)] = list(
                np.mean(ys[camera][ccd], axis=0)
            )  # / pixel_scale

            # This needs to be the crpix of the camera
            xcent[int(camera)][int(ccd)] = np.mean(xcent[camera][ccd])  # / pixel_scale
            ycent[int(camera)][int(ccd)] = np.mean(ycent[camera][ccd])  # / pixel_scale
    log.debug("Building SIP data.")
    sip = {}
    for attr in ["a", "b", "ap", "bp"]:
        sip[str(attr)] = {
            int(camera): {
                int(ccd): np.mean(
                    [
                        _wcs_dicts[int(sector)][int(camera)][int(ccd)][f"sip_{attr}"]
                        for sector in _wcs_dicts.keys()
                        if _wcs_dicts[sector][camera][ccd]["sip_a_order"] == 4
                    ],
                    axis=0,
                ).tolist()
                for ccd in np.arange(1, 5)
            }
            for camera in np.arange(1, 5)
        }
    M = {1: {}, 2: {}, 3: {}, 4: {}}
    for camera in np.arange(1, 5):
        for ccd in np.arange(1, 5):
            truth = np.asarray(
                [
                    xs[camera][ccd] - xcent[camera][ccd],
                    ys[camera][ccd] - ycent[camera][ccd],
                ]
            )
            approx = ((pixel_corners - [1045, 1001])).T
            M[int(camera)][int(ccd)] = get_M(truth.T, approx.T).tolist()

    for var in ["xs", "ys", "xcent", "ycent", "sip", "M"]:
        filename = f"{PACKAGEDIR}/data/TESS_wcs_{var}.json.bz2"
        json_data = json.dumps(locals()[var])
        with bz2.open(filename, "wt", encoding="utf-8") as f:
            f.write(json_data)
    return


def _fix_keys(dict):
    return {
        int(key): np.asarray(item) if isinstance(item, list) else item
        for key, item in dict.items()
    }


def _load_support_dicts():
    log.debug("Loading support dictionaries from file")

    support_dicts = []

    for var in ["xs", "ys", "xcent", "ycent", "M"]:
        filename = f"{PACKAGEDIR}/data/TESS_wcs_{var}.json.bz2"
        with bz2.open(filename, "rt", encoding="utf-8") as f:
            support_dict = json.load(f)
        support_dict = {int(key): _fix_keys(dict) for key, dict in support_dict.items()}
        support_dicts.append(support_dict)

    filename = f"{PACKAGEDIR}/data/TESS_wcs_sip.json.bz2"
    with bz2.open(filename, "rt", encoding="utf-8") as f:
        sip = json.load(f)
    sip_dict = {
        key: _fix_keys({int(key): _fix_keys(dict) for key, dict in sip_dict.items()})
        for key, sip_dict in sip.items()
    }
    support_dicts.append(sip_dict)
    return support_dicts


def _load_warp_matrices():
    log.debug("Loading warp matrices dictionary from file")
    filename = f"{PACKAGEDIR}/data/TESS_wcs_Ms.json.bz2"
    if not os.path.isfile(filename):
        log.warn(
            "No warp matrices found. Either download them or fit them using `tesswcs.tesswcs._build_warp_matrices`"
        )
        return None, None
    with bz2.open(filename, "rt", encoding="utf-8") as f:
        Ms = json.load(f)
    Ms = {int(key): _fix_keys(dict) for key, dict in Ms.items()}
    filename = f"{PACKAGEDIR}/data/TESS_wcs_offset_weights.json.bz2"
    if not os.path.isfile(filename):
        log.warn(
            "No warp matrices found. Either download them or fit them using `tesswcs.tesswcs._build_warp_matrices`"
        )
        return None, None
    with bz2.open(filename, "rt", encoding="utf-8") as f:
        offset_weights = json.load(f)
    offset_weights = {int(key): _fix_keys(dict) for key, dict in offset_weights.items()}
    return Ms, offset_weights


def plot_geometry(ax=None):
    # This lets us reorganize four coordinate corners to plot a square
    s = [0, 1, 3, 2, 0]

    xs, ys, xcent, ycent, M, _ = _load_support_dicts()
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(5, 5))
    for camera in np.arange(1, 5):
        for ccd in np.arange(1, 5):
            ax.plot(
                ys[camera][ccd][s] * pixel_scale,
                xs[camera][ccd][s] * pixel_scale,
                c=f"C{camera-1}",
            )
            ax.scatter(
                ys[camera][ccd].mean() * pixel_scale,
                xs[camera][ccd].mean() * pixel_scale,
                c=f"C{camera-1}",
            )
    lim = 60 * pixel_scale
    ax.set(
        xlim=(-lim, lim),
        ylim=(-lim, lim),
        xlabel="Column [pixels]",
        ylabel="Row [pixels]",
        title="TESS Observatory Geometry",
    )
    return ax
