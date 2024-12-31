"""Functions for building the database of inputs"""

import bz2
import json
import os

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import Angle, SkyCoord
from botocore.exceptions import ClientError
from tqdm import tqdm

from . import PACKAGEDIR, log, pixel_corners, pointings, rcolumns, rrows
from .utils import _load_wcs_database, get_M


def update_database():
    _build_wcs_database()
    _build_warp_matrices()


def _build_support_dicts():
    """This function creates dictionaries of the best fitting CCD coordinates and SIP polynomials"""
    wcs_dicts = _load_wcs_database()
    # Contains the corners of each CCD w.r.t the boresight
    xs, ys = (
        {
            1: {1: [], 2: [], 3: [], 4: []},
            2: {1: [], 2: [], 3: [], 4: []},
            3: {1: [], 2: [], 3: [], 4: []},
            4: {1: [], 2: [], 3: [], 4: []},
        },
        {
            1: {1: [], 2: [], 3: [], 4: []},
            2: {1: [], 2: [], 3: [], 4: []},
            3: {1: [], 2: [], 3: [], 4: []},
            4: {1: [], 2: [], 3: [], 4: []},
        },
    )

    # Contains the center of each CCD w.r.t the boresight
    xcent, ycent = (
        {
            1: {1: [], 2: [], 3: [], 4: []},
            2: {1: [], 2: [], 3: [], 4: []},
            3: {1: [], 2: [], 3: [], 4: []},
            4: {1: [], 2: [], 3: [], 4: []},
        },
        {
            1: {1: [], 2: [], 3: [], 4: []},
            2: {1: [], 2: [], 3: [], 4: []},
            3: {1: [], 2: [], 3: [], 4: []},
            4: {1: [], 2: [], 3: [], 4: []},
        },
    )
    log.debug("Building camera/CCD data.")
    for camera in np.arange(1, 5):
        for ccd in np.arange(1, 5):
            for sector in wcs_dicts.keys():
                phi, r = (
                    wcs_dicts[sector][camera][ccd]["corner_pa"]
                    + wcs_dicts[sector]["roll"] * u.deg,
                    wcs_dicts[sector][camera][ccd]["corner_sep"].deg,
                )
                x, y = r * np.cos(phi.rad), r * np.sin(phi.rad)

                xs[camera][ccd].append(x)
                ys[camera][ccd].append(y)

                phi, r = (
                    wcs_dicts[sector][camera][ccd]["pa"]
                    + wcs_dicts[sector]["roll"] * u.deg,
                    wcs_dicts[sector][camera][ccd]["sep"].deg,
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
                        wcs_dicts[int(sector)][int(camera)][int(ccd)][f"sip_{attr}"]
                        for sector in wcs_dicts.keys()
                        if wcs_dicts[sector][camera][ccd]["sip_a_order"] == 4
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
            approx = (pixel_corners - [1045, 1001]).T
            M[int(camera)][int(ccd)] = get_M(truth.T, approx.T).tolist()

    for var in ["xs", "ys", "xcent", "ycent", "sip", "M"]:
        filename = f"{PACKAGEDIR}/data/TESS_wcs_{var}.json.bz2"
        json_data = json.dumps(locals()[var])
        with bz2.open(filename, "wt", encoding="utf-8") as f:
            f.write(json_data)
    return


def _add_sector_data(sector, c0):
    try:
        from tesscube import TESSCube
    except ImportError:
        raise ImportError(
            "Must have `tesscube` installed to build the wcs dictionary. "
            "If you installed `tesswcs` with pip make sure you include the database extras "
            "`pip install tesswcs[database]`."
        )
    sector_dict = {}
    for cam in np.arange(1, 5):
        cam = int(cam)
        sector_dict[cam] = {}
        wcss = []
        for ccd in np.arange(1, 5):
            try:
                wcss.append(TESSCube(sector, cam, ccd).wcs)
            except ClientError:
                # Data does not exist in cloud for this sector.
                return None
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

            sector_dict[cam][idx + 1] = {
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
                    SkyCoord(*wcss[idx].wcs_pix2world([[rrows, 0]], 0)[0], unit="deg"),
                    SkyCoord(
                        *wcss[idx].wcs_pix2world([[rrows, rcolumns]], 0)[0],
                        unit="deg",
                    ),
                ]
            )
            sector_dict[cam][idx + 1]["corner_pa"] = list(
                c0.position_angle(corners).rad
            )
            sector_dict[cam][idx + 1]["corner_sep"] = list(c0.separation(corners).deg)
            sector_dict[cam][idx + 1]["corner_ra"] = list(corners.ra.deg)
            sector_dict[cam][idx + 1]["corner_dec"] = list(corners.dec.deg)

            corners = SkyCoord(
                [
                    SkyCoord(*wcss[idx].all_pix2world([[0, 0]], 0)[0], unit="deg"),
                    SkyCoord(
                        *wcss[idx].all_pix2world([[0, rcolumns]], 0)[0], unit="deg"
                    ),
                    SkyCoord(*wcss[idx].all_pix2world([[rrows, 0]], 0)[0], unit="deg"),
                    SkyCoord(
                        *wcss[idx].all_pix2world([[rrows, rcolumns]], 0)[0],
                        unit="deg",
                    ),
                ]
            )
            sector_dict[cam][idx + 1]["corner_sip_pa"] = list(
                c0.position_angle(corners).rad
            )
            sector_dict[cam][idx + 1]["corner_sip_sep"] = list(
                c0.separation(corners).deg
            )
            sector_dict[cam][idx + 1]["corner_sip_ra"] = list(corners.ra.deg)
            sector_dict[cam][idx + 1]["corner_sip_dec"] = list(corners.dec.deg)
    return sector_dict


def _build_wcs_database():
    filename = f"{PACKAGEDIR}/data/TESS_wcs_data.json.bz2"
    if os.path.isfile(filename):
        with bz2.open(filename, "rt", encoding="utf-8") as f:
            sep_dicts = json.load(f)
    else:
        sep_dicts = {}
    for cycle, sector, ra, dec, roll, start, end in tqdm(
        pointings.iterrows(),
        desc="Loading TESS Sector",
        total=len(pointings),
        leave=True,
        position=0,
    ):
        if f"{sector}" in sep_dicts.keys():
            continue
        else:
            sector = int(sector)
            c0 = SkyCoord(ra, dec, unit="deg", frame="icrs")
            sep_dict = _add_sector_data(sector, c0)
            if sep_dict is not None:
                sep_dicts[sector] = sep_dict
                sep_dicts[sector]["ra"] = ra
                sep_dicts[sector]["dec"] = dec
                sep_dicts[sector]["roll"] = roll
            else:
                # No more sector data exists
                break
    json_data = json.dumps(sep_dicts)
    with bz2.open(filename, "wt", encoding="utf-8") as f:
        f.write(json_data)
    _build_support_dicts()


def _build_warp_matrices(plot=False):
    """Find the best fitting warp matrices

    This function uses all the archival TESS data to build a small warp matrix and updated
    CRPIX value for the WCS for each camera and CCD.

    This function works by building the true, archival WCS, and then comparing with the
    predicted WCS. We find the sky projected points in each case, and then find the transformation
    matrix that makes the predicted sky coordinates match the true sky coordinates. We then translate
    those warped coordinates back to pixel space using the predicted WCS. We can fit the warp between the
    new pixel space and the true pixel grid, to find the small change that needs to be applied to our PC matrix.

    We also find the change in CRPIX (offset term) that is best fit for each pointing.

    This term is a strong function of the roll angle w.r.t the ecliptic, so we specify a delta CRPIX as a function of the ecliptic roll.
    """
    from .tesswcs import WCS

    R, C = np.meshgrid(
        np.arange(0, rrows, 10), np.arange(0, rcolumns, 10), indexing="ij"
    )

    Ms = {
        1: {1: [], 2: [], 3: [], 4: []},
        2: {1: [], 2: [], 3: [], 4: []},
        3: {1: [], 2: [], 3: [], 4: []},
        4: {1: [], 2: [], 3: [], 4: []},
    }

    wcs_dicts = _load_wcs_database()

    for sector in tqdm(wcs_dicts.keys()):
        if plot:
            fig, ax = plt.subplots(4, 4, figsize=(10, 10), sharex=True, sharey=True)
        for camera in np.arange(1, 5):
            for ccd in np.arange(1, 5):
                wcs_t = WCS.from_sector(sector=sector, camera=camera, ccd=ccd)
                truth = wcs_t.wcs_pix2world(np.asarray([R.ravel(), C.ravel()]).T, 0)

                wcs_c = WCS.predict(
                    ra=wcs_t.ra,
                    dec=wcs_t.dec,
                    roll=wcs_t.roll,
                    camera=wcs_t.camera,
                    ccd=wcs_t.ccd,
                    warp=False,
                )
                prediction = wcs_c.wcs_pix2world(
                    np.asarray([R.ravel(), C.ravel()]).T, 0
                )

                if (wcs_t.wcs.crval[0] > 90) & (wcs_t.wcs.crval[0] < 270):
                    wrap = 360 * u.deg
                else:
                    wrap = 180 * u.deg
                truth = np.asarray(
                    [Angle(truth[:, 0] * u.deg).wrap_at(wrap).deg, truth[:, 1]]
                ).T
                prediction = np.asarray(
                    [
                        Angle(prediction[:, 0] * u.deg).wrap_at(wrap).deg,
                        prediction[:, 1],
                    ]
                ).T

                M = get_M(truth, prediction)
                transformation = M.dot(
                    np.asarray(
                        [prediction[:, 0], prediction[:, 1], prediction[:, 0] ** 0]
                    )
                )[:2].T
                pix = wcs_c.wcs_world2pix(transformation, 0)
                M2 = np.linalg.inv(get_M(np.asarray([R.ravel(), C.ravel()]).T, pix))

                wcs_c.wcs.pc = wcs_c.wcs.pc.dot(M2[:2, :2])
                prediction = wcs_c.wcs_pix2world(
                    np.asarray([R.ravel(), C.ravel()]).T, 0
                )
                prediction = np.asarray(
                    [
                        Angle(prediction[:, 0] * u.deg).wrap_at(wrap).deg,
                        prediction[:, 1],
                    ]
                ).T
                M2[:2, 2] = (
                    wcs_c.wcs_world2pix([prediction.mean(axis=0)], 0)[0]
                    - wcs_c.wcs_world2pix([truth.mean(axis=0)], 0)[0]
                )

                wcs_c.wcs.crpix += M2[:2, 2]

                prediction = wcs_c.wcs_pix2world(
                    np.asarray([R.ravel(), C.ravel()]).T, 0
                )
                prediction = np.asarray(
                    [
                        Angle(prediction[:, 0] * u.deg).wrap_at(wrap).deg,
                        prediction[:, 1],
                    ]
                ).T
                if plot:
                    sep = np.hypot(*(truth - prediction).T).reshape(R.shape)
                    ax[camera - 1][ccd - 1].pcolormesh(
                        R, C, sep, vmin=0, vmax=3 * 21 / 3600
                    )
                    ax[camera - 1][ccd - 1].set(title=f"Camera {camera}, CCD {ccd}")
                Ms[int(camera)][int(ccd)].append(M2.copy().tolist())

    ra = np.asarray([wcs_dicts[sdx]["ra"] for sdx in wcs_dicts.keys()])
    dec = np.asarray([wcs_dicts[sdx]["dec"] for sdx in wcs_dicts.keys()])
    roll = np.asarray([wcs_dicts[sdx]["roll"] for sdx in wcs_dicts.keys()])

    c = SkyCoord(ra, dec, unit="deg")
    ecl = c.transform_to("geocentricmeanecliptic")

    # The roll angle w.r.t the ecliptic plane
    ecl_roll = ecl.position_angle(
        c.directional_offset_by(-roll * u.deg, 1 * u.arcsecond).transform_to(
            "geocentricmeanecliptic"
        )
    )

    X = np.asarray([np.sin(ecl_roll), np.cos(ecl_roll), roll**0]).T

    offset_weights = {
        int(camera): {
            int(ccd): np.linalg.solve(
                X.T.dot(X), X.T.dot(np.asarray(Ms[camera][ccd])[:, :2, 2])
            ).tolist()
            for ccd in np.arange(1, 5)
        }
        for camera in np.arange(1, 5)
    }

    for camera in np.arange(1, 5):
        for ccd in np.arange(1, 5):
            Ms[int(camera)][int(ccd)] = np.median(Ms[camera][ccd], axis=0)[
                :2, :2
            ].tolist()

    filename = f"{PACKAGEDIR}/data/TESS_wcs_Ms.json.bz2"
    json_data = json.dumps(Ms)
    with bz2.open(filename, "wt", encoding="utf-8") as f:
        f.write(json_data)

    filename = f"{PACKAGEDIR}/data/TESS_wcs_offset_weights.json.bz2"
    json_data = json.dumps(offset_weights)
    with bz2.open(filename, "wt", encoding="utf-8") as f:
        f.write(json_data)
