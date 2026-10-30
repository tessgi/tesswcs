"""Work with WCS objects for TESS"""

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS as astropyWCS
from astropy.wcs import Sip

from . import log  # noqa: E402
from . import CRPIX_EPOCHS, get_epoch, pixel_corners, pointings, rcolumns, rrows
from .utils import (
    _load_support_dicts,
    _load_warp_matrices,
    _load_wcs_database,
    angle_to_matrix,
    deprecated,
    get_M,
)


class WCS(astropyWCS):
    """A special subclass of astropy.wcs.WCS

    This class allows us to add attributes to the WCS header and add class methods
    """

    wcs_dicts = _load_wcs_database()

    (
        xs,
        ys,
        xcent,
        ycent,
        M,
        sip_dict,
    ) = _load_support_dicts()

    Ms, offset_weights = _load_warp_matrices()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    @deprecated("from_sector")
    def from_archive(cls, sector: int, camera: int, ccd: int):
        """Load a WCS from archival TESS data"""
        return cls.from_sector(sector=sector, camera=camera, ccd=ccd)

    @classmethod
    def from_sector(cls, sector: int, camera: int, ccd: int):
        """Load a WCS from archival TESS data"""
        if sector in list(cls.wcs_dicts.keys()):
            epoch = get_epoch(sector)
            wcs = cls(naxis=2)
            wcs.sector, wcs.camera, wcs.ccd = sector, camera, ccd
            wcs.ra = cls.wcs_dicts[sector]["ra"]
            wcs.dec = cls.wcs_dicts[sector]["dec"]
            wcs.roll = cls.wcs_dicts[sector]["roll"]
            wcs.pixel_shape = (rcolumns, rrows)
            wcs.wcs.ctype = ["RA---TAN-SIP", "DEC--TAN-SIP"]
            wcs.wcs.cunit = ["deg", "deg"]
            wcs.wcs.radesys = "ICRS"
            wcs.wcs.crpix = cls.wcs_dicts[sector][camera][ccd]["crpix_hdr"]
            wcs.wcs.crval = cls.wcs_dicts[sector][camera][ccd]["crval_hdr"]
            wcs.wcs.cdelt = [1, 1]
            wcs.wcs.pc = cls.wcs_dicts[sector][camera][ccd]["cd"]
            wcs.sip = Sip(
                *[
                    cls.sip_dict[attr][epoch][camera][ccd]
                    for attr in ["a", "b", "ap", "bp"]
                ],
                cls.wcs_dicts[sector][camera][ccd]["crpix_hdr"],
            )
            return wcs
        else:
            if sector not in pointings["Sector"]:
                raise ValueError(
                    f"Sector {sector} does not yet have predicted pointing information."
                )
            log.warning(
                f"Data for Sector {sector} has not been archived yet. This function will return the predicted WCS for Sector {sector}."
            )
            ra, dec, roll = [
                i
                for i in pointings[pointings["Sector"] == sector][
                    ["RA", "Dec", "Roll"]
                ][0].values()
            ]
            wcs = cls.predict(
                ra=ra,
                dec=dec,
                roll=roll,
                camera=camera,
                ccd=ccd,
                sector=sector,
                warp=True,
            )
            wcs.sector = sector
            return wcs

    @classmethod
    def predict(
        cls,
        ra: float,
        dec: float,
        roll: float,
        camera: int,
        ccd: int,
        sector: int = None,
        epoch: int = None,
        warp: bool = True,
    ):
        """Predict a WCS based on RA, Dec and roll

        Parameters
        ----------
        ra : float
            Spacecraft boresight RA in degrees.
        dec : float
            Spacecraft boresight Dec in degrees.
        roll : float
            Spacecraft roll angle in degrees.
        camera : int
            Camera number (1 to 4).
        ccd : int
            CCD number (1 to 4).
        sector : int, optional
            TESS sector number. Used to determine the epoch and nominal reference pixel
            (CRPIX = [1045, 1001] for sectors <= 90; [1069, 1025] for sectors >= 91).
        epoch : int, optional
            Explicit epoch index (1 or 2). Overrides sector if provided. Defaults to 1
            if neither sector nor epoch is specified.
        warp : bool
            Whether to warp by the best fit TESS Camera/CCD warp. This is recommended,
            otherwise solution will be off by up to 10 pixels. This must be switched off
            if fitting for the warp.
        """
        if epoch is None:
            epoch = get_epoch(sector) if sector is not None else 1

        crpix = CRPIX_EPOCHS[epoch]

        if (ccd == 1) | (ccd == 3):
            C = np.asarray([[1, -1], [-1, 1]])
        elif (ccd == 2) | (ccd == 4):
            C = np.asarray([[1, 1], [1, 1]])

        c_boresight = SkyCoord(ra, dec, unit="deg")
        R = angle_to_matrix(-roll)

        # initialize object
        wcs = cls(naxis=2)
        wcs.ra, wcs.dec, wcs.roll, wcs.camera, wcs.ccd = ra, dec, roll, camera, ccd
        wcs.pixel_shape = (rcolumns, rrows)

        # center of the CCD
        xc, yc = cls.xcent[epoch][camera][ccd], cls.ycent[epoch][camera][ccd]
        r, phi = np.hypot(xc, yc) * u.deg, np.arctan2(yc, xc) * u.rad
        c_ccd = c_boresight.directional_offset_by(
            position_angle=phi - roll * u.deg, separation=r
        )

        wcs.wcs.ctype = ["RA---TAN-SIP", "DEC--TAN-SIP"]
        wcs.wcs.cunit = ["deg", "deg"]
        wcs.wcs.radesys = "ICRS"
        wcs.wcs.crpix = np.asarray(crpix, dtype=float)
        wcs.wcs.cdelt = [-1, 1]
        # center is the ccd center in world space
        wcs.wcs.crval = [c_ccd.ra.deg, c_ccd.dec.deg]
        wcs.sip = Sip(
            *[
                cls.sip_dict[attr][epoch][camera][ccd]
                for attr in ["a", "b", "ap", "bp"]
            ],
            crpix,
        )
        # rotate the corners of the CCD to the correct roll
        x, y = R.dot(np.vstack([cls.xs[camera][ccd], cls.ys[camera][ccd]]))
        r, phi = np.hypot(x, y) * u.deg, np.arctan2(y, x) * u.rad

        # positions of the corners in world space
        c = c_boresight.directional_offset_by(position_angle=phi, separation=r)
        r, phi = c_ccd.separation(c), c_ccd.position_angle(c)

        # position of the corners in detector space
        truth = np.asarray([r * np.cos(phi), r * np.sin(phi)])

        # Our estimate of where the corners should in world space
        c = wcs.pixel_to_world(*pixel_corners.T)
        r, phi = c_ccd.separation(c), c_ccd.position_angle(c)
        # Our estimate of where the corners should in detector space
        approx = np.asarray([r * np.cos(phi), r * np.sin(phi)])

        # fit a transformation matrix
        matrix = get_M(truth.T, approx.T)[:2, :2]
        wcs.wcs.pc = matrix.dot(angle_to_matrix(180)).dot(
            np.asarray([[0.0593726, 0.00123252], [0.00123252, 0.0593726]]) * (C)
        )
        if warp:
            Ms = cls.Ms.get(epoch, cls.Ms.get(1))
            offset_weights = cls.offset_weights.get(epoch, cls.offset_weights.get(1))
            wcs.wcs.pc = wcs.wcs.pc.dot(Ms[camera][ccd])
            ecl_roll = wcs.ecl_roll.rad
            X = np.asarray([np.sin(ecl_roll), np.cos(ecl_roll), 1])
            wcs.wcs.crpix += X.dot(offset_weights[camera][ccd])
        return wcs

    def to_header(self, key=None):
        hdr = super().to_header(relax=True, key=key)
        self._update_header(hdr)
        return hdr

    def to_fits(self, **kwargs):
        hdu = fits.PrimaryHDU(header=self.to_header())
        hdulist = fits.HDUList([hdu])
        return hdulist

    def _update_header(self, hdr):
        hdr["AUTHOR"] = "tesswcs"
        for attr in ["ra", "dec", "roll", "sector", "camera", "ccd"]:
            if hasattr(self, attr):
                hdr[attr] = getattr(self, attr)
        return hdr

    @property
    def ecl_roll(self):
        """returns the roll angle of the pointing with respect to the ecliptic plane.
        This is handy for figuring out which are northern hemisphere pointings, and helping with warping.
        """
        c = SkyCoord(self.ra, self.dec, unit="deg")
        ecl = c.transform_to("geocentricmeanecliptic")

        # The roll angle w.r.t the ecliptic plane
        ecl_roll = ecl.position_angle(
            c.directional_offset_by(-self.roll * u.deg, 1 * u.arcsecond).transform_to(
                "geocentricmeanecliptic"
            )
        )
        return ecl_roll
