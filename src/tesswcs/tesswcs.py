"""Work with WCS objects for TESS"""

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS as astropyWCS
from astropy.wcs import Sip

from . import log, pixel_corners, pointings, rcolumns, rrows  # noqa: E402
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
            wcs = cls(naxis=2)
            wcs.sector, wcs.camera, wcs.ccd = sector, camera, ccd
            wcs.ra = cls.wcs_dicts[sector]["ra"]
            wcs.dec = cls.wcs_dicts[sector]["dec"]
            wcs.roll = cls.wcs_dicts[sector]["roll"]
            # This looks wrong because wcs is column major
            wcs.pixel_shape = (rcolumns, rrows)
            wcs.wcs.ctype = ["RA---TAN-SIP", "DEC--TAN-SIP"]
            wcs.wcs.cunit = ["deg", "deg"]
            wcs.wcs.radesys = "ICRS"
            wcs.wcs.crpix = cls.wcs_dicts[sector][camera][ccd]["crpix0"]
            wcs.wcs.crval = cls.wcs_dicts[sector][camera][ccd]["crval0"]
            wcs.wcs.cdelt = [1, 1]
            wcs.wcs.pc = cls.wcs_dicts[sector][camera][ccd]["cd"]
            wcs.sip = Sip(
                *[cls.sip_dict[attr][camera][ccd] for attr in ["a", "b", "ap", "bp"]],
                [1045, 1001],
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
                ra=ra, dec=dec, roll=roll, camera=camera, ccd=ccd, warp=True
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
        warp: bool = True,
    ):
        """Predict a WCS based on RA, Dec and roll


        warp : bool
            Whether to warp by the best fit TESS Camera/CCD warp. This is recommended,
            otherwise solution will be off by up to 10 pixels. This must be switched off
            if fitting for the warp.
        """
        if (ccd == 1) | (ccd == 3):
            C = np.asarray([[1, -1], [-1, 1]])
        elif (ccd == 2) | (ccd == 4):
            C = np.asarray([[1, 1], [1, 1]])

        c_boresight = SkyCoord(ra, dec, unit="deg")
        R = angle_to_matrix(-roll)

        # initialize object
        wcs = cls(naxis=2)
        wcs.ra, wcs.dec, wcs.roll, wcs.camera, wcs.ccd = ra, dec, roll, camera, ccd
        wcs.pixel_shape = (rrows, rcolumns)

        # center of the CCD
        xc, yc = cls.xcent[camera][ccd], cls.ycent[camera][ccd]
        r, phi = np.hypot(xc, yc) * u.deg, np.arctan2(yc, xc) * u.rad
        c_ccd = c_boresight.directional_offset_by(
            position_angle=phi - roll * u.deg, separation=r
        )

        wcs.wcs.ctype = ["RA---TAN-SIP", "DEC--TAN-SIP"]
        wcs.wcs.cunit = ["deg", "deg"]
        wcs.wcs.radesys = "ICRS"
        #        wcs.wcs.crpix = np.asarray([rrows / 2, rcolumns / 2])
        wcs.wcs.crpix = np.asarray([1045, 1001])
        wcs.wcs.cdelt = [-1, 1]
        # center is the ccd center in world space
        wcs.wcs.crval = [c_ccd.ra.deg, c_ccd.dec.deg]
        wcs.sip = Sip(
            *[cls.sip_dict[attr][camera][ccd] for attr in ["a", "b", "ap", "bp"]],
            [1045, 1001],
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
            wcs.wcs.pc = wcs.wcs.pc.dot(cls.Ms[camera][ccd])
            ecl_roll = wcs.ecl_roll.rad
            X = np.asarray([np.sin(ecl_roll), np.cos(ecl_roll), 1])
            wcs.wcs.crpix += X.dot(cls.offset_weights[camera][ccd])
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
