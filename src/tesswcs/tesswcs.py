"""Work with WCS objects for TESS"""
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS as astropyWCS
from astropy.wcs import Sip
from astropy.io import fits

from . import pixel_corners, rcolumns, rrows
from .utils import _load_support_dicts, _load_wcs_data, angle_to_matrix, get_M

# from . import pointings  # noqa: E402


class WCS(astropyWCS):
    """A special subclass of astropy.wcs.WCS

    This class allows us to add attributes to the WCS header and add class methods
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _load_support_dicts(self):
        (
            self._xs,
            self._ys,
            self._xcent,
            self._ycent,
            self._sip_dict,
        ) = _load_support_dicts()

    def _load_wcs_data(self):
        self._wcs_dicts = _load_wcs_data()

    @property
    def wcs_dicts(self):
        if not hasattr(self, "_wcs_dicts"):
            self._load_wcs_data()
        return self._wcs_dicts

    @property
    def xs(self):
        if not hasattr(self, "_xs"):
            self._load_support_dicts()
        return self._xs

    @property
    def ys(self):
        if not hasattr(self, "_ys"):
            self._load_support_dicts()
        return self._ys

    @property
    def xcent(self):
        if not hasattr(self, "_xcent"):
            self._load_support_dicts()
        return self._xcent

    @property
    def ycent(self):
        if not hasattr(self, "_ycent"):
            self._load_support_dicts()
        return self._ycent

    @property
    def sip_dict(self):
        if not hasattr(self, "_sip"):
            self._load_support_dicts()
        return self._sip_dict

    @classmethod
    def from_archive(cls, sector: int, camera: int, ccd: int):
        """Load a WCS from archival TESS data"""
        wcs = cls(naxis=2)
        wcs.sector, wcs.camera, wcs.ccd = sector, camera, ccd
        wcs.ra = wcs.wcs_dicts[sector]["ra"]
        wcs.dec = wcs.wcs_dicts[sector]["dec"]
        wcs.roll = wcs.wcs_dicts[sector]["roll"]

        wcs.wcs.ctype = ["RA---TAN-SIP", "DEC--TAN-SIP"]
        wcs.wcs.cunit = ["deg", "deg"]
        wcs.wcs.radesys = "ICRS"
        wcs.wcs.crpix = wcs.wcs_dicts[sector][camera][ccd]["crpix0"]
        wcs.wcs.crval = wcs.wcs_dicts[sector][camera][ccd]["crval0"]
        wcs.wcs.cdelt = [1, 1]
        wcs.wcs.pc = wcs.wcs_dicts[sector][camera][ccd]["cd"]
        wcs.sip = Sip(
            *[wcs.sip_dict[attr][camera][ccd] for attr in ["a", "b", "ap", "bp"]],
            [1045, 1001]
        )
        return wcs

    @classmethod
    def predict(cls, ra: float, dec: float, roll: float, camera: int, ccd: int):
        """Predict a WCS based on RA, Dec and roll"""
        if (ccd == 1) | (ccd == 3):
            C = np.asarray([[1, -1], [-1, 1]])
        elif (ccd == 2) | (ccd == 4):
            C = np.asarray([[1, 1], [1, 1]])

        c_boresight = SkyCoord(ra, dec, unit="deg")
        R = angle_to_matrix(-roll)

        # initialize object
        wcs = cls(naxis=2)
        wcs.ra, wcs.dec, wcs.roll, wcs.camera, wcs.ccd = ra, dec, roll, camera, ccd

        # center of the CCD
        xc, yc = wcs.xcent[camera][ccd], wcs.ycent[camera][ccd]
        r, phi = np.hypot(xc, yc) * u.deg, np.arctan2(yc, xc) * u.rad
        c_ccd = c_boresight.directional_offset_by(
            position_angle=phi - roll * u.deg, separation=r
        )

        wcs.wcs.ctype = ["RA---TAN-SIP", "DEC--TAN-SIP"]
        wcs.wcs.cunit = ["deg", "deg"]
        wcs.wcs.radesys = "ICRS"
        wcs.wcs.crpix = np.asarray([rrows / 2, rcolumns / 2])
        wcs.wcs.cdelt = [-1, 1]
        # center is the ccd center in world space
        wcs.wcs.crval = [c_ccd.ra.deg, c_ccd.dec.deg]
        wcs.sip = Sip(
            *[wcs.sip_dict[attr][camera][ccd] for attr in ["a", "b", "ap", "bp"]],
            [1045, 1001]
        )
        # rotate the corners of the CCD to the correct roll
        x, y = R.dot(np.vstack([wcs.xs[camera][ccd], wcs.ys[camera][ccd]]))
        r, phi = np.hypot(x, y) * u.deg, np.arctan2(y, x) * u.rad

        # positions of the corners in world space
        c = c_boresight.directional_offset_by(position_angle=phi, separation=r)
        r, phi = c_ccd.separation(c), c_ccd.position_angle(c)

        # position of the corners in detector space
        truth = np.asarray([r * np.cos(phi), r * np.sin(phi)])

        # Our estimate of where the corners should in world space
        c = wcs.pixel_to_world(
            *(pixel_corners + np.asarray([rrows / 2, rcolumns / 2])[:, None])
        )
        r, phi = c_ccd.separation(c), c_ccd.position_angle(c)
        # Our estimate of where the corners should in detector space
        approx = np.asarray([r * np.cos(phi), r * np.sin(phi)])

        # fit a transformation matrix
        matrix = get_M(truth.T, approx.T)[:2, :2]
        wcs.wcs.pc = matrix.dot(angle_to_matrix(180)).dot(
            np.asarray([[0.0593726, 0.00123252], [0.00123252, 0.0593726]]) * (C)
        )

        if (ccd == 1) | (ccd == 3):
            C = np.asarray([[1, -1], [-1, 1]])
        elif (ccd == 2) | (ccd == 4):
            C = np.asarray([[1, 1], [1, 1]])

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
