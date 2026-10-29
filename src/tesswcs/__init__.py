# Standard library
import os  # noqa

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

from importlib.metadata import PackageNotFoundError, version  # noqa


def get_version():
    try:
        return version("tesswcs")
    except PackageNotFoundError:
        return "unknown"


__version__ = get_version()


# Standard library
import logging  # noqa: E402

log = logging.getLogger("tesswcs")

import numpy as np  # noqa: E402
from astropy.table import Table  # noqa: E402

pointings = Table.read(f"{PACKAGEDIR}/data/pointings.csv", encoding="utf-8-sig")

# Real rows and columns in a CCD
rrows, rcolumns = (2078, 2136)

# pixel size in mm
pixel_size = 0.015

# Reference pixels per epoch
# Epoch 1: Sectors 1-90
# Epoch 2: Sectors 91+
CRPIX_EPOCHS = {
    1: [1045, 1001],
    2: [1069, 1025],
}
CRPIX = CRPIX_EPOCHS[1]


def get_epoch(sector: int) -> int:
    """Return the epoch index (1 for sectors <= 90, 2 for sectors >= 91)."""
    return 1 if int(sector) <= 90 else 2


def get_crpix(sector: int) -> list[int]:
    """Return nominal reference pixel for a given sector."""
    return CRPIX_EPOCHS[get_epoch(sector)]



def get_epoch(sector: int) -> int:
    """Return the epoch index (1 for sectors <= 90, 2 for sectors >= 91)."""
    return 1 if int(sector) <= 90 else 2


def get_crpix(sector: int) -> list[int]:
    """Return nominal reference pixel for a given sector."""
    return CRPIX_EPOCHS[get_epoch(sector)]

# pixel scale in degrees
pixel_scale = 21 / 3600

# This sets the orientations of the ccds within the cameras
cdelt = {1: [-1, 1], 2: [-1, 1], 3: [1, -1], 4: [1, -1]}

# This sets up the corners of each CCD
pixel_corners = np.array(
    [
        [0, 0],  # Bottom left
        [0, 1],  # Top left
        [1, 0],  # Bottom right
        [1, 1],  # Top right
    ]
) * np.asarray((rrows, rcolumns))

from .tesswcs import WCS  # noqa: E402, F401
from .utils import _load_wcs_database  # noqa: E402, F401

latest_sector = np.max(list(_load_wcs_database().keys()))
