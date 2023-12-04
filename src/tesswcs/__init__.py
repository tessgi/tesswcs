# Standard library
import os  # noqa

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

version = "1.1.2"

# Standard library
import logging  # noqa: E402

# This library lets us have log messages with syntax highlighting
from rich.logging import RichHandler  # noqa: E402

log = logging.getLogger("tesswcs")
log.addHandler(RichHandler(markup=True))

import numpy as np  # noqa: E402
from astropy.table import Table  # noqa: E402

pointings = Table.read(f"{PACKAGEDIR}/data/pointings.csv", encoding="utf-8-sig")

# Real rows and columns in a CCD
rrows, rcolumns = (2078, 2136)

# pixel size in mm
pixel_size = 0.015

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

from .utils import _load_support_dicts  # noqa: E402, F401
from .utils import _load_warp_matrices  # noqa: E402, F401
from .utils import _load_wcs_data  # noqa: E402, F401

wcs_dicts = _load_wcs_data()

(
    xs,
    ys,
    xcent,
    ycent,
    M,
    sip_dict,
) = _load_support_dicts()

Ms, offset_weights = _load_warp_matrices()

from .tesswcs import WCS  # noqa: E402, F401
