import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time

from tesswcs.locate import check_observability, get_pixel_locations


def test_locate():
    t = Time("2023-11-06T18:13:23Z")
    c = SkyCoord(113.4482, 29.2245, unit="degree")

    observable = check_observability(c, sector=1)
    assert isinstance(observable, Table)
    assert len(observable) == 16
    assert len(observable.keys()) == 4
    observable = check_observability([c], sector=1)
    assert len(observable.keys()) == 4
    observable = check_observability([c, c, c, c], sector=1)
    assert len(observable) == 16
    assert len(observable.keys()) == 7
    observable = check_observability(SkyCoord([c, c, c, c]), sector=1)
    assert len(observable.keys()) == 7
    observable = check_observability(SkyCoord([c, c, c, c]), sector=90)
    assert len(observable.keys()) == 7
    observable = check_observability(c, time=t)
    sector, camera, ccd, _ = np.asarray(observable[observable["targ_0001"]])[0]
    assert sector == 71
    assert camera == 4
    assert ccd == 1

    pixel_locations = get_pixel_locations(c, cycle=6)
    assert len(pixel_locations) == 2
    assert len(pixel_locations.keys()) == 6
    assert pixel_locations["Sector"][0] == 71
    pixel_locations = get_pixel_locations(c, time=t)
    assert len(pixel_locations) == 1


def test_pixel_location():
    # Time and coordinate correspond to asteroid 1998 YT6, as queried on JPL Horizons.
    t = Time(2458489.8075233004, format="jd")
    c = SkyCoord("05 45 56.64 -00 16 15.3", frame="icrs", unit=(u.hourangle, u.deg))

    pixel_locations = get_pixel_locations(c, time=t)

    assert len(pixel_locations) == 1
    assert pixel_locations["Sector"][0] == 6
    assert pixel_locations["Camera"][0] == 1
    assert pixel_locations["CCD"][0] == 1

    # Expected col, row calculated by:
    # 1. using tessrip to get average WCS from sector/camera/ccd
    # 2. converting coord to col, row with wcs_world2pix()
    assert np.round(pixel_locations["Row"][0], 1) == 1107.6
    assert np.round(pixel_locations["Column"][0], 1) == 1087.8


def test_pixel_location_skycoord_array():
    # Time is in Sector 6, targets are TOI-700, L 98-59, and SN 1987 A
    # TOI-700 and SN 1987 A were observed in Sector 6, L 98-59 was not
    t = Time(2458489.8075233004, format="jd")
    coordinates = ["97.09679 -65.57931",
                   "124.53176 -68.313",
                   "83.86658 -69.2696"]
    c = SkyCoord(coordinates, frame="icrs", unit=u.deg)

    pixel_locations = get_pixel_locations(c, time=t)

    assert len(pixel_locations) == 2
    assert pixel_locations["Sector"][0] == 6
    assert pixel_locations["Camera"][0] == 4
    assert pixel_locations["CCD"][0] == 1
    assert pixel_locations["Target Index"][0] == 0

    assert pixel_locations["Sector"][1] == 6
    assert pixel_locations["Camera"][1] == 4
    assert pixel_locations["CCD"][1] == 3
    assert pixel_locations["Target Index"][1] == 2


