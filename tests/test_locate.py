import numpy as np
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
