<a href="https://github.com/tessgi/tesswcs/actions/workflows/tests.yml"><img src="https://github.com/tessgi/tesswcs/workflows/pytest/badge.svg" alt="Test status"/></a> [![Generic badge](https://img.shields.io/badge/documentation-live-blue.svg)](https://tessgi.github.io/tesswcs/)
[![PyPI version](https://badge.fury.io/py/tesswcs.svg)](https://badge.fury.io/py/tesswcs)

<p align="center">
  <img src="https://github.com/tessgi/tesswcs/blob/main/docs/images/logo.png?raw=true" width="350" alt="tesswcs logo">
</p>

# tesswcs

This package will enable you to create an [`astropy` World Coordinate System](https://docs.astropy.org/en/stable/wcs/) for any pointing of the TESS telescope. You can access both the true WCS from archival data, and predict the WCS for a given RA, Dec, and spacecraft roll.

## Installation

You can install `tesswcs` using `pip`

```
pip install --upgrade tesswcs
```

## Usage

Below is an example of how to obtain a WCS for a given sector

```python
import tesswcs

wcs = tesswcs.WCS.from_sector(sector=1, camera=1, ccd=1)
wcs
```

```
WCS Keywords

Number of WCS axes: 2
CTYPE : 'RA---TAN-SIP' 'DEC--TAN-SIP' 
CRVAL : 13.594227829971574 -20.630003889665105 
CRPIX : 1046.0046042531758 1002.004711738613 
PC1_1 PC1_2  : -0.0052945863349177745 0.0022643466111789898 
PC2_1 PC2_2  : 0.002130969047581938 0.005238488386839341 
CDELT : -1.0 1.0 
NAXIS : 2078  2136
```

This returns an `astropy` World Coordinate System object so you can use all the features of `astropy`. For example, using `astropy`'s WCS interface you can now either work with this object, for example you can obtain the sky position of row and column positions

```python
wcs.pixel_to_world(row, column)
```

or obtain row and column positions from sky positions (using the `astropy.coordinates.SkyCoord` object).

```python
from astropy.coordinates import SkyCoord
wcs.world_to_pixel(SkyCoord.from_name("HD 209458"))
```

You can also save these objects to fits files

```python
wcs.to_fits('wcs.fits')
```

If you have a RA, Dec and roll in degrees you can also predict a WCS.

```python
import tesswcs

wcs = tesswcs.WCS.predict(ra=0, dec=0, roll=0, camera=1, ccd=1)
```

You can use tesswcs to better understand what sources will be obervable on TESS pixels, or to make figures like the ones below! [Check out the tutorials](https://tessgi.github.io/tesswcs/tutorial3/) for more information on how to produce these figures.

![Figure showing the predicted TESS WCS](figures/tess_1_predict.png)

If you have a list of targets and you want to see if they are observable you can use

```python
from tesswcs.locate import check_observability

check_observability(coords, sector=sector)
```

Which will check if the input `astropy.coordinates.SkyCoord` object is observable. This will be faster if you pass in a sector to check against .

If you are interested in what interesting targets are available in a given TESS sector you can use

```python
from tesswcs.locate import get_interesting_targets

get_interesting_targets(sector=sector)
```

which will return a `pd.DataFrame` object containing astronomical targets with high numbers of paper references.

## What's new?

In version 1.2 and higher `tesswcs` now includes the expected pointing parameters for TESS EM3 (sectors 97-134). These are **expected** pointings and are subject to change. As TESS takes data `tesswcs` will be updated to replace the predicted WCS for future sectors with the measured WCS when data is archived.

## Limits on accuracy

TESS observes large parts of the sky with a rapid cadence during ~27 day sectors. Since the spacecraft orbits the Earth, which in turn orbits the Sun, the spacecraft moving at high velocity. Because of its large field of view, this imparts a differential velocity aberration (i.e. targets at different parts of the detector have a different apparent motion.) `tesswcs` provides a static World Coordinate System solution at the mean time of the observation, which can not account for this differential velocity aberration. This can cause small (usually subpixel) apparent inaccuracies of the pointing over time. As such, `tesswcs` alone can not be used for astrometric measurements with TESS.

## Updating `tesswcs`

To update `tesswcs` with, for example, new pointing information, follow these steps:

1. `git clone` this repository. Ensure that all your files are up to date with a `git pull` if you already have the repository.
2. Update the required information (for example, update the pointings table in `src/tesswcs/data/pointings.csv')
3. Locally, run the tests using the make file. You can run `make` or `make black`, `make isort`, `make pytest`. Check that everything passes.
4. Update the version number. The **major** version number should be saved for significant API changes. The **patch** number is saved for automatic updates of the database when new TESS data is available. You should update the **minor** version number if you update the pointings table.
5. Open a pull request against this repository, ensure tests pass online, and merge.
6. Release a new version on github
7. On your local machine, release a new version to pypi using `poetry build` and then `poetry publish`.
8. Check the documentation using `make serve`. Update the documentation by running `make deploy`.

## Changelog

- v1.5.0 removes a bug identified by @ben-cassese and @jgagneastro where predicted sectors were mistakenly hard coded to only return Sector 3. Updating minor version number to reflect severity of this bug!
- v.1.6.1 removes a bug identified by @altuson where the shape of the detector was row major not column major
