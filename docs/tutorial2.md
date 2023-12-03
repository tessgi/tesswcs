# Using `tesswcs` to find observable targets

`tesswcs` has some convenience functions to help you find out which targets are observable with TESS. The `get_pixel_locations` functon will accept sky coordinates, and return a table of which sector, camera, CCD, row, and column those targets fall on (if any). 

Let's take a look at how to use it. 

## Example 1: Finding where transiting planets fall on TESS pixels


```python
from tesswcs.locate import get_pixel_locations
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.table import Table
```

First, we need a target list to check the observability of. For this purpose, I will use the list of confirmed exoplanets. Below I use astroquery to get the coordinates of all the confirmed, transiting exoplanets from the [Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/index.html).


```python
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
confirmed_exoplanets_table = NasaExoplanetArchive.query_criteria(table="ps", select="pl_name,ra,dec,sy_tmag", where="default_flag=1 AND tran_flag=1")
confirmed_exoplanets_table.sort('sy_tmag')
```


```python
confirmed_exoplanets_table
```


```python
confirmed_exoplanets_table['sky_coord']
```

I can now pass these `SkyCoord`s to `get_pixel_location` to get a table of the regions those targets fall on. Any targets that do not fall on a pixel will be omitted.


```python
observable = get_pixel_locations(confirmed_exoplanets_table['sky_coord'])
```

We can use pandas to learn more about this table. For example, how many planets are observable for the first time during cycle 7?


```python
g = observable.to_pandas().groupby('Target Index').min()
g
```

This groupby operation takes the `min` across all the columns for cases of the same `'Target Index'`


```python
cycle7_mask = (g['Sector'] >= 84) & (g['Sector'] <= 96)
target_index = list(g[cycle7_mask].index)
print(len(target_index))
```

There are 92 confirmed transiting planets that will be observed for the first time in Cycle 7!


```python
confirmed_exoplanets_table[target_index]
```

## Example 2: Finding whether a transient will be observable with TESS

Sometimes we might find an interesting transient, and we want to know whether it will fall on a TESS pixel so that we can recover the optical time-series when the data is downlinked. For example, [GRB 231106A](https://gcn.nasa.gov/circulars/34956) occured in November 2023 and happened to be observable during the TESS survey. 

Here we show how we would find out if a transient has been (or will be) observable with TESS. 

To do this we will need two things, the coordinate of the object, and the time at which we would like to observe it. Here we'll use the location and detection time of GRB 231106A.


```python
from astropy.time import Time
from astropy.coordinates import SkyCoord

t = Time('2023-11-06T18:13:23Z')
c = SkyCoord(113.4482, 29.2245, unit='degree')
```


```python
get_pixel_locations(coords=c, time=t)
```

We can see that this returned a hit, we are able to detect the target during Sector 71, on Camera 4, CCD 1! 

In instances where the time passed is recent (or in the future!) and there is no archival TESS WCS to use, `tesswcs` will predict the WCS for future TESS pointings.
