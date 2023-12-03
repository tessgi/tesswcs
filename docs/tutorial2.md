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




<div><i>QTable masked=True length=4176</i>
<table id="table4546443248" class="table-striped table-bordered table-condensed">
<thead><tr><th>pl_name</th><th>ra</th><th>dec</th><th>sy_tmag</th><th>sky_coord</th></tr></thead>
<thead><tr><th></th><th>deg</th><th>deg</th><th>mag</th><th>deg,deg</th></tr></thead>
<thead><tr><th>str19</th><th>float64</th><th>float64</th><th>float64</th><th>SkyCoord</th></tr></thead>
<tr><td>HD 219134 c</td><td>348.3372026</td><td>57.1696255</td><td>4.6278</td><td>348.3372026,57.1696255</td></tr>
<tr><td>HD 219134 b</td><td>348.3372026</td><td>57.1696255</td><td>4.6278</td><td>348.3372026,57.1696255</td></tr>
<tr><td>HD 136352 b</td><td>230.4401147</td><td>-48.3188174</td><td>5.0494</td><td>230.4401147,-48.3188174</td></tr>
<tr><td>HD 136352 c</td><td>230.4401147</td><td>-48.3188174</td><td>5.0494</td><td>230.4401147,-48.3188174</td></tr>
<tr><td>HD 136352 d</td><td>230.4401147</td><td>-48.3188174</td><td>5.0494</td><td>230.4401147,-48.3188174</td></tr>
<tr><td>pi Men c</td><td>84.2992798</td><td>-80.4646041</td><td>5.1054</td><td>84.2992798,-80.4646041</td></tr>
<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
<tr><td>WD 1856+534 b</td><td>284.415675</td><td>53.5090244</td><td>16.338</td><td>284.415675,53.5090244</td></tr>
<tr><td>Lupus-TR-3 b</td><td>232.5774696</td><td>-42.9799461</td><td>17.7596</td><td>232.5774696,-42.9799461</td></tr>
<tr><td>HAT-P-63 b</td><td>279.7967696</td><td>-20.0242684</td><td>19.2571</td><td>279.7967696,-20.0242684</td></tr>
<tr><td>OGLE-TR-56 b</td><td>269.1479583</td><td>-29.5392222</td><td>nan</td><td>269.1479583,-29.5392222</td></tr>
<tr><td>SWEEPS-11 b</td><td>269.7583333</td><td>-29.1981944</td><td>nan</td><td>269.7583333,-29.1981944</td></tr>
<tr><td>SWEEPS-4 b</td><td>269.7246667</td><td>-29.1890556</td><td>nan</td><td>269.7246667,-29.1890556</td></tr>
</table></div>




```python
confirmed_exoplanets_table['sky_coord']
```




    <SkyCoord (ICRS): (ra, dec) in deg
        [(348.3372026,  57.1696255), (348.3372026,  57.1696255),
         (230.4401147, -48.3188174), ..., (269.1479583, -29.5392222),
         (269.7583333, -29.1981944), (269.7246667, -29.1890556)]>



I can now pass these `SkyCoord`s to `get_pixel_location` to get a table of the regions those targets fall on. Any targets that do not fall on a pixel will be omitted.


```python
observable = get_pixel_locations(confirmed_exoplanets_table['sky_coord'])
```

We can use pandas to learn more about this table. For example, how many planets are observable for the first time during cycle 7?


```python
g = observable.to_pandas().groupby('Target Index').min()
g
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Sector</th>
      <th>Camera</th>
      <th>CCD</th>
      <th>Row</th>
      <th>Column</th>
    </tr>
    <tr>
      <th>Target Index</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>17</td>
      <td>2</td>
      <td>1</td>
      <td>384.5</td>
      <td>105.2</td>
    </tr>
    <tr>
      <th>1</th>
      <td>17</td>
      <td>2</td>
      <td>1</td>
      <td>384.5</td>
      <td>105.2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>12</td>
      <td>1</td>
      <td>1</td>
      <td>204.0</td>
      <td>107.7</td>
    </tr>
    <tr>
      <th>3</th>
      <td>12</td>
      <td>1</td>
      <td>1</td>
      <td>204.0</td>
      <td>107.7</td>
    </tr>
    <tr>
      <th>4</th>
      <td>12</td>
      <td>1</td>
      <td>1</td>
      <td>204.0</td>
      <td>107.7</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>4171</th>
      <td>12</td>
      <td>1</td>
      <td>1</td>
      <td>147.9</td>
      <td>1114.4</td>
    </tr>
    <tr>
      <th>4172</th>
      <td>92</td>
      <td>1</td>
      <td>2</td>
      <td>418.6</td>
      <td>1154.5</td>
    </tr>
    <tr>
      <th>4173</th>
      <td>39</td>
      <td>1</td>
      <td>1</td>
      <td>356.6</td>
      <td>27.1</td>
    </tr>
    <tr>
      <th>4174</th>
      <td>65</td>
      <td>1</td>
      <td>1</td>
      <td>1057.1</td>
      <td>1334.6</td>
    </tr>
    <tr>
      <th>4175</th>
      <td>65</td>
      <td>1</td>
      <td>1</td>
      <td>1055.5</td>
      <td>1339.7</td>
    </tr>
  </tbody>
</table>
<p>4169 rows × 5 columns</p>
</div>



This groupby operation takes the `min` across all the columns for cases of the same `'Target Index'`


```python
cycle7_mask = (g['Sector'] >= 84) & (g['Sector'] <= 96)
target_index = list(g[cycle7_mask].index)
print(len(target_index))
```

    92


There are 92 confirmed transiting planets that will be observed for the first time in Cycle 7!


```python
confirmed_exoplanets_table[target_index]
```




<div><i>QTable masked=True length=92</i>
<table id="table5011763008" class="table-striped table-bordered table-condensed">
<thead><tr><th>pl_name</th><th>ra</th><th>dec</th><th>sy_tmag</th><th>sky_coord</th></tr></thead>
<thead><tr><th></th><th>deg</th><th>deg</th><th>mag</th><th>deg,deg</th></tr></thead>
<thead><tr><th>str19</th><th>float64</th><th>float64</th><th>float64</th><th>SkyCoord</th></tr></thead>
<tr><td>HD 106315 c</td><td>183.4724742</td><td>-0.3934357</td><td>8.56338</td><td>183.4724742,-0.3934357</td></tr>
<tr><td>HD 106315 b</td><td>183.4724742</td><td>-0.3934357</td><td>8.56338</td><td>183.4724742,-0.3934357</td></tr>
<tr><td>Wolf 503 b</td><td>206.8461979</td><td>-6.1393369</td><td>9.2493</td><td>206.8461979,-6.1393369</td></tr>
<tr><td>HD 137496 b</td><td>231.7420543</td><td>-16.5090014</td><td>9.29961</td><td>231.7420543,-16.5090014</td></tr>
<tr><td>K2-292 b</td><td>205.3758786</td><td>-9.9460704</td><td>9.31975</td><td>205.3758786,-9.9460704</td></tr>
<tr><td>K2-233 b</td><td>230.4799011</td><td>-20.2317975</td><td>9.8125</td><td>230.4799011,-20.2317975</td></tr>
<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
<tr><td>K2-316 c</td><td>233.6247336</td><td>-23.2594663</td><td>14.2731</td><td>233.6247336,-23.2594663</td></tr>
<tr><td>K2-316 b</td><td>233.6247336</td><td>-23.2594663</td><td>14.2731</td><td>233.6247336,-23.2594663</td></tr>
<tr><td>K2-52 b</td><td>246.6114157</td><td>-24.9699462</td><td>14.2904</td><td>246.6114157,-24.9699462</td></tr>
<tr><td>K2-315 b</td><td>228.0210932</td><td>-20.1081634</td><td>14.3273</td><td>228.0210932,-20.1081634</td></tr>
<tr><td>K2-378 b</td><td>205.1596569</td><td>-11.1172332</td><td>14.8943</td><td>205.1596569,-11.1172332</td></tr>
<tr><td>HAT-P-63 b</td><td>279.7967696</td><td>-20.0242684</td><td>19.2571</td><td>279.7967696,-20.0242684</td></tr>
</table></div>



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




<div><i>Table length=1</i>
<table id="table4429965392" class="table-striped table-bordered table-condensed">
<thead><tr><th>Target Index</th><th>Sector</th><th>Camera</th><th>CCD</th><th>Row</th><th>Column</th></tr></thead>
<thead><tr><th>int64</th><th>int64</th><th>int64</th><th>int64</th><th>float64</th><th>float64</th></tr></thead>
<tr><td>0</td><td>71</td><td>4</td><td>1</td><td>764.6</td><td>1451.4</td></tr>
</table></div>



We can see that this returned a hit, we are able to detect the target during Sector 71, on Camera 4, CCD 1! 

In instances where the time passed is recent (or in the future!) and there is no archival TESS WCS to use, `tesswcs` will predict the WCS for future TESS pointings.
