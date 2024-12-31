# How to show a TESS Observability Map

`tesswcs` is a tool to provide you with the World Coordinate System (WCS) for all past TESS observations, and predict the WCS of future observations. This notebook shows how to use `tesswcs` to plot how much of the sky has been observed by TESS by the end of Cycle 7 observations. 


```python
import tesswcs
from tesswcs.locate import get_observability_mask
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from tqdm import tqdm
```

Make sure you're using the most recent version, (1.2.0 or higher)


```python
tesswcs.__version__
```




    '1.2.0'




```python
# All the TESS pointings
pointings = tesswcs.pointings[['RA', "Dec", "Roll"]].to_pandas().values
```

To plot observability we will need a grid of RAs and Decs over the sky to check. I've just done a uniform grid in RA and Dec, which is not ideal, you could switch this to a different method.


```python
# Grid of RA and Dec to check
# 2000 points in RA
# 1201 points in Dec

# If you increase the resolution this will take longer to calculate.
RA, Dec = np.mgrid[:360:2000j, -90:90:1201j]
```

Below we loop through all the pointings, cameras, and CCDs and calculate which points in the RA and Dec grid fall on a camera. We are going to show only data that will be observed up to and including Cycle 7. This can take a few minutes depending on the resolution of the grid. 


```python
# Array to accumulate number of observations
nobs = np.zeros(RA.shape, dtype=float)
k = tesswcs.pointings['Cycle'] <= 7
# Loop through all the ra, dec and roll of the pointings
for ra, dec, roll in tqdm(pointings[k], desc='Pointing', leave=True, position=0):
    # Loop through each camera
    for camera in np.arange(1, 5):
        # Loop through each CCD
        for ccd in np.arange(1, 5):
            wcs = tesswcs.WCS.predict(ra, dec, roll, camera, ccd)    
            mask = get_observability_mask(wcs, SkyCoord(RA, Dec, unit='deg')).astype(int)
            nobs += mask
```

    Pointing: 100%|█████████████████████████████████| 96/96 [04:22<00:00,  2.73s/it]


Now we have finished, `nobs` is the number of times TESS is able to observe that point in the sky. We can plot this map;


```python
fig, ax = plt.subplots(dpi=250)
cmap = plt.get_cmap('viridis')
cmap.set_extremes(under='r')
im = ax.pcolormesh(RA, Dec, nobs, cmap=cmap, vmin=1, vmax=30, shading='nearest')
cbar = plt.colorbar(im, ax=ax)
ax.set(xlabel='RA [deg]', ylabel='Dec [dec]', title=f"TESS Sectors 1-{tesswcs.pointings['Sector'].max()}")
cbar.set_label("Number of Observations")
```


    
![png](tutorial3_files/tutorial3_10_0.png)
    


If we sum up the number of points that are not zero, we get the sky coverage.


```python
print(100*(nobs!=0).sum()/np.prod(nobs.shape), "% of the sky observed by the end of Cycle 7")
```

    98.64458784346378 % of the sky observed by the end of Cycle 7

