# About `tesswcs`

`tesswcs` is available to help make understanding what sources will fall on which pixels in the TESS focal plane. This two works by taking the archival WCS of TESS and using them to create a general WCS solution. `tesswcs` then can apply this general solution to any RA, Dec, and Roll of the boresight of the telescope, to predict the WCS is any given pointing. 

`tesswcs` relies on `tessrip` to access the TESSCut datacubes stored by MAST at AWS. These cubes contain within them the WCS solution for every FFI taken with TESS. `tessrip` downloads and averages the WCS solutions for every FFI in a given sector, camera, and CCD and averages them. `tesswcs` then uses these averages to calculate a general solution which works for any pointing.

