# Recovery Project

Analysis of large scale seasonal co2 flux dynamics and it's drivers in the northern hemisphere

### Notes

Consider doing reduction in Earth Engine first, then downloading results at lower resolution.
This strategy did not work when trying to reproject MODIS 500m resolution data. Earth Engine returns errors. Data size too big.

#### Geometries 

Geometry objects are used to mask or filter the data, e.g. continents, above 50N, ecoregions, test region.

#### Loading raster data

Check units. Could be that they are per pixel and not per area (Carboscope inversion seems to be so).

#### Applying same crs to all data

Using the common WGS84: EPSG 4326. This should be set to all rasters used. If the crs was WGS84 but the property was not set, then use:  
`mydata.rio.write_crs("epsg:4326", inplace=True)`  
Else, change the crs with   
`mydata.rio.reproject("EPSG:4326")`  
When adding more datasets, these can be adjusted to the first using:  
`mydata2 = mydata2.rio.reproject_match(mydata)`  

#### Check if missing data value is set
`mydata.rio.nodata` or `mydata.rio.encoded_nodata` will show the fill value if it is set
`mydata.rio.set_nodata(-9999, inplace=True)` # will set the nadata attrribute without modifying the data
`mydata.rio.write_nodata(-9999, inplace=True)` # will write to the array (I guess replacing the existing missing data value?) Need to test.

Note that the reproject_match method from above will modify the nodata value of mydata2 to match that of mydata.  

Use the following to mask the missing data:  
```
nodata = raster.rio.nodata
raster = raster.where(raster != nodata)
raster.rio.write_nodata(nodata, encoded=True, inplace=True)
```