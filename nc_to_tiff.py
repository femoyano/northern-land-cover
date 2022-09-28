# packages 
import xarray as xr 
import rioxarray as rio

# Open data with Gdal and get the CRS information
file_path = '/Users/moyanofe/BigData/GeoSpatial/Carboscope/Inversions/s85oc_v2022_daily.nc' 
# Open data with xarray
nc = xr.open_dataset(file_path)
print(nc)

co2 = nc['co2flux_land'] # extract the co2flux_land variable
co2 = co2.reindex(lat=co2.lat[::-1])
co2.rio.set_spatial_dims(x_dim='lon', y_dim='lat', inplace=True)

print(co2.rio.crs) # check if there is already a crs defined
co2.rio.write_crs("epsg:4326", inplace=True) # Add a crs if not. Global data mostly in WGS84=epsg4326
print(co2.rio.crs)
co2.rio.to_raster(r"/Users/moyanofe/BigData/GeoSpatial/Carboscope/Inversions/s85oc_v2022_daily.tiff")
