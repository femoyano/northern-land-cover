# %% [markdown]
# ### Creating new time series rasters for each land cover
# 
# **Steps**
# - Define arrays with LC indeces and LC names according to the MODIS LC type used. (see the data documentation)
# - Define paths to files, create sorted list of file names
# - Create array of years covered by the data.
# - Loop over LC types. For each:
#     - Loop over files. For each:
#         - Create the file path
#         - Read the data for the LC type as an xarray
#         - Add a 'time' dimension
#         - Add to already read file data by concatenating on the time dimension
#     - Rename array and attributes as necessary and set a time coordinate with years values
#     - Create an xarray array using np.array and lat, lon, time dims
#     - Add array to xr dataset and save to file

# %%
# Import packages
import os
import xarray as xr
import rioxarray as rio
import numpy as np
import rasterio
import numpy as np

# %% [markdown]
# ### General Settings

# %%
# Default data properties
nodata = -9999.0
epsg_crs = 4326

# Chose the inversion data version. Later (shorter) time series are driven by data from more stations.
inversion = 's99'  # Inversion starting year. Options: s85, s99, s06, s10

input_dir = '../data_input/'
output_dir = '../data_output/'

lc_path_in = '/Users/moyanofe/BigData/GeoSpatial/LandCover/LandCover_MODIS_MCD12/MCD12C1'
lc_path_out = '/Users/moyanofe/BigData/GeoSpatial/LandCover/LandCover_MODIS_MCD12/MCD12C1_proc'

# Index for MCD12C1 Land_Cover_Type_1_Percent: IGBP land cover types
lcIndex = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
lcNames = ['Water', 'ENForest', 'EBForest', 'DNForest', 'DBForest', 
    'MixForest', 'ClosedShrub', 'OpenShrub', 'WoodySavanna',
    'Savanna', 'Grassland', 'PermWetland', 'Cropland',
    'Urban', 'CropNatMosiac', 'PermSnowIce', 'Barren']

# %% [markdown]
# ### Processing LC data

# %%
# Extract individual lc data from files, concatenate into single xarrays and save.

files_in = os.listdir(lc_path_in)
files_in.sort() # Sort to order files by year (which is part of the file name)
files_in = files_in[0:2] # shorten for testing only
years = np.arange(stop=len(files_in), step=1) + 2001

for i in range(len(lcIndex)): # [0]: #
    print(lcNames[i])
    # Take each LC from the array and create a new array for each by combining all yearly file
    for j in range(len(files_in)): # 
        print(j)
        file = files_in[j]
        print(file)
        path = os.path.join(lc_path_in, file)
        ds_file = rio.open_rasterio(path, masked=True, variable='Land_Cover_Type_1_Percent').squeeze()
        ar = ds_file.Land_Cover_Type_1_Percent[i]
        ar = ar.expand_dims('time')
        # print(ar)
        if j == 0:
            stacked = ar
        else:
            stacked = xr.concat([stacked, ar], dim='time')
    stacked = stacked.rename(lcNames[i])
    stacked.attrs['long_name'] = lcNames[i]
    stacked.attrs['units'] = 'percent in integers'
    stacked['time'] = ('time', years)
    ds_modislc_igbp = xr.Dataset()
    ds_modislc_igbp[lcNames[i]] = stacked

    # Save each land caover time series to file
    file_out = 'MCD12C1.A2001-2021.061.LCtype1.' + lcNames[i] + '.nc'
    filepath_out = os.path.join(lc_path_out, file_out)
    ds_modislc_igbp.to_netcdf(filepath_out)


# %%
# Save at lower resolution: same as nee data

# Load the nee data for reproject_match ----
file_neeAmpStats = os.path.join(output_dir, 'neeAmpStats_2010-2021_Inv2010.nc')
neeamp = rio.open_rasterio(file_neeAmpStats)

# Path to files
ds = xr.Dataset()

# Loop over files, reproject and save
for i in range(len(lcIndex)): # [0]: #
    print(lcNames[i])
    file_in = 'MCD12C1.A2001-2021.061.LCtype1.' + lcNames[i] + '.nc'
    filepath_in = os.path.join(lc_path_out, file_in)
    ds_in = rio.open_rasterio(filepath_in, masked=True, variable=lcNames[i])
    da = ds_in[lcNames[i]]
    da.rio.write_crs("epsg:4326", inplace=True)
    da = da.rio.reproject_match(neeamp, resampling=rasterio.enums.Resampling.average)
    ds[lcNames[i]] = da

file_out = 'MCD12C1.A2001-2021.061.LCtype1.All.lr.nc'
filepath_out = os.path.join(lc_path_out, file_out)
ds.to_netcdf(filepath_out)


