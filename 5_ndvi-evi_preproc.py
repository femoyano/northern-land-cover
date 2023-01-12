# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     custom_cell_magics: kql
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: landcover
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## Process NDVI or EVI

# %%
# Magic commands
# %load_ext autoreload
# %autoreload 2

# Packages
# Import packages
import os
import xarray as xr
import rioxarray as rio
import numpy as np
import rasterio
import numpy as np
import pandas as pd

# %%
# Uncomment lines to choose processing NDVI or EVI

# var_raw = 'CMG 0.05 Deg Monthly NDVI'
# var_name = 'ndvi'
# var_longname = "0.05Deg Monthly NDVI"

var_raw = 'CMG 0.05 Deg Monthly EVI'
var_name = 'evi'
var_longname = "0.05Deg Monthly EVI"

# %%
mod13c2_path_in = '/Users/moyanofe/BigData/GeoSpatial/MODIS/MOD13C2/v061_0.05deg'
mod13c2_path_out = '/Users/moyanofe/BigData/GeoSpatial/MODIS/MOD13C2/v061_0.05deg_'+var_name

files_all = os.listdir(mod13c2_path_in)

years = np.arange(stop=23, step=1) + 2000

for year in years:
    files_in = [k for k in files_all if 'MOD13C2.A'+str(year) in k]
    files_in.sort() # Sort to order files by year (which is part of the file name)

    if(year==2000):
        dates = pd.date_range("2000-02", periods=len(files_in), freq='M')
    else:
        dates = pd.date_range(str(year)+"-01", periods=len(files_in), freq='M')

    for j in range(len(files_in)):
        print(j)
        file = files_in[j]
        print(file)
        path = os.path.join(mod13c2_path_in, file)
        xrds = rio.open_rasterio(path, masked=True, variable=var_raw).squeeze()
        ar = xrds[var_raw]
        ar = ar.expand_dims('time')
        # print(ar)
        if j == 0:
            stacked = ar
        else:
            stacked = xr.concat([stacked, ar], dim='time')
    stacked = stacked.rename(var_name)
    stacked.attrs['long_name'] = var_longname
    stacked.attrs['units'] = ''
    stacked['time'] = ('time', dates)
    ds_mod13c2 = xr.Dataset()
    ds_mod13c2[var_name] = stacked

    # Save to file
    file_out = 'MOD13C2.A_' + var_name + str(year) + '.061.nc'
    filepath_out = os.path.join(mod13c2_path_out, file_out)
    ds_mod13c2.to_netcdf(filepath_out)

# %%
# Read yearly files, calculate yearly means, save to a single file

ndvi_path_in = '/Users/moyanofe/BigData/GeoSpatial/MODIS/MOD13C2/v061_0.05deg_'+var_name
ndvi_path_out = '/Users/moyanofe/BigData/GeoSpatial/MODIS/MOD13C2'

files_all = os.listdir(ndvi_path_in)
years = np.arange(stop=23, step=1) + 2000

for year in years:
    file = 'MOD13C2.A_' + var_name + str(year) + '.061.nc'
    print(file)
    path = os.path.join(ndvi_path_in, file)
    xrds = rio.open_rasterio(path, masked=True, variable=var_name).squeeze()
    ar = xrds[var_name]
    ar = ar.mean(['time'])
    ar = ar.expand_dims('time')
    # print(ar)
    if year == 2000:
        stacked = ar
    else:
        stacked = xr.concat([stacked, ar], dim='time')

stacked = stacked.rename(var_name)
stacked.attrs['long_name'] = var_longname
stacked.attrs['units'] = ''
stacked['time'] = ('time', years)
ds_ndvi = xr.Dataset()
ds_ndvi[var_name] = stacked

# Save to file
file_out = 'MOD13C2.A_'+ var_name + '2000-2022'+'.061.nc'
filepath_out = os.path.join(ndvi_path_out, file_out)
ds_ndvi.to_netcdf(filepath_out)
