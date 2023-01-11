# -*- coding: utf-8 -*-
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
# ## NEE flux analysis
#
# This script processes the raw inversion data and saves it in standardized netcdf files.
#
# Both CarboScope and CAMS ivnersion data are processed and saved.
#
# Steps include (not necessarily complete or in order):
# - load inversion NEE data
# - adjust units
# - set nodata
# - create or adjust the time variable
# - mask using a continents roi
#
# **Note on inversion datasets:**
#
# Carboscope inversions differ in the amount of stations data used to create them. Versions starting earlier have less stations with the excpetion fo sEXT, which uses a relatiosnhip between NEE adn temperature.
#
# CAMS inversions have two versions: surface and satellite. Surface is similar to Carboscope, i.e. driven by surface CO2 data, but has non zero priors (I thinkn model based). Satellite uses satellite column CO2 data. This data starts in mid-2009, but the file used here starts in 2010-01 in order to have whole years.
#
# Currently (2023-01) I'm using CAMS monthly inversion data, as downloaded from the net. The daily Carboscope data is then averaged from daily to monthly to have equal frequency. But there is also the option to download 'instantaneous' CAMS data, which is 3-hourly and could be averaged to daily, and then analyze everything with daily frequency.
# %%
# Magic commands to updated functions when python module is modified
# %load_ext autoreload
# %autoreload 2

# Get all required packages

# Packages
import os
import rioxarray as rio
import xarray as xr
import geopandas as gpd
import pandas as pd
import numpy as np

# %% [markdown]
# ### General Settings

# %%
# Directories
input_dir = '../data_input/'
output_dir = '../data_output/'
geospat_dir = '/Users/moyanofe/BigData/GeoSpatial/'

# Default data properties
epsg_crs = 4326

# Region files
conts_file = os.path.join(input_dir, 'continents.geojson')
# north50_file = os.path.join(input_dir, 'north50.geojson')
# north30_file = os.path.join(input_dir, 'north30.geojson')
# north0_file = os.path.join(input_dir, 'north0.geojson')
# testrect_file = os.path.join(input_dir, 'testrect.geojson')

# Set the region to use for the analysis
roi_file = conts_file

# Get the Carboscope projection to reproject other data
CS_file = geospat_dir + 'Carboscope/Inversions/' + 's10oc_v2022_daily.nc'
CS_inv = rio.open_rasterio(CS_file, variable='co2flux_land')
CS_inv = CS_inv['co2flux_land']
CS_inv.rio.write_crs(4326, inplace=True)
# %%
CS_inv.rio.transform()

# %%

# CarboScope inversions names and starting years
cs_inversions = ['sEXT', 's76', 's81', 's85', 's93', 's99', 's06', 's10']  # Inversion versions and starting year

# CAMS inversion names and starting years
cams_inversions = ['sur', 'sat'] # surface and satellite inversions

# %%
# Define functions to process data

def proc_CSinv(neeInv_file, roi_file = None, nodata_val = None):
    # Function to process CarboScope inversion files
    # Importantly, the data is averaged from daily to monthly in order to make it same as CAMS

    # All three file reading methods read in the attributes differently.
    # neeInvRaw = rio.open_rasterio(file_path)
    # neeInvRaw = xr.open_dataset(file_path, engine="rasterio", variable='co2flux_land')
    neeInvRaw_List = rio.open_rasterio(neeInv_file, variable='co2flux_land')
    neeInvRaw = neeInvRaw_List['co2flux_land']  # Extracting the dataarray from the list
    neeInvRaw = neeInvRaw.rename({'mtime': 'time'})  # Rename the time dimension and coordinate to something standard.
    neeInvRaw.attrs['units'] = 'PgC/yr'  # For some reason the units is a tuple of repeating values. Replacing with single value.
    neeInvRaw.rio.write_crs(4326, inplace=True)  # Set the crs: in this case the data is in epsg4326. This creates the spatial_ref coordinate.
    
    # Setting nodata
    if(nodata_val is not None):
        neeInvRaw = neeInvRaw.rio.write_nodata(nodata_val)  # Define what value represents nodata
        neeInvRaw = neeInvRaw.rio.write_nodata(neeInvRaw.rio.nodata, encoded=True, inplace=True)  # Set as encoded. This changes the values to NaN
    
    # Get the monthly mean
    neeInv = neeInvRaw.resample(time="1M").mean()

    # Get the pixel area and change the units to tC/ha/y
    pixarea_List = rio.open_rasterio(neeInv_file, variable='area')
    pixarea = pixarea_List['area']
    landarea = pixarea[0]
    landarea.attrs['units'] = 'm2'
    landarea.coords.__delitem__('rt')

    # Convert the NEE units: pixel to sqm, then to ha; and petagram to tons, so PgC y-1 -> tC ha-1 y-1
    neeInv = neeInv / landarea * 10000 * 10**9 # to ha, to tons.

    # Add some important attributes
    neeInv.name = 'co2flux_land'
    neeInv.attrs['units'] = 'tC/ha/y'
    neeInv.attrs['long_name'] = 'Net Ecosystem Exchange'

    if(roi_file is not None):
        roi = gpd.read_file(roi_file)
        neeInvS = neeInv.rio.clip(roi.geometry, roi.crs, all_touched=True)

    return(neeInvS)


def proc_CAMSinv(neeInv_file, inv, roi_file = None, nodata_val = None):
    # Function to process CAMS inversion files
    # Among other things, the time coordinate is set for each inversion separately

    # All three file reading methods read in the attributes differently
    # neeInvRaw = rio.open_rasterio(file_path)
    # neeInvRaw = xr.open_dataset(file_path, engine="rasterio", variable='co2flux_land')
    neeInvRaw_List = rio.open_rasterio(neeInv_file, variable='flux_apos_bio')
    neeInvRaw = neeInvRaw_List['flux_apos_bio']  # Extracting the dataarray from the list
    neeInvRaw.attrs['units'] = 'kgC/m2/month'  # For reference only, define the raw units
    neeInvRaw.rio.write_crs(4326, inplace=True)  # Set the crs: in this case the data is in epsg4326. This creates the spatial_ref coordinate.

    # Reproject to match resolution of Carboscope data
    neeInvRaw = neeInvRaw.rio.reproject_match(CS_inv)

    # Setting nodata
    if(nodata_val is not None):
        neeInvRaw = neeInvRaw.rio.write_nodata(nodata_val)  # Define what value represents nodata
        neeInvRaw = neeInvRaw.rio.write_nodata(neeInvRaw.rio.nodata, encoded=True, inplace=True)  # Set as encoded. This changes the values to NaN
    
    # Create a time coordinate variable
    neeInvRaw = neeInvRaw.rename({'band': 'time'})  # Rename the time dimension and coordinate to something standard.

    # Convert the NEE units: sqm to ha, kilogram to tons, month to year
    neeInv = neeInvRaw * 10000 / 1000 * 12 # to ha, to tons.
    neeInv.attrs['units'] = 'tC/ha/y'

    # Add a time coordinate
    if inv == 'sur':
        months = pd.date_range("1979-01", periods=504, freq='M')
    if inv == 'sat':
        months = pd.date_range("2010-01", periods=144, freq='M')
    neeInv["time"] = ('time', months)

    if(roi_file is not None):
        roi = gpd.read_file(roi_file)
        neeInvS = neeInv.rio.clip(roi.geometry, roi.crs, all_touched=True)
    
    # Add some important attributes
    neeInvS.attrs['long_name'] = 'Net Ecosystem Exchange'
    neeInvS.name = 'co2flux_land' # Make the co2 flux variable name equal to carboscope 

    return(neeInvS)


# %%
# Loop over inversion versions and analysis periods

# Process CarboScope inversions to obtain amplitudes
for inv in cs_inversions:

    # Chose the inversion data version
    fname = inv + 'oc_v2022_daily.nc'
    neeInv_file = geospat_dir + 'Carboscope/Inversions/' + fname
    
    print('Analysing inversion file: ', fname)

    # Call function for inversion NEE preprocessing
    neeProc = proc_CSinv(neeInv_file=neeInv_file, roi_file=roi_file)

    # Save to file
    file_out = os.path.join(output_dir, 'nee_proc/neeProc_CS'+inv+'.nc')
    neeProc.to_netcdf(file_out)


# %%
# Loop over inversions and analysis periods
    
# Process CAMS inversions to get amplitudes
for inv in cams_inversions:

    if inv == 'sur':
        neeInv_file = geospat_dir + 'CAMS/cams73_v20r2_co2_flux_surface_mm_1979-2020.nc'
    if inv == 'sat':
        neeInv_file = geospat_dir + 'CAMS/cams73_v21r2_co2_flux_satellite_mm_2010-2020.nc'
    
    print('Analysing inversion file: ', neeInv_file)

    neeProc = proc_CAMSinv(neeInv_file=neeInv_file, inv=inv, roi_file=roi_file)

    # Save to file
    file_out = os.path.join(output_dir, 'nee_proc/neeProc_CAMS'+inv+'.nc')
    neeProc.to_netcdf(file_out)

# %%
# -----
# Save the trend data to GeoTiff (useful if reading into Earth Engine is required)
# Saving the database (not as an array) keeps the band name info. (But seems like Earth Engßine does not import band names)
# file_out = os.path.join(output_dir, 'neeAmpStats_Inv' + ident + '.tif')
# neeAmpStats.rio.to_raster(file_out)

# Save the NEE flux seasonal difference data to GeoTiff
# file_out = os.path.join(output_dir, 'neeAmp_Inv' + ident + '.tif')
# da_fluxamp = neeAmp['nee_yearlyamp']
# da_fluxamp.rio.to_raster(file_out)
