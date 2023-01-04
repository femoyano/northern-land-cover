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
# ## NEE yearly flux amplitude analysis
#
# This script gets the processed inversion data and calculates yearly amplitudes and related stats 
#
# Steps include:
# - Load inversion NEE data
# - Create a grouping coordinate that divides growing and non-growing seasons.   
#     - This can be done manually, but consider fitting a harmonic to divide seasons.  
# - For each year and season, calculate the NEE flux sum.
# - For each year, get the difference between growing and non-growing season.
#     - This is the NEE flux amplitude that can then be related to atmospheric CO2 and driving variables. 
#

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
from xarrayutils.utils import linear_trend

# Local modules
from mod_ra1 import *

# %% [markdown]
# ### General Settings

# %%
# Directories
input_dir = '../data_input/'
output_dir = '../data_output/'
proc_nee_dir = '../data_output/nee_proc/'

# Default data properties
epsg_crs = 4326

# Region files available
conts_file = os.path.join(input_dir, 'continents.geojson')
# north50_file = os.path.join(input_dir, 'north50.geojson')
# north30_file = os.path.join(input_dir, 'north30.geojson')
# north0_file = os.path.join(input_dir, 'north0.geojson')
# testrect_file = os.path.join(input_dir, 'testrect.geojson')

# Set the region to use for the analysis
roi_file = conts_file

# %%
# Inversion names and starting years
inversions = {'CSsEXT':1957, 'CSs76':1976, 'CSs81':1981, 'CSs85':1985, 'CSs93':1993, 'CSs99':1999,
            'CSs06':2006, 'CSs10':2010, 'CAMSsur':1979, 'CAMSsat': 2010}
# cs_inversions = {'CSsEXT':1957}  # To run selected only

# Periods for temporal analysis
# periods = [(1957,2021), (1976,2021), (1981,2021), (1985,2021), (1993,2021), (1999,2021), (2006,2021), (2010,2021)] # based on CarboScope periods
periods = [(1961,2021), (1981,2021), (2001,2021), (2010,2021)]
# periods = [(1993,2021)]

# %%
# Loop over inversions and analysis periods

# Process inversions to obtain amplitudes
for inv in inversions:
    inv = 'CAMSsur'
    # Chose the inversion data version
    fname = 'neeProc_' + inv + '.nc'
    neeInv_file = proc_nee_dir + fname
    
    inv_start = inversions[inv]

    print('Analysing inversion file: ', fname)

    neeProc = xr.open_dataset(neeInv_file)

    # Call function to get the yearly seasonal NEE flux amplitude
    neeAmp = get_neeSeasAmp(neeProc)

    # Save amps to file
    file_out_amp = os.path.join(output_dir, 'neeAmp_CS_'+inv_start+'.nc')
    neeAmp.to_netcdf(file_out_amp)


# %%
# Calculate amplitude stats

for inv in inversions:

    inv_start = inversions[inv]

    # read amps file
    file_in_amp = os.path.join(output_dir, 'neeAmp_Inv'+str(inv_start)+'.nc')
    neeAmp = rio.open_rasterio(file_in_amp)
    
    # Loop over periods to get temporal statistics for each:

    for period in periods:

        print('inv_start: ', inv_start)
        print('period: ', period)
    
        if(inv_start > period[0]):
            print('Period not within inversion run. Skipping ...')
            continue

        print('Executing analysis with: ', inv, ', period: ', period)

        # Call function to calculate temporal trends, mean, variance, etc in seasonal amplitudes
        neeAmpStats = get_neeAmpStats(neeAmp, period=period)

        # Save amp stats
        ident = str(period[0])+'-'+str(period[1])+'_Inv'+str(inv_start)
        file_out_stats = os.path.join(output_dir, 'neeAmpStats_'+ident+'.nc')
        neeAmpStats.to_netcdf(file_out_stats)


