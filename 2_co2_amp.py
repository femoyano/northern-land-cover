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

# Inversion names and starting years
inversions = {'CSsEXT':1957, 'CSs76':1976, 'CSs81':1981, 'CSs85':1985, 'CSs93':1993, 'CSs99':1999,
            'CSs06':2006, 'CSs10':2010, 'CAMSsur':1979, 'CAMSsat': 2010}
# cs_inversions = {'CSsEXT':1957}  # To run selected only

# Periods for temporal analysis
# periods = [(1957,2021), (1976,2021), (1981,2021), (1985,2021), (1993,2021), (1999,2021), (2006,2021), (2010,2021)] # based on CarboScope periods
periods = [(1961,2021), (1981,2021), (2001,2021), (2010,2021)]
# periods = [(1993,2021)]

# %%
# Define functions to calculate flux amplitude and stats

def get_neeSeasAmp(nee):
    """
    This function calculates the yearly aplitude between winter and summer months C fluxes.
    nee must be an xarray array with monthly data.
    
    # Calculating seasonal amplitudes and temporal amplitude trends
    1. First create a 'groups' coordinate from the combination of years and seasons (winter and summer: Oct-Mar, Apr-Sep).
    2. Calculate the mean flux rate for each year-season group.
        - Why the mean and not sum? Because these are rates, not totals. Could calculate the total seasonal flux by multiplying the mean rate by the 6 month time period.
    3. Calculate the winter-summer flux difference (amplitudes).
    4. Calculate the temporal trends in the amplitudes.
    """

    # Extract the month values
    months = nee.time.dt.month.data
    # Split year in two periods of 6 months each (two seasons) and code as: 0 = Oct-Mar, 1 = Apr-Sep
    season = np.where((months < 4) | (months > 9), 0, 1)
    # Create year-season groups for grouped calculations 
    groups = nee.time.dt.year.data + (season/10)
    # Set as coordinate
    nee['groups'] = ('time', groups)

    # Calculate the mean for each year-season group
    # neeInvTG = invnee.groupby(['time.year', 'season']).mean()  # Does not work! Grouping by two variables is not supported in xarray!
    neeSeasRate = nee.groupby('groups').mean()  # Get the mean values of surface C fluxes by year-season
    # Calculate the total seasonal flux by multiplying the mean by the 6 month time period (i.e. 0.5 years)
    neeSeas = neeSeasRate * 0.5
    # Add again the coordinate representing years
    years = np.round(neeSeas.coords['groups'].data)
    neeSeas.coords['year'] = ('groups', years)

    # Check results:
    # print(neeInvTG[0:,0,0])
    # print(neeInvTG.coords['groups'])
    # print(neeInvTG.coords['year'])
    # neeInvTG.where(neeInvTG.year==1985, drop=True)

    # Calculate the difference between winter (positive) and summer (negative) mean fluxes. For the southern hemisphere the sign is corrected afterwards.
    def diff(x):
        return(x[0] - x[1])  # Oct-Mar average flux minus Apr-Sep average flux. Signs are correct for NH. SH sign is corrected below.
    neeSeasDiff = neeSeas.groupby('year').map(diff)  # Apply the difference function for each year
    neeSeasDiff = neeSeasDiff.where(neeSeasDiff.y > 0, neeSeasDiff * (-1))  # Correct the sign in the Southern Hemisphere

    # Rename the xarray data array and setting the crs
    neeSeasDiff = neeSeasDiff.rename('nee_yearlyamp')
    neeSeasDiff.rio.write_crs(4326, inplace=True)
    neeSeasDiff.attrs['long_name'] = 'NEE seasonal amplitude'
    neeSeasDiff.attrs['units'] = 'tC/ha/y'

    return(neeSeasDiff)


def get_neeAmpStats(neeAmp, period=None):

    # Selecting a time period
    if(period is not None):
        neeAmp = neeAmp.sel(year=slice(period[0], period[1]))

    fluxamp_stats = linear_trend(neeAmp, 'year')
    fluxamp_stats.slope.attrs['units'] = 'tC/ha/y'
    fluxamp_stats.slope.attrs['long_name'] = 'Trend in seasonal NEE amplitude'
    fluxamp_stats.rio.write_crs(4326, inplace=True)

    # Variance
    fluxamp_var = neeAmp.var(dim='year')
    fluxamp_var.attrs['units'] = 'tC/ha/y'
    fluxamp_var.attrs['long_name'] = 'Temporal variance in seasonal NEE amplitude'
    fluxamp_var.rio.write_crs(4326, inplace=True)
    fluxamp_var = fluxamp_var.rename('neeamp_variance')

    # Mean
    fluxamp_mean = neeAmp.mean(dim='year')
    fluxamp_mean.attrs['units'] = 'tC/ha/y'
    fluxamp_mean.attrs['long_name'] = 'Temporal mean in seasonal NEE amplitude'
    fluxamp_mean.rio.write_crs(4326, inplace=True)
    fluxamp_mean = fluxamp_mean.rename('neeamp_mean')

    fluxamp_stats['mean'] = fluxamp_mean
    fluxamp_stats['variance'] = fluxamp_var

    return(fluxamp_stats)



# %%
# Loop over inversions and analysis periods

# Process inversions to obtain amplitudes
for inv in inversions:

    # Chose the inversion data version
    fname = 'neeProc_' + inv + '.nc'
    neeInv_file = proc_nee_dir + fname
    
    inv_start = inversions[inv]

    print('Analysing inversion file: ', fname)

    neeProcDS = xr.open_dataset(neeInv_file)
    neeProc = neeProcDS['co2flux_land']

    # Call function to get the yearly seasonal NEE flux amplitude
    neeAmp = get_neeSeasAmp(neeProc)

    # Save amps to file
    file_out_amp = os.path.join(output_dir, 'neeAmp_'+ inv +'.nc')
    neeAmp.to_netcdf(file_out_amp)


# %%
# Calculate amplitude stats

for inv in inversions:

    inv_start = inversions[inv]
    print('------------\nAnalysing inversion: ', inv)
    print('Inversion start: ', inv_start)

    # read amps file
    file_in_amp = os.path.join(output_dir, 'neeAmp_'+inv+'.nc')
    neeAmp = rio.open_rasterio(file_in_amp)
    
    # Loop over periods to get temporal statistics for each:

    for period in periods:

        if(inv_start > period[0]):
            print('Period ', period, ' not within inversion run. Skipping ...')
            continue

        print('Analysing period: ', period)

        # Call function to calculate temporal trends, mean, variance, etc in seasonal amplitudes
        neeAmpStats = get_neeAmpStats(neeAmp, period=period)

        # Save amp stats
        ident = str(period[0]) + '-' + str(period[1]) + '_' + inv
        file_out_stats = os.path.join(output_dir, 'neeAmpStats_'+ident+'.nc')
        neeAmpStats.to_netcdf(file_out_stats)

