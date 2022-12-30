# %% [markdown]
# ## NEE flux analysis
# 
# This script gets the raw inversion data and calculates 
# 
# Later (shorter) time series are driven by data from more stations.
# 
# ### Methods:
# 
# 1. Load inversion NEE data and preprocess:
#     - select roi
#     - adjust units
#     - set nodata
#     - ...
# 3. Create a grouping coordinate that divides growing and non-growing seasons.   
#     - This can be done manually, but consider fitting a harmonic to divide seasons.  
# 4. For each year and season, calculate the NEE flux sum.
# 5. For each year, get the difference between growing and non-growing season.
#     - This is the NEE flux amplitude that can then be related to atmospheric CO2 and driving variables. 
# 

#%%
# Magic commands to updated functions when python module is modified
# %load_ext autoreload
# %autoreload 2

# Get all required packages

# Packages
import os
import rioxarray as rio
import xarray as xr

# Local modules
from mod_ra1 import *

# %% [markdown]
# ### Options

# %%
# Choose the list of inversion data version to analyse. Later (shorter) time series are driven by data from more stations.
inversions = {'sEXT':1957, 's76':1976, 's81':1981, 's85':1985, 's93':1993, 's99':1999, 's06':2006, 's10':2010}  # Inversion versions and starting year
# inversions = {'s76':1976}  # To run selected only

# Choose periods for temporal stats
periods = [(1957,2021), (1976,2021), (1981,2021), (1985,2021), (1993,2021), (1999,2021), (2006,2021), (2010,2021)]
# periods = [(1993,2021)]

# Choose what analysis step to make
do_getamp = True
do_getstats = True

# %% [markdown]
# ### General Settings

# %%
# Directories
input_dir = '../data_input/'
output_dir = '../data_output/'

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
# Loop over inversion versions and analysis periods
if do_getamp:
    
    for inv in inversions:

        # Chose the inversion data version. Later (shorter) time series are driven by data from more stations.
        fname = inv + 'oc_v2022_daily.nc'
        neeInv_file = '/Users/moyanofe/BigData/GeoSpatial/Carboscope/Inversions/' + fname
        
        inv_start = inversions[inv]

        print('Analysing inversion file: ', fname)

        # Call function for inversion NEE preprocessing
        neeProc = proc_neeInv(neeInv_file=neeInv_file, roi_file=roi_file)

        # Call function to get the yearly seasonal NEE flux amplitude
        neeAmp = get_neeSeasAmp(neeProc)

        # Save amps to file
        file_out_amp = os.path.join(output_dir, 'neeAmp_Inv'+str(inv_start)+'.nc')
        neeAmp.to_netcdf(file_out_amp)


# -----
# Save the trend data to GeoTiff (useful if reading into Earth Engine is required)
# Saving the database (not as an array) keeps the band name info. (But seems like Earth Engine does not import band names)
# file_out = os.path.join(output_dir, 'neeAmpStats_Inv' + ident + '.tif')
# neeAmpStats.rio.to_raster(file_out)

# Save the NEE flux seasonal difference data to GeoTiff
# file_out = os.path.join(output_dir, 'neeAmp_Inv' + ident + '.tif')
# da_fluxamp = neeAmp['nee_yearlyamp']
# da_fluxamp.rio.to_raster(file_out)

# %%
if do_getstats:

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


