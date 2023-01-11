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
# ## Plotting and Regional Analysis

# %%
# Magic commands
# %load_ext autoreload
# %autoreload 2

# Packages
import os
import rioxarray as rio
import pandas as pd
import geopandas as gpd
import shapely as sh
from xarrayutils.utils import linear_trend
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# Local modules
from mod_ra1_plots import *

# %%
# Directories
input_dir = '../data_input/'
output_dir = '../data_output/'
figs_dir = '../figures/'

inversions = ['CSsEXT', 'CSs76', 'CSs81', 'CSs85', 'CSs93', 'CSs99', 'CSs06', 'CSs10', 'CAMSsur', 'CAMSsat']

# Choose the list of inversion data version to analyse. Later (shorter) time series are driven by data from more stations.
inv_startyear = {'CSsEXT':1957, 'CSs76':1976, 'CSs81':1981, 'CSs85':1985, 'CSs93':1993, 'CSs99':1999,
            'CSs06':2006, 'CSs10':2010, 'CAMSsur':1979, 'CAMSsat': 2010}  # Inversion names and starting year
# inv_startyear = {'CSs99':1999}

# Number of stations of each inversion version
inv_nstations = {'CSsEXT': 169,'CSs76':9, 'CSs81':14, 'CSs85':21, 'CSs93':35, 'CSs99':49, 'CSs06':59, 'CSs10':78, 'CAMSsur':-999, 'CAMSsat':0}

# List of start and end years to define analysis periods
periods = [(1961,2021), (1981,2021), (2001,2021), (2010,2021)]
# periods = [(1957,2021), (1976,2021), (1981,2021), (1985,2021), (1993,2021), (1999,2021), (2006,2021), (2010,2021)] # These corrspond to carboscope inv start years
# periods = [(2001,2021)]

# %%
# Functions for plotting maps

def plotOrthographic(xr_val, ident, fig_size=[16,6], cmap='RdYlGn', reverse_cmap=True):
    xr_val.attrs['units'] = 'tC/ha/y'
    var = xr_val.name
    fig = plt.figure(figsize=fig_size)
    ax1 = fig.add_subplot(121, projection = ccrs.Orthographic(-70, 20), facecolor="white")
    ax2 = fig.add_subplot(111, projection = ccrs.Orthographic(70, 20), facecolor="white")
    c_map=plt.cm.get_cmap(cmap)
    if reverse_cmap:
        c_map = c_map.reversed()
    p1 = xr_val.plot(ax = ax1, transform=ccrs.PlateCarree(), cmap=c_map, add_colorbar=False)
    p2 = xr_val.plot(ax = ax2, transform=ccrs.PlateCarree(), cmap=c_map, add_colorbar=False)
    ax1.set_global()
    ax1.coastlines()
    ax1.set_title(xr_val.attrs['long_name'])
    ax2.set_global()
    ax2.coastlines()
    ax2.set_title(xr_val.attrs['long_name'])
    plt.colorbar(p1, label=xr_val.attrs['units'])
    # plt.colorbar(p2, label=co2InvT.isel(time=0).attrs['units'])
    # plt.show()
    fname =  os.path.join(figs_dir, 'co2_landfluxes', var+'_ortho_'+ident+'.png')
    plt.savefig(fname)


def plotLambert(xr_val, ident, fig_size=[12,5], cmap='RdYlGn'):
    fig = plt.figure(figsize=fig_size)
    var = xr_val.name
    ax = fig.add_subplot(
        111,
        projection = ccrs.LambertConformal(  # projection=ccrs.Orthographic(0, 90)
            central_longitude=10.0,
            central_latitude=39.0,
            false_easting=0.0,
            false_northing=0.0,
            standard_parallels=(33, 45),
            globe=None, cutoff=30),
        facecolor="white"
    )
    p = xr_val.plot(
        ax = ax,
        transform=ccrs.PlateCarree(),
        cmap=cmap,
        vmin=-0.08, vmax=0.08,
        add_colorbar=False
    )
    # ax.set_global() # good for what?
    ax.coastlines()
    ax.set_title(xr_val.attrs['long_name']+'\n'+ident)
    plt.colorbar(p, label=xr_val.attrs['units'])
    # plt.show()
    fname =  os.path.join(figs_dir, 'co2_landfluxes', var+'_lambert_'+ident+'.png')
    plt.savefig(fname)

def plotOrthoSig(xr_val, ident, p_lim=None, xr_pval=None, vmin=None, vmax=None, fig_size=[16,6], cmap='RdYlGn'):
    # Map plot of flux trends globally
    # Plot the region of interest
    var = xr_val.name
    fig = plt.figure(figsize=fig_size)
    ax1 = fig.add_subplot(121, projection = ccrs.Orthographic(-70, 20), facecolor="white")
    ax2 = fig.add_subplot(111, projection = ccrs.Orthographic(70, 20), facecolor="white")
    if(xr_pval is not None):
        xr_val = xr_val.where(xr_pval < p_lim) # Filter for trend significance
    p1 = xr_val.plot(ax = ax1, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax, add_colorbar=False)
    p2 = xr_val.plot(ax = ax2, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax, add_colorbar=False)
    ax1.coastlines()
    ax2.coastlines()
    ax1.set_title('')
    ax2.set_title('')
    fig.suptitle(xr_val.attrs['long_name']+'\n'+ident)
    plt.colorbar(p1, label=xr_val.attrs['units'])
    # plt.colorbar(p2, label=co2InvT.isel(time=0).attrs['units'])
    # plt.show()
    fname =  os.path.join(figs_dir, 'co2_landfluxes', var+'_'+'ortho_sig'+str(p_lim)+'_'+ident+'.png')
    plt.savefig(fname)



# %%
# 
for inv in inv_startyear:

    inv_start = inv_startyear[inv]

    for period in periods:

        print('inv_start: ', inv_start)
        print('period: ', period)
    
        if(inv_start > period[0]):
            print('Continuing ...')
            continue

        print('Reading data: ', inv, ', period: ', period)

        ident = str(period[0])+'-'+str(period[1])+'_'+inv
        file_in_stats = os.path.join(output_dir, 'neeAmpStats_'+ident+'.nc')
        neeAmpStats = rio.open_rasterio(file_in_stats)

        # Map plot of flux trends
        # plotLambert(xr_val=neeAmpStats.slope, ident=ident)

        # Plot slope with significance filter = 0.1
        # plotOrthoSig(xr_val=neeAmpStats.slope, ident=ident, p_lim=0.1, xr_pval=neeAmpStats.p_value, vmin=-0.08, vmax=0.08)

        # Plot slope without significance filter
        plotOrthoSig(xr_val=neeAmpStats.slope, ident=ident, p_lim=1, xr_pval=neeAmpStats.p_value, vmin=-0.08, vmax=0.08)

        # Plot variance without significance filter
        plotOrthoSig(xr_val=neeAmpStats.variance, ident=ident, p_lim=1, xr_pval=neeAmpStats.p_value, vmin=0, vmax=0.2, cmap='Oranges')
            

# %% [markdown]
# ### Regional Analysis
#
# Specific regions that coincide with NEE-amp trend or variance hotspots are selected for analysis.
# A mask for selection is created by either:
#
# 1. clipping the data to a continent and selecting the KG climate region that overlaps with the hotspot, or
# 2. clipping to any geometry, such as latitudinal bands

# %%
# Plot functions for time series

def plotInvMeans(inv_names, regions):
    # This function will create one plot per inversion with lines showing the temporal mean for different region
    for inv in inv_startyear:
        file_in_amp = os.path.join(output_dir, 'neeAmp_'+ inv +'.nc')
        xr_val = rio.open_rasterio(file_in_amp)
        var = xr_val.name
        fig = plt.figure()
        ax1 = fig.add_subplot()
        for roi_name, roi in regions.items():
            # Calculate the spatial mean and plot.
            try:
                regDiffMean = xr_val.where(roi).mean(['x','y'])
            except:
                regDiffMean = xr_val.rio.clip([roi]).mean(['x','y'])
            label = roi_name
            del regDiffMean['spatial_ref']  # Remove coordinate so it doesn't show as title.
            regDiffMean.attrs['units'] = 'tC/ha'  # Add the new units
            regDiffMean.attrs['long_name'] = "Land CO2 flux winter-summer amplitude"
            regDiffMean.plot(ax=ax1, label=label)
        plt.legend()
        ax1.set_title(regDiffMean.attrs['long_name']+'\n'+inv)
        fname =  os.path.join(figs_dir, 'time_series', var+'_'+inv+'.png')
        # plt.show
        plt.savefig(fname)

def plotRegionMeans(regions, inv_names):
    # This function will create one plot per region with lines showing the temporal mean of different inversions
    for roi_name, roi in regions.items():
        tsDict = dict()
        fig = plt.figure(figsize=[10,5])
        ax1 = fig.add_subplot()
        for inv in inv_names:
            label = inv
            file_in_amp = os.path.join(output_dir, 'neeAmp_'+ inv +'.nc')
            xr_val = rio.open_rasterio(file_in_amp)
            # Calculate the spatial mean and plot.
            try:
                regDiffMean = xr_val.where(roi).mean(['x','y'])
            except:
                regDiffMean = xr_val.rio.clip([roi]).mean(['x','y'])
            regDiffMean.attrs['units'] = 'tC/ha'  # Add the new units
            regDiffMean.attrs['long_name'] = "Land CO2 flux winter-summer amplitude"
            regDiffMean.plot(ax=ax1, label=label)
        plt.legend()
        ax1.set_title(regDiffMean.attrs['long_name']+'\n'+roi_name)
        var = xr_val.name
        fname =  os.path.join(figs_dir, 'time_series', var+'_'+roi_name+'.png')
        plt.savefig(fname)


# %%
# ---- Create masks and geometries for filtering and plotting data ----

# Load the nee data to reproject_match other data
file_neeAmpStats = os.path.join(output_dir, 'neeAmpStats_2010-2021_CSsEXT.nc')
neeAmp = rio.open_rasterio(file_neeAmpStats)

# Load continent data to use as filter 
conts = gpd.read_file('../data_input/continents.geojson')
cont_bounds = {'Asia': [24, 0, 190, 81], 'Europe': [-31, 35, 69, 81], 'North America': [-178, 0, -15, 84]} # This avoids plotting the entire globe

# Create masks using Koeppen-Geiger climate regions ---

# Load the Koeppen-Geiger data
kopp = rio.open_rasterio('../data_input/VU-VIENA/KG_1986-2010.grd', mask_and_scale=False)
kopp.rio.write_nodata(32, inplace=True) # Sets water as missing data
kopp = kopp.where(kopp != kopp.rio.nodata) # Not sure this is necessary
kopp_rp = kopp.rio.reproject_match(neeAmp) # Reproject to match inversion

# Define regions for further filtering
# geom_eastasia_tempwarm = sh.geometry.box(100, 20, 125, 35)
# geom_eastasia_tempcold = sh.geometry.box(115, 35, 135, 50)
geom_eastasia = sh.geometry.box(100, 0, 179, 90)
# geom_europe_temp = sh.geometry.box(-10, 40, 40, 55)
# geom_northamer_tempcold = sh.geometry.box(-130, 40. -50, 60)

# Create regional masks using Koeppen-Geiger climate regions ---
kopp_mask = dict()

cont = conts[conts['CONTINENT'] == 'Europe']
kopp_sel = kopp_rp.rio.clip(cont.geometry, cont.crs, all_touched=True, drop=False).squeeze()
# mask = kopp_sel.where(kopp_sel.isin([10])) # This would select the values instead of making a 0,1 mask
kopp_mask['europe_temp'] = kopp_sel.isin([10])

cont = conts[conts['CONTINENT'] == 'Asia']
kopp_sel = kopp_rp.rio.clip([geom_eastasia], drop=False).rio.clip(cont.geometry, cont.crs, all_touched=True, drop=False).squeeze()
kopp_mask['asia_tempwarm'] = kopp_sel.isin([9, 15])
kopp_mask['asia_tempcold'] = kopp_sel.isin([6, 26, 27])

cont = conts[conts['CONTINENT'] == 'North America']
kopp_sel = kopp_rp.rio.clip(cont.geometry, cont.crs, all_touched=True, drop=False).squeeze()
kopp_mask['northamer_tempcold'] = kopp_sel.isin([18, 19])


# Create geometries of latitudinal bands ---

geom_lat10_30 = sh.geometry.box(-179, 10, 179, 30)
geom_lat30_50 = sh.geometry.box(-179, 30, 179, 50)
geom_lat50_70 = sh.geometry.box(-179, 50, 179, 70)
geom_lat70_90 = sh.geometry.box(-179, 70, 179, 90)

lat_regions = {'lat10-30': geom_lat10_30, 'lat30-50': geom_lat30_50, 'lat50-70': geom_lat50_70, 'lat70-90': geom_lat70_90}


# %%

# List of inversions to plot.
inv_names = ['CSsEXT', 'CSs76', 'CSs81', 'CSs85', 'CSs93', 'CSs99', 'CSs06', 'CSs10', 'CAMSsur', 'CAMSsat']
# inv_names = ['CSs76', 'CSs81', 'CSs85', 'CSs93', 'CSs99', 'CSs06', 'CSs10']


plotInvMeans(inv_names=inv_names, regions=kopp_mask)
plotInvMeans(inv_names=inv_names, regions=lat_regions)

plotRegionMeans(lat_regions, inv_names)
plotRegionMeans(kopp_mask, inv_names)
