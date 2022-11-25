# mod_ra1.py 
# Python module file: Recovery Analysis 1
#
# Functions to analyse CO2 fluxes and other land surface variables

# Imports
import rioxarray as rio
import geopandas as gpd
import numpy as np
from xarrayutils.utils import linear_trend
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


def procInvCO2(invco2_file, roi_file = None, nodata_val = None):

    # All three file reading methods read in the attributes differently.
    # co2InvRaw = rio.open_rasterio(file_path)
    # co2InvRaw = xr.open_dataset(file_path, engine="rasterio", variable='co2flux_land')
    co2InvRaw_List = rio.open_rasterio(invco2_file, variable='co2flux_land')
    co2InvRaw = co2InvRaw_List['co2flux_land']  # Extracting the dataarray from the list
    co2InvRaw = co2InvRaw.rename({'mtime': 'time'})  # Rename the time dimension and coordinate to something standard.
    co2InvRaw = co2InvRaw.rename({'x': 'longitude'})  # Rename the time dimension and coordinate to something standard.
    co2InvRaw = co2InvRaw.rename({'y': 'latitude'})  # Rename the time dimension and coordinate to something standard.
    co2InvRaw.attrs['units'] = 'PgC/yr'  # For some reason the units is a tuple of repeating values. Replacing with single value.
    co2InvRaw.rio.write_crs(4326, inplace=True)  # Set the crs (I assume this is the correct one). This creates the spatial_ref coordinate.

    # Setting nodata
    if(nodata_val is not None):
        co2InvRaw = co2InvRaw.rio.write_nodata(nodata_val)  # Define what value represents nodata
        co2InvRaw = co2InvRaw.rio.write_nodata(co2InvRaw.rio.nodata, encoded=True, inplace=True)  # Set as encoded. This changes the values to NaN
    
    # Get the pixel area and change the units to tC/ha/y
    pixarea_List = rio.open_rasterio(invco2_file, variable='area')
    pixarea = pixarea_List['area']
    pixarea = pixarea.rename({'x': 'longitude'})  # Rename the time dimension and coordinate to something standard.
    pixarea = pixarea.rename({'y': 'latitude'})  # Rename the time dimension and coordinate to something standard.
    landarea = pixarea[0]
    landarea.attrs['units'] = 'm2'
    landarea.coords.__delitem__('rt')

    # Convert the co2 flux units: pixel to sqm, then to ha; and petagram to tons, so PgC y-1 -> tC ha-2 y-1
    co2Inv = co2InvRaw / landarea * 10000 * 10**9 # to ha, to tons.

    # Add some important attributes
    co2Inv.attrs['units'] = 'tC/ha/y'
    co2Inv.attrs['long_name'] = 'Land-Atmosphere CO2 flux'

    if(roi_file is not None):
        roi = gpd.read_file(roi_file)
        co2InvS = co2Inv.rio.clip(roi.geometry, roi.crs, all_touched=True)

    return(co2InvS)

def get_CO2fluxSeasAmp(co2flux):
    """
    This function calculates the yearly aplitude between winter and summer months co2 fluxes.
    co2flux must be an xarray array with monthly data.
    
    # Calculating seasonal amplitudes and temporal amplitude trends
    1. First create a 'groups' coordinate from the combination of years and seasons (winter and summer: Oct-Mar, Apr-Sep).
    2. Calculate the mean flux rate for each year-season group.
        - Why the mean and not sum? Because these are rates, not totals. Could calculate the total seasonal flux by multiplying the mean rate by the 6 month time period.
    3. Calculate the winter-summer flux difference (amplitudes).
    4. Calculate the temporal trends in the amplitudes.
    """

    # Extract the month values
    months = co2flux.time.dt.month.data
    # Split year in two periods of 6 months each (two seasons) and code as: 0 = Oct-Mar, 1 = Apr-Sep
    season = np.where((months < 4) | (months > 9), 0, 1)
    # Create year-season groups for grouped calculations 
    groups = co2flux.time.dt.year.data + (season/10)
    # Set as coordinate
    co2flux['groups'] = ('time', groups)

    # Calculate the mean for each year-season group
    # co2InvTG = invco2flux.groupby(['time.year', 'season']).mean()  # Does not work! Grouping by two variables is not supported in xarray!
    co2fluxSeasRate = co2flux.groupby('groups').mean()  # Get the mean values of surface co2 fluxes by year-season
    # Calculate the total seasonal flux by multiplying the mean by the 6 month time period (i.e. 0.5 years)
    co2fluxSeas = co2fluxSeasRate * 0.5
    # Add again the coordinate representing years
    years = np.round(co2fluxSeas.coords['groups'].data)
    co2fluxSeas['year'] = ('groups', years)

    # Check results:
    # print(co2InvTG[0:,0,0])
    # print(co2InvTG.coords['groups'])
    # print(co2InvTG.coords['year'])
    # co2InvTG.where(co2InvTG.year==1985, drop=True)

    # Calculate the difference between winter (positive) and summer (negative) mean fluxes. For the southern hemisphere the sign is corrected afterwards.
    def diff(x):
        return(x[0] - x[1])  # Oct-Mar average flux minus Apr-Sep average flux. Signs are correct for NH. SH sign is corrected below.
    co2fluxSeasDiff = co2fluxSeas.groupby('year').map(diff)  # Apply the difference function for each year
    co2fluxSeasDiff = co2fluxSeasDiff.where(co2fluxSeasDiff.latitude > 0, co2fluxSeasDiff * (-1))  # Correct the sign in the Southern Hemisphere

    # Rename the xarray data array and setting the crs
    co2fluxSeasDiff = co2fluxSeasDiff.rename('co2flux_yearlyamp')
    co2fluxSeasDiff.rio.write_crs(4326, inplace=True)

    return(co2fluxSeasDiff)


def get_co2fluxampstats(co2fluxAmp):

    fluxamp_stats = linear_trend(co2fluxAmp, 'year')
    fluxamp_stats.slope.attrs['units'] = 'tC/ha/y'
    fluxamp_stats.slope.attrs['long_name'] = 'Trend in seasonal CO2 flux amplitude'
    fluxamp_stats.rio.write_crs(4326, inplace=True)

    # Variance
    fluxamp_var = co2fluxAmp.var(dim='year')
    fluxamp_var.attrs['units'] = 'tC/ha/y'
    fluxamp_var.attrs['long_name'] = 'Temporal variance in seasonal CO2 flux amplitude'
    fluxamp_var.rio.write_crs(4326, inplace=True)
    fluxamp_var = fluxamp_var.rename('co2fluxamp_variance')

    # Mean
    fluxamp_mean = co2fluxAmp.mean(dim='year')
    fluxamp_mean.attrs['units'] = 'tC/ha/y'
    fluxamp_mean.attrs['long_name'] = 'Temporal mean in seasonal CO2 flux amplitude'
    fluxamp_mean.rio.write_crs(4326, inplace=True)
    fluxamp_mean = fluxamp_mean.rename('co2fluxamp_mean')

    fluxamp_stats['mean'] = fluxamp_mean
    fluxamp_stats['variance'] = fluxamp_var

    return(fluxamp_stats)

def plotOrthographic(co2fluxProc, fig_size=[16,6], cmap='RdYlGn', reverse_cmap=True):
    fig = plt.figure(figsize=fig_size)
    ax1 = fig.add_subplot(121, projection = ccrs.Orthographic(-70, 20), facecolor="white")
    ax2 = fig.add_subplot(111, projection = ccrs.Orthographic(70, 20), facecolor="white")
    c_map=plt.cm.get_cmap(cmap)
    if reverse_cmap:
        c_map = c_map.reversed()
    p1 = co2fluxProc.isel(time=0).plot(ax = ax1, transform=ccrs.PlateCarree(), cmap=c_map, add_colorbar=False)
    p2 = co2fluxProc.isel(time=0).plot(ax = ax2, transform=ccrs.PlateCarree(), cmap=c_map, add_colorbar=False)
    ax1.set_global()
    ax1.coastlines()
    ax1.set_title(co2fluxProc.isel(time=0).attrs['long_name'])
    ax2.set_global()
    ax2.coastlines()
    ax2.set_title(co2fluxProc.isel(time=0).attrs['long_name'])
    plt.colorbar(p1, label=co2fluxProc.isel(time=0).attrs['units'])
    # plt.colorbar(p2, label=co2InvT.isel(time=0).attrs['units'])
    plt.show()

def plotRegMeans(co2fluxAmp, lats):
    # Calculate the spatial mean and plot.

    for lat in lats:
        regDiffMean = co2fluxAmp.where(co2fluxAmp.latitude > lat[0]).where(co2fluxAmp.latitude < lat[1]).mean(['longitude','latitude'])
        label = str(lat[0])+'N to '+str(lat[1])+'N'
        del regDiffMean['spatial_ref']  # Remove coordinate so it doesn't show as title.
        regDiffMean.attrs['units'] = 'tC/ha'  # Add the new units
        regDiffMean.attrs['long_name'] = " Land CO2 flux winter-summer amplitude"
        regDiffMean.plot(label=label)
    plt.legend()

def plotLambert(xarray, fig_size=[12,5], cmap='RdYlGn'):
    fig = plt.figure(figsize=fig_size)
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
    p = xarray.plot(
        ax = ax,
        transform=ccrs.PlateCarree(),
        cmap=cmap,
        add_colorbar=False
    )
    # ax.set_global() # good for what?
    ax.coastlines()
    ax.set_title(xarray.attrs['long_name'])
    plt.colorbar(p, label=xarray.attrs['units'])
    plt.show()

def plotOrthoSig(xr_val, xr_pval, p_lim, fig_size=[16,6], cmap='RdYlGn'):
    # Map plot of flux trends globally
    # Plot the region of interest
    fig = plt.figure(figsize=fig_size)
    ax1 = fig.add_subplot(121, projection = ccrs.Orthographic(-70, 20), facecolor="white")
    ax2 = fig.add_subplot(111, projection = ccrs.Orthographic(70, 20), facecolor="white")
    xr_val_sig = xr_val.where(xr_pval < p_lim) # Change the p+value limit to filter for trend significance
    p1 = xr_val_sig.plot(ax = ax1, transform=ccrs.PlateCarree(), cmap=cmap, add_colorbar=False)
    p2 = xr_val_sig.plot(ax = ax2, transform=ccrs.PlateCarree(), cmap=cmap, add_colorbar=False)
    # ax1.set_global()
    ax1.coastlines()
    ax1.set_title(xr_val.attrs['long_name'])
    # ax2.set_global()
    ax2.coastlines()
    ax2.set_title(xr_val.attrs['long_name'])
    plt.colorbar(p1, label=xr_val.attrs['units'])
    # plt.colorbar(p2, label=co2InvT.isel(time=0).attrs['units'])
    plt.show()
