# mod_ra1.py 
# Python module file: Recovery Analysis 1
#
# Functions to analyse C fluxes and other land surface variables

# Imports
import rioxarray as rio
import geopandas as gpd
import numpy as np
from xarrayutils.utils import linear_trend


def proc_neeInv(neeInv_file, roi_file = None, nodata_val = None):

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
    
    # Get the pixel area and change the units to tC/ha/y
    pixarea_List = rio.open_rasterio(neeInv_file, variable='area')
    pixarea = pixarea_List['area']
    landarea = pixarea[0]
    landarea.attrs['units'] = 'm2'
    landarea.coords.__delitem__('rt')

    # Convert the NEE units: pixel to sqm, then to ha; and petagram to tons, so PgC y-1 -> tC ha-1 y-1
    neeInv = neeInvRaw / landarea * 10000 * 10**9 # to ha, to tons.

    # Add some important attributes
    neeInv.attrs['units'] = 'tC/ha/y'
    neeInv.attrs['long_name'] = 'Net Ecosystem Exchange'

    if(roi_file is not None):
        roi = gpd.read_file(roi_file)
        neeInvS = neeInv.rio.clip(roi.geometry, roi.crs, all_touched=True)

    return(neeInvS)

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
    neeSeas['year'] = ('groups', years)

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
