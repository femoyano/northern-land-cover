# mod_ra1_plots.py 
# Python module file: Plot fucntions for Recovery Analysis 1

# Imports
from xarrayutils.utils import linear_trend
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


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
