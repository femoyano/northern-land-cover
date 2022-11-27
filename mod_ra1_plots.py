# mod_ra1_plots.py 
# Python module file: Plot fucntions for Recovery Analysis 1

# Imports
import os
from xarrayutils.utils import linear_trend
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# Directories
output_dir = '../figures/co2_landfluxes'


def plotOrthographic(xr_val, ident, fig_size=[16,6], cmap='RdYlGn', reverse_cmap=True):
    xr_val.attrs['units'] = 'tC/ha/y'
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
    fname =  os.path.join(output_dir, ident+'_ortho'+'.png')
    plt.savefig(fname)


def plotRegMeans(xr_val, lats, ident):
    # Calculate the spatial mean and plot.
    lats_str = str(lats)
    fig = plt.figure()
    ax1 = fig.add_subplot()
    for lat in lats:
        regDiffMean = xr_val.where(xr_val.y > lat[0]).where(xr_val.y < lat[1]).mean(['x','y'])
        label = str(lat[0])+'N to '+str(lat[1])+'N'
        del regDiffMean['spatial_ref']  # Remove coordinate so it doesn't show as title.
        regDiffMean.attrs['units'] = 'tC/ha'  # Add the new units
        regDiffMean.attrs['long_name'] = " Land CO2 flux winter-summer amplitude"
        regDiffMean.plot(ax=ax1, label=label)
    plt.legend()
    fname =  os.path.join(output_dir, 'regmean_'+lats_str+ident+'.png')
    # plt.show
    plt.savefig(fname)

def plotLambert(xarray, ident, fig_size=[12,5], cmap='RdYlGn'):
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
        vmin=-0.08, vmax=0.08,
        add_colorbar=False
    )
    # ax.set_global() # good for what?
    ax.coastlines()
    ax.set_title(xarray.attrs['long_name'])
    plt.colorbar(p, label=xarray.attrs['units'])
    # plt.show()
    fname =  os.path.join(output_dir, 'lambert_'+ident+'.png')
    plt.savefig(fname)

def plotOrthoSig(xr_val, xr_pval, p_lim, ident, fig_size=[16,6], cmap='RdYlGn'):
    # Map plot of flux trends globally
    # Plot the region of interest
    fig = plt.figure(figsize=fig_size)
    ax1 = fig.add_subplot(121, projection = ccrs.Orthographic(-70, 20), facecolor="white")
    ax2 = fig.add_subplot(111, projection = ccrs.Orthographic(70, 20), facecolor="white")
    xr_val_sig = xr_val.where(xr_pval < p_lim) # Change the p+value limit to filter for trend significance
    p1 = xr_val_sig.plot(ax = ax1, transform=ccrs.PlateCarree(), cmap=cmap, vmin=-0.08, vmax=0.08, add_colorbar=False)
    p2 = xr_val_sig.plot(ax = ax2, transform=ccrs.PlateCarree(), cmap=cmap, vmin=-0.08, vmax=0.08, add_colorbar=False)
    # ax1.set_global()
    ax1.coastlines()
    ax1.set_title(xr_val.attrs['long_name'])
    # ax2.set_global()
    ax2.coastlines()
    ax2.set_title(xr_val.attrs['long_name'])
    plt.colorbar(p1, label=xr_val.attrs['units'])
    # plt.colorbar(p2, label=co2InvT.isel(time=0).attrs['units'])
    # plt.show()
    fname =  os.path.join(output_dir, 'ortho_sig'+str(p_lim)+'_'+ident+'.png')
    plt.savefig(fname)
