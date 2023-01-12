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
# ## Regional Analysis
#
# ### Steps
#
# #### Step 1: Selection of regions of interest
#
# In this step we rank regions by their NEE flux seasonal amplitude variance over time.
# For each region, the temporal NEE amp variance and total area are calculated.
# Region to analyze further can then be selcted by ranking them according to variance or weighted variance
#
# 1. Load regions data
# 2. Loop over regions. For each:
#     - Calculate the variance in NEE flux seasonal amplitude
#     - Calculate the total area
# 3. Calculate global sum of variance*area
# 4. Replace raster region indices with variance values and plot maps
#
# ### Notes:
#
# For regions, currently using the Koppen-Geiger data from : http://koeppen-geiger.vu-wien.ac.at/
#

# %%
import os
import pandas as pd
import rioxarray as rio
import xarray as xr
import geopandas as gpd
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import shapely as sh

# %%
# Settings

figs_dir = '../figures/regions/'
output_dir = '../data_output/'
conts = gpd.read_file('../data_input/continents.geojson')

# %%

# Reading in regions and setting nodata values

# kopp = rio.open_rasterio('../data_input/KoeppenGeiger3_KG3_CRUTS32_Hist_7100.tif', mask_and_scale=False)
kopp = rio.open_rasterio(
    '../data_input/VU-VIENA/KG_1986-2010.grd', mask_and_scale=False)
kopp.rio.write_nodata(32, inplace=True)
kopp = kopp.where(kopp != kopp.rio.nodata)
kopp.rio.write_nodata(32, encoded=True, inplace=True)
kopp = kopp.rio.reproject("EPSG:4326")
kopp.name = 'KG_climate_class'
# print(f"nodata: {kopp.rio.nodata}")
# print(f"encoded_nodata: {kopp.rio.encoded_nodata}")

# %%

# Define required variables

# cont_bounds are used so that map figures are restricted to these areas (otherwise Asia and North America figures cover the entire globe)
cont_bounds = {
    'Asia': sh.geometry.box(24,0,190,81), 
    'Europe': sh.geometry.box(-31, 35, 69, 81), 
    'North America': sh.geometry.box(-178, 0, -15, 84)
    }

kopp_ind = np.unique(kopp.data[~np.isnan(kopp.data)]).astype(int)

# Color list and labels (including oceans as layer 32)
kopp_cols = np.array(["#960000", "#FF0000", "#FF6E6E", "#FFCCCC", "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64", "#007800", "#005000",
          "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00", "#550055", "#820082", "#C800C8",
          "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF",
          "#64FFFF", "#F5FFFF"])

kopp_labels = np.array(['Af', 'Am', 'As', 'Aw', 'BSh', 'BSk', 'BWh', 'BWk', 'Cfa', 'Cfb',
            'Cfc', 'Csa', 'Csb', 'Csc', 'Cwa', 'Cwb', 'Cwc', 'Dfa', 'Dfb', 'Dfc',
            'Dfd', 'Dsa', 'Dsb', 'Dsc', 'Dsd', 'Dwa', 'Dwb', 'Dwc', 'Dwd', 'EF',
            'ET', 'Ocean'])
            
# To set ocean layer 32 as missing, do:
kopp_cols = kopp_cols[0:-1]
kopp_labels = kopp_labels[0:-1]


# %%
# Plotting function: map plot for Koeppen-Geiger data

def plot_KG(rast, cols, ticks, labels, fig_size, ident):
    # cols is either a lost of colors or the name of an inbuilt matplotlib colormap
    
    var = rast.name

    if isinstance(cols, str):
        cmap = mpl.cm.get_cmap(cols)
    else:
        cmap = ListedColormap(cols)

    norm = mpl.colors.BoundaryNorm(np.append(ticks, ticks[-1]+1), ncolors=cmap.N)
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot(
        111,
        projection=ccrs.PlateCarree(),
        facecolor="white"
    )
    p = rast.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmap,
        norm=norm,
        add_colorbar=False,
        # vmin=1., # activate vmin and vmax if not using norm
        # vmax=31., 
    )
    # ax.set_global()
    ax.coastlines()
    ax.set_title('Koeppen-Geiger Regions')
    cbar = fig.colorbar(p, label=None, ticks=ticks+np.append(np.diff(ticks)/2, 0.5)) # This sets the tick positions at the middle of the colorbar boxes
    cbar = cbar.ax.set_yticklabels(labels)
    fname =  os.path.join(figs_dir, var+'_'+ident+'.png')
    plt.savefig(fname)


# %%
# Plot regions to check all looks good.
plot_KG(kopp, kopp_cols, kopp_ind, kopp_labels, [17,7], ident='global')

# %% [markdown]
# ## Regional analysis

# %%
# Prepare variables ----

# Load any of the inversion files to get area and do reproject_match
neeAmp = rio.open_rasterio(os.path.join(output_dir, 'neeAmp_CSs10.nc'))

# Create area grid for neeAmp grid cells
latcos = np.cos(np.deg2rad(neeAmp.y))
latcos.name = "weights"
area_invres = neeAmp[1].squeeze()
# Aproximate grid area at Equator
area_invres = area_invres.where(area_invres.isnull(), other=(2 * 111) * (2.5 * 111))

# Reproject kopp to inversion resolution
kopp_rp = kopp.rio.reproject_match(neeAmp)

# Create a template dataframe to save results
df_regions = pd.DataFrame({'kopp_ind': kopp_ind, 'regname': kopp_labels,
                        'regcolor': kopp_cols, 'area': np.nan, 'neeamp_var': np.nan})

# %%
# Function to get region area and variance for a continent

def reg_analysis(cont_name):
    
    # cont_name is the continent name string as found in the conts geopandas object
    cont = conts[conts['CONTINENT'] == cont_name]
    
    # Clip to the continent of interest (continent bounds make sure the mapped area only covers the continent)
    kopp_sel = kopp_rp.rio.clip(cont.geometry, cont.crs, all_touched=True).rio.clip(
        [cont_bounds[cont_name]], crs="EPSG:4326").squeeze()  # squeeze should remove dimensions of length 1.

    # Remove rid of the 'band' dimension
    del kopp_sel['band']
    
    # Make copy of DF to store results
    df_results = df_regions.copy()

    # Create list of region indices
    regind_sel = np.unique(kopp_sel.data[~np.isnan(kopp_sel.data)]).astype(int)

    plot_KG(kopp_sel, kopp_cols[(regind_sel-1)], regind_sel, kopp_labels[(regind_sel-1)], [20,7], ident=cont_name)

    # Loop over region indices and calculate statistics of nee fluxes
    for ri in regind_sel:  # [15]: #
        reg_neeamp = neeAmp.where(kopp_sel == ri)
        reg_area = area_invres.where(kopp_sel == ri).weighted(latcos)
        totarea = reg_area.sum().data
        mean_var = reg_neeamp.var(dim='year').mean().data
        df_results.loc[df_results['kopp_ind'] ==
                       ri, 'ampvar'] = np.round(mean_var, 3)
        df_results.loc[df_results['kopp_ind'] == ri, 'area'] = np.round(totarea)

    df_results = df_results.loc[~df_results['area'].isna(),:]
    df_results['ampvararea'] = np.round(
        df_results['area'] * df_results['ampvar'])
    df_results = df_results.sort_values('ampvararea', ascending=False)
    df_results['ampvararea_sum_cont'] = df_results['ampvararea'].cumsum()
    amparea_sum = df_results['ampvararea'].sum()

    print('Total amparea sum = : ', amparea_sum, '\n')
    print('80% of amparea sum = : ', amparea_sum * 0.8, '\n')

    # Explore data
    df_results.sort_values('area', ascending=False).plot.bar('regname', 'area')
    df_results.sort_values('ampvar', ascending=False).plot.bar(
        'regname', 'ampvar')
    df_results.sort_values('ampvararea', ascending=False).plot.bar(
        'regname', 'ampvararea')
    
    df_results['continent'] = cont_name

    return (df_results)


# %%

# Choose inversion data version to analyse.
inversions = ['CSsEXT', 'CSs76', 'CSs81', 'CSs85', 'CSs93', 'CSs99', 'CSs06', 'CSs10', 'CAMSsur', 'CAMSsat']

for inv in inversions:
    
    file_in_amp = os.path.join(output_dir, 'neeAmp_'+inv+'.nc')
    neeAmp = rio.open_rasterio(file_in_amp)
    ident=inv

    # Analyze each continent
    df_nee_as = reg_analysis("Asia")
    df_nee_eu = reg_analysis("Europe")
    df_nee_na = reg_analysis("North America")

    # Concatenate data frames
    df_nee_glob = pd.concat([df_nee_as, df_nee_eu, df_nee_na], axis=0)
    # Calculate the cumulative sum of ampvararea
    df_nee_glob = df_nee_glob.sort_values('ampvararea', ascending=False)
    df_nee_glob['ampvararea_sum_glob'] = df_nee_glob['ampvararea'].cumsum()

    df_nee_glob.to_csv(os.path.join(output_dir, 'KG-neeAmpStats_'+ident+'.csv'))

# %% [markdown]
# #### Recode rasters to create maps showing statistics by regions

# %%
# Function to create rasters with ampvar and area*ampvar values instead of region indices
def recodeConts(df_nee_glob):

    # Create an empty xr dataset and dictionaries to hold results
    ds_kopp_nee = xr.Dataset({'region_index': kopp_rp})
    ampvararea_out = dict()
    ampvar_out = dict()

    # Names of continents to analyse
    cont_names = ['Asia', 'Europe', 'North America']

    # Loop over continents and replace values
    for cont_name in cont_names:
        cont = conts[conts['CONTINENT'] == cont_name]

        # Clip to the selected continent
        kopp_sel = kopp_rp.rio.clip(cont.geometry, cont.crs, all_touched=True).rio.clip(
            [cont_bounds[cont_name]], crs="EPSG:4326").squeeze()
        del kopp_sel['band']

        # Create list of region indices
        regind_sel = np.unique(kopp_sel.data[~np.isnan(kopp_sel.data)]).astype(int)
        kopp_ampvararea = kopp_sel
        kopp_ampvar = kopp_sel

        # Replace region indices with their respective amp variance values (or ampvararea)
        for ri in regind_sel:
            ampvararea = df_nee_glob.loc[(df_nee_glob['continent'] == cont_name) & (df_nee_glob['kopp_ind'] == ri), 'ampvararea'].to_numpy()
            ampvar = df_nee_glob.loc[(df_nee_glob['continent'] == cont_name) & (df_nee_glob['kopp_ind'] == ri), 'ampvar'].to_numpy()
            kopp_ampvararea = xr.where(kopp_ampvararea == ri, ampvararea, kopp_ampvararea)
            kopp_ampvar = xr.where(kopp_ampvar == ri, ampvar, kopp_ampvar)
        
        # Save resulting xr arrays in dictionaries
        ampvararea_out[cont_name] = kopp_ampvararea
        ampvar_out[cont_name] = kopp_ampvar
    
    # Combine all continent xr arrays into a single array (like a mosaic)
    kopp_ampvararea = ampvararea_out['Asia'].combine_first(ampvararea_out['Europe']).combine_first(ampvararea_out['North America'])
    kopp_ampvar = ampvar_out['Asia'].combine_first(ampvar_out['Europe']).combine_first(ampvar_out['North America'])

    # Save results into the xr dataset
    ds_kopp_nee['ampvararea'] = kopp_ampvararea
    ds_kopp_nee['ampvar'] = kopp_ampvar

    return(ds_kopp_nee)

# Function to plot values
def plotGradMap(rast, cols, title, ident, vmax=None, vmin=None):
    var = rast.name
    cmap = mpl.cm.get_cmap(cols)
    fig = plt.figure(figsize=[14,6])
    ax = fig.add_subplot(
        111,
        projection=ccrs.PlateCarree(),
        facecolor="white"
    )
    p = rast.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        # add_colorbar=False,
    )
    ax.coastlines()
    ax.set_title(title+'\n'+ident)
    fname =  os.path.join(figs_dir, var+'_'+ident+'.png')
    plt.savefig(fname)


# %%
# Choose inversion data version to analyse.
inversions = ['CSsEXT', 'CSs76', 'CSs81', 'CSs85', 'CSs93', 'CSs99', 'CSs06', 'CSs10', 'CAMSsur', 'CAMSsat']

for inv in inversions:
    
    file_df_nee_glob = os.path.join(output_dir, 'KG-neeAmpStats_'+inv+'.csv')
    df_nee_glob = pd.read_csv(file_df_nee_glob)
    ident=inv
    
    ds_kopp_nee = recodeConts(df_nee_glob)

    # Make map plots of (area weighted) seasonal NEE flux amp variance
    ds_kopp_nee['ampvararea'].attrs['long_name'] = 'NEE Flux Amplitude Variance * Region Area (sqkm)'
    title = 'Regional Area-Weighted Temporal Variance in NEE Seasonal Amplitude'
    plotGradMap(ds_kopp_nee['ampvararea'], 'Oranges', title, ident=ident, vmax=400000, vmin=0)

    ds_kopp_nee['ampvar'].attrs['long_name'] = 'NEE Flux Amplitude Variance'
    title = 'Regional Temporal Variance in NEE Seasonal Amplitude'
    plotGradMap(ds_kopp_nee['ampvar'], 'Oranges', title, ident=ident, vmax=0.16, vmin=0)


