# %% [markdown]
# ## LCC analysis
# 
# ### Methods:
# 
# 1. Load geometries, co2flux amp trends, and LCC data and preprocess as necessary.
# 2. Clip to ROI (region of interest). Use a small region for testing puposes.
# 3. Reduce land cover data to the resolution of the co2flux data, creating new images with bands representing the fraction of cover per pixel per land type
# 4. Calculate 

# %%
import rioxarray as rio
import xarray as xr
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import folium

# %%
import ee
ee.Initialize()

co2amp = ee.Image("projects/augs-geo-recovery/assets/CO2InvSeasAmpTrend_CarboScope_s85ocv2022")
modis_lc = ee.ImageCollection('MODIS/006/MCD12Q1')
# Initial date of interest (inclusive).
i_date = '2001-01-01'
# Final date of interest (exclusive).
f_date = '2022-01-01'
# Selection of appropriate bands and dates for LST.
igbp_lc = modis_lc.select('LC_Type1').filterDate(i_date, f_date)

# %%

lcIndex = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
lcNames = ['ENForest', 'EBForest', 'DNForest', 'DBForest', 
    'MixForest', 'ClosedShrub', 'OpenShrub', 'WoodySavanna',
    'Savanna', 'Grassland', 'PermWetland', 'Cropland',
    'Urban', 'CropNatMosiac', 'PermSnowIce', 'Barren']
def masklc(img):
    tmp = img.select('LC_Type1').eq(ind).rename(name)
    return img.addBands(tmp)

for i in range(len(lcIndex)):
    ind = ee.Number(lcIndex[i])
    name = ee.String(lcNames[i])
    igbp_lc = igbp_lc.map(masklc)

# %%
# Code for creating an interactive map

def add_ee_layer(self, ee_image_object, vis_params, name):
    """Adds a method for displaying Earth Engine image tiles to folium map."""
    map_id_dict = ee.Image(ee_image_object).getMapId(vis_params)
    folium.raster_layers.TileLayer(
        tiles=map_id_dict['tile_fetcher'].url_format,
        attr='Map Data &copy; <a href="https://earthengine.google.com/">Google Earth Engine</a>',
        name=name,
        overlay=True,
        control=True
    ).add_to(self)

# Add Earth Engine drawing method to folium.
folium.Map.add_ee_layer = add_ee_layer

# Set visualization parameters for land cover.
lc_vis_params = {
    'min': 0,'max': 1,
    'palette': ['000000', '05450a']
}

lc = 'MixForest'
igbp_enf_2001 = igbp_lc.first()
lat, lon = 45.77, 4.855
my_map = folium.Map(location=[lat, lon], zoom_start=2)
mask = igbp_enf_2001.select(lc)
plotimg = igbp_enf_2001.select(lc).updateMask(mask)
my_map.add_ee_layer(plotimg, lc_vis_params, 'land cover')
my_map.add_child(folium.LayerControl())
display(my_map)


