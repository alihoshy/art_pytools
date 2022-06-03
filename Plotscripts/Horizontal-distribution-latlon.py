'''
-------------
 Description:
 Python script to visualize the horizontal distribution of ICON data on lat/lon grid.
 Valid data types: NetCDF, GRIB
 There may be warnings when reading GRIB data, which can be ignored.

 31.05.2022  Anika Rohde, KIT
-------------
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib as mpl
import cartopy.feature as cfeature
import math, logging

#--  variable
var              = ''                   #-- variable name

#-- data file, grid file, plot saving path, and height info
icon_data        = ''                   #-- ICON data file
plotpath         = ''                   #-- directory where plot should be saved
plotname         = ''                   #-- plot file name

#-- plot settings
#varmin           =                     #-- Min value on colorbar (optional)
#varmax           =                     #-- Max value on colorbar (optional)
nlevs            = 25                   #-- number of levels on colorbar
plottitle        = ''                   #-- plot title
fs               = 10                   #-- fontsize
dpi              = 200                  #-- resolution of plot
n_lons           = 9                    #-- number of ticks along longitudes
n_lats           = 7                    #-- number of ticks along latitudes
cbar_col         = 'Spectral_r'         #-- colorbar
crs_data         = ccrs.PlateCarree()   #-- projection
lons             = [-180, 180]          #-- min and max lon
lats             = [-90,90]             #-- min and max lat


#-- load icon data (netcdf)
ds_icon                    = xr.open_dataset(icon_data)                      #-- open dataset as netcdf

#-- alternatively (grib)
#logging.disable(logging.ERROR)
#ds_icon                    = xr.open_dataset(icon_data, engine="cfgrib", indexpath='')  #-- open dataset as grib
#logging.disable(logging.NOTSET)

#-- select time step (adjust if necessary)
if 'time' in ds_icon.dims:                        #-- check if 'time' is a dimension in the dataset
  ds_icon = ds_icon.isel(time=0)                  #-- if so, chose first time
if 'step' in ds_icon.dims:                        #-- check if 'step' is a dimension in the dataset
  ds_icon = ds_icon.isel(step=0)                  #-- if so, chose first time step


data = ds_icon[var]                               #-- select data of variable out of the dataset
unit = data.units                                 #-- get unit info
dims = data.dims                                  #-- get the dimensions

#-- 3D variables:
if len(dims) == 3:
  #ilev = 0
  #data = data[ilev,:]                            #-- select level ilev

  # or process data: example - sum over all levels
  data = data.sum(dims[0])


#-- Colorbar for plot

#-- set min and max if not defined
if not 'varmin' in locals():
  varmin = np.nanmin(data)
if not 'varmax' in locals():
  varmax = np.nanmax(data)

#-- check if all data is in the range of varmin and var max -> adjust look of colorbar
extendarrow = 'neither'
if (varmax < np.nanmax(data)) & (varmin > np.nanmin(data)):
  extendarrow = 'both'
elif varmax < np.nanmax(data):
  extendarrow = 'max'
elif varmin > np.nanmin(data):
  extendarrow = 'min'


collev       = np.linspace(varmin,varmax, nlevs)               #-- set nlevs segments on colorbar
levticks     = collev[0::4]                                    #-- set colorbar labels (every 4th)
cmap         = plt.get_cmap(cbar_col, len(collev)-1)
colors       = np.zeros([len(data),4],np.float32)              #-- assign color array for triangles
ncol         = cmap.N                                          #-- number of colors
norm         = mpl.colors.BoundaryNorm(collev, ncol)           #-- normalize for color bar

extendarrow = 'neither'
if (varmax < np.nanmax(data)) & (varmin > np.nanmin(data)):
  extendarrow = 'both'
elif varmax < np.nanmax(data):
  extendarrow = 'max'
elif varmin > np.nanmin(data):
  extendarrow = 'min'


#-- start figure
fig = plt.figure()
ax  = fig.add_axes([0.09,0.0,0.71,0.85], projection=crs_data)                    #-- margin of the map (values: left corner x, y + width and height of map) + projection
ax.set_extent([np.min(lons), np.max(lons), np.min(lats), np.max(lats)], crs=crs_data)    # -- set extend of the map (here min and max lon/lat of data)

ax.coastlines(resolution='110m', lw=0.7)                                         #-- coastlines
#ax.add_feature(cfeature.BORDERS, linestyle=':', lw=1, edgecolor='white')        #-- country boarders in white (for dark colors)
ax.add_feature(cfeature.BORDERS, lw=0.3)                                         #-- country boarders in black (for bright colors)
im =data.plot(norm=norm, cmap=cmap, add_colorbar=False)                          #-- plot data

#-- add lat lon ticks
lon_ticks   = np.linspace(np.min(lons), np.max(lons), n_lons)                    #-- create lon ticks array
lat_ticks   = np.linspace(np.min(lats), np.max(lats), n_lats)                    #-- create lat ticks array
lon_labels  = np.array(['{0:.{1}f}'.format(x,1)+'$^\circ$' for x in lon_ticks])  #-- create string array, round to one decimal (lon labels)
lat_labels  = np.array(['{0:.{1}f}'.format(x,1)+'$^\circ$' for x in lat_ticks])  #-- create string array, round to one decimal (lat labels)
lon_labels[range(0,n_lons,2)] = ''                                               #-- empty every second value (cleaner look)
lat_labels[range(0,n_lats,2)] = ''                                               #-- empty every second value (cleaner look)
ax.set_xticks(lon_ticks, crs=crs_data)                                           #-- add longitudes ticks
ax.set_yticks(lat_ticks, crs=crs_data)                                           #-- add latitudes ticks
ax.set_xticklabels(lon_labels, fontsize=fs)                                      #-- add lon tick labels
ax.set_yticklabels(lat_labels, fontsize=fs)                                      #-- add lat tick labels
ax.tick_params(direction='out', right=True, top=True, pad=4)                     #-- ticks inside figure and space (pad) between numbers and axis
ax.set_xlabel('Longitude', fontsize=fs)                                          #-- add x-axis label 
ax.set_ylabel('Latitude', fontsize=fs)                                           #-- add y-axis label
ax.set_title(plottitle, pad=10)                                                  #-- set title



#-- add colorbar
ax   = fig.add_axes([0.83, 0.17,0.025,0.51])                                     #-- margin of colorbar (values: left corner x, y + width and height)
cbar = plt.colorbar(mappable=im, cax=ax,  extend=extendarrow) 

cbar.ax.tick_params(labelsize=fs)
cbar.set_label(unit, fontsize=fs)                                                #-- label colorbar

plt.savefig(plotpath + '/' + plotname, dpi=dpi, bbox_inches='tight')             #-- save figure

