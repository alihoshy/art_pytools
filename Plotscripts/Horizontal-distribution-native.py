'''
-------------
 Description: 
 Python script to visualize the horizontal distribution of ICON data on the native grid.
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
from   matplotlib.collections import PolyCollection
import cartopy.feature as cfeature
import math, logging 

#--  variable
var              = ''                   #-- variable name

#-- data file, grid file, plot saving path, and height info
icon_data        = ''                   #-- ICON data file
icon_grid        = ''                   #-- respective grid file
plotpath         = ''                   #-- directory where plot should be saved
plotname         = ''                   #-- plot file name

#-- plot settings
#varmin           =                     #-- Min value on colorbar (optional)
#varmax           =                     #-- Max value on colorbar (optional)
nlevs            = 25                   #-- number of levels on colorbar
plottitle        = ''                   #-- plot title
fs               = 10                   #-- fontsize
dpi              = 200                  #-- resolution of plot
n_ticks          = 6                    #-- number of latitude/longitude ticks in plot
cbar_col         = 'Spectral_r'         #-- colorbar
crs_data         = ccrs.PlateCarree()   #-- projection


#-- load grid information
grid                       = xr.open_dataset(icon_grid)[['clon_vertices','clat_vertices']].rename({'cell': 'ncells'})
rad2deg                    = 180.0/np.pi                                     #-- factor to calculate degeree from radian 
vlon                       = grid['clon_vertices']*rad2deg                   #-- convert grid icon info to degrees (lats/lons of triangle vertices)
vlat                       = grid['clat_vertices']*rad2deg 
ncells, nv                 = vlon.shape[0], vlon.shape[1]                    #-- number of cells and edges
vlon, vlat                 = np.array(vlon), np.array(vlat)                  #-- convert to numpy arrays


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
print(ds_icon.dims)

data = ds_icon[var]                               #-- select data of variable out of the dataset
unit = data.units                                 #-- get unit info
dims = data.dims                                  #-- get the dimensions

#-- 3D variables:
if len(dims) == 2:        
  #ilev = 0                       
  #data = data[ilev,:].values                     #-- select level ilev

  # or process data: example - sum over all levels
  data = data.sum(dims[0]).values      

#-- 2D variables:
elif len(dims) == 1:
  data = data.values



#-- workaround (approximation) for global data
#-- prevents distortions on the map caused by triangles crossing -180 / 180 longitude
#-- solution: cut triangles at 180/-180 and copy the other part to other side
tri_idx         = np.where(np.any(vlon <-100, axis=1) & np.any(vlon >100, axis=1))[0]     #-- find triangles with vertices on both sides

if len(tri_idx >0):
  copy_lon = vlon[tri_idx]                                    #-- copy lons of vertices crossing -180/180

  for i in range(len(tri_idx)):
    copy_lon[i]     [vlon[tri_idx[i]]<0]  = 180               #-- recreate triangles on the other side, set negative lons to 180
    vlon[tri_idx[i]][vlon[tri_idx[i]]>0]  = -180              #-- modify/cut original triangles at -180 (set positive lons to -180)

  data = np.append(data,  data[tri_idx])                      #-- add data of additional triangles (not modified - copy)
  vlat = np.vstack((vlat, vlat[tri_idx]))                     #-- add lats of additional triangles (not modified - copy)
  vlon = np.vstack((vlon, copy_lon))                          #-- add lons of additional triangles (modified)


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
    
#-- set color index of all cells in between levels
for m in range(0,nlevs-1):
  colors[np.where((data    >= collev[m]) & (data < collev[m+1])),:] = cmap(m)

colors[np.where(data < varmin),:]    = cmap(0)                 #-- set color below levels
colors[np.where(data >= varmax),:]   = cmap(ncol-1)            #-- set color index for cells greater than levels[nlevs-1]
 
#-- create list of vertices
triangles = [np.column_stack((vlon[i,:],vlat[i,:])) for i in range(0,len(vlon))]
coll = PolyCollection(triangles, facecolor=colors,edgecolors="none", closed=True, transform=crs_data)



#-- start figure
fig = plt.figure()
ax  = fig.add_axes([0.09,0.0,0.71,0.85], projection=crs_data)                    #-- margin of the map (values: left corner x, y + width and height of map) + projection
ax.set_extent([np.min(vlon), np.max(vlon), np.min(vlat), np.max(vlat)], crs=crs_data)    # -- set extend of the map (here min and max lon/lat of data)

ax.coastlines(resolution='110m', lw=0.7)                                         #-- coastlines
#ax.add_feature(cfeature.BORDERS, linestyle=':', lw=1, edgecolor='white')        #-- country boarders in white (for dark colors)
ax.add_feature(cfeature.BORDERS, lw=0.3)                                         #-- country boarders in black (for bright colors)
ax.add_collection(coll)                                                          #-- add triangles (data) 

#-- add lat lon ticks
n_lons = 9                                                                       #-- number of ticks along longitudes 
n_lats = 7                                                                       #-- number of ticks along latitudes
lon_ticks   = np.linspace(np.min(vlon), np.max(vlon), n_lons)                    #-- create lon ticks array 
lat_ticks   = np.linspace(np.min(vlat), np.max(vlat), n_lats)                    #-- create lat ticks array
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
ax  = fig.add_axes([0.83330.17,0.025,0.51])                                      #-- margin of colorbar (values: left corner x, y + width and height) 
cbar = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, spacing='proportional',   
                              ticks=levticks,  orientation='vertical', extend=extendarrow)

cbar.ax.tick_params(labelsize=fs)
cbar.set_label(unit, fontsize=fs)                                                #-- label colorbar   

plt.savefig(plotpath + '/' + plotname, dpi=dpi, bbox_inches='tight')             #-- save figure






