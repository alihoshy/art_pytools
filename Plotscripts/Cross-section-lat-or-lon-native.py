'''
-------------
 Description:
 Python script to visualize a cross section along a latitude or longitude.

 Valid data types: NetCDF, GRIB
 There may be warnings when reading GRIB data, which can be ignored.
 
 01.04.2022  Anika Rohde, KIT
-------------
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr
import logging 


#--  variable, and location
var              = ''                    #-- variable to be plotted               
#cross_lat        = 48                   #-- Latitude  
#lon_bounds       = [50, 70]             #-- Min and Max latitude to be plotted (optional)
                                         #-- or
cross_lon        = 8.4                   #-- Longitude      
#lat_bounds       = [0, 90]               #-- Min and Max latitude to be plotted (optional)

#-- data file, grid file, plot saving path, and height info
icon_data        = ''                    #-- ICON data file
icon_grid        = ''                    #-- respective grid file
plotpath         = ''                    #-- directory where plot should be saved
plotname         = 'Cross-section.png'   #-- plot file name

alt_var          = ''                    #-- height info var e.g. 'z_ifc' or 'HHL'
#-- icon_data shoud include altitude info like z_ifc or HHL, otherwise levels are shown on y-axis
#-- you can load altitude info from a different file:
#alt_data         = ''

#-- ICON topography from extpar (optional) 
#extpar           = '/work/bb1070/b380982/Ali/icon_extpar_0051_R02B07_N02_20181011_tiles.nc'
#-- if icon data includes different time steps
#i_step           = 0



#-- plot settings
#varMin           =                  #-- Min value on colorbar (optional)
#varMax           =                  #-- Max value on colorbar (optional) 
#alt_plot         = [0.,20.]          #-- Min and Max altitude to be plotted (km) (optional)
nlevs            = 40                #-- number of levels on colorbar
plottitle        = ''                #-- plot title
fs               = 10                #-- fontsize
dpi              = 400               #-- resolution of plot
n_ticks          = 6                 #-- number of latitude/longitude ticks in plot
cbar_col         = 'Spectral_r'      #-- colorbar


def distance(origin, destination):
    """
    This function computes the distance of a location (lon & lat, here: origin)
    to another location (destination) or a set of lats & lons.
    Unit: km.
    """
    lon1, lat1 = origin
    lon2, lat2 = destination
    radius = 6371 # km

    dlat = np.radians(lat2-lat1)
    dlon = np.radians(lon2-lon1)
    a = np.sin(dlat/2) * np.sin(dlat/2) + np.cos(np.radians(lat1)) \
      * np.cos(np.radians(lat2)) * np.sin(dlon/2) * np.sin(dlon/2)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    d = radius * c

    return d



#-- Load ICON data
try:
  ds_icon      = xr.open_dataset(icon_data)                                           #-- open icon dataset
except ValueError:                                                                    #-- except not readable as netcdf
  logging.disable(logging.ERROR)                                                      #-- suppress error messages
  ds_icon      = xr.open_dataset(icon_data, engine="cfgrib", indexpath='')            #-- open dataset as grib
  logging.disable(logging.NOTSET)                                                     #-- enable error messages


if 'time' in ds_icon.coords:                                                          #-- check if 'time' is a dimension in the dataset 
  ds_icon = ds_icon.isel(time=0)                                                      #-- if so, remove the dimension
if 'step' in ds_icon.coords:                                                          #-- check if 'step' is a dimension in the dataset 
  ds_icon = ds_icon.isel(step=i_step)                                                 #-- if so, chose time step

data      = ds_icon[var]
units     = data.attrs['units']                                                       #-- get units info
lev_num   = data.coords[data.dims[0]].values                                          #-- icon level numbers 
data      = data.values                                                               #-- extract values
lat_icon  = np.array(xr.open_dataset(icon_grid)["clat"]) * 180.0/np.pi                #-- get latidues from gridfile
lon_icon  = np.array(xr.open_dataset(icon_grid)["clon"]) * 180.0/np.pi                #-- get longitudes from gridfile
try:
  topo      = xr.open_dataset(extpar)[['topography_c']].topography_c.values/1000      #-- read in topography in km
except:
  print("Note: no extpar file provided for topography.")


#-- get height info
try:
  altitudes    = np.array(ds_icon[alt_var])                                           #-- check if altitude variable' is included in the dataset
  heightinfo   = True                                                                 #-- indicates height info is available
except  KeyError:
  if 'alt_data' in locals():                                                          #-- check if a supplementary file with height info was specified
    try:
      altitudes    = xr.open_dataset(alt_data)                                        #-- read in the supplementary data containing height info
    except ValueError:
      logging.disable(logging.ERROR)                                                  #-- suppress error messages
      altitudes    = xr.open_dataset(alt_data, engine="cfgrib", indexpath='')         #-- open dataset as grib
      logging.disable(logging.NOTSET)                                                 #-- enable error messages
    heightinfo   = True                                                               #-- indicates height info is available
    if 'time' in altitudes.coords:                                                    #-- check if supplementary data has a time dimension
      altitudes = altitudes.isel(time=0)                                              #-- remove time dimension
    altitudes   = altitudes[alt_var].values                                           #-- get values of height info
  else:
    print('  Warning: altitude info is missing!')                                     #-- if no height info is available
    altitudes  = np.arange(len(data),0, -1)                                           #-- enunmerate layers instead
    heightinfo = None                                                                 #-- indicates height info is NOT available

if not heightinfo == None:
  if altitudes.shape[0] == len(data)+1:
    altitudes   = altitudes[1:,]                                                      #-- if variable has less levels than height info
  altitudes = altitudes/1000                                                          #-- convert m to km


#-- auto assign alt_plot, lat_bounds and lon_bounds if not assigned
if not 'alt_plot' in locals():
  if not heightinfo == None:
    alt_plot = [np.min(altitudes), np.max(altitudes)]
if not 'lat_bounds' in locals():
  lat_bounds = [np.min(lat_icon), np.max(lat_icon)]
if not 'lon_bounds' in locals():
  lon_bounds = [np.min(lon_icon), np.max(lon_icon)]

#-- creat lat and lon array describing a path along entered latitude or longitude
if 'cross_lat' in locals():
  y_lats = np.repeat(cross_lat, 2000)
  x_lons = np.linspace(max([-180, min(lon_bounds)]),min([180, max(lon_bounds)]), 2000)
elif 'cross_lon' in locals():
  y_lats = np.linspace(min([90, max(lat_bounds)]), max(-90, min(lat_bounds)), 2000)
  x_lons = np.repeat(cross_lon, 2000)


#-- prepare empty arrays
icon_data_on_path    = np.zeros([len(lev_num), len(y_lats)])                            #-- empty array with shape (levels, datapoints on path xy)
icon_lev_on_path     = np.zeros([len(lev_num), len(y_lats)])                            #-- empty array with shape (levels, datapoints on path xy)
icon_data_xyz        = np.zeros(len(y_lats))                                            #-- empty array with shape (datapoints on path xyz) (to be filled with icon data)
topo_icon            = np.zeros(len(y_lats))                                            #-- empty array with shape (datapoints on path xy) (to be filled with topography info, optional) 

#-- loop over all lon & lat points on path xy to find closest icon cell
for i in range(len(y_lats)):
    distances                  = distance([x_lons[i], y_lats[i]], [lon_icon, lat_icon]) #-- compute distances of all icon cells to individual point xy
    icon_loc_index             = np.argmin(distances)                                   #-- chose icon cell index closest to individual point xy
    icon_data_on_path[:, i]    = data[:, icon_loc_index]                                #-- fill icon data (all levels) closest to xy
    icon_lev_on_path[:, i]     = altitudes[:, icon_loc_index]                           #-- fill icon level heights (all) closest to xy
    if 'topo' in locals():
      topo_icon[i]             = topo[icon_loc_index]                                   #-- fill topography info      

#-- arange data for plotting
plot_data = icon_data_on_path                                                         

if 'cross_lat' in locals():
  plot_x  = np.tile(x_lons, (len(lev_num),1))
elif 'cross_lon' in locals():
   plot_x  = np.tile(y_lats, (len(lev_num),1))
 
#-- plot settings
cmap        = plt.get_cmap(cbar_col, nlevs-1)                                            #-- read the color map

if not 'varMin' in locals():                                                             #-- check if Min and Max values for colorbar are assigned
  varMin    = np.nanmin(plot_data)                                                       #-- otherwise set Min to Min of the selcted data
  varMax    = np.nanmax(plot_data)                                                       #-- and set Max to Max of the selcted data

levels      = np.linspace(varMin, varMax, nlevs)                                         #-- create levels on colorbar
ticklevels  = np.linspace(varMin, varMax, int(nlevs/4))                                  #-- create ticks on colorbar
norm        = mpl.colors.BoundaryNorm(levels, cmap.N)                                    #-- normalize for colorbar


#-- create figure
pfig    = plt.figure()
ax1     = pfig.add_subplot(1,1,1)                                                        #-- first x and y-axis (latitudes)

c = ax1.contourf(plot_x,icon_lev_on_path, plot_data, levels=levels,cmap=cmap,norm=norm, extend='both')        #-- plot data


#-- Longitudes on x-axis
if 'cross_lat' in locals():
  lab_xticks_lon = np.round(np.linspace(x_lons[0], x_lons[x_lons.size-1], n_ticks),2)    #-- create lon-ticks on x axis, (amount = n_ticks)
  lon_labels     = ['{0:.{1}f}'.format(x,2) for x in lab_xticks_lon]                     #-- round, 2 decimals (lon labels) and convert to string
  ax1.set_xlim(lon_bounds)
  ax1.set_xticks(lab_xticks_lon)
  ax1.set_xticklabels(lon_labels, fontsize=fs)
  ax1.set_xlabel('Longitudes ($^\circ$)',  fontsize=fs)
  if 'topo' in locals():
    ax1.plot(x_lons, topo_icon, color='red',lw=0.7)                                      #-- draw surface

#-- Latitudes on x-axis
elif 'cross_lon' in locals():
  lab_xticks_lat = np.round(np.linspace(y_lats[0], y_lats[y_lats.size-1], n_ticks),2)    #-- create lat-ticks on x axis, (amount = n_ticks) 
  lat_labels     = ['{0:.{1}f}'.format(x,2) for x in lab_xticks_lat]                     #-- round, 2 decimals (lat labels) and convert to string
  ax1.set_xlim(lat_bounds)                                                               #-- set map limits (Latitudes)
  ax1.set_xticks(lab_xticks_lat)                                                         #-- set xticks
  ax1.set_xticklabels(lat_labels, fontsize=fs)                                           #-- add latitude labels
  ax1.set_xlabel('Latitudes ($^\circ$)',  fontsize=fs)
  if 'topo' in locals():
    ax1.plot(y_lats, topo_icon, color='red',lw=0.7)                                      #-- draw surface


#-- adjust yticks and ylabels (altitude)
if not heightinfo == None:
  yticks = np.arange(alt_plot[0], alt_plot[1]+5, 5)                                      #-- arange position of altitude yticks
  ax1.set_yticks(yticks)                                                                 #-- set a altitude yticks
  ytick_lab =  ['{0:.{1}f}'.format(x,1) for x in yticks]
  ax1.set_yticklabels(ytick_lab, fontsize=fs)          #-- add altitude labels 
  ax1.set_ylabel('Altitude (km)', fontsize=fs)                                           #-- add a-axis label 
  ax1.set_ylim(alt_plot)                                                                 #-- set limits of y-axis

#-- title
ax1.set_title(plottitle, pad=10, fontsize=fs+fs/3)                                       #-- add plot title

#-- colorbar
cbar = pfig.colorbar(c, ticks=ticklevels, pad=0.01)                                      #-- add a colorbar
cbar.ax.tick_params(labelsize=fs)                                                        #-- set fontsize of ticklabels on colorbar   
cbar_label = units                                                              
cbar.set_label(cbar_label, fontsize=fs)                                                  #-- add units as colorbar label

plt.savefig(plotpath + '/' + plotname, dpi=dpi, bbox_inches='tight')                     #-- save figure



