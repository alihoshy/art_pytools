'''
-------------
 Description:
 Python script to visualize a cross section along a defined xy path.
 This script focuses on ICON data on native grid.
 A path xyz can be vizalized as well as external data on xyz.   

 Valid data types: NetCDF, GRIB
 There may be warnings when reading GRIB data, which can be ignored.
 
 28.02.2022  Anika Rohde, KIT
-------------
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr
import logging 



#--  variable, data file and grid 
var              = ''               #-- variable name that should be plotted
icon_data        = ''               #-- ICON data file
icon_grid        = ''               #-- respective grid file
plotpath         = ''               #-- directory where plot should be saved
#-- icon_data must include z_ifc - otherwise load z_ifc from a different file here:
#z_ifc_data       = ''
#-- ICON topography from extpar (optional) 
#extpar           = ''
#-- if icon data includes different time steps
i_step           = 0


#-- plot settings
#varMin           =                  #-- Min value on colorbar (optional)
#varMax           =                  #-- Max value on colorbar (optional) 
alt_plot         = [0.,15.]          #-- Min and Max altitude to be plotted (km)
lat_bounds       = [50, 70]          #-- Min and Max latitude to be plotted
nlevs            = 40                #-- number of levels on colorbar
icon_plottitle   = ''                #-- plot title
icon_plotname    = 'ICON.png'        #-- plot file name
fs               = 10                #-- fontsize
dpi              = 400               #-- resolution of plot
n_ticks          = 6                 #-- number of latitude/longitude ticks in plot
cbar_col         = 'Spectral_r'      #-- colorbar


#-- x/y cross section
y_lats           =                   #-- insert lat array here 
x_lons           =                   #-- insert lon array here

#-- plot xyz path or plot external data along xyz (optional) 
plot_xyz         = False                              #-- add flight path xyz -> set True
z_height         =                                    #-- insert z array in km here 
#data_xyz         =                                   #-- insert external data array here (optional)



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
  ds_icon      = xr.open_dataset(icon_data)                                               #-- open icon dataset
except ValueError:                                                                        #-- except not readable as netcdf
  logging.disable(logging.ERROR)                                                          #-- suppress error messages
  ds_icon      = xr.open_dataset(icon_data, engine="cfgrib", indexpath='')                #-- open dataset as grib
  logging.disable(logging.NOTSET)                                                         #-- enable error messages

#-- get z_ifc info
try:
  z_ifc        = np.array(ds_icon['z_ifc'])                                                        #-- check if 'z_ifc' is included in the dataset
except  KeyError:
  try:
    z_ifc        = np.array(xr.open_dataset(z_ifc_data)['z_ifc'])                                  #-- otherwise read in the additional data containing z_ifc
  except ValueError:
    logging.disable(logging.ERROR)                                                                 #-- suppress error messages
    z_ifc        = np.array(xr.open_dataset(z_ifc_data, engine="cfgrib", indexpath='')['z_ifc'])   #-- open dataset as grib
    logging.disable(logging.NOTSET)                                                                #-- enable error messages

if 'time' in ds_icon.coords:                                                                       #-- check if 'time' is a dimension in the dataset 
  ds_icon = ds_icon.isel(time=0)                                                                   #-- if so, remove the dimension
if 'step' in ds_icon.coords:                                                                       #-- check if 'step' is a dimension in the dataset 
  ds_icon = ds_icon.isel(step=i_step)                                                              #-- if so, chose time step



ds_icon   = ds_icon[var]                                                              #-- select variable from dataset
lev_num   = ds_icon.coords[ds_icon.dims[0]].values                                    #-- icon level numbers 
data      = ds_icon.values                                                            #-- extract values
lat_icon  = np.array(xr.open_dataset(icon_grid)["clat"]) * 180.0/np.pi                #-- get latidues from gridfile
lon_icon  = np.array(xr.open_dataset(icon_grid)["clon"]) * 180.0/np.pi                #-- get longitudes from gridfile
try:
  topo      = xr.open_dataset(extpar)[['topography_c']].topography_c.values/1000      #-- read in topography in km
except:
  print("Note: no extpar file provided for topography.")

if z_ifc.shape[0] == len(lev_num)+1:    
  z_ifc   = z_ifc[1:,]                                                                #-- if variable has less levels than z_ifc
z_ifc = z_ifc/1000                                                                    #-- convert z_ifc from m to km




#-- select data along xy

#-- prepare empty arrays
icon_data_on_path    = np.zeros([len(lev_num), len(y_lats)])                            #-- empty array with shape (levels, datapoints on path xy)
icon_lev_on_path     = np.zeros([len(lev_num), len(y_lats)])                            #-- empty array with shape (levels, datapoints on path xy)
icon_data_xyz        = np.zeros(len(y_lats))                                            #-- empty array with shape (x,y,z) (to be filled with icon data)
icon_height          = np.zeros(len(y_lats))
topo_icon            = np.zeros(len(y_lats)) 

#-- loop over all lon & lat points on path xy to find closest icon cell
for i in range(len(y_lats)):
    distances                  = distance([x_lons[i], y_lats[i]], [lon_icon, lat_icon]) #-- compute distances of all icon cells to individual point xy
    icon_loc_index             = np.argmin(distances)                                   #-- chose icon cell index closest to individual point xy
    icon_data_on_path[:, i]    = data[:, icon_loc_index]                                #-- fill icon data (all levels) closest to xy
    icon_lev_on_path[:, i]     = z_ifc[:, icon_loc_index]                               #-- fill icon level heights (all) closest to xy
    if 'topo' in locals():
      topo_icon[i]             = topo[icon_loc_index]                                    

#-- loop over all lon & lat points on path xy to find closest level to z
for j in range(icon_data_on_path.shape[1]):
    icon_lev_ind     = np.argmin(abs(icon_lev_on_path[:, j] - z_height[j]))             #-- find index of smallest difference between z_ifc and z
    icon_data_xyz[j] = icon_data_on_path[icon_lev_ind, j]                               #-- fill icon data on path xyz
    icon_height[j]   = icon_lev_on_path[icon_lev_ind, j]                                #-- get respective z_ifc of nearest icon level
   
#-- arange data for plotting
plot_data = icon_data_on_path                             
plot_lats = np.tile(y_lats, (len(lev_num),1))
 
#-- plot settings
cmap        = plt.get_cmap(cbar_col, nlevs-1)                                            #-- read the color map
units       = ds_icon.attrs['units']                                                     #-- get units info

if not 'varMin' in locals():                                                             #-- check if Min and Max values for colorbar are assigned
  varMin    = np.nanmin(plot_data)                                                       #-- otherwise set Min to Min of the selcted data
  varMax    = np.nanmax(plot_data)                                                       #-- and set Max to Max of the selcted data

levels      = np.linspace(varMin, varMax, nlevs)                                         #-- create levels on colorbar
ticklevels  = np.linspace(varMin, varMax, int(nlevs/4))                                  #-- create ticks on colorbar
norm        = mpl.colors.BoundaryNorm(levels, cmap.N)                                    #-- normalize for colorbar

#-- adjust xticks and xlabels
lab_xticks_lat = np.round(np.linspace(y_lats[0], y_lats[y_lats.size-1], n_ticks),2)      #-- create lat-ticks on x axis, (amount = n_ticks) 
lab_xticks_lon = np.round(np.linspace(x_lons[0], x_lons[x_lons.size-1], n_ticks),2)      #-- create lon-ticks on x axis, (amount = n_ticks)
lat_labels     = ['{0:.{1}f}'.format(x,2) for x in lab_xticks_lat]                       #-- round, 2 decimals (lat labels) and convert to string
lon_labels     = ['{0:.{1}f}'.format(x,2) for x in lab_xticks_lon]                       #-- round, 2 decimals (lon labels) and convert to string
lat_labels[0]  = 'Lat   '+ lat_labels[0]                                                 #-- add 'Lat' at first tick
lon_labels[0]  = 'Lon   ' + lon_labels[0]                                                #-- add 'Lon' at first tick
pos_xticks     = lab_xticks_lat                                                          #-- position of ticks on the x-axis (ax1)




#-- create figure
pfig    = plt.figure()
ax1     = pfig.add_subplot(1,1,1)                                                        #-- first x and y-axis (latitudes)
ax2     = ax1.twiny()                                                                    #-- second x-axis (longitudes)

c = ax1.contourf(plot_lats,icon_lev_on_path, plot_data, levels=levels,cmap=cmap,norm=norm, extend='both')        #-- plot data
if 'topo' in locals():
  ax1.plot(y_lats, topo_icon, color='red')                                               #-- draw surface

#-- plot external or xyz path with icon data
if plot_xyz:                                                                             #-- if xyz is shoud be drawn 
  if 'data_xyz' in locals():                                                             #-- if provided, plot external data  
    ax1.scatter(y_lats, z_height, c=data_xyz,  cmap=cmap,norm=norm, s=10, edgecolor='black', linewidth=0.2)
  else:                                                                                  #-- otherwise mark icon data on xyz
    ax1.scatter(y_lats, z_height, c=icon_data_xyz,  cmap=cmap,norm=norm, s=10, edgecolor='black', linewidth=0.2) 

#-- Latitudes on x-axis
ax1.set_xlim(lat_bounds)                                                                 #-- set map limits (Latitudes)
ax1.set_xticks(pos_xticks)                                                               #-- set xticks
ax1.set_xticklabels(lat_labels, fontsize=fs, ha='right')                                 #-- add latitude labels

#-- Longitudes on x-axis
x_axis = ax2.get_xaxis()
ax2.set_xticks(pos_xticks)                                                               #-- set xticks  
ax2.set_xlim(lat_bounds)
x_axis.set_visible(True)                                                                 #-- show second x-axis (longitudes)
ax2.set_xticklabels(lon_labels, fontsize=fs, ha='right')                                 #-- add longitude labels
ax2.xaxis.set_ticks_position("bottom")                                                   #-- show ticks of second x-axis on the bottom
ax2.xaxis.set_label_position("bottom")                                                   #-- show labels of second x-axis (longitudes) on the bottom
ax2.tick_params(axis='x', which='major', pad=13)                                         #-- offset second x-axis

#-- adjust yticks and ylabels (altitude)
pos_yticks = np.arange(alt_plot[0], int(alt_plot[1])+5, 5).astype(int)                   #-- arange position of altitude yticks
ax1.set_yticks(pos_yticks)                                                               #-- set a altitude yticks
ax1.set_yticklabels(np.arange(alt_plot[0], int(alt_plot[1])+5, 5), fontsize=fs)          #-- add altitude labels 
ax1.set_ylabel('Altitude (km)', fontsize=fs)                                             #-- add a-axis label 
ax1.set_ylim(alt_plot)                                                                   #-- set limits of y-axis

#-- title
ax1.set_title(icon_plottitle, pad=10, fontsize=fs+fs/3)                                  #-- add plot title

#-- colorbar
cbar = pfig.colorbar(c, ticks=ticklevels, pad=0.01)                                      #-- add a colorbar
cbar.ax.tick_params(labelsize=fs)                                                        #-- set fontsize of ticklabels on colorbar   
cbar_label = units                                                              
cbar.set_label(cbar_label, fontsize=fs)                                                  #-- add units as colorbar label

plt.savefig(plotpath + '/' + icon_plotname, dpi=dpi, bbox_inches='tight')                #-- save figure



