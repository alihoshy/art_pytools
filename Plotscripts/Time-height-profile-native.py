'''
-------------
 Description:
 Python script to create a time height profile at a lat/lon position.

 Valid data types: NetCDF, GRIB
 There may be warnings when reading GRIB data, which can be ignored.
 
 01.04.2022  Anika Rohde, KIT
-------------
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr
import logging, os 
from datetime import datetime


#--  variable, and location 
var              = ''                          #-- variable to be plotted
lat              = 49.0068901                  #-- latitude
lon              = 8.4036527                   #-- longitude


#-- data file, grid file, plot saving path, and height info
folder           = ''                                                                                 #-- path to folder containing data files
prefix           = ''                                                                                 #-- prefix of data to sort out other data in folder
icon_grid        = ''                                                                                 #-- respective grid file
plotpath         = ''                                                                                 #-- directory where plot should be saved
plotname         = 'Profile.png'                                                                      #-- plot file name

alt_var          = 'z_ifc'                                                                            #-- height info var e.g. 'z_ifc' or 'HHL' 
#-- icon_data shoud include altitude info like z_ifc or HHL, otherwise levels are shown on y-axis 
#-- you can load altitude info from a different file:
alt_data         = ''                                                   


#-- plot settings
#varMin           =                  #-- Min value on colorbar (optional)
#varMax           =                  #-- Max value on colorbar (optional) 
plotincr         = 6                 #-- ticks in plot ('hourly' or 'daily' or number in h)
#alt_plot         = [0.,20.]          #-- Min and Max altitude to be plotted (km) (optional)
nlevs            = 40                #-- number of levels on colorbar
plottitle        = ''                #-- plot title
fs               = 10                #-- fontsize
dpi              = 400               #-- resolution of plot
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


#-- find gridcell
lat_icon      = np.array(xr.open_dataset(icon_grid)["clat"]) * 180.0/np.pi                #-- get latidues from gridfile
lon_icon      = np.array(xr.open_dataset(icon_grid)["clon"]) * 180.0/np.pi                #-- get longitudes from gridfile
distances     = distance([lon, lat], [lon_icon, lat_icon])                                #-- distances between every grid cell and defined location
cell_idx      = np.argmin(distances)                                                      #-- cell index (smallests distance)


#-- find files in folder
folderfiles       = sorted(np.array(os.listdir(folder)))                                  #-- all files in folder
usedata           = [prefix in x for x in folderfiles]                                    #-- look for prefix
filenames         = np.array(folderfiles)[usedata]                                        #-- choose only data with prefix

time              = []                                                                    #-- generate time array
time_interv       = np.zeros(len(filenames))                                              #-- generate a time difference array (difference to start date)

#-- load icon data in a loop (create a column of data for every file)
for i in range(len(filenames)):          
  #-- Load ICON data
  try:
    ds_icon      = xr.open_dataset(folder + '/' + filenames[i])                                        #-- open icon dataset
  except ValueError:                                                                                   #-- except not readable as netcdf
    logging.disable(logging.ERROR)                                                                     #-- suppress error messages
    ds_icon      = xr.open_dataset(folder + '/' + filenames[i], engine="cfgrib", indexpath='')         #-- open dataset as grib
    logging.disable(logging.NOTSET)                                                                    #-- enable error messages
  
  if 'time' in ds_icon.coords:                                                                         #-- check if 'time' is a dimension in the dataset 
    ds_icon = ds_icon.isel(time=0)                                                                     #-- if so, remove the dimension

  data = ds_icon[var][:, cell_idx].values                                                              #-- get data at cell index
  time.append(str(ds_icon.coords['time'].values))                                                      #-- get time 

  #-- get zinfo and var attributes
  if i == 0:                                                                                           #-- execute only when first file
  
    start_date = datetime.strptime(str(ds_icon.coords['time'].values),'%Y-%m-%dT%H:%M:00.000000000')   #-- set start time and date
    units   = ds_icon[var].attrs['units']                                                              #-- get units info
    n_lev   = len(data)                                                                                #-- number of icon levels
    profile = np.zeros([n_lev, len(filenames)])                                                        #-- create array for plotting
    height  = np.zeros([n_lev, len(filenames)])                                                        #-- create array with height info

    try:
      altitudes    = np.array(ds_icon[alt_var])                                                        #-- check if altitude variable' is included in the dataset
      heightinfo   = True                                                                              #-- indicates height info is available
    except  KeyError:
      if 'alt_data' in locals():                                                                       #-- check if a supplementary file with height info was specified 
        try:
          altitudes    = xr.open_dataset(alt_data)                                                     #-- read in the supplementary data containing height info
        except ValueError:
          logging.disable(logging.ERROR)                                                               #-- suppress error messages
          altitudes    = xr.open_dataset(alt_data, engine="cfgrib", indexpath='')                      #-- open dataset as grib
          logging.disable(logging.NOTSET)                                                              #-- enable error messages
        heightinfo   = True                                                                            #-- indicates height info is available
        if 'time' in altitudes.coords:                                                                 #-- check if supplementary data has a time dimension
          altitudes = altitudes.isel(time=0)                                                           #-- remove time dimension  
        altitudes   = altitudes[alt_var].values                                                        #-- get values of height info
      else: 
        print('  Warning: altitude info is missing!')                                                  #-- if no height info is available
        altitudes  = np.arange(n_lev,0, -1)                                                            #-- enunmerate layers instead  
        heightinfo = None                                                                              #-- indicates height info is NOT available  
   
    if not heightinfo == None: 
      if altitudes.shape[0] == n_lev+1:
        altitudes   = altitudes[1:,]                                                                   #-- if variable has less levels than height info
      altitudes = altitudes[:, cell_idx]/1000                                                          #-- convert m to km 
    

  this_time = datetime.strptime(str(ds_icon.coords['time'].values),'%Y-%m-%dT%H:%M:00.000000000')      #-- get the time and date of this file
  time_interv[i]  = (this_time-start_date).total_seconds()/60                                          #-- compute time difference to start date
  profile[:,i] = np.flip(data)                                                                         #-- extract data at cell index (reverse levels)
  height[:,i]  = np.flip(altitudes)                                                                    #-- extract z_if at cell index (reverse levels)
 
#-- plot settings
cmap        = plt.get_cmap(cbar_col, nlevs-1)                                          #-- read the color map

if not 'varMin' in locals():                                                           #-- check if Min and Max values for colorbar are assigned
  varMin    = np.nanmin(profile)                                                       #-- otherwise set Min to Min of the selcted data
  varMax    = np.nanmax(profile)                                                       #-- and set Max to Max of the selcted data

if not 'alt_plot' in locals():                                                         #-- check if altitude boundaries are set
  if not heightinfo == None:
    alt_plot = [np.min(altitudes), np.max(altitudes)]                                  #-- otherwise set it to Min and Max z_ifc

time       = np.array([datetime.strptime(str(x),'%Y-%m-%dT%H:%M:00.000000000') for x in time])   #-- convert date and time to datetime format
minutes    = np.array([float(x.strftime('%H'))*60 + float(x.strftime('%M')) for x in time])      #-- convert date abd time into minutes of the day
time_index = np.arange(len(time))                                                                #-- create array for time index

if plotincr == 'daily':                                                                #-- index beginning of day 
  time_index = time_index[minutes == 0]                                         
elif plotincr == 'hourly':                                                             #-- index beginning of hour
  time_index = time_index[minutes % 60. == 0]
else:                                                                                  #-- index defined time increment
  time_index = time_index[time_interv  % (60.*plotincr)  == 0.]  

dates = np.array([x.strftime('%d-%m-%Y %H:%M') for x in time])                         #-- covert time for plot  (remove seconds)

levels      = np.linspace(varMin, varMax, nlevs)                                       #-- create levels on colorbar
ticklevels  = np.linspace(varMin, varMax, int(nlevs/4))                                #-- create ticks on colorbar
norm        = mpl.colors.BoundaryNorm(levels, cmap.N)                                  #-- normalize for colorbar

plot_x  = np.tile(time_interv, (n_lev,1))                                              #-- create time/height array

#-- create figure
pfig    = plt.figure()                                                                 #-- create figure
ax      = pfig.add_subplot(1,1,1)                                                      #-- create axis 

c = ax.contourf(plot_x, height, profile, levels=levels,cmap=cmap,norm=norm, extend='both')        #-- plot data

#-- adjust yticks and ylabels (altitude)
if not heightinfo == None:
  yticks = np.arange(alt_plot[0], alt_plot[1]+5, 5)                                    #-- arange position of altitude yticks
  ax.set_yticks(yticks)                                                                #-- set a altitude yticks
  ytick_lab =  ['{0:.{1}f}'.format(x,1) for x in yticks]
  ax.set_yticklabels(ytick_lab, fontsize=fs)                                           #-- add altitude labels
  ax.set_ylabel('Altitude (km)', fontsize=fs)                                          #-- add a-axis label
  ax.set_ylim(alt_plot)                                                                #-- set limits of y-axis

ax.set_xticks(time_interv[time_index])                                                 #-- set position of ticks (timeintervals)
ax.set_xticklabels(dates[time_index], fontsize=fs, rotation=60, ha='right')            #-- add date to ticks

#-- title
ax.set_title(plottitle, pad=10, fontsize=fs+fs/3)                                      #-- add plot title

#-- colorbar
cbar = pfig.colorbar(c, ticks=ticklevels, pad=0.01)                                    #-- add a colorbar
cbar.ax.tick_params(labelsize=fs)                                                      #-- set fontsize of ticklabels on colorbar   
cbar_label = units                                                              
cbar.set_label(cbar_label, fontsize=fs)                                                #-- add units as colorbar label

plt.savefig(plotpath + '/' + plotname, dpi=dpi, bbox_inches='tight')                   #-- save figure



