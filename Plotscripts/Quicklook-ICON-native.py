'''
-------------
 Description: 
 Python script to visualize ICON data on unstructured grid.
 Valid data types: NetCDF, GRIB
 There may be warnings when reading GRIB data, which can be ignored.
 
 03.02.2022  Anika Rohde, KIT
-------------
'''

import numpy as np
import xarray as xr 
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.colors as mcolors
from   matplotlib.collections import PolyCollection
import cartopy.feature as cfeature
from termcolor import colored
import math, logging, sys

sys.tracebacklimit=0                                                                                    #-- reduce error messages (Tracebacks)

print('')
datafilepath               = str(input(colored('DATA - enter path and datafile name:   ', 'blue')))     #-- read path to dataset
try:
  raw_data                 = xr.open_dataset(datafilepath)                                              #-- try open dataset as netcdf 
except ValueError:                                                                                      #-- except not readable as netcdf 
  logging.disable(logging.ERROR)                                                                        #-- supress error messages
  raw_data                   = xr.open_dataset(datafilepath, engine="cfgrib", indexpath='')             #-- open dataset as grib
  logging.disable(logging.NOTSET)                                                                       #-- enable error messages

grid                       = xr.open_dataset(str(input(colored('GRID - enter path and gridfile name:   ', 'blue'))))[['clon_vertices','clat_vertices']].rename({'cell': 'ncells'})
rad2deg                    = 180.0/np.pi                                                                #-- factor to calculate degeree from radian 

print('')
grid['clon_vertices'], grid['clat_vertices'] = grid['clon_vertices']*rad2deg, grid['clat_vertices']*rad2deg  # convert grid icon info to degrees (lats/lons of triangle vertices)

#-- check if dataset and grid matches in length:
try:
  samesize = grid.sizes['ncells'] == raw_data.sizes['ncells']        #-- netcdf  
except KeyError:
   samesize = grid.sizes['ncells'] == raw_data.sizes['values']       #-- grib
if not samesize:
  print(colored('Error: Dataset and ICON grid do not match!', 'red')); print('') 
  exit()                                                             #-- stop running if diffrent lengths

#-- start processing data variables
continueloop = 'y'                                                   #-- start loop to examine variables 

showvar = str(input('Show available data? (y/n):            '))                                         

if showvar == 'y':
  print(raw_data.data_vars)                                          #-- show variables                          
  print(raw_data.coords)                                             #-- show coordinates

while continueloop == 'y':
  print('')
  data     =  raw_data                                               #-- reset data (for loop) 

  var_dims = None                                                    #-- empty array "var_dims" for variable dimensions
  while var_dims is None:                                            #-- as long as var_dims is empty try to select a variable
    try:
      var      =  str(input('Enter variable:                        '))
      var_dims =  data[var].dims
    except KeyError:                                                 #-- as long as entered variable is part of the dataset var_dims cannot be assigned, loop continues
      print('Error: ' + colored(var, 'red') + ' is not a variable in this dataset!')
      print('')
      pass

  try:                                                               #-- try to get information about unit of variable 
    units    =  data[var].units
  except AttributeError:                                             #-- information is empty for dimensionless variables
    units = ''                                                       #-- in this case set unit ''
  
  if var_dims[0]== 'time':                                           #-- check if time is a dimension in dataset
    if len(data.coords['time']) == 1:                                #-- if time is a dimension with one element
      data  = data.isel(time=0)                                      #-- select first element
      var_dims = data[var].dims                                      #-- assign residual dimensions to 'var_dims'
    else:                                                            #-- if time is a dimension with multiple elements
      print('')
      print('Available time:      ')                                 #-- print available dates
      for i in range(len(data.coords['time'])):
        print(str(data.coords['time'][i].values) +  colored('  ->  '+ str(i), 'red'))
      print('')
      timeint = None
      while timeint is None:                                         #-- loop as long as no valid date was selected
        try:
          data  = data.isel(time=int(input('Enter number:                          '))) #-- Select a date
          timeint = 1
        except (IndexError, ValueError):                             #-- continue loop if number was wrong
          print('Error: wrong number!')
          print('')
          pass
      var_dims = data[var].dims                                      #-- assign residual dimensions to 'var_dims'

  if var_dims[0]== 'step':                                           #-- check if step is a dimension in dataset
    if len(data.coords['step']) == 1:                                #-- if step is a dimension with one element
      data  = data.isel(step=0)                                      #-- select first element 
      var_dims = data[var].dims                                      #-- assign residual dimensions to 'var_dims'
    else:                                                            #-- if step is a dimension with multiple elements
      print('')
      print('Available timesteps: ')                                 #-- print available timesteps 
      for i in range(len(data.coords['step'])):
        print(str(data.coords['step'][i].values) +  colored('  ->  '+ str(i), 'red'))
      print('')
      timeint2 = None
      while timeint2 is None:                                        #-- loop as long no valid timestep was selected
        try:
          data  = data.isel(step=int(input('Enter number:                          '))) #-- Select a timestep
          timeint2 = 1
        except (IndexError, ValueError):                             #-- continue loop if number was wrong
          print('Error: wrong number!')
          print('')
          pass
      var_dims = data[var].dims                                      #-- assign residual dimensions to 'var_dims'

  data  = data[var]                                                  #-- pick variable from dataset
  
  if len(var_dims) == 2:                                             #-- check if there are two residual dimensions (icon levels & number of cells)
    if data.shape[0] ==1:                                            #-- if there are 2 dimensions, but the levels dimension has only one element (e.g. t2m)  
      z        = data[0,:].values                                    #-- assign the values of the first level to z
    else:                                                            #-- if there are multiple levels
      z = None
      while z is None:                                               #-- loop until data of one level is assigned to z
        try:
          lev_index = int(input('Enter level:                           '))        #-- select a level
          z         = data[lev_index-1,:].values
        except (IndexError, ValueError):
          print('Error: There are only ' + colored(str(data.shape[0]) +' levels', 'red')+'!')
          print('')
          pass
  
  elif len(var_dims) == 1:                                            #-- if there is only one dimension left (number of cells)
    z        = data.values                                            #-- values of variable is directly assigned to z
 
  else:                                                               #-- Exception/Error: there are more than 2 dimensions left? Why?
    print('Error: too many dimensions!')
    z = np.array([0,0])                                               #-- assign zeros to z just to continue the loop
 
  varmin, varmax  = np.nanmin(z), np.nanmax(z)                        #-- get min and max of the selected data
  
  if (varmin == varmax):                                                
    print('Error: Variable has the same Min and Max!')                #-- if data is only one value dont create a plot + print value
    print(colored('       ' + var + ' = ' + str(varmin) + ' ' + units, 'red'))

  else:                                                  
    print('')
    print('Mean :' + str(np.nanmean(z)))                              #-- print mean
    print('Min  :' + str(np.nanmin(z)))                               #-- print Min
    print('Max  :' + str(np.nanmax(z)))                               #-- print Max
    
    if math.isnan(np.mean(z)):                                        #-- check for NaNs
      print('')
      print(colored('Warning:  NaNs in ' + var + '!!!', 'red' ))
    print('')
    
    levadjust = str(input('Adjust Min and Max of Colorbar? (y/n): ')) #-- let user change min and max on colorbar
    if levadjust == 'y':
      print('')
      varmin  = float(input('Enter Colorbar Min:                    '))
      varmax  = float(input('Enter Colorbar Max:                    '))

    if (varmin < 0) and (varmax >0):

      #-- make colorbar symmetric to zero
      section        = (varmax -varmin)/33                               #-- get the induvidual size of 33 sections for the colorbar
      neg_sec        = np.floor((abs(varmin)-0.5*section)/section)       #-- amount of full negative sections on colorbar
      pos_sec        = np.floor((abs(varmax)-0.5*section)/section)       #-- amount of full positive sections of colorbar
      neg_part       = (neg_sec+1)/(neg_sec+pos_sec+2)                   #-- negative fraction of colorbar 
      pos_part       = (pos_sec+1)/(neg_sec+pos_sec+2)                   #-- positive fraction of colorbar
      lower          = max(0.5- neg_part/max(neg_part,pos_part)/2., 0.)  #-- norm lower limit between 0 and 0.5 to get associated color 
      upper          = min(0.5+ pos_part/max(neg_part,pos_part)/2., 1.)  #-- norm upper limit between 0.5 and 1. to get associated color
           
      #-- create new continuos colormap with adjusted colors with 'lower' and 'upper' (neg -> blue, pos -> red)
      cmap_complete  = plt.get_cmap('seismic')
      cmap_new       = mcolors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap_complete, a=lower, b=upper),
                                                                  cmap_complete(np.linspace(lower, upper, 265)))
      plt.register_cmap(cmap=cmap_new, name='my_cmap')                   #-- temporarily register color shifted cmap as 'my_cmap'

      #-- create actual colorbar with colors from 'my_cmap'
      full_min       = neg_sec*-section                                  #-- get lower limit of the lowest full section on colorbar)
      full_max       = pos_sec* section                                  #-- get upper limit of the top full section on colorbar)
      collev         = np.sort(np.concatenate([np.array([varmin]),np.arange(0.5*section, full_min-section, -section), np.arange(0.5*section, full_max+section, section), np.array([varmax])]))
      levticks       = np.sort(np.concatenate([np.arange(1.5*section, full_min-3*section,  -3*section), np.arange(1.5*section, full_max+3*section, section*3)]))
      cmap           = plt.get_cmap('my_cmap', len(collev)-1)

    else:
      collev       = np.linspace(varmin,varmax, 41)               #-- set 40 segments on colorbar
      levticks     = np.linspace(varmin,varmax, 11)               #-- set 11 colorbar labels
      cmap  = plt.get_cmap('Spectral_r', len(collev)-1)

    crs_data       = ccrs.PlateCarree()                           #-- projection
    nlevs          = len(collev)                                  #-- number of levels
    
    #-- define the x-, y-values and the polygon points
    vlon, vlat     = grid.clon_vertices, grid.clat_vertices       #-- lon and lat in degrees
    ncells, nv     = vlon.shape[0], vlon.shape[1]                 #-- number of cells and edges
    vlon, vlat     = np.array(vlon), np.array(vlat)
   
    #-- workaround (approximation) for global data
    #-- prevents distortions on the map caused by triangles crossing -180 / 180 longitude
    #-- solution: cut triangles at 180/-180 and copy the other part to other side
    tri_idx         = np.where(np.any(vlon <-100, axis=1) & np.any(vlon >100, axis=1))[0]     #-- find triangles with vertices on both sides 

    if len(tri_idx >0):                     
      copy_lon = vlon[tri_idx]                                    #-- copy lons of vertices crossing -180/180 
     
      for i in range(len(tri_idx)):
        copy_lon[i]     [vlon[tri_idx[i]]<0]  = 180               #-- recreate triangles on the other side, set negative lons to 180  
        vlon[tri_idx[i]][vlon[tri_idx[i]]>0]  = -180              #-- modify/cut original triangles at -180 (set positive lons to -180) 
     
      z    = np.append(z,     z   [tri_idx])                      #-- add data of additional triangles (not modified - copy)          
      vlat = np.vstack((vlat, vlat[tri_idx]))                     #-- add lats of additional triangles (not modified - copy)
      vlon = np.vstack((vlon, copy_lon))                          #-- add lons of additional triangles (modified)

    #-- define color map
    colors   = np.zeros([len(z),4],np.float32)                    #-- assign color array for triangles
    ncol     = cmap.N                                             #-- number of colors
    norm     = mpl.colors.BoundaryNorm(collev, ncol)              #-- normalize for color bar
    
    #-- set color index of all cells in between levels
    for m in range(0,nlevs-1):
      colors[np.where((z    >= collev[m]) & (z < collev[m+1])),:] = cmap(m)
    
    colors[np.where(z < varmin),:]    = cmap(0)                   #-- set color below levels
    colors[np.where(z >= varmax),:]   = cmap(ncol-1)              #-- set color index for cells greater than levels[nlevs-1]
     
    #-- create list of vertices
    triangles = [np.column_stack((vlon[i,:],vlat[i,:])) for i in range(0,len(vlon))]
    coll = PolyCollection(triangles, facecolor=colors,edgecolors="none", closed=True, transform=crs_data)
  
    print('')
    print('Creating figure...')
    print('')
    fig = plt.figure()
    ax  = fig.add_axes([0.09,0.0,0.71,0.85], projection=crs_data)                    #-- margin of the map (values: left corner x, y + width and height of map) + projection
    ax.set_extent([np.min(vlon), np.max(vlon), np.min(vlat), np.max(vlat)], crs=crs_data)    # -- set extend of the map (here min and max lon/lat of data)
  
    ax.coastlines(resolution='50m')                                                  #-- coastlines
    ax.add_feature(cfeature.BORDERS, linestyle=':', lw=1, edgecolor='white')         #-- country boarders in white (for dark colors)
    ax.add_feature(cfeature.BORDERS, linestyle=':', lw=0.8)                          #-- country boarders in black (for bright colors)
    ax.add_collection(coll)                                                          #-- add triangles (data) 
    ax.set_title(var)                                                                #-- set title

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
    ax.set_xticklabels(lon_labels)                                                   #-- add lon labels
    ax.set_yticklabels(lat_labels)                                                   #-- add lat_labels
    ax.tick_params(direction='in', right=True, top=True, pad=4)                      #-- ticks inside figure and space (pad) between numbers and axis  
  
    #-- add colorbar
    ax  = fig.add_axes([0.82, 0.05,0.025,0.8])                                       #-- margin of colorbar (values: left corner x, y + width and height) 
    cbar = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, spacing='proportional',   
                                  ticks=levticks,  orientation='vertical',drawedges=True)
   
    cbar.set_label(units)                                                            #-- label colorbar   
    cbar.dividers.set_color('white')                                                 #-- create white lines between colors
    cbar.dividers.set_linewidth(0.3)                                                 #-- reduce linewidth of the lines between colors
    plt.show()
    plt.close()

    print('Figure closed.')
    print('----------------------------------------')
    print('')

    continueloop = str(input('Visualize other variables? (y/n):      '))
    print('')

print('Bye-bye!')
print('')


