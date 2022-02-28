'''
-------------
 Description:
 Python script to visualize CALIPSO and ICON data on unstructured grid.
 Valid data types: NetCDF, GRIB
 There may be warnings when reading GRIB data, which can be ignored.
 
 Colors based on  ccplot/cmap/calipso-backscatter.cmap 
 Regriding and averaging based on NASA-DEVELOP/VOCAL/calipso/plot 

 24.02.2022  Anika Rohde, KIT
-------------
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from   datetime import datetime, timedelta
from scipy import interpolate
from numpy import ma
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
from pyhdf.SD import SD, SDC
from pyhdf import VS 
from pyhdf.HDF import HDF
import logging 



#--  define path, grid and variable
calipso_var      = 'Total_Attenuated_Backscatter_532'
icon_var         = 'bsc_dust_532nm'
calipso_data     = ''
icon_data        = ''
icon_grid        = '' 
plotpath         = ''
#-- if icon data includes different time steps
i_step           = 0


#-- if z_ifc should be loaded from a different file (optional)
z_ifc_data       = ''

#-- set maximum of altitude to be plotted and amount of CALIOP columns to be averaged horizontally 
max_alt          = 15.     #-- km
horiz_avg        = 121     #-- columns  (adjust according to the horizontal resolution of your icon data)

#-- set boundaries for the plots (optional)
lat_bounds       = [30, 70]   # or
#lon_bounds       = [-10, 70] 


#-- plot CALIPSO flight track (optional)
plot_track       = True
map_lat          = [30, 70]
map_lon          = [-10, 70]

#-- settings of plots
calipso_plottitle = 'Total Attenuated Backscatter 532nm (km$^{-1}$ sr$^{-1}$)'
calipso_plotname  = 'Calipso.png'
icon_plottitle    = 'Dust Attenuated Backscatter 532nm (km$^{-1}$ sr$^{-1}$)'
icon_plotname     = 'ICON.png'
fs                = 10                                   #-- fontsize
fs_cbar           = fs/2                                 #-- fontsize of colorbar
dpi               = 400                                  #-- resolution of plot
n_ticks           = 8                                    #-- number of latitude/longitude ticks in plot


#--- start processing CALIOP data
hdf_interface     = HDF(calipso_data)                    #-- open calipso file to retrieve metadata
vs_interface      = hdf_interface.vstart()               #-- necessary to get height information 
meta              = vs_interface.attach("metadata")   
field_infos       = meta.fieldinfo()               
all_data          = meta.read(meta._nrecs)[0]      
meta.detach()

data_dictionary     = {}
field_name_index    = 0
for field_info, data in zip(field_infos, all_data):
  data_dictionary[field_info[field_name_index]] = data       
lidar_altitudes     = np.array(data_dictionary['Lidar_Data_Altitudes'])                          #-- height info


datafile    = SD(calipso_data, SDC.READ)                                                         #-- open calipso file to get data 
data_cal    = datafile.select(calipso_var).get()                                                 #-- select data


lat_cal     = datafile.select('Latitude').get()[:,0]                                             #-- Latitudes
lon_cal     = datafile.select('Longitude').get()[:,0]                                            #-- Longitudes
time_cal    = datafile.select('Profile_UTC_Time').get()[:,0]                                     #-- CALIPSO time format: yymmdd.ffffffff    'ffffffff' is the fractional part of the day
surface     = datafile.select('Surface_Elevation').get()[:,0]                                    #-- surface info 
frac_day    = np.array([float(str(x)[6:]) for x in time_cal])                                    #-- ectract fractional part of day
date        = np.array([datetime.strptime(str(x)[:6], '%y%m%d') for x in time_cal])              #-- extract date
date_cal    = np.array([date[x] +  timedelta(hours=frac_day[x]*24) for x in range(date.size)])   #-- add fractional part of days in hours to the date

datafile.end()                                                                                   #-- close data



#-- required functions for plotting

#-- CALIPSO colors
def calipso_cmap():

    from matplotlib.colors import ListedColormap

    ### BOUNDS
    tmp1 = np.arange(0.0001,0.001,0.0001)
    tmp2 = np.arange(0.001,0.0085,0.0005)
    tmp3 = np.arange(0.01,0.11,0.01)
    bounds = np.concatenate((tmp1,tmp2,tmp3))

    ### UNDER_OVER_BAD_COLORS
    under, over, bad  = [0, 42, 127],[255, 255, 255],[0, 42, 127]

    ### COLORS
    colors = [[ 18,  52, 140],[ 35,  61, 153],[ 62, 119, 185],[ 67, 161, 217],
              [ 65, 200, 240],[108, 202, 220],[109, 199, 182],[110, 195, 146],
              [ 12, 128, 128],[ 25, 170,  86],[242, 234,  26],[241, 235,  26],
              [250, 212,   4],[250, 168,  25],[244, 127,  31],[240,  85,  36],
              [235,  32,  35],[237,  47,  90],[238,  84, 126],[241, 127, 169],
              [ 70,  70,  70],[100, 100, 100],[130, 130, 130],[155, 155, 155],
              [180, 180, 180],[200, 200, 200],[225, 225, 225],[235, 235, 235],
              [240, 240, 240],[242, 242, 242],[245, 245, 245],[249, 249, 249],
              [253, 253, 253]]

    ## transform rgb values
    colors = [ [ c[rgb]/255. for rgb in range(3) ] for c in colors ]
    under  = [ rgb/255. for rgb in under ]
    over   = [ rgb/255. for rgb in over ]
    bad    = [ rgb/255. for rgb in bad ]

    ### Define colorbar
    cmap = ListedColormap(colors)
    cmap.set_under(under)
    cmap.set_over(over)
    cmap.set_bad(bad)

    return cmap,bounds


def uniform_alt_2(max_altitude, old_altitude_array):
    """
    Description (8/2013):  Builds a uniformly spaced altitude grid above region 2
    of the CALIOP lidar data.  The spacing is 30 km.  From what I have been told
    the idea here is to have a 30 km spacing instead of 60 km.  Note that the altitude
    is stored in descending order.

    Parameters:
    max_altitude        - [in] maximum altitude for grid (should be above region 2)
    old_altitude_array  - [in] region 2 altitudes

    new_alt             - [out] output array with region 2 and above

    """
    D_ALT = 0.03 # spacing is 30 km
    MID_RES_TOP = 288
    MID_RES_BOT = 576

    # Altitude indices for high res altitude region (region 2):
    # 288:576

    alt2 = old_altitude_array[MID_RES_TOP:MID_RES_BOT]

    new_num_bins = int(np.ceil((max_altitude-alt2[0])/D_ALT))

    new_length = int(new_num_bins + len(alt2))

    new_alt = np.zeros(int(new_length))
    new_alt[int(new_num_bins):int(new_length)] = alt2

    upper_altitudes =  (np.arange(new_num_bins) + 1.)*D_ALT
    new_alt[:int(new_num_bins)] = new_alt[int(new_num_bins)] + upper_altitudes[::-1]

    return new_alt

def regrid_lidar(alt, inMatrix, new_alt, method = 'linear'):
    """
    This function will regrid the matrix inMatrix defined by the (Nx1) vector 'alt'
    onto the new grid (Jx1) 'new_alt'.
    The assumption is that the horizontal dimension changes column by
    column, and altitude is stored row by row (e.g. row x col == alt x (dist
    or time).

    Note that all values outside of bounds are returned as NaN's
    For interp1d to work, the ordinate array has to be monotonically increasing
    This is why the altitude and inMatrix arrays have been reversed in their
    common dimension.
    """

    regrided = interpolate.interp1d(alt[::-1], np.transpose(inMatrix[::-1,:]), kind=method,
                                      axis=0, bounds_error=False)(new_alt)

    return np.transpose(regrided)


def avg_horz_data(data, N):
    """
    This function will average lidar data for N profiles.
    Inputs:
        data - the lidar data to average. Profiles are assumed to be stored in columns and to be a masked array.
        N    - the number of profile to average
    Outputs:
        out - the averaged data array.
    """
    nAlts = data.shape[1]
    nProfiles = data.shape[0]


    nOutProfiles = np.floor(nProfiles/N)
    out = np.zeros((int(nAlts), int(nOutProfiles)))

    for i in range(0, int(nOutProfiles) - 1):
        out[:, i] = ma.mean(np.transpose(data)[:, i*N:(i+1)*N - 1], axis=1)

    return np.transpose(out)


def cut_data(data, start, end, lons, lats, cuttype):
    """
    This function will cut out data along latitude or longitude boundaries
    Inputs:
        data    - CALIPSO data
        start   - minimum latitude or longitude
        end     - maximum latitude or longitude
        lons    - array of longitudes according to CALIPSO data
        lats    - array of latitudes according to CALIPSO data
        cuttype - string 'lat' or 'lon' indicating which to be cut
    Outputs:
        new_data    - adjusted CALIPSO data
        lons[index] - new array of longitudes according to adjusted CALIPSO data
        lats[index] - new array of latitudes according to adjusted CALIPSO data
        index       - index of adjusted data to reproduce cut in other arrays
    """
    if cuttype == 'lat':
       index   = np.where((lats>=start) & (lats<=end))[0]
    elif cuttype == 'lon':
       index   = np.where((lons>=start) & (lons<=end))[0]

    new_data = data[index]

    return new_data, lons[index], lats[index], index

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


#-- prepare CALIOP data
if max_alt > max(lidar_altitudes): max_alt = max(lidar_altitudes)        #-- check if max_alt is bigger than max CALIPSO altitude        

#-- cut data (optional) -> lat_bounds or lon_bounds need to be defined at the top
if 'lat_bounds' in locals():
    data_cal, lon_cal, lat_cal, cut_index =  cut_data(data_cal, min(lat_bounds), max(lat_bounds), lon_cal, lat_cal, cuttype='lat')
    surface = surface[cut_index]
elif 'lon_bounds' in locals():
    data_cal, lon_cal, lat_cal, cut_index =  cut_data(data_cal, min(lon_bounds), max(lon_bounds), lon_cal, lat_cal, cuttype='lon')
    surface = surface[cut_index]


data_cal_flip   = np.flip(data_cal)                                      #-- flip data 
data_cal_flip   = avg_horz_data(data_cal_flip, horiz_avg)                #-- horizontal average by columns
new_alt         = uniform_alt_2(max_alt, lidar_altitudes)                #  
regrid_data     = regrid_lidar(lidar_altitudes,data_cal_flip,new_alt)    #-- interpolate data


#-- plot CALIPSO data
plot_data = np.transpose(regrid_data)

cmap, levels = calipso_cmap()                                            #-- get CALIPSO colormap
norm         = mpl.colors.BoundaryNorm(levels, cmap.N)                   

pfig    = plt.figure(figsize=(regrid_data.shape[0]*(4*horiz_avg/7)/regrid_data.shape[1],4))  #-- start figure
ax1     = pfig.add_subplot(1,1,1)                                                            #-- latitudes x-axis 
ax2     = ax1.twiny()                                                                        #-- longitudes x-axis
ax3     = ax1.twinx()                                                                        #-- second y-axis for surface
im      = ax1.imshow(plot_data,    cmap=cmap, aspect='auto',                                 #-- plot data
                     norm=norm,    interpolation='nearest')

#-- adjust yticks and ylabels (ax1)
pos_yticks = [int(np.where(min(abs(new_alt-x)) == abs(new_alt-x))[0]) for x in range(0, int(max_alt)+5, 5)] 
ax1.set_yticks(pos_yticks)
ax1.set_yticklabels(range(0, int(max_alt)+5, 5), fontsize=fs)
ax1.set_ylabel('Altitude (km)', fontsize=fs)
ax1.set_ylim([regrid_data.shape[1],0])


#-- adjust xticks and xlabels (ax1 & ax2)
lab_xticks_lat = np.round(np.linspace(lat_cal[0], lat_cal[len(lat_cal)-1], n_ticks),2)
lab_xticks_lon = np.round(np.linspace(lon_cal[0], lon_cal[len(lon_cal)-1], n_ticks),2)
lat_labels     = ['{0:.{1}f}'.format(x,2) for x in lab_xticks_lat]                     #-- round, 2 decimals (lat labels)
lon_labels     = ['{0:.{1}f}'.format(x,2) for x in lab_xticks_lon]                     #-- round, 2 decimals (lon labels)
lat_labels[0]  = 'Lat   '+ lat_labels[0]
lon_labels[0]  = 'Lon   ' + lon_labels[0]
pos_xticks     = np.linspace(0, regrid_data.shape[0], n_ticks)

#-- Latiitudes on x-axis
ax1.set_xlim([0, regrid_data.shape[0]])
ax1.set_xticks(pos_xticks)
ax1.set_xticklabels(lat_labels, fontsize=fs, ha='right')

#-- Longitudes on x-axis
x_axis = ax2.get_xaxis()
x_axis.set_visible(True)
ax2.set_xticks(pos_xticks)
ax2.set_xticklabels(lon_labels, fontsize=fs, ha='right')
ax2.xaxis.set_ticks_position("bottom")                                                 #-- show ticks of second x-axis on the bottom
ax2.xaxis.set_label_position("bottom")                                                 #-- show labels of second x-axis (longitudes) on the bottom
ax2.tick_params(axis='x', which='major', pad=13)                                       #-- offset second x-axis  


#-- title
ax1.set_title(calipso_plottitle, pad=10, fontsize=fs+fs/3)                             #-- plot title

#-- draw surface (optional)
ax3.set_ylim([min(new_alt), max(new_alt)])                                             #-- set min and max of second y-axis
ax3.plot(np.linspace(0,regrid_data.shape[0], surface.size), surface, color='red')      #-- draw surface
y_axis = ax3.get_yaxis()
y_axis.set_visible(False)                                                              #-- dont show second y-axis

#-- colorbar
tick_label = ['1.0x10$^{-4}$','2.0','3.0','4.0','5.0','6.0','7.0','8.0','9.0',
              '1.0x10$^{-3}$','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0',
              '5.5','6.0','6.5','7.0','7.5','8.0',
              '1.0x10$^{-2}$','2.0','3.0','4.0','5.0','6.0','7.0','8.0','9.0',
              '1.0x10$^{-1}$']

cbar = pfig.colorbar(im, ticks=levels, pad=0.01)
cbar.set_ticklabels(tick_label)
cbar.ax.tick_params(labelsize=fs_cbar)
#-- colorbar label:
#cbar_label = 'Total Attenuated Backscatter 532nm (km$^{-1}$ sr$^{-1}$)'    
#cbar.set_label(cbar_label, fontsize=fs)                                   

plt.savefig(plotpath + calipso_plotname, dpi=dpi, bbox_inches='tight')                  #-- save figure





#-- Prepare ICON data
try:
  ds_icon      = xr.open_dataset(icon_data)                                               #-- open icon dataset
except ValueError:                                                                        #-- except not readable as netcdf
  logging.disable(logging.ERROR)                                                          #-- supress error messages
  ds_icon      = xr.open_dataset(datafilepath, engine="cfgrib", indexpath='')             #-- open dataset as grib
  logging.disable(logging.NOTSET)                                                         #-- enable error messages


try:
  z_ifc        = np.array(ds_icon['z_ifc'])                                             #-- check if 'z_ifc' is included in the dataset
except  KeyError:
  try:
    z_ifc        = np.array(xr.open_dataset(z_ifc_data)['z_ifc'])                                  #-- otherwise read in the additional data containing z_ifc
  except ValueError:
    logging.disable(logging.ERROR)                                                                 #-- supress error messages
    z_ifc        = np.array(xr.open_dataset(z_ifc_data, engine="cfgrib", indexpath='')['z_ifc'])   #-- open dataset as grib
    logging.disable(logging.NOTSET)                                                                #-- enable error messages


ds_icon  =  ds_icon[icon_var]  
  
if 'time' in ds_icon.coords:                                                            #-- check if 'time' is a dimension in the dataset 
  ds_icon = ds_icon.isel(time=0)                                                        #-- if so, remove the dimension
if 'step' in ds_icon.coords:                                                            #-- check if 'step' is a dimension in the dataset
  ds_icon = ds_icon.isel(step=i_step)                                                   #-- if so, chose time step


lat_icon     = np.array(xr.open_dataset(icon_grid)["clat"]) * 180.0/np.pi               #-- get latidues from gridfile
lon_icon     = np.array(xr.open_dataset(icon_grid)["clon"]) * 180.0/np.pi               #-- get longitudes from gridfile

#-- reduce icon data (saves computing time)
icon_ind     = np.unique(np.where((lat_icon >= min(lat_cal)) & (lat_icon <= max(lat_cal)) & (lon_icon >= min(lon_cal)) & (lon_icon <= max(lon_cal)))[0] )
ds_icon      = ds_icon [:,icon_ind]
lat_icon     = lat_icon[icon_ind]
lon_icon     = lon_icon[icon_ind]
z_ifc        = z_ifc   [:,icon_ind]

#-- find the respetive icon data to the Calipso data (closest cells in an array icon_index)
icon_data   = np.zeros([data_cal.shape[0], ds_icon.shape[0]])
icon_height = np.zeros([data_cal.shape[0], ds_icon.shape[0]]) 
icon_index  = np.zeros(data_cal.shape[0])

for j in range(0,len(lat_cal)):

    dist = distance([lon_cal[j],lat_cal[j]],[lon_icon,lat_icon])      
    location = np.argmin(dist)
    
    icon_data[j,:]      = ds_icon[:, location]
    icon_height[j,:]    = z_ifc[:z_ifc.shape[0]-1, location]
    icon_index[j]       = location
    


#-- plot ICON data
plot_data = icon_data*1000.                               #-- convert backscatter:   CALIPSO in km-1 sr-1, ICON in m-1 sr-1
lat    = np.tile(lat_cal, (len(icon_height[0,:]),1))
lat    = np.swapaxes(lat,0,1)
height = icon_height/1000                                 #-- convert m in km
 
cmap, levels = calipso_cmap()
norm         = mpl.colors.BoundaryNorm(levels, cmap.N)

pfig    = plt.figure(figsize=(regrid_data.shape[0]*(4*horiz_avg/7)/regrid_data.shape[1],4))
ax1     = pfig.add_subplot(1,1,1)                                                        #-- first x and y-axis (latitudes)
ax2     = ax1.twiny()                                                                    #-- second x-axis (longitude)

c = ax1.contourf(lat,height,plot_data, levels,cmap=cmap,norm=norm, extend='both')        #-- plot data
ax1.plot(lat_cal, height[:,89], color='red')                                             #-- draw surface


#-- adjust xticks and xlabels
lab_xticks_lat = np.round(np.linspace(lat_cal[0], lat_cal[len(lat_cal)-1], n_ticks),2)
lab_xticks_lon = np.round(np.linspace(lon_cal[0], lon_cal[len(lon_cal)-1], n_ticks),2)
lat_labels  = ['{0:.{1}f}'.format(x,2) for x in lab_xticks_lat]                          #-- round, 2 decimals (lat labels)
lon_labels  = ['{0:.{1}f}'.format(x,2) for x in lab_xticks_lon]                          #-- round, 2 decimals (lon labels)
lat_labels[0] = 'Lat   '+ lat_labels[0]
lon_labels[0] = 'Lon   ' + lon_labels[0]
pos_xticks     = np.linspace(0, regrid_data.shape[0], n_ticks)

#-- Latitudes on x-axis
ax1.set_xlim([0, regrid_data.shape[0]])
ax1.set_xticks(pos_xticks)
ax1.set_xticklabels(lat_labels, fontsize=fs, ha='right')

#-- Longitudes on x-axis
x_axis = ax2.get_xaxis()
x_axis.set_visible(True)
ax2.set_xticks(pos_xticks)
ax2.set_xticklabels(lon_labels, fontsize=fs, ha='right')
ax2.xaxis.set_ticks_position("bottom")                                                   #-- show ticks of second x-axis on the bottom
ax2.xaxis.set_label_position("bottom")                                                   #-- show labels of second x-axis (longitudes) on the bottom
ax2.tick_params(axis='x', which='major', pad=13)                                         #-- offset second x-axis

#-- adjust yticks and ylabels
pos_yticks = np.arange(0, int(max_alt)+5, 5).astype(int)
ax1.set_yticks(pos_yticks)
ax1.set_yticklabels(np.arange(0, int(max_alt)+5, 5), fontsize=fs)
ax1.set_ylabel('Altitude (km)', fontsize=fs)
ax1.set_ylim([min(new_alt),max_alt])


#-- title
ax1.set_title(icon_plottitle, pad=10, fontsize=fs+fs/3)

#-- colorbar
tick_label = ['1.0x10$^{-4}$','2.0','3.0','4.0','5.0','6.0','7.0','8.0','9.0',
              '1.0x10$^{-3}$','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0',
              '5.5','6.0','6.5','7.0','7.5','8.0',
              '1.0x10$^{-2}$','2.0','3.0','4.0','5.0','6.0','7.0','8.0','9.0',
              '1.0x10$^{-1}$']

cbar = pfig.colorbar(im, ticks=levels, pad=0.01)
cbar.set_ticklabels(tick_label)
cbar.ax.tick_params(labelsize=fs_cbar)
#-- colorbar label:
# cbar_label = 'Total Attenuated Backscatter 532nm (km$^{-1}$ sr$^{-1}$)'    
# cbar.set_label(cbar_label, fontsize=fs) 

plt.savefig(plotpath + icon_plotname, dpi=dpi, bbox_inches='tight')                     #-- save figure






#-- flight track  (optional)

if plot_track:

  fs          = 10                                                   #-- fontsize
  crs_data    = ccrs.PlateCarree()                                   #-- projection  
  
  pfig    = plt.figure()                                             #-- start figure
  fig     = pfig.add_subplot(1,1,1,projection=crs_data)              #-- create map
  
  #-- background
  fig.coastlines(resolution='10m', lw=0.5)                           #-- draw coastlines 
  fig.stock_img()                                                    #-- draw background color
  fig.add_feature(cfeature.BORDERS, linestyle='--', lw=0.2)          #-- draw country boarders
  fig.add_feature(cfeature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', '10m'), facecolor='None', edgecolor='lightskyblue', alpha=0.5, lw=0.3)
  fig.add_feature(cfeature.NaturalEarthFeature('physical', 'lakes', '10m'), alpha=0.5)
  
  #-- CALIPSO flight track
  fig.plot(lon_cal, lat_cal, markersize=130, transform=ccrs.PlateCarree(), c="magenta", markeredgecolor="magenta", zorder=12)  #-- draw flight track 
  
  if not 'map_lat' in locals():
    map_lat = [min(lat_cal), max(lat_cal)]                           #-- adjust map margin if map_lat is defined 

  if not 'map_lon' in locals():
    map_lon = [min(lon_cal), max(lon_cal)]                           #-- adjust map margin if map_lon is defined

  #-- xticks (longitude)
  xticks = np.linspace(-180, 180, 19)                                                          
  xticks = xticks[np.where((xticks >= min(map_lon)) & (xticks <= max(map_lon)))[0]]
  fig.set_xticks(xticks, crs=ccrs.PlateCarree())
  fig.set_xlim(map_lon)
  fig.set_xlabel('Longitude ($^\circ$)', fontsize= fs)  
  
  #-- yticks (latitude)
  yticks = np.linspace(90, -90, 19)
  yticks = yticks[np.where((yticks >= min(map_lat)) & (yticks <= max(map_lat)))[0]]
  fig.set_ylim(map_lat)
  fig.set_yticks(yticks, crs=ccrs.PlateCarree())
  fig.set_ylabel('Latitude ($^\circ$)', fontsize= fs)     
  
  #-- ticks around map and grid
  fig.tick_params(direction='in', right=True, top=True, left=True, bottom=True )
  fig.grid(ls='-', lw=0.4, color='white', alpha=0.8)                        #-- add a grid
  
  #-- legend
  fig.legend(['CALIPSO'])                                                   #-- add a legend
  
  #-- maximize and save the PNG file
  plt.savefig(plotpath + 'Flighttrack.png', dpi=dpi, bbox_inches='tight')   #-- save figure
