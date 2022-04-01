'''
-------------
 Description: 
 Python quicklook script to visualize ICON meteogram 3D and 2D data. 
 File "meteogram_varlims.txt" is required for min and max values.

 Valid data types: NetCDF
 
 01.04.2022  Anika Rohde, KIT
-------------
'''

import numpy as np
import xarray as xr 
import matplotlib.pyplot as plt 
import matplotlib as mpl
from termcolor import colored
import math,  sys, csv
from datetime import datetime
from matplotlib.ticker import MultipleLocator, NullFormatter, LogLocator                        


import warnings
warnings.filterwarnings("ignore")

fs               = 10                                                                                 #-- fontsize
cbar_col         = 'Spectral_r'                                                                       #-- colorbar
cbar_div         = 'RdBu_r'                                                                           #-- colorbar
plotincr         = 12                                                                                 #-- ticks in plot (h)
minorticks       = 6                                                                                  #-- minor ticks in plot (x-axis)  
nlevs            = 41                                                                                 #-- number of levels on colorbar
alt_plot         = [0.,15.]                                                                           #-- Min and Max altitude to be plotted (km) (optional)
yticks           = np.arange(alt_plot[0], alt_plot[1]+5, 5)                                           #-- arange position of altitude yticks
ytick_lab        = np.array(['{0:.{1}f}'.format(x,1) for x in yticks])


print('')
datafilepath             = str(input(colored('Enter File:                            ', 'blue')))     #-- request input dataset
dataset                  = xr.open_dataset(datafilepath)                                              #-- open dataset  

nstations                = dataset.nstations.values                                                   #-- numbers of stations in dataset  
stationnames             = dataset.station_name.values                                                #-- station names in dataset
stations                 = np.array([x.decode("utf-8") for x in stationnames])                        #-- convert station names to strings
stationsup               = np.array([x.upper() for x in stations])                                    #-- station names in upper case
station_height           = dataset.station_hsurf.values/1000                                          #-- station height in km

#-- 3D Vars
nvars                    = np.array(dataset.nvars)                                                    #-- numbers of 3D variables in dataset
var_name                 = dataset.var_name.values                                                    #-- 3D variable names in dataset
variables                = np.array([x.decode("utf-8").upper() for x in var_name])                    #-- convert variable names to strings
var_unit                 = dataset.var_unit.values                                                    #-- units of 3D variables in dataset 
units                    = np.array([x.decode("utf-8").replace('[', '').replace(']', '') for x in var_unit]) #-- convert units of 3D variables to strings
units                    = np.array([x.replace('**', '').replace('^', '') for x in units])                   #-- remove [, ], ^, **,
units                    = np.array([x.replace('/s2', ' s-2').replace('/m', ' m-').replace('/kg', ' kg-1').replace('/s', ' s-1') for x in units]) 
var_group_id             = dataset.var_group_id.values                                                #-- get var group IDs
heights                  = dataset.heights                                                            #-- altitudes of 3D variables in dataset

#-- 2D Vars
nsfcvars                 = np.array(dataset.nsfcvars)                                                 #-- numbers of 2D variables in dataset
sfcvar_name              = dataset.sfcvar_name.values                                                 #-- 2D variable names in dataset
sfcvar                   = np.array([x.decode("utf-8").upper() for x in sfcvar_name])                 #-- 2D variable names in dataset
sfcvar_unit              = dataset.sfcvar_unit.values                                                 #-- units of 2D variables in dataset
sfc_units                = np.array([x.decode("utf-8").replace('[', '').replace(']', '') for x in sfcvar_unit])  #-- remove [, ], ^, **
sfc_units                = np.array([x.replace('**', '').replace('^', '') for x in sfc_units]) 
sfc_units                = np.array([x.replace('/s2', ' s-2').replace('/m', ' m-').replace('/kg', ' kg-1').replace('/s', ' s-1') for x in sfc_units]) 

atmvars                  = dataset['values']                                                          #-- values of 3D variables 
sfcvars                  = dataset['sfcvalues']                                                       #-- values of 2D variables

vardates                 = dataset.date.values                                                                  #-- get date and time  
dates                    = np.array([datetime.strptime(x.decode("utf-8"), '%Y%m%dT%H%M%SZ') for x in vardates]) #-- convert date and time into datetime object
dates_out                = np.array([datetime.strftime(x, '%d-%m-%Y %H:%M') for x  in dates])                   #-- convert date to string for printing
time_interv              = np.array([(x - dates[0]).total_seconds()/60 for x in dates])                         #-- compute time difference to start date and time
minutes                  = np.array([float(x.strftime('%H'))*60 + float(x.strftime('%M')) for x in dates])      #-- convert date abd time into minutes of the day
time_index               = np.arange(len(dates))
time_index               = time_index[time_interv  % (60.*plotincr)  == 0.]                                     #-- compute wich elements of time and date have the time increment plotincr 
plotdates                = np.array([x.strftime('%d.%m.  %H:00') for x in dates])                               #-- covert date to string for plot  (remove seconds)
fcrun                    = datetime.strftime(dates[0], '%Y%m%d%H')                                              #-- get first date for the plotname (fcst-run)  
#days                     = np.array([x.strftime('%d.%m.  %H:00') for x in dates[minutes == 0]])
#plotdates                = np.array([x.strftime('        %H:00') for x in dates])                         #-- covert time for plot  (remove seconds)
#plotdates[minutes == 0]  = days


#-- Reading Min and Max of all variables from external file
varlimname, varlimmin, varlimmax, whitevals, varspace = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
with open('Meteogram_varlims.txt', newline='', mode='r') as f:    

  spamreader = csv.reader(f, delimiter=',')
  next(spamreader)
  for row in spamreader: 
    varlimname = np.hstack([varlimname, row[0].replace(' ' , '').upper()])
    varlimmin  = np.hstack([varlimmin , float(row[1].replace(' ' , '')) ])
    varlimmax  = np.hstack([varlimmax , float(row[2].replace(' ' , '')) ])
    whitevals  = np.hstack([whitevals , float(row[3].replace(' ' , '')) ])
    varspace   = np.hstack([varspace  , row[4].replace(' ' , '').upper()])


print('')

continueloop = 'y'                                                   #-- start loop to examine variables 
showvar = str(input('Show available data? (y/n):            '))                                         

if showvar == 'y':
  print('')
  print('Stations:')
  print(stations)                                                    #-- show variables                          
  print('')
  print('3D Vars:')
  print(variables)                                                   #-- show coordinates
  print('')
  print('2D Vars:')
  print(sfcvar)                                                      #-- show coordinates
  print('')
  print('time:')
  print(dates_out)                                                   #-- show date/time
  print('')

while continueloop == 'y':
  
  varassign = False 
  compsum        = False 

  while varassign == False:
    var       = str(input('Enter variable:                        '))       #-- input variable
    varup      = var.upper()                                                #-- variable in upper case    

    warningmessage = 'Error: ' + colored(var, 'red') + ' is not a variable in this dataset!' 
   
    if varup in variables:                                         #-- 3D vars                    
      atmovar   = True                                             #-- indicate 3D var                                    
      varassign = True                                             #-- indicate variable found in dataset
      varidx    = int(nvars[variables==varup])                     #-- var index      
      varunit   = units[varidx]                                    #-- var unit         
      vargroup  = var_group_id[varidx]                             #-- group ID        
   
    elif varup in sfcvar:                                          #-- 2D vars 
      atmovar   = False                                            #-- indicate 2D var 
      varassign = True                                             #-- indicate variable found in dataset
      varidx    = int(nsfcvars[sfcvar==varup])                     #-- var index 
      varunit   = sfc_units[varidx]                                #-- var unit     
    
    elif varup + 'A' in np.hstack([variables, sfcvar]):            #-- input is ending with dust (for dusta, dustb and dustc)
      varassign = True                                             #-- indicate variable found in dataset
      dustvars  = [varup + 'A',   varup + 'B', varup + 'C']        #-- assemble dusta, dustb, and dustc variables
      compsum   = True                                             #-- set "compute sum" True
      if varup + 'A' in variables:                                 #-- 3D vars
        varidx    = int(nvars[variables==varup + 'A'])               #-- var index
        atmovar   = True                                             #-- indicate 3D var
        varunit   = units[varidx]                                    #-- var unit
        vargroup  = var_group_id[varidx]                             #-- group ID
      elif varup + 'A' in sfcvar:                                  #-- 2D vars 
        varidx    = int(nsfcvars[sfcvar==varup + 'A'])               #-- var index
        atmovar   = False                                            #-- indicate 2D var
        varunit   = sfc_units[varidx]                                #-- var unit
      else: 
        print(warningmessage)
  
    elif varup[:-1] + 'A0' in np.hstack([variables, sfcvar]):      #-- input is ending with dust0 (for dusta0, dustb0 and dustc0)
      varassign = True                                             #-- indicate variable found in dataset
      dustvars  = [varup[:-1] + 'A0',   varup[:-1] + 'B0', varup[:-1] + 'C0']     #-- assemble dusta0, dustb0, and dustc0 variables
      compsum   = True                                             #-- set "compute sum" True
      if varup[:-1] + 'A0' in variables:                           #-- 3D vars
        varidx    = int(nvars[variables==varup[:-1] + 'A0'])         #-- var index
        atmovar   = True                                             #-- indicate 3D var
        varunit   = units[varidx]                                    #-- var unit
        vargroup  = var_group_id[varidx]                             #-- group ID
      elif varup[:-1] + 'A0'  in sfcvar:                           #-- 2D vars
        varidx    = int(nsfcvars[sfcvar==varup[:-1] + 'A0'])         #-- var index
        atmovar   = False                                            #-- indicate 2D var
        varunit   = sfc_units[varidx]                                #-- var unit
      else:
        print(warningmessage)
  
    else:
      print('Error: ' + colored(var, 'red') + ' is not a variable in this dataset!')
      print('')


  #-- Create colorlevels
  varMin  = varlimmin[varlimname==varup][0]      
  varMax  = varlimmax[varlimname==varup][0]
  white   = whitevals[varlimname==varup][0]
  spacing = varspace[varlimname==varup][0]
 
  if varassign:
    stationassign  = False
        
    while stationassign  == False:                                            #-- repeat until station input is correct
      station_enter  = str(input('Enter Station:                         ')) 
      try:
        stindex = int(nstations[stationsup == station_enter.upper()][0])      #-- find station index
        stationassign = True                                                  #-- if station is found, indicate
      except (IndexError, TypeError):      
        stationassign = False
        print('Error: ' + colored(station_enter, 'red')+' is not a station in this dataset!')       
        print('')
        
    st_atmvars =  atmvars.isel(nstations=stindex)                             #-- select 3D vars of this station from dataset 
    st_sfcvars =  sfcvars.isel(nstations=stindex)                             #-- select 2D vars of this station from dataset

    #-- 3D variables
    if atmovar:                                                           
      #-- create colorbar (linear, increase, or log)  
      if spacing == 'LINEAR':

        barspacing  ='proportional'
        levels      = np.linspace(varMin, varMax, nlevs)
        levels[0]   = np.max([white, varMin])
        ticklevels  = np.arange(varMin, varMax, (levels[2]-levels[1])*4) +  (levels[2]-levels[1])*4

        if white > varMin:
          ticklevels = np.hstack([white,ticklevels])
        else:
          ticklevels = np.hstack([varMin, ticklevels])
        contours    = ticklevels
    #    ticklabels  = np.array(['{0:.2}'.format(x) for x in ticklevels])

      elif spacing == 'LOG': 

        barspacing  ='uniform'
        levels      = np.logspace(np.log10(varMin), np.log10(varMax), nlevs)
        levels[0]   = np.max([white, varMin])
        ticklevels  = levels[range(0,nlevs, 4)]
        contours    = ticklevels
    #    ticklabels  = np.array(['{:.1e}'.format(x) for x in ticklevels])

      elif spacing == 'INCREASE':

        barspacing  ='proportional'
        levels      = np.arange(varMin,(varMax-varMin)/5., (varMax-varMin)/50.)
        levels[0]   = np.max([white, varMin])
        levels      = np.hstack([levels, np.arange((varMax-varMin)/5., (varMax-varMin)/5.*3., (varMax-varMin)/20.) ])
        levels      = np.unique(np.hstack([levels, np.arange((varMax-varMin)/5.*3, varMax + (varMax-varMin)/10., (varMax-varMin)/10.) ]))
        ticklevels  = np.arange(varMin, varMax, (varMax-varMin)/10.)+(varMax-varMin)/10.

        if white > varMin:
          ticklevels = np.hstack([white,ticklevels])
        else:
          ticklevels = np.hstack([varMin, ticklevels])
        contours     = ticklevels
    #    ticklabels   = np.array(['{0:.2}'.format(x) for x in ticklevels])
 
      if -varMin == varMax:                                            #-- if negative and pos values with same extreme vallues, colormap is set blue to red
        cmap       = plt.get_cmap(cbar_div, len(levels)-1)
      else:                                                            #-- otherwise colormap is Spectral
        cmap       = plt.get_cmap(cbar_col, len(levels)-1)                                 
      norm       = mpl.colors.BoundaryNorm(levels, cmap.N)             #-- normalize for colorbar    

      #-- compute sum of dusta, b and c (or dusta0, b0 and c0)
      if compsum:
        varidx       = int(nvars[variables==dustvars[0]])
        vardata      = st_atmvars.isel(nvars=varidx)
 
        varidx       = int(nvars[variables==dustvars[1]])
        vardatab     = st_atmvars.isel(nvars=varidx)
  
        varidx       = int(nvars[variables==dustvars[2]])
        vardatac     = st_atmvars.isel(nvars=varidx)

        vardata      = np.add(np.array(vardata), np.array(vardatab), np.array(vardatac))
        xr.Dataset.close(vardatab)
        xr.Dataset.close(vardatac)

      else:
        vardata      = st_atmvars.isel(nvars=varidx).values                          #-- get data
  
      varheights   = heights.isel(nstations=stindex).isel(nvars=varidx).values/1000  #-- get respective altitudes and convert m to km

      #-- convert kg-1 to m-3                
      if ('DUST' in varup) & (varunit[-4:] == 'kg-1'):
        rho_data = st_atmvars.isel(nvars=int(nvars[variables=='RHO']))                
        vardata  = vardata * rho_data.values                                         #-- multiply data with air density "RHO"
        varunit  = varunit[:-4] +  'm-3'                                             #-- overwrite unit for plot 
        xr.Dataset.close(rho_data) 
      elif varunit[-4:] == 'kg-1':                                                   #-- if not dust but unit kg-1 ask if data should be converted
        convert = 'n'                                            
        convert = str(input('Convert ' + colored(varunit, 'red') + ' to ' + colored(varunit[:-4] +  'm-3', 'red') + '? (y/n):    '))
        if convert == 'y':
          rho_data = st_atmvars.isel(nvars=int(nvars[variables=='RHO']))              
          vardata  = vardata * rho_data.values                                       #-- multiply data with air density "RHO"
          varunit  = varunit[:-4] +  'm-3'                                           #-- overwrite unit for plot
          xr.Dataset.close(rho_data)

      #-- remove first level (ground level)
      if vargroup == 1: 
        varheights = varheights[:varheights.size-1]
        vardata    = vardata[:,:vardata.shape[1]-1]
        
      #-- prepare plot data 
      plot_data    = np.flip(vardata, axis=1).T                                      #-- flip y-axis  
      plot_height  = np.tile(np.flip(varheights), [vardata.shape[0],1]).T            #-- generate altitude arrays for plot
      pltunit      = varunit.replace('mu', '\u03BC')                                 #-- replace mu with greek letter
      plot_x  = np.tile(time_interv, (varheights.size,1))                            #-- create time/height array  

#      maxval       = np.nanmax(vardata)
#      if abs(maxval) > 10:
#        maxval       = '{:.1f}'.format(maxval)
#      elif abs(maxval) > 0.1:
#        maxval       = '{:.2f}'.format(maxval)
#      else:
#        maxval       = '{:.{n}f}'.format(maxval, n=int(np.ceil(np.log10(abs(maxval)))*-1+2))
#      plottitle    = plottitle +  ', max: ' + maxval + ' ' + pltunit

      print('')
      print('Mean :' + str(np.nanmean(plot_data)) + ' ' + varunit)                              #-- print mean
      print('Min  :' + str(np.nanmin(plot_data)) + ' ' + varunit)                               #-- print Min
      print('Max  :' + str(np.nanmax(plot_data)) + ' ' + varunit)                               #-- print Max
  
      if math.isnan(np.mean(plot_data)):                                                        #-- check for NaNs
        print('')
        print(colored('Warning:  NaNs in ' + var + '!!!', 'red' ))
      print('')

      print('Creating figure...')
      print('')


      pfig    = plt.figure()                                                                           #-- create figure
      ax      = pfig.add_subplot(1,1,1)                                                                #-- create axis
      pfig.subplots_adjust(left=0.18, bottom=0.26)                                                     #-- adjust space on the left, bottom and top

      c = ax.contourf(plot_x, plot_height, plot_data, levels=levels,cmap=cmap,norm=norm, extend='max') #-- plot data
      c.cmap.set_over('black')                                                                         #-- set values above max color level black 
      c.cmap.set_under('white')                                                                        #-- set values below lowest color level / white val white
      plt.contour(plot_x, plot_height, plot_data, levels=contours, colors='black', linewidths=0.4)     #-- add contours


      #-- adjust yticks and ylabels (altitude)
      ax.set_yticks(yticks)                                                                            #-- set a altitude yticks
      ax.set_yticklabels(ytick_lab, fontsize=fs)                                                       #-- add altitude labels
      ax.set_ylim(station_height[stindex],alt_plot[1])                                                 #-- set limits of y-axis
      ax.set_ylabel('AMSL [km]', fontsize=fs)                                                          #-- add a-axis label
      
      ax.set_xticks(time_interv[time_index])                                                           #-- set position of ticks (timeintervals)
      ax.set_xticklabels(plotdates[time_index], fontsize=fs, rotation=60, ha='right', rotation_mode='anchor')            #-- add date to ticks
      ax.set_xlabel('[UTC]', fontsize=fs)
      ax.tick_params(which='major', length=5)                                                          #-- increase length of major ticks  
      ax.minorticks_on()                                                                               #-- add minor ticks
      ax.xaxis.set_minor_locator(MultipleLocator(minorticks*60))                                       #-- set distance between minor ticks
 
      #-- title
      ax.set_title('fcst-run: ' + fcrun + ', ' + station_enter , pad=18, fontsize=fs+fs/3)             #-- add plot title
      #-- colorbar
      cbar = pfig.colorbar(c, ticks=ticklevels, pad=0.01,spacing=barspacing)                           #-- add a colorbar
      cbar.ax.tick_params(labelsize=fs)                                                                #-- set fontsize of ticklabels on colorbar
      cbar.set_label(var + ' [' + pltunit + ']', fontsize=fs)                                          #-- add label (variable and unit)    
      plt.show()

    #-- 2D variables
    else:

      #-- compute sum of dusta, b and c (or dusta0, b0 and c0)
      if compsum:
        varidx       = int(nsfcvars[sfcvar==dustvars[0]])
        vardata      = st_sfcvars.isel(nsfcvars=varidx).values
        for i in range(1,3):
          varidx       = int(nsfcvars[sfcvar==dustvars[i]])
          vardata      = vardata +  st_sfcvars.isel(nsfcvars=varidx).values
      else:
        vardata      = st_sfcvars.isel(nsfcvars=varidx).values


      pltunit   = varunit.replace('mu', '\u03BC')                                            #-- replace mu with greek letter
#      maxval       = np.nanmax(vardata)
#      if abs(maxval) > 10:
#        maxval       = '{:.1f}'.format(maxval)
#      elif abs(maxval) > 0.1:
#        maxval       = '{:.2f}'.format(maxval)
#      else:
#        maxval       = '{:.{n}f}'.format(maxval, n=int(np.ceil(np.log10(abs(maxval)))*-1+2))

      print('')
      print('Mean :' + str(np.nanmean(vardata)) + ' ' + varunit)                             #-- print mean
      print('Min  :' + str(np.nanmin(vardata)) + ' ' + varunit)                              #-- print Min
      print('Max  :' + str(np.nanmax(vardata)) + ' ' + varunit)                              #-- print Max

      if math.isnan(np.mean(vardata)):                                                       #-- check for NaNs
        print('')
        print(colored('Warning:  NaNs in ' + var + '!!!', 'red' ))
      print('')

      print('Creating figure...')
      print('')

      pfig    = plt.figure()                                                                 #-- create figure
      pfig.subplots_adjust(left=0.15, bottom=0.26)                                           #-- adjust space on the left and bottom
      ax      = pfig.add_subplot(1,1,1) 
      
      ax.set_title('fcst-run: ' + fcrun + ', ' + station_enter , pad=18, fontsize=fs+fs/3)   #-- add plot title
      ax.set_xticks(time_interv[time_index])                                                 #-- set position of ticks (timeintervals)
      ax.set_xlabel('[UTC]', fontsize=fs)                                                    #-- add xlabel
      ax.set_xticklabels(plotdates[time_index], fontsize=fs, rotation=60, ha='right', rotation_mode='anchor')            #-- add date to ticks
      ax.tick_params(which='major', length=5)                                                #-- adjust length of major ticks 
      ax.minorticks_on()                                                                     #-- turn on minor ticks
      ax.xaxis.set_minor_locator(MultipleLocator(minorticks*60))                             #-- set distance between minor ticks
#      ax.set_xlim([min(time_interv), max(time_interv)])                                     #-- set limits of y-axis
      plt.grid(linestyle=":")                                                                #-- add grid lines (and set linestyle)
      if spacing == 'LOG':                                                                   #-- check if axis should be logarithmic
        ax.set_yscale('log')
        y_minor = LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
        ax.yaxis.set_minor_locator(y_minor)
        ax.yaxis.set_minor_formatter(NullFormatter())
      ax.set_ylabel(var + ' [' + pltunit + ']', fontsize=fs)                                 #-- add ylabel
      ax.set_ylim([varMin, varMax])                                                          #-- set y-axis limits

      plt.plot(time_interv, vardata)                                                         #-- plot data
      plt.show()                                                                             #-- show figure
    plt.close()                                                                              #-- close figure
      

    xr.Dataset.close(st_sfcvars)
    xr.Dataset.close(st_atmvars)

    continueloop = str(input('Visualize other variables? (y/n):      '))
    print('')

xr.Dataset.close(dataset)
print('Bye-bye!')
