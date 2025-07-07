# -*- coding: iso-8859-1 -*-
'''
field_area_mean_timeseries_CMIP6 

Date: 19 March 2023
Purpose: Build timeseries of area-means for CMIP6 inputs  
Author: Mike Bell   

Change history
 
   23 Oct 2024:   Revised to read one field at a time 

'''

def field_area_mean_timeseries_CMIP6(c_model, c_expt_type):

   import sys        
   import netCDF4
   import numpy as np

   import matplotlib.pyplot as plt
   import matplotlib.ticker as mticker
   import cartopy
   import cartopy.crs as ccrs
   import cartopy.img_transform

   from utilities import period_string, rd_ncdf_var_check_one
   from utilities_CMIP6 import cmip6_file_names, lat_lon_area_mask_CMIP6

###----------------------------------------------------------------------------------
# Start of user inputs 

   dirout = '../Time_series/'
   
   std_out = '../std_out/field_area_mean_timeseries_CMIP6.txt' 
   
   cvar = 'hfds' # 'hfds', 'tos', 'tauu' 
   
   nmonths_ts = 6000   # 23 Oct 2024    number of months in 500 year output timeseries

   c_dir, c_descrip, c_period, c_end = cmip6_file_names(c_model, c_expt_type, cvar)
   print ( 'c_model, c_descrip, c_period, cvar = ', c_model, c_descrip, c_period, cvar ) # , file = std_out )
   
   c_type = 'integral'   #  'integral', 'mean', 'int_pos' (integral of positive values), 'int_neg' ;   NOT 'anom' - that is created by another script  
   
   lat_lon_typ = True # false option removed
   
# define the region over which the mean or integral is calculated    ; 2.0 or 3.0 for 'tauuo' ; 5.0 or 10.0 for 'hfds' 
   if lat_lon_typ : 
      lat_min =    -90.0          # -3.0 or 5.0 or 27.0 or -40.0; -90.0 for global
      lat_max =     90.0          #  3.0 or 5.0 or 43.0 or -10.0;  90.0
   
#  Code assumes that -180. <= lon <= 180.  
      lon_min =           -180.0      # 130.0 for west Pacific ;  -150.0 for 'east'   ; -60.0 for eq Atlantic; 35.0 for Indian eq
      lon_max =            180.0      # should be -70.0 for east Pacific ; -90.0 for Nino3; 30.0 for eq Atlantic; 105.0 for Indian eq
      areaname= str(lat_min)+'-lat-'+str(lat_max) + '_' + str(lon_min)+'-lon-'+str(lon_max)  
                              # 135.0 137.0 west Pacific strip ; -84.0 -82.0 east Pacific strip ; -47.0 to -45.0 west Atl strip ; 40.0 to 35.0 west indian strip
                              # -150.0 to -90.0 Nino 3 ; -170.0 to  -120.0 Nino3.4; 160.0 to -150.0 Nino4 
                              # -90.0, -40.0 for West Atlantic Gulf Stream heat loss; 120.0, 180.0 for Kuroshio
			      # -81.0 to -80.0 east eq Pac;   8.0 to 9.0 east eq Atl 
			      # 40.0 to 100.0  Indian Ocean  
   
# define the output file for the time means
   outfile=open(dirout + c_model +'/' + c_model +'_' + c_descrip +'_' + c_period+'_'+ cvar +'_' + areaname + '_'+ c_type+'.txt','w')   

# if plot_mask: plots a geographical map of outmask 
   plot_mask = True    
   outmask='../std_out/'+ c_model  +'_' + c_descrip +'_' + cvar + '_' + areaname + '.png'

#----------------------------------------------------------------------------------------------
# read in static fields 

   cvar_static = cvar 
   if cvar == 'hfds' : cvar_static = 'tos' # because the land-sea mask is in the tos file 

   c_dir, c_descrip, c_period, c_end = cmip6_file_names(c_model, 'piControl', cvar_static)
   nav_lat, nav_lon, nav_area, land_mask = lat_lon_area_mask_CMIP6(cvar, c_model, c_descrip, c_period)
   c_dir, c_descrip, c_period, c_end = cmip6_file_names(c_model, c_expt_type, cvar)

###----------------------------------------------------------------------------------
# determine the area mask (in_area) and if c_type == 'mean' the total area (total_area)   

# 05 Jan 2023: ensure the dataset longitudes are in range -180 to 180; no longer left to the user to remember to do this correctly     
   nav_lon = np.where( nav_lon >  180., nav_lon - 360., nav_lon )  # 
   nav_lon = np.where( nav_lon < -180., nav_lon + 360., nav_lon )  #  added for GFDL output which is in the range -300 to 60 

   sea_point = land_mask > 0.5 
   
   if lat_lon_typ : 
      if lon_min > 180.0 or lon_max > 180. : 
         print ( ' user error set -180 < lon_min, lon_max < 180: lon_min, lon_max = ', lon_min, lon_max ) 
         sys.exit() 

      in_lat = np.logical_and( nav_lat > lat_min, nav_lat < lat_max) 

      if lon_max > lon_min : 
         in_lon = np.logical_and( nav_lon > lon_min, nav_lon < lon_max)
      else : 
         in_lon = np.logical_or( nav_lon > lon_min, nav_lon < lon_max)
      
      in_area = np.logical_and( in_lat, in_lon ) 
      in_area = np.logical_and( in_area, sea_point )  

   print ('shape(in_area) = ', np.shape(in_area) )
      
   in_area = np.where (in_area, 1.0, 0.0)    # probably not necessary ! convert in_area from a logical to a float 

   grd_pt_area  =       nav_area * in_area 
   total_area = np.sum(grd_pt_area) 
   print ( ' total_area = ', total_area )
   if total_area == 0.0 : 
      print (' total_area = 0.0: user error specifying area '  ) 
      sys.exit()

   print ( ' areaname = ', areaname, '; cvar = ', cvar, '; c_type = ', c_type, ' total_area = ', total_area, file = outfile) 

###----------------------------------------------------------------------------------
#  Read in the fields for the chosen variable

   file_var = c_dir + c_model +'_' + c_descrip +'_' + c_period  + c_end + '.nc'
   g_var = netCDF4.Dataset(file_var,'r')

# 23 Oct 2024 all fields previously read in at once outside the loop over imth. This required a lot of memory
#   field = g_var.variables[cvar]
#   nmonths, npj, npi = np.shape(field) 
#   print ( ' nmonths, npj, npi = ', nmonths, npj, npi ) 
      
   npj_lat, npi_lat = np.shape(nav_lat)  
   print ( ' c_model, cvar, npj_lat, npi_lat = ', c_model, cvar, npj_lat, npi_lat )

# 23 Oct 2024 this check moved inside the loop 
#   if npj_lat != npj or npi_lat != npi : 
#      print ( ' stopping because sizes do not match: npj_lat, npi_lat = ', npj_lat, npi_lat ) 
#      sys.exit()

###----------------------------------------------------------------------------------
# loop through the time-varying fields calculating and outputting the area mean or integral 
      
   for imth in range ( nmonths_ts ) :   

        field = g_var.variables[cvar][imth]
        npj, npi = np.shape(field) 

        if imth == 0 : 
          if npj_lat != npj or npi_lat != npi : 
             print ( ' stopping because sizes do not match: npj_lat, npi_lat = ', npj_lat, npi_lat ) 
             sys.exit()

        weighted_var = field * nav_area * in_area

        if c_type == 'int_pos': 
           weighted_var = np.where(weighted_var>0.0, weighted_var, 0.0) 
        elif c_type == 'int_neg': 
           weighted_var = np.where(weighted_var<0.0, weighted_var, 0.0) 

        var_output = np.sum(weighted_var)

        if c_type == 'mean':  
           var_output = var_output / total_area

        if imth%100 == 1 : print (' imth, var_output = ', imth, var_output ) 

        print ( imth, var_output , file = outfile )     # changed from iyr+1 on 4 Oct 2020 

###----------------------------------------------------------------------------------
#  Plot out  grd_pt_area

   if plot_mask : 
   
      projection=ccrs.PlateCarree()

# Set the figure size
      plt.figure(figsize=(8, 5))

# Declare the projection of the plot
      plt.gcf().add_subplot(111, projection=projection)

# Determine number of pixels in the subplot
      bbox = plt.gca().get_window_extent().transformed(plt.gcf().dpi_scale_trans.inverted())
      nx = int ( bbox.width * plt.gcf().get_dpi() )
      ny = int ( bbox.height * plt.gcf().get_dpi() )

      print ('nx, ny = ', nx, ny)
 
# Reproject the data onto a regular grid (with dimensions set by the number of pixels in the subplot, as above)
      x_extent = plt.gca().get_extent()[:2]
      y_extent = plt.gca().get_extent()[-2:]
      x, y = cartopy.img_transform.mesh_projection(projection, nx, ny, x_extents=x_extent, y_extents=y_extent)[:2]
      grd_pt_area_rgd = cartopy.img_transform.regrid(grd_pt_area, nav_lon, nav_lat, ccrs.PlateCarree(), projection, x, y)

#      colorlevels = np.linspace( -1, 1, 3)
      c = plt.contourf(x, y, grd_pt_area_rgd, cmap='coolwarm')  # levels = colorlevels ,

# Add coastlines, contour labels and a colour bar
      plt.gca().coastlines()
#     plt.clabel(c, inline=False, colors='k')
#     plt.colorbar(c, orientation='horizontal', extend='both')

# Add latitude and longitude labels and then matching gridlines
      plt.gca().set_xticks(np.linspace(-180, 180, 7), crs=projection)
      plt.gca().set_yticks(np.linspace(-90, 90, 7), crs=projection) 
      g1 = plt.gca().gridlines(color='k')
      g1.xlocator = mticker.FixedLocator( np.linspace(-180, 180, 7) )
      g1.ylocator = mticker.FixedLocator( np.linspace(-90,90, 7) )

      print ( 'mask plotted in file : ', outmask )   
      plt.savefig(outmask)
      plt.close()   
      
   
def run_field_area_mean_timeseries_CMIP6():

   c_model_list = ['ACCESS-CM2', 'CanESM5', 'GFDL-CM4', 'IPSL-CM6A-LR', 'MIROC6', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM' ]  
#   c_model_list = ['IPSL-CM6A-LR'] # , 'MIROC6']
   c_expt_type = 'piControl' # '4xCO2' or 'piControl'

   for c_model in c_model_list : 
      field_area_mean_timeseries_CMIP6(c_model, c_expt_type) 
      

run_field_area_mean_timeseries_CMIP6()
