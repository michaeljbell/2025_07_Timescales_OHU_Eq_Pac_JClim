# -*- coding: iso-8859-1 -*-
'''
scatter_plot_CMIP6_decadal_means  

Date: 10 March 2025
Purpose: scatter plots of area means of CMIP6 decadal means for expt minus control ; 
         the user is able to choose the variable, area, expt and control for each of the x-axis and y-axis  
	 and can aggregrate the decadal means - or plot each of them too. 
	 This version allows the user to choose to calculate the ratio of each choice with the global SST change   

Author: Mike Bell   

'''

def area_mean_CMIP6(mean_type, nav_lat, nav_lon, field, cellarea, lat_min, lat_max, lon_min, lon_max, file_std) : 

   import numpy as np

   ll_sea_point = np.abs(field) < 1.E30    # missing values are actually 1.E36
   
   in_lat = np.logical_and( nav_lat > lat_min, nav_lat < lat_max) 

   if lon_max > lon_min : 
      in_lon = np.logical_and( nav_lon > lon_min, nav_lon < lon_max)
   else : 
      in_lon = np.logical_or( nav_lon > lon_min, nav_lon < lon_max)
      
   in_area = np.logical_and( in_lat, in_lon ) 
   in_area = np.logical_and( in_area, ll_sea_point )  

   grd_pt_area  = cellarea * in_area 
   total_area   = np.sum(grd_pt_area) 

   weighted_var = field * cellarea * in_area

   var_output = np.sum(weighted_var)

   if mean_type == 'mean':  
      var_output = var_output / total_area
      
   print ( ' total_area, var_output = ', total_area, var_output, file = file_std)  

   return var_output, in_area

###----------------------------------------------------------------------------------
#  Plot out  grd_pt_area

def plot_mask ( nav_lon, nav_lat, grd_pt_area, file_out) :

# nav_lat           latitudes of grid (2D numpy array)
# nav_lon           longitudes
# grid_pt_area      0 outside area; 1 inside area 
# file_out          name of output file 

   import matplotlib.pyplot as plt
   import matplotlib.ticker as mticker
   import cartopy
   import cartopy.crs as ccrs
   import cartopy.img_transform
   import numpy as np

   level_diag = 10    
   
# projection used to plot the field    
   projection=ccrs.PlateCarree()

# Set the figure size
   plt.figure(figsize=(8, 5))

# Declare the projection of the plot
   plt.gcf().add_subplot(111, projection=projection)

# Determine number of pixels in the subplot
   bbox = plt.gca().get_window_extent().transformed(plt.gcf().dpi_scale_trans.inverted())
   nx = int ( bbox.width * plt.gcf().get_dpi() )
   ny = int ( bbox.height * plt.gcf().get_dpi() )

   if level_diag > 5 :   print ('nx, ny = ', nx, ny)
 
# Reproject the data onto a regular grid (with dimensions set by the number of pixels in the subplot, as above)
   x_extent = plt.gca().get_extent()[:2]
   y_extent = plt.gca().get_extent()[-2:]
   x, y = cartopy.img_transform.mesh_projection(projection, nx, ny, x_extents=x_extent, y_extents=y_extent)[:2]
   grd_pt_area_rgd = cartopy.img_transform.regrid(grd_pt_area, nav_lon, nav_lat, ccrs.PlateCarree(), projection, x, y)

#      colorlevels = np.linspace( -1, 1, 3)
   c = plt.contourf(x, y, grd_pt_area_rgd, cmap='coolwarm')  # levels = colorlevels ,

# Add coastlines, contour labels and a colour bar
   plt.gca().coastlines()
#  plt.clabel(c, inline=False, colors='k')
#  plt.colorbar(c, orientation='horizontal', extend='both')

# Add latitude and longitude labels and then matching gridlines
   plt.gca().set_xticks(np.linspace(-180, 180, 7), crs=projection)
   plt.gca().set_yticks(np.linspace(-90, 90, 7), crs=projection) 
   g1 = plt.gca().gridlines(color='k')
   g1.xlocator = mticker.FixedLocator( np.linspace(-180, 180, 7) )
   g1.ylocator = mticker.FixedLocator( np.linspace(-90,90, 7) )

   if level_diag > 5 : print ( 'mask plotted in file : ', file_out )   
   plt.savefig(file_out)
   plt.close()   
   return
   
#-------------------------------------------------------------------------------------------------------   

def scatter_plot_CMIP6_decadal_means():

   import sys
   import math         
   import numpy as np
   import matplotlib.pyplot as plt
   import numpy.linalg as linalg
   import netCDF4
   from utilities import rd_ncdf_var_check_one
   from utilities_CMIP6 import cmip6_convert_int_to_str4

# 1. User choices 

   infolder='/data/users/mike.bell/CMIP6/'
   outfolder = '../scatter_plots/'

# 1A. Specify which models to use. The data available depend on the model and variable  

# list for historical and ssp245 not including tauuo
   CMIP6_list = ['ACCESS-CM2' ] #  ['CanESM5', 'GFDL-CM4', 'HadGEM3-GC31-LL', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'ACCESS-CM2' ]  # 'MIROC6',   

# list for ssp245 including tauuo
#   CMIP6_list = ['CanESM5', 'GFDL-CM4', 'HadGEM3-GC31-LL', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'NorESM2-LM' ]  # 'MIROC6',   

# list for abrupt-4xCO2   
#   CMIP6_list = ['GFDL-CM4', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'ACCESS-CM2', 'CanESM5' ]  # 'MIROC6',  'NorESM2-LM' - first year is missing from piControl  

# list for ssp585 not including tauuo
#   CMIP6_list = ['GFDL-CM4', 'HadGEM3-GC31-LL', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'ACCESS-CM2', 'CanESM5', 'HadGEM3-GC31-MM' ] # 'MIROC6',   

# list for ssp585 including tauuo
#   CMIP6_list = ['GFDL-CM4', 'HadGEM3-GC31-LL', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'HadGEM3-GC31-MM' ] #  'NorESM2-LM', 'ACCESS-CM2', 'CanESM5', 'MIROC6',   


# 1B. Further definition of what is to be plotted
   
# we can combine decades into periods of 20, 30, 50 years etc. 
   number_periods = 1   # typically between 1 and 10 
   
#  Define the area means to be plotted; item 0 is plotted on x-axis; item 1 on y-axis 

   dec_start_list =   [ 190,          60          ]     # 60 for 1pctCO2, abrupt-4xCO2 ;  190 for ssp245; 200 for ssp585  ; 0 for 'abrupt-4xCO2' ; 140 for 'historical'
   expt_list      =   [ 'ssp245',     '1pctCO2'   ]
   cntl_list      =   [ 'piControl',  'piControl' ] 
   mean_period_list = [ 40,           40          ]     # must be a multiple of 10 
   
# Set colour cycle and symbol cycle 
   color_cycle=['tab:green', 'tab:purple', 'tab:olive', 'tab:red', 'tab:cyan', 'tab:brown', 'tab:orange' ]
   color_cycle=['tab:orange']
   ncolors = 7
   marker_cycle=['o', 'v', 'x']
   nmarkers= 3
   icycle = -1      # initialising 

# Plotting choices
   plot_legend = True #  legend of CMIP6 
   ll_plot_area = False   #   plots area of first field read in in '../std_out/'+filename_out+'.png'

   ll_publication = False   # True for figure suitable for publication ; False includes axes and legend
   cc_markersize = 150      # size of colour markers on the plot 
   label_size = 16  # size of numbers on axes (20)  except 16 for wind stress 

# 1C. Choose variables and areas to be plotted

# choose what will be plotted. The x-axis is specified by the first element of each list; the y-axis by the second element
# One can use two quantities on each axis. This is controled by two_variable_list which can be 'no', 'ratio' or 'difference'
# two_variable_list can be used for example to construct sst Nino4 - sst Nino3 or to normalise a variable by the change in global SST.  
#
# The second quantity "2" is specified first. Then the first quantity. This makes it easy to normalise by the global SST 
# 

   cc_user_plot_type = 'global'   # 'global', 'Nino3', 'heat_heat_eq_Pac', 'wind_wind_stress', 'heat_wind_stress', 
#                                       # 'heat_tos_grad', 'tos_tos_grad', 'heat_global_eqPac', 'heat_heat_equatorial', 'heat_equatorial_global'

#   mean_type_1_list    = [ 'mean', 'mean' ]     # options are 'mean' 'integral' 
   mean_type_1_list    = [ 'integral', 'integral' ]     # options are 'mean', 'integral' 

   ll_ratio_global_mean_SST_change_list = [ False, False ]   # can choose either x-axis or y-axis to be True or False  
#   ll_ratio_global_mean_SST_change_list = [ True, True ]   # can choose either x-axis or y-axis to be True or False  
   
   ll_subtract_control_list = [ True, True ]  #  can choose either x-axis or y-axis to be True or False

# Best to set ll_plot_x_eq_y_line to False first to find suitable max and min 
   ll_plot_x_eq_y_line = True
   x_min =               0.5 #  -0.4  -0.001 # -4.0  
   x_max =               0.9 #   0.4 #   0.008 #  4.0

#   multiplication_factor = 1.0
   multiplication_factor = 1.0E-15

# 2.  The user should not normally need to alter the code beyond this point 

# 2A. Specify the second quantity 
# User is free to choose any of the options in the next line 

   two_variable_list  = [ 'no',   'no'   ]  #  ['no', 'no'] or ['difference', 'no'] or ['no, 'difference'] or ['difference', 'difference']

# the following elements are read but not used if the corresponding element in two_variable_list = 'no'. 
# The values are set for the Nino3 area 
   variable_2_list   = [  'tos',   'tos' ]        
   mean_type_2_list  = [ 'mean', 'mean'  ]  # options are 'mean', 'integral' 
   lon2_min_list     = [ -150.0, -150.0  ]    
   lon2_max_list     = [  -90.0,  -90.0  ] 
   lat2_min_list     = [   -5.0,   -5.0  ] 
   lat2_max_list     = [    5.0,    5.0  ] 

# 2B. Specify the first quantity - 

   if cc_user_plot_type == 'global' :
     cc_user_global_var = 'hfds'           # options are 'hfds' 'tos' 'tauuo'
     variable_1_list   = [ cc_user_global_var, cc_user_global_var ]     
     lon1_min_list = [  -180.0, -180.0 ] 
     lon1_max_list = [   180.0,  180.0 ] 
     lat1_min_list = [ -90.0, -90.0 ] 
     lat1_max_list = [  90.0,  90.0 ] 

   elif cc_user_plot_type == 'Nino3' :
     cc_user_global_var = 'tos'           # only used if cc_user_plot_type = 'global'
     variable_1_list   = [ cc_user_global_var, cc_user_global_var ]     # options are 'hfds' 'tos' 'tauuo' 
     lon1_min_list = [  -150.0, -150.0 ] 
     lon1_max_list = [   -90.0,  -90.0 ] 
     lat1_min_list = [    -5.0,  -5.0 ] 
     lat1_max_list = [     5.0,   5.0 ] 

   elif cc_user_plot_type == 'heat_heat_eq_Pac' : 
     variable_1_list   = [ 'hfds', 'hfds' ]     # options are 'hfds' 'tos' 'tauuo' 
     lon1_min_list = [  130.0, 130.0 ] 
     lon1_max_list = [  -70.0, -70.0 ] 
     lat1_min_list = [   -5.0,  -5.0 ] 
     lat1_max_list = [    5.0,   5.0 ] 

   elif cc_user_plot_type == 'wind_wind_stress' : 
     variable_1_list   = [ 'tauuo', 'tauuo' ]     # options are 'hfds' 'tos' 'tauuo' 
     lon1_min_list = [  130.0, 130.0 ] 
     lon1_max_list = [  -70.0, -70.0 ] 
     lat1_min_list = [   -3.0,  -3.0 ] 
     lat1_max_list = [    3.0,   3.0 ] 

   elif cc_user_plot_type == 'heat_wind_stress' : 
     variable_1_list   = [ 'tauuo', 'hfds' ]     # options are 'hfds' 'tos' 'tauuo' 
     variable_2_list   = [ 'tauuo', 'hfds' ]     # options are 'hfds' 'tos' 'tauuo'  
     lon1_min_list = [  130.0, 130.0 ] 
     lon1_max_list = [  -70.0, -70.0 ] 
     lat1_min_list = [   -3.0,  -5.0 ] 
     lat1_max_list = [    3.0,   5.0 ] 

   elif cc_user_plot_type == 'heat_tos_grad' : 
     variable_1_list   = [ 'tos',        'hfds']     # options are 'hfds' 'tos' 'tauuo' 
     lon1_min_list = [  160.0,  130.0 ] 
     lon1_max_list = [ -150.0,  -70.0 ] 
     lat1_min_list = [   -5.0,   -5.0 ] 
     lat1_max_list = [    5.0,    5.0 ] 

   elif cc_user_plot_type == 'tos_tos_grad' : 
     variable_1_list   = [ 'tos',        'tos']     # options are 'hfds' 'tos' 'tauuo' 
     lon1_min_list = [  160.0,  160.0 ] 
     lon1_max_list = [ -150.0, -150.0 ] 
     lat1_min_list = [   -5.0,   -5.0 ] 
     lat1_max_list = [    5.0,    5.0 ] 

   elif cc_user_plot_type == 'heat_global_eqPac' : 
     variable_1_list   = [ 'hfds', 'hfds' ]     # options are 'hfds' 'tos' 'tauuo' 
     lon1_min_list = [  -180.0, 130.0 ] 
     lon1_max_list = [   180.0, -70.0 ] 
     lat1_min_list = [   -90.0,  -5.0 ] 
     lat1_max_list = [    90.0,   5.0 ] 

   elif cc_user_plot_type == 'heat_heat_equatorial' : 
     variable_1_list   = [ 'hfds', 'hfds' ]     # options are 'hfds' 'tos' 'tauuo' 
     lon1_min_list = [  -180.0, -180.0 ] 
     lon1_max_list = [   180.0,  180.0 ] 
     lat1_min_list = [    -15.0,   -15.0 ] 
     lat1_max_list = [     15.0,    15.0 ] 

   elif cc_user_plot_type == 'heat_eqPac_equatorial' : 
     variable_1_list   = [ 'hfds', 'hfds' ]     # options are 'hfds' 'tos' 'tauuo' 
     lon1_min_list = [   130.0,   130.0 ] 
     lon1_max_list = [   -70.0,   -70.0 ] 
     lat1_min_list = [    -5.0,   -15.0 ] 
     lat1_max_list = [     5.0,    15.0 ] 

   elif cc_user_plot_type == 'heat_equatorial_global' : 
     variable_1_list   = [ 'hfds', 'hfds' ]     # options are 'hfds' 'tos' 'tauuo' 
     lon1_min_list = [  -180.0, -180.0 ] 
     lon1_max_list = [   180.0,  180.0 ] 
     lat1_min_list = [    -5.0,  -90.0 ] 
     lat1_max_list = [     5.0,   90.0 ] 

   else : 
      print ( 'user error: cc_user_plot_type = ', cc_user_plot_type, ' is not supported' ) 
      sys.exit()

# End of user inputs

# 3. Set output filename and axis labels

   title_axes = ['','']
   filename_out = ''
   for item in [0, 1] : 
      lon1_min = lon1_min_list[item]
      lon1_max = lon1_max_list[item]
      if lon1_max > 180.0 or lon1_max < -180. or lon1_min > 180.0 or lon1_min < -180. : 
         print ( ' user error set -180 < lon1_min, lon1_max < 180: item, lon1_min, lon1_max = ', item, lon1_min, lon1_max ) 
         sys.exit() 
      lon2_min = lon2_min_list[item]
      lon2_max = lon2_max_list[item]
      if lon2_max > 180.0 or lon2_max < -180. or lon2_min > 180.0 or lon2_min < -180. : 
         print ( ' user error set -180 < lon2_min, lon2_max < 180: item, lon2_min, lon2_max = ', item, lon2_min, lon2_max ) 
         sys.exit() 

      title_axes[item] = expt_list[item] 
      if ll_subtract_control_list[item] : 
         title_axes[item] = title_axes[item] +'_'+ cntl_list[item] 	 
      title_axes[item] = title_axes[item]+'_'+ str(dec_start_list[item]) +'_'+ str(mean_period_list[item])

      title_axes[item] = title_axes[item] +'_'+ variable_1_list[item] +'_'+ mean_type_1_list[item] 
      title_axes[item] = title_axes[item] +'_'+ str(lon1_min_list[item]) +'_'+ str(lon1_max_list[item]) 
      title_axes[item] = title_axes[item] +'_'+ str(lat1_min_list[item]) +'_'+ str(lat1_max_list[item])

      if two_variable_list[item] != 'no' : 
         title_axes[item] = title_axes[item] +'_'+ two_variable_list[item] +'_'+ variable_2_list[item] +'_'+ mean_type_2_list[item]
         title_axes[item] = title_axes[item] +'_'+ str(lon2_min_list[item]) +'_'+ str(lon2_max_list[item]) 
         title_axes[item] = title_axes[item] +'_'+ str(lat2_min_list[item]) +'_'+ str(lat2_max_list[item])

      if ll_ratio_global_mean_SST_change_list[item] : 
         title_axes[item] = title_axes[item] +'_ratio_glob_mean_SST_chng' 
      
      filename_out = filename_out + title_axes[item]
      if item == 0 : 
         filename_out = filename_out + '_'
         
   title_periods = str(number_periods) + '_' + str( multiplication_factor ) #    str(mean_period) +'_'+ str(number_periods)
   filename_out = filename_out +'_'+ title_periods  

   if item == 1 and ll_publication : 
      filename_out = filename_out + '_publication' 

   file_std = open('../std_out/'+filename_out+'.txt','w') 

# 4. Start the scatter plot 
   fig, ax = plt.subplots(figsize=(6,5))

   if ll_plot_x_eq_y_line : 
      print ( ' plotting x = y line, x_min, x_max = ', x_min, x_max ) 
      pt_min = [x_min, x_max]
      pt_max = [x_min, x_max]
      plt.plot( pt_min, pt_max, color = 'black', lw=1, ls='--' ) 

# 5. Start loop over the models read in lats, lons and cell areas  
   
   for c_model in CMIP6_list :   

      icycle = icycle + 1 
      if icycle > ncolors*nmarkers :
         print ( ' stopping because not enough markers or colours for the number of models plotted ' )
         print ( ' icycle, ncolors, nmarkers = ', icycle, ncolors, nmarkers )  
         sys.exit()
      cc_color  = color_cycle  [      icycle % ncolors    ] 
      cc_marker = marker_cycle [ int( icycle / ncolors )  ]

      file_in_areacello = infolder + c_model +'/dec/areacello_'+ c_model +'_piControl.nc'      
      g_cellarea = netCDF4.Dataset(file_in_areacello,'r')
      cellarea = rd_ncdf_var_check_one(g_cellarea, 'areacello')
      nav_lat  = rd_ncdf_var_check_one(g_cellarea, 'nav_lat')
      nav_lon  = rd_ncdf_var_check_one(g_cellarea, 'nav_lon')
      print (' c_model, np.shape(cellarea), np.shape(nav_lat), np.shape(nav_lon) = ', c_model, np.shape(cellarea), np.shape(nav_lat), np.shape(nav_lon) )
      g_cellarea.close()

# ensure the dataset longitudes are in range -180 to 180; left to the user to remember to do this correctly     
      nav_lon = np.where( nav_lon >  180., nav_lon - 360., nav_lon )  # 
      nav_lon = np.where( nav_lon < -180., nav_lon + 360., nav_lon )  #  added for GFDL output which is in the range -300 to 60 

# 6. Start to loop over the points that will be plotted for each model; item = 0 is for x-axis; item = 1 for y-axis

      pts_to_plot_x = []  
      pts_to_plot_y = []

      ll_first_field = True   # only used to decide whether to plot out the area selected 

      for iper in range( number_periods )  : 

         for item in [0, 1] :
            ll_ratio_global_mean_SST_change = ll_ratio_global_mean_SST_change_list[item]
            ll_subtract_control      = ll_subtract_control_list[item] 
            dec_start                = dec_start_list[item]
            mean_period              = mean_period_list[item]

            expt      = expt_list[item]
            cntl      = cntl_list[item]
            var1      = variable_1_list[item] 

            mean_typ1 = mean_type_1_list[item] 
            lon1_min   = lon1_min_list[item] 
            lon1_max   = lon1_max_list[item] 
            lat1_min   = lat1_min_list[item] 
            lat1_max   = lat1_max_list[item] 

            c_var2     = two_variable_list[item]
            var2       = variable_2_list[item]    
            mean_typ2 = mean_type_2_list[item] 
            lon2_min   = lon2_min_list[item] 
            lon2_max   = lon2_max_list[item] 
            lat2_min   = lat2_min_list[item] 
            lat2_max   = lat2_max_list[item] 

            if iper == 0 : 	    
               print ( 'expt, cntl, var1, mean_typ1, lon1_min, lon1_max,lat1_min, lat1_max = ', expt, cntl, var1, mean_typ1, lon1_min, lon1_max,lat1_min, lat1_max, file = file_std)
               print ( 'c_var2, var2, mean_typ2, lon2_min, lon2_max,lat2_min, lat2_max = ', c_var2, var2, mean_typ2, lon2_min, lon2_max,lat2_min, lat2_max, file = file_std ) 

# 7. Loop over the decadal means within this period; read in the fields for expt and cntl  

            dec_first = dec_start + iper * mean_period
            num_sub_periods = 0 
            for iyr in range ( dec_first, dec_first + mean_period, 10) : 
               num_sub_periods = num_sub_periods + 1	       
               cyr = cmip6_convert_int_to_str4(iyr) 
               cdate = cyr+'01-'+cyr[0:3]+'912'    # date string for the decadal mean file 
               print ('cdate = ', cdate, file = file_std ) 
               print ( 'c_model, item, expt, var1, iyr, cdate = ', c_model, item, expt, var1, iyr, cdate )  
               print ( 'c_model, item, expt, var1, iyr, cdate = ', c_model, item, expt, var1, iyr, cdate, file = file_std )  
               file_in_stub = infolder + c_model +'/dec/'+ var1 +'_'+ c_model +'_'

               file_in_expt = file_in_stub + expt +'_'+ cdate + '.nc'
               g_expt = netCDF4.Dataset(file_in_expt,'r')
               field_main  = rd_ncdf_var_check_one(g_expt, var1)
               g_expt.close()
               print (' np.shape(field_main) = ', np.shape(field_main), file = file_std )

               if ll_subtract_control : 
                  file_in_cntl = file_in_stub + cntl +'_'+ cdate + '.nc'
                  g_cntl = netCDF4.Dataset(file_in_cntl,'r')
                  field_cntl  = rd_ncdf_var_check_one(g_cntl, var1)       	       
                  field_main = field_main - field_cntl
                  g_cntl.close()
                  print (' np.shape(field_cntl) = ', np.shape(field_cntl), file = file_std )

               if c_var2 != 'no' : 
                  file_in_stub = infolder + c_model +'/dec/'+ var2 +'_'+ c_model +'_'

                  file_in_expt = file_in_stub + expt +'_'+ cdate + '.nc'
                  g_expt = netCDF4.Dataset(file_in_expt,'r')
                  field2_main  = rd_ncdf_var_check_one(g_expt, var2)
                  g_expt.close()

                  if ll_subtract_control : 
                     file_in_cntl = file_in_stub + cntl +'_'+ cdate + '.nc'
                     g_cntl = netCDF4.Dataset(file_in_cntl,'r')
                     field2_cntl  = rd_ncdf_var_check_one(g_cntl, var2)  
                     field2_main = field2_main - field2_cntl
                     g_cntl.close()

               if ll_ratio_global_mean_SST_change :
                  file_in_stub = infolder + c_model +'/dec/'+ 'tos' +'_'+ c_model +'_'

                  file_in_expt = file_in_stub + expt +'_'+ cdate + '.nc'
                  g_expt = netCDF4.Dataset(file_in_expt,'r')
                  tos_expt  = rd_ncdf_var_check_one(g_expt, 'tos')

                  file_in_cntl = file_in_stub + cntl +'_'+ cdate + '.nc'
                  g_cntl = netCDF4.Dataset(file_in_cntl,'r')
                  tos_cntl  = rd_ncdf_var_check_one(g_cntl, 'tos')
                  tos_change = tos_expt - tos_cntl

                  g_expt.close()
                  g_cntl.close()
                  print (' np.shape(tos_change) = ', np.shape(tos_change), file = file_std )
	        

# 8. Calculate the area means (or integrals), accumulate them and at the end of the loop divide by the number of sub-periods

               mean_main, in_area = area_mean_CMIP6(mean_typ1, nav_lat, nav_lon, field_main, cellarea, lat1_min, lat1_max, lon1_min, lon1_max, file_std) 
               print ( ' var1, mean_main = ', var1, mean_main, file = file_std )  

               if ll_plot_area and ll_first_field : 
                  file_mask = '../std_out/mask_'+c_model+'_'+title_axes[0] +'.png'
                  plot_mask ( nav_lon, nav_lat, in_area, file_mask) 
                  ll_first_field = False

               if iyr == dec_first : 
                  total_change = mean_main
               else : 
                  total_change = total_change + mean_main

               if c_var2 != 'no' : 
                  mean_main, in_area = area_mean_CMIP6(mean_typ2, nav_lat, nav_lon, field2_main, cellarea, lat2_min, lat2_max, lon2_min, lon2_max, file_std) 
                  print ( ' var2, mean_main = ', var2, mean_main, file = file_std )  

                  if iyr == dec_first : 
                     total_change2 = mean_main
                  else : 
                     total_change2 = total_change2 + mean_main

               if ll_ratio_global_mean_SST_change : 
                  mean_tos_change, in_area = area_mean_CMIP6('mean', nav_lat, nav_lon, tos_change, cellarea, -90.0, 90.0, -180.0, 180.0, file_std) 
                  if iyr == dec_first : 
                     total_change_tos = mean_tos_change
                  else : 
                     total_change_tos = total_change_tos + mean_tos_change
                  print ( ' tos: mean_tos_change, total_change_tos = ', mean_tos_change, total_change_tos, file = file_std )  
     	             
	              
# Note that following line should NOT depend on mean_type - because that is an area-mean type 
            out_change = total_change / float(num_sub_periods) 
            print ( 'item, num_sub_periods, total_change, out_change = ', item, num_sub_periods, total_change, out_change, file = file_std )

            if c_var2 != 'no' : 
               out_change2 = total_change2 / float(num_sub_periods) 
               print ( 'item, num_sub_periods, total_change2, out_change2 = ', item, num_sub_periods, total_change2, out_change2, file = file_std )

               if c_var2 == 'ratio' : 
                  out_change = out_change / out_change2      # should not be using variables for which out_change2 is near zero
               elif c_var2 == 'difference' : 
                  out_change = out_change - out_change2      
               print ( 'c_var2, out_change = ', c_var2, out_change, file = file_std ) 

            if ll_ratio_global_mean_SST_change : 
               out_change_tos = total_change_tos / float(num_sub_periods) 
               out_change = out_change / out_change_tos 
               print ( 'tos_change, out_change = ', tos_change, out_change, file = file_std ) 

# 9. Store the points (from the loop started in 4.) 

            out_change = out_change * multiplication_factor 
            if item == 0 :
               pts_to_plot_x.append(out_change) 
            else :
               pts_to_plot_y.append(out_change) 

# 10. print out the points for each c_model 

      print (  'pts_to_plot_x = ', pts_to_plot_x, file = file_std )    
      print (  'pts_to_plot_y = ', pts_to_plot_y, file = file_std )    
      ax.scatter(pts_to_plot_x,pts_to_plot_y, s = cc_markersize, color= cc_color, marker = cc_marker, label=c_model)


# 11. complete the plots 

   if plot_legend and ( not ll_publication ) : 
      plt.legend()
   else : 
      filename_out = filename_out + '_nolegend'

   plt.tick_params(axis='x', labelsize=label_size )      
   plt.tick_params(axis='y', labelsize=label_size )

   if not ll_publication : 
      plt.xlabel(title_axes[0], fontsize=6)
      plt.ylabel(title_axes[1], fontsize=6)
      plt.title(title_periods, fontsize=10)
   plt.savefig(outfolder+filename_out+'.png') 
   plt.close('all') 
     
scatter_plot_CMIP6_decadal_means()
