# -*- coding: iso-8859-1 -*-
'''
field_area_mean_timeseries  

Date: 22 April 2016
Purpose: Calculate the area mean or the area integral of surface fields   
Author: Mike Bell   

'''

# 04 Oct 2020:  now writes iyr not iyr+1 as the year in each line of the time-series. This means all timeseries 
#               written earlier than 04 Oct 2020 should not be used with those written on or after 04 Oct 2020 

# 6 March 2023: extended script so that monthly time-series can start from the calendar month set by mth_cal_st
#               timeseries_lsq_fit.py is able to start reading from any point in a timeseries - so this extension was not needed

def field_area_mean_timeseries(c_model):

   import sys        
   import netCDF4
   import numpy as np
   import matplotlib.pyplot as plt
   import matplotlib.ticker as mticker
   import cartopy
   import cartopy.crs as ccrs
   import cartopy.img_transform
   from utilities import period_string, rd_ncdf_var_check_one
   from utilities_CMIP6 import cmip6_convert_int_to_str4, cmip6_int_to_str_mth

###----------------------------------------------------------------------------------
# Start of user inputs 

   data_type = 'CMIP6' # ocean_model, EN4, atmos_model, era5, DEEPC or GloSea5 or ORAS5 or CMIP6 

   expt_id = 'ar766'     # ak306 ar766 am120 bg555 bm758

   grid_type = 'T' # 'T' or 'U' or 'F'
   datadir = '/data/users/mike.bell/CMIP6HighResMIP' #     CMIP6HighResMIP or CMIP6AMIP or GC5
   varname='tos'  # tos, sos, hfds, empmr, tauuo, C, G, t20d so20chgt
   config_res = 'eORCA1-GO6'   # eORCA1-GO6 or eORCA025_GC5 or N96
   
   ll_read_mask = True

   mth_cal_st = '06'    #  using '03' for seas and  '01' for mth ;  needing to use '03' for aj368 and ay490
   cc_ann_seas_mth_mean = 'ann'    # 'ann' or 'seas' or 'mth' or 'one_seas'
   
   if cc_ann_seas_mth_mean == 'ann' : 
      y_m     = 'y'
      y_m_dir_in  = 'y_' + mth_cal_st
      y_m_dir_out = 'y_' + mth_cal_st
   elif cc_ann_seas_mth_mean == 'seas' : 
      y_m = 's' 
      y_m_dir_in  = 's'
      y_m_dir_out = 's_' + mth_cal_st 
   elif cc_ann_seas_mth_mean == 'mth'  : 
      y_m = 'm' 
      y_m_dir_in  = 'm'
      y_m_dir_out = '1m' #   was   'm_' + mth_cal_st 
   elif cc_ann_seas_mth_mean == 'one_seas' : 
      y_m = 's' 
      y_m_dir_in  = 's'
      y_m_dir_out = '1s_' + mth_cal_st 
   else : 
      print ( ' stopping: user error specifying cc_ann_seas_mth_mean = ', cc_ann_seas_mth_mean ) 
      sys.exit()

              
   if data_type == 'ocean_model' : 
      iyr_st  = 1850    #    1950 or 1850  1985 1979
      iyr_end = 2348    #    2048 or 2348  2016 2077
      inmesh =datadir+'/MESH_MASK/mesh_mask_'+config_res+'.nc'  
      folder=datadir+'/u-'+expt_id+'/on'+y_m_dir_in+'/'     
      infolder=folder+'surface_'+grid_type+'/'
      infilestub = 'nemo_'+expt_id+'o_1'+y_m+'_'
      outfilestub = expt_id+'_1'+y_m_dir_out+'_'+str(iyr_st)+'_'+str(iyr_end)+'_'
      outdir='../Time_series/u-'+expt_id+'/u-'
      grid_type_in = '_grid-' + grid_type

   elif data_type == 'atmos_model' :   # this option is only half written - because these timeseries not required yet 
      iyr_st  = 1960    #    1950 or 1850  1985
      iyr_end = 2176    #    2048 or 2348  2016
      inmesh = datadir+'/MESH_MASK/mesh_mask_'+config_res+'_'+grid_type+'.nc'  
      infolder=datadir+'/u-'+expt_id+'/an'+y_m_dir_in+'/surface_'+grid_type+'/'     
      infilestub = expt_id+'_1'+y_m+'_'
      outfilestub = expt_id+'_1'+y_m_dir_out+'_'+str(iyr_st)+'_'+str(iyr_end)+'_'
      outdir='../Time_series/u-'+expt_id+'/u-'
      grid_type_in = '' # '' or '_grid-' + grid_type

   elif data_type == 'era5' :
      iyr_st  = 1960    #    1960 or 1985
      iyr_end = 2018    #    2018 or 2016
      inmesh = '/data/users/mike.bell/ERA5_pm90/mesh_mask.nc'  
      infolder='/data/users/mike.bell/ERA5_pm90/on'+y_m_dir_in+'/' #    'on_'+y_m_dir_in+'/'    is optional
      infilestub='era5_1'+y_m + '_'  #  + varname + '_'  is optional 
      grid_type_in='_pm90'                         # _pm30 is optional
      outfilestub = 'era5_'+y_m_dir_out+'_'+str(iyr_st)+'_'+str(iyr_end)+'_'
      outdir='../Time_series/era5/'

   elif data_type == 'DEEPC' :
      iyr_st  = 1985    #      1985
      iyr_end = 2016    #      2016
      inmesh  = '/data/users/mike.bell/DEEPC_v05'+'/mesh_mask.nc'   
      infolder= '/data/users/mike.bell/DEEPC_v05/on'+y_m_dir_in+'/'    
      infilestub='DEEPC_v05_1'+y_m+'_'   # 
      grid_type_in=''                         # _pm30 is optional
      outfilestub = 'DEEPC_v05_'+y_m_dir_out+'_'+str(iyr_st)+'_'+str(iyr_end)+'_'
      outdir='../Time_series/DEEPC_v05/'

   elif data_type == 'EN4' : 
      iyr_st  = 1985    #    1985
      iyr_end = 2013    #    2016   or 2017 when building dhddt
      inmesh = '/data/users/mike.bell/EN4/mesh_mask.nc'  
      infolder='/data/users/mike.bell/EN4/on'+y_m_dir_in+'/2D/'     
      infilestub = 'EN.4.2.2.f.analysis.g10.' + '' #   + 't20d.' or + ''
      outfilestub = 'en4_'+y_m_dir_out+'_'+str(iyr_st)+'_'+str(iyr_end)+'_'
      outdir='../Time_series/en4/'
      grid_type_in = '' 

   elif data_type == 'GloSea5' :
      iyr_st  = 1993    #    1993
      iyr_end = 2013    #    2016
      inmesh  = '/data/users/mike.bell/GloSea5/mesh_mask.nc'   
      infolder= '/data/users/mike.bell/GloSea5/on'+y_m_dir_in+'/2D/'    
      if y_m == 'm' :
         infilestub='monthmean_'+varname+'.'   # 
      elif y_m == 'y' :
         infilestub='GloSea5_1y_'   # 
      grid_type_in=''                         # _pm30 is optional
      outfilestub = 'GloSea5_'+y_m_dir_out+'_'+str(iyr_st)+'_'+str(iyr_end)+'_'
      outdir='../Time_series/GloSea5/'

   elif data_type == 'ORAS5' :
      iyr_st  = 1985    #    1985
      iyr_end = 2013    #    2014
      inmesh  = '/data/users/mike.bell/GloSea5/mesh_mask.nc'            # this is correct - same grid as GloSea5
      infolder= '/data/users/mike.bell/ORAS5/on'+y_m_dir_in+'/surface_T/'    
      if y_m == 'm' :
         infilestub=varname+'_control_monthly_highres_2D_'   #  surface 20oC depth 
      elif y_m == 'y' :
         infilestub='ORAS5_1y_'   #  surface 20oC depth 
      grid_type_in='_CONS_v0.1'                         # _pm30 is optional
      outfilestub = 'ORAS5_'+y_m_dir_out+'_'+str(iyr_st)+'_'+str(iyr_end)+'_'
      outdir='../Time_series/ORAS5/'

   elif data_type == 'CMIP6' :
#   c_model provided as argument to the function
      c_expt   = 'abrupt-4xCO2'    #             piControl ssp245 1pctCO2  abrupt-4xCO2 ssp585 historical
      iyr_st  =    0    #    1993   0   165      10                            0  
      iyr_end =  100   #    2016  80   249     200                          100
      infolder='/data/users/mike.bell/CMIP6/'+c_model+'/1'+y_m_dir_in+'/'
      infilestub=varname+'_'+c_model+'_'+c_expt+'_'
      inmesh='/data/users/mike.bell/CMIP6/'+c_model+'/dec/'+'areacello_'+c_model+'_piControl.nc'
      outfilestub = c_model+'_'+c_expt+'_1'+y_m_dir_out+'_'+str(iyr_st)+'_'+str(iyr_end)+'_'
      outdir='../Time_series/CMIP6/'
      grid_type_in=''

   c_type = 'mean'   #  'integral', 'mean', 'int_pos' (integral of positive values), 'int_neg' ;   NOT 'anom' - that is created by another script
   
   lat_lon_typ = True # false option re-instated!

# define the region over which the mean or integral is calculated    ; 2.0 or 3.0 for 'tauuo' ; 5.0 or 10.0 for 'hfds' 
   if lat_lon_typ : 
      lat_min =     -5.0          # -3.0 or 5.0 or 27.0 or -40.0; -90.0 for global
      lat_max =      5.0          #  3.0 or 5.0 or 43.0 or -10.0;  90.0
   
#  Code assumes that -180. <= lon <= 180.  
      lon_min =    130.0      # 130.0 for west Pacific ;  -150.0 for 'east'   ; -60.0 for eq Atlantic; 35.0 for Indian eq
      lon_max =    -70.0      # should be -70.0 for east Pacific ; -90.0 for Nino3; 30.0 for eq Atlantic; 105.0 for Indian eq
      areaname= str(lat_min)+'-lat-'+str(lat_max) + '_' + str(lon_min)+'-lon-'+str(lon_max)  
                              # 135.0 137.0 west Pacific strip ; -84.0 -82.0 east Pacific strip ; -47.0 to -45.0 west Atl strip ; 40.0 to 35.0 west indian strip
                              # -150.0 to -90.0 Nino 3 ; -170.0 to  -120.0 Nino3.4; 160.0 to -150.0 Nino4 
                              # -90.0, -40.0 for West Atlantic Gulf Stream heat loss; 120.0, 180.0 for Kuroshio
			      # -81.0 to -80.0 east eq Pac;   8.0 to 9.0 east eq Atl 
			      # 40.0 to 100.0  Indian Ocean  

   else : 
 #    ii_max and jj_max ARE included in the area (they are not in the boxes) 
 #                       ne_atl_50 ;  NW Pac bdy t20d ; NW Atl bdy t20d ; SW Pac bdy ; SW Atl bdy;  SW Atl meet; N Guinea tail; Philippines
      ii_min = 80      # 245       ; 52               ; 208             ; 82         ; 240       ; 248         ; 79           ; 54
      ii_max = 150      # 294       ; 53               ; 209             ; 83         ; 241       ; 248         ; 80           ; 55
      jj_min = 236     # 260       ; 226              ; 232             ; 140        ; 140       ; 149         ; 160          ; 204
      jj_max = 248     # 309       ; 228              ; 234             ; 143        ; 143       ; 150         ; 163          ; 206
      areaname= str(jj_min)+'-jj-'+str(jj_max) + '_' + str(ii_min)+'-ii-'+str(ii_max)  
   
# define the output file for the time means
   outfilename = outdir+outfilestub+varname+'_'+areaname + '_'+ c_type+'.txt'
   outfile=open(outfilename,'w')   

# if plot_mask: plots a geographical map of outmask 
   plot_mask = False    
   outmask='../std_out/'+ infilestub + grid_type + '_' + areaname + '.png'

###----------------------------------------------------------------------------------
# End of user inputs 

#----------------------------------------------------------------------------------------------
# read in static fields 
   
   g = netCDF4.Dataset(inmesh,'r')
   
   nav_lat = rd_ncdf_var_check_one(g, 'nav_lat')
   nav_lon = rd_ncdf_var_check_one(g, 'nav_lon')

# 05 Jan 2023: ensure the dataset longitudes are in range -180 to 180; no longer left to the user to remember to do this correctly     
   nav_lon = np.where( nav_lon > 180., nav_lon - 360., nav_lon )  # 
   
   lowercg = grid_type.lower()
   ce1 = 'e1' + lowercg  
   ce2 = 'e2' + lowercg
   if data_type == 'CMIP6' :
      area = rd_ncdf_var_check_one(g, 'areacello')
   else : 
      e1 = rd_ncdf_var_check_one(g, ce1)
      e2 = rd_ncdf_var_check_one(g, ce2)

   if ll_read_mask and data_type != 'CMIP6' :
      cmask =  lowercg + 'mask' 
      mask = rd_ncdf_var_check_one(g, cmask)
      if mask.ndim == 3 :
         mask = mask[0]      # surface mask 
   else :
      mask = np.ones_like ( nav_lat ) 
   
   print ('shape(nav_lat) = ', np.shape(nav_lat) )
   print ('shape(nav_lon) = ', np.shape(nav_lon) )
   if data_type == 'CMIP6' :
      print ('shape(area)     = ', np.shape(area) )
   else : 
      print ('shape(e1)     = ', np.shape(e1) )
      print ('shape(e2)     = ', np.shape(e2) )
   print ('shape(mask)   = ', np.shape(mask) )

   
###----------------------------------------------------------------------------------
# determine the area mask (in_area) and if c_type == 'mean' the total area (total_area)   

   sea_point = mask > 0.5 
   
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

   else : 

      in_area = np.zeros_like( nav_lat ) 
      for jj in range( jj_min, jj_max + 1 ) : 
        if  ii_min < ii_max :
          for ii in range( ii_min, ii_max + 1 ) : 
            in_area[jj,ii] = sea_point[jj,ii]	  
        else : 
          for ii in range ( ii_min, npi ) : 
            in_area[jj,ii] = sea_point[jj,ii]	  
          for ii in range ( 0, ii_max ) : 
            in_area[jj,ii] = sea_point[jj,ii]	  

   print ('shape(in_area) = ', np.shape(in_area) )
      
   in_area = np.where (in_area, 1.0, 0.0)    # probably not necessary ! convert in_area from a logical to a float 

   if data_type != 'CMIP6' : area = e1 * e2

   grd_pt_area  =  area * in_area 
   total_area = np.sum(grd_pt_area) 
   print ( ' total_area = ', total_area )
   if total_area == 0.0 : 
      print (' total_area = 0.0: user error specifying area '  ) 
      sys.exit()

   print ( 'infilestub = ' , infilestub, '; areaname = ', areaname, '; varname = ', varname, '; c_type = ', c_type, ' total_area = ', total_area, file = outfile) 

###----------------------------------------------------------------------------------
# loop through the time-varying fields calculating and outputting the area mean or integral 

   for iyr in range ( iyr_st, iyr_end + 1  ) :   

     i_mth_cal_st = int(mth_cal_st)
     print ( 'iyr, i_mth_cal_st = ', iyr, i_mth_cal_st ) 
     n_month_end = 12 
     if cc_ann_seas_mth_mean == 'ann':
        nmonths = 12
     elif cc_ann_seas_mth_mean == 'seas' : 
        nmonths = 3
     elif cc_ann_seas_mth_mean == 'mth' : 
        nmonths = 1
     elif cc_ann_seas_mth_mean == 'one_seas' : 
       n_month_end = 3 
       nmonths = 3
      
     for imth in range ( i_mth_cal_st, i_mth_cal_st+n_month_end, nmonths ) :  

        date_str = period_string( iyr, imth, nmonths )
        if data_type == 'CMIP6' :
           start_month = int(mth_cal_st)
           yr_str = cmip6_convert_int_to_str4(iyr)
           if start_month == 1 : 
              yr_end = iyr
           else : 
              yr_end = iyr+1
           yr_end_str = cmip6_convert_int_to_str4(yr_end)

           end_month = (start_month + 10) % 12 + 1       # start_month =1 has end_month = 12 
           c_mth_end   = cmip6_int_to_str_mth( end_month)         
           date_str = yr_str+mth_cal_st+'-'+yr_end_str+c_mth_end

        if cc_ann_seas_mth_mean == 'mth' :
          if data_type == 'ORAS5' : 
            date_str = date_str[0:6] 
	
        infile=infolder+infilestub+date_str+grid_type_in+'.nc'
        print ('infile = ', infile)
	
        g = netCDF4.Dataset(infile,'r')
        var = rd_ncdf_var_check_one(g, varname)   
        g.close()

        print ('shape(var) = ', np.shape(var) )

        if data_type == 'CMIP6' :  
           ll_sea_point = np.abs(var) < 1.E30
           var = var * ll_sea_point

        weighted_var = var * area * in_area

        if c_type == 'int_pos': 
           weighted_var = np.where(weighted_var>0.0, weighted_var, 0.0) 
        elif c_type == 'int_neg': 
           weighted_var = np.where(weighted_var<0.0, weighted_var, 0.0) 

        var_output = np.sum(weighted_var)

        if c_type == 'mean':  
           var_output = var_output / total_area

        if cc_ann_seas_mth_mean == 'ann' : 
           print ( iyr, var_output , file = outfile )     # changed from iyr+1 on 4 Oct 2020 
        elif cc_ann_seas_mth_mean == 'seas' : 
           print ( 4*(iyr-iyr_st) + int(imth/3), var_output , file = outfile )
        elif cc_ann_seas_mth_mean == 'mth' : 
           print ( 12*(iyr-iyr_st) + imth, var_output , file = outfile )
        elif cc_ann_seas_mth_mean == 'one_seas' : 
           print ( iyr, var_output , file = outfile )

   print ( ' outfilename = ', outfilename ) 
###----------------------------------------------------------------------------------
#  Plot out in_mask 

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
      in_area_rgd = cartopy.img_transform.regrid(in_area, nav_lon, nav_lat, ccrs.PlateCarree(), projection, x, y)

      colorlevels = np.linspace( -1, 1, 3)
      c = plt.contourf(x, y, in_area_rgd, levels = colorlevels , cmap='coolwarm')

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

def run_field_area_mean_timeseries() :

# list for ssp585 including tauuo
#   c_model_list = [ 'GFDL-CM4', 'HadGEM3-GC31-LL', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'HadGEM3-GC31-MM' ]
# list for ssp585 not including tauuo
   c_model_list = ['GFDL-CM4', 'HadGEM3-GC31-LL', 'IPSL-CM6A-LR', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'ACCESS-CM2', 'CanESM5', 'HadGEM3-GC31-MM', 'MIROC6' ]   
#   c_model_list = ['HadGEM3-GC31-LL']

   for c_model in c_model_list :
      field_area_mean_timeseries(c_model)

run_field_area_mean_timeseries()      
