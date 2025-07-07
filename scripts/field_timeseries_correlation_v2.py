# -*- coding: iso-8859-1 -*-
'''
field_timeseries_correlation 

Date: 22 June 2020
Purpose: Calculate fields of correlations with a single time-series. Uses time-mean and std devn fields as inputs (in addition to the time-series of the fields)     
Author: Mike Bell   

'''

# 11 Nov 2020 addition of 1 to the year removed (consistent with field_area_mean_timeseries.py) 

# 23 Oct 2022 revised so that running mean of the single timeseries is calculated within this function 

def field_timeseries_correlation():

   import netCDF4
   import numpy as np
   import math
   from utilities import period_string, rd_ncdf_var_check_one, create_2D_file
   from timeseries_read_write import timeseries_read
   from timeseries_running_mean import timeseries_running_mean
   
   datadir = '/data/users/frmk/CMIP6HighResMIP' #     CMIP6HighResMIP or CMIP6AMIP 
   data_type = 'model_expt' # 'model_expt' or 'era5' or 'DEEPC_v05' or 'en4'
   grid_type = 'T'  # 'T' or 'U' or 'F'
   varname_fld='tos'  # tos, sos, hfds, empmr, tauuo, G 

   mth_cal_st = '06'    #  using '01' for both seas and mth ;  needing to use '03' for aj368 and ay490 
   cc_ann_seas_mth_mean = 'y'   # 'y'      #   's' or 'm' are not supported yet
   
   y_m_dir = cc_ann_seas_mth_mean + '_' + mth_cal_st

   if data_type == 'model_expt' : 
      expt_id = 'ak306'   #   am120 or ak306
      config_res = 'eORCA1'
      inmesh = datadir+'/MESH_MASK/mesh_mask_'+config_res+'-GO6.nc'  # _GC5 or -GO6
      infolder=datadir+'/u-'+expt_id                   
      filestub='/surface_'+grid_type+'/nemo_'+expt_id+'o_1'+ cc_ann_seas_mth_mean + '_'
      fileend='.nc'
   elif data_type == 'era5' :
      infolder='/data/users/frmk/ERA5_pm90'    # _pm30 is optional
      inmesh = infolder+'/mesh_mask.nc'  
      filestub='/era5_1y_' #  +varname_fld +'_'
      fileend='_pm90.nc'                       # _pm30 is optional
   elif data_type == 'DEEPC_v05' :
      infolder='/data/users/frmk/DEEPC_v05'   
      inmesh = infolder+'/mesh_mask.nc'  
      filestub='/DEEPC_v05_1y_' #   +varname_fld +'_'
      fileend='.nc'                       # _pm30 is optional
   elif data_type == 'en4' :
      infolder='/data/users/frmk/EN4'    # _pm30 is optional
      inmesh = infolder+'/mesh_mask.nc'  
      filestub='/2D/EN.4.2.2.f.analysis.g10.' #    t20d.'
      fileend='.nc'                       # 
   
   ll_running_mean = False
   len_time_window = 1 

   if not ll_running_mean : len_time_window = 1   # to ensure this is the case   

   data_type_ts = 'model_expt' # 'model_expt' or 'era5' or 'DEEPC_v05' or 'en4'
   varname_ts = 'tos' 
   areaname= '35.0-lat-45.0_-45.0-lon--30.0'
   
   iyr_st  = 1950   # usually 1960 for era5  1985 for DEEP-C;   1950 or 1980 for model_expt
   iyr_end = 2048   # the last year used IS iyr_end ; usually 2016 for DEEP-C; 2018 for era5; 2048 for model_expt 

# calculation of running mean start and end period - this must be consistent wtih field_running_mean 
   iyr_st_fields  = iyr_st  + int(len_time_window/2)
   iyr_end_fields = iyr_end + int(len_time_window/2) - ( len_time_window - 1) 

   if data_type in [ 'model_expt', 'DEEPC_v05', 'era5', 'en4' ] : 
      period_str = str(iyr_st_fields)+mth_cal_st+'01-'+str(iyr_end_fields+1)+mth_cal_st+'01'    # changed iyr_end+1 to iyr_end 23 Oct 2022 
                                                                                 # need iyr_end for running mean fields?? 
   else :
      period_str = str(iyr_st_fields)+'_'+str(iyr_end_fields)

   ll_omit_coords = False 

#----------------------------------------------------------------------------------------------

   ll_ann_means = True  # working with annual means 
   imth = int(mth_cal_st) 
   nmonths = 12 

   if not ll_running_mean :
      onyrf = '/on' + y_m_dir  
      str_run_var = ''
      str_tw = ''
   else : 
      onyrf = '/on' + y_m_dir + '_rf' 
      str_run_var = varname_fld + '_' + str(len_time_window) + '_'
      str_tw = '_tw_'+ str(len_time_window) 

   infilestub = infolder + onyrf + filestub
   name_out = 'corrln_' + varname_fld + str_tw + '_' + varname_ts  + '_' + areaname + '_' + period_str
   out_file = infilestub + name_out + fileend
   print ( 'out_file = ', out_file ) 
   file_std = open( '../std_out/' + name_out + '_.txt','w') 

   if data_type_ts == 'model_expt' :
      infile_ts = '../Time_series/u-'+expt_id+'/u-'+expt_id+'_1y_' + mth_cal_st + '_' + str(iyr_st)+'_'+str(iyr_end)+'_'+varname_ts+'_'+areaname + '_mean.txt'

   else: 
      infile_ts = '../Time_series/'+data_type_ts+'/'+data_type_ts+'_1y_' + mth_cal_st + '_' + str(iyr_st)+'_'+str(iyr_end)+'_'+varname_ts+'_'+areaname + '_mean.txt'
      

#----------------------------------------------------------------------------------------------

# Mesh-mask file used 
   g = netCDF4.Dataset(inmesh,'r')
   cvar_dep = 'gdept_1d' 
   if data_type == 'en4' :  cvar_dep = 'depth'  
   depth   = g.variables[cvar_dep][0] 
   nav_lat = rd_ncdf_var_check_one(g, 'nav_lat')
   nav_lon = rd_ncdf_var_check_one(g, 'nav_lon')
   g.close()
   (npj, npi) = np.shape(nav_lat)
   print ('shape(nav_lat) = ', np.shape(nav_lat) )

   if data_type in [ 'model_expt', 'en4' ] :
     depth = depth[0]

# Read in the mean and std deviation fields 

   file_mean = infilestub + str_run_var + 'mean_'+ varname_fld + '_' + period_str + '_grid-' +grid_type+ fileend
   g = netCDF4.Dataset(file_mean,'r') 
   mean_fld = rd_ncdf_var_check_one(g, 'mean_'+varname_fld) 
   g.close()

   file_stdevn = infilestub + str_run_var + 'std_' + varname_fld +'_' +  period_str + '_grid-' +grid_type+ fileend 
   g = netCDF4.Dataset(file_stdevn,'r')
   std_fld = rd_ncdf_var_check_one(g,'std_'+varname_fld) 
   g.close()

   print ('shape(mean_fld) = ', np.shape(mean_fld) )

# Read in the timeseries for correlation. 

# Change 11 Nov 2020: The year written in the time-series for fields value on year iyr is now iyr . 
# Change 11 Nov 2020: So we DO NOT add 1 to iyr_st and iyr_end as we pass them into timeseries_read  

   full_ts  = timeseries_read(infile_ts,  iyr_st, iyr_end, file_std ) 

   if ll_running_mean : 
      full_ts = timeseries_running_mean(full_ts, len_time_window)

   mean_ts = np.mean(full_ts) 
   variance_ts  =  np.var(full_ts)             
   std_ts  = math.sqrt(variance_ts)
   print (' std_ts = ', std_ts ) 

   anom_ts = full_ts - mean_ts
   
   print ('dimensions of anom_ts = ', np.shape(anom_ts) )  

   covariance = np.zeros( (npj, npi) ) 

# Loop over input fields calculating the covariance 

   for iyr in range( iyr_st_fields, iyr_end_fields+1 ) :             
      if data_type in [ 'model_expt', 'DEEPC_v05', 'era5', 'en4' ] : 
         print ( ' iyr, imth = ', iyr, imth ) 
         date_str = period_string(iyr, imth, nmonths ) 
      else :
         date_str = str(iyr)
   
      if not ll_running_mean : 
         infile=infilestub+date_str
      else : 
         infile=infilestub+str_run_var+date_str 

      if data_type == 'model_expt' : 
         infile = infile + '_grid-'+grid_type+'.nc' 
      else :
         infile = infile + fileend

      if iyr == iyr_st_fields: print ('infile = ', infile)
	    
      g = netCDF4.Dataset(infile,'r')
      full_field = rd_ncdf_var_check_one(g,varname_fld) 
      g.close()

      if iyr == iyr_st_fields: print ( ' np.shape(full_field) = ',  np.shape(full_field) ) 
          
      anom_fld = full_field - mean_fld

      covariance = covariance + anom_fld * anom_ts[iyr - iyr_st_fields] 

   denom =  float ( iyr_end_fields + 1 - iyr_st_fields) * std_ts * std_fld 

   correlation = np.where( denom > 0.0, covariance / denom, 0.0 )   

   create_2D_file ( nav_lat, nav_lon, depth, correlation, 'corrln_'+varname_fld+'_'+varname_ts, out_file, ll_omit_coords )  

   print ( 'out_file = ', out_file ) 
         
field_timeseries_correlation()
