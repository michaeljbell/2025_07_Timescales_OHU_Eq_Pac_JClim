# -*- coding: iso-8859-1 -*-
'''
field_timeseries_correlation 

Date: 22 June 2020
Purpose: Calculate fields of correlations with a single time-series. Uses time-mean and std devn fields as inputs (in addition to the time-series of the fields)     
Author: Mike Bell   

'''

# 11 Nov 2020 addition of 1 to the year removed (consistent with field_area_mean_timeseries.py) 

# 23 Oct 2022 revised so that running mean of the single timeseries is calculated within this function 

def field_timeseries_correlation(c_model, c_descrip, c_period, cvar, c_mth_start):

   import sys
   import netCDF4
   import numpy as np
   import math
   from utilities import create_2D_file, rd_ncdf_var_check_one
   from timeseries_read_write import timeseries_read
   from timeseries_running_mean import timeseries_running_mean

   dirin = '/data/users/frmk/CMIP6/'   # CMIP6AMIP or CMIP6HighResMIP
      
   ll_running_mean = True
   len_time_window = 5 

   cvar_ts     = 'hfds' 
   areaname_ts = '-5.0-lat-5.0_130.0-lon--70.0'
   meantype_ts = '_integral_ann_mean_'   # '_mean_ann_mean_' '_integral_ann_mean_'
   
   cyr_st  = '00'     # usually 1960 for era5  ;   1950 or 1980 for model_expt
   iyr_st = int(cyr_st) 
   iyr_end = 498   # the last year used IS iyr_end ; usually 2019 for era5; 2049 for model_expt 

# calculation of running mean start and end period - this must be consistent wtih field_running_mean 
#   iyr_st_fields  = iyr_st  + int(len_time_window/2)
#   iyr_end_fields = iyr_end + int(len_time_window/2) - ( len_time_window - 1)    
#   print (' c_model, iyr_st_fields, iyr_end_fields = ',  c_model, iyr_st_fields, iyr_end_fields )

   ll_omit_coords = False 

#----------------------------------------------------------------------------------------------

   if cvar == 'tauu' and c_model == 'GFDL-CM4'   : 
      c_descrip_fld = 'piControl_r1i1p1f1_gr1'    
   else :
      c_descrip_fld = c_descrip    

   if not ll_running_mean :
      len_time_window = 1   # to ensure this is the case 
      str_tw = ''
   else : 
      str_tw = '_rf_'+ str(len_time_window) 
   
   file_stub = c_model +'_' + c_descrip_fld +'_' + c_period 
   c_years = '_'+cyr_st+c_mth_start+'01-'+str(iyr_end+1)+c_mth_start+'01'

   in_file_var    = dirin+c_model+'/'     + cvar+'_'+ file_stub + str_tw + c_years +'.nc'
   in_file_mean   = dirin+c_model+'/mean_'+ cvar+'_'+ file_stub + str_tw + c_years +'.nc'
   in_file_stdevn = dirin+c_model+'/std_' + cvar+'_'+ file_stub + str_tw + c_years +'.nc'

   file_stub_ts = c_model +'_' + c_descrip +'_' + c_period 
   infile_ts = '../Time_series/'+c_model+'/' + file_stub_ts + '_' + cvar_ts + '_' + areaname_ts + meantype_ts + c_mth_start + '.txt'

   out_file = dirin+c_model+'/corrln_' + cvar + '_' + file_stub + str_tw + c_years + '_' + cvar_ts  + '_' + areaname_ts + '.nc' 

   print ( ' out_file  = ', out_file ) 
   file_std = open( '../std_out/' + file_stub + '_.txt','w') 
      
   ll_ann_means = True  # working with annual means 
   imth = 0 

#----------------------------------------------------------------------------------------------

# Read in the fields 

   g = netCDF4.Dataset(in_file_mean,'r')
   mean_fld = rd_ncdf_var_check_one (g, 'mean_'+cvar)         
   g.close()

   g = netCDF4.Dataset(in_file_stdevn,'r')
   std_fld = rd_ncdf_var_check_one (g, 'std_'+cvar)         
   g.close()

   g = netCDF4.Dataset(in_file_var,'r')
   nav_lat = rd_ncdf_var_check_one (g, 'nav_lat') 
   nav_lon = rd_ncdf_var_check_one (g, 'nav_lon') 
   fields  = rd_ncdf_var_check_one (g,cvar) 
   g.close()
   
   print ('np.shape(mean_fld), np.shape(fields) = ', np.shape(mean_fld), np.shape(fields) )
   npt, npj, npi = np.shape(fields)

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
   
   if npt != len(anom_ts) : 
     'stopping because lengths of input fields and timeseries do not match' 
     sys.exit() 

   covariance = np.zeros( (npj, npi) ) 

# Loop over input fields calculating the covariance 

   for iyr in range( npt ) :             
    
      anom_fld = fields[iyr] - mean_fld

      covariance = covariance + anom_fld * anom_ts[iyr] 

   denom =  float(npt) * std_ts * std_fld 

   correlation = np.where( denom > 0.0, covariance / denom, 0.0 )   

   depth = 0.0
   create_2D_file ( nav_lat, nav_lon, depth, correlation, 'corrln_'+cvar+'_'+cvar_ts, out_file, ll_omit_coords )  
         

def run_field_timeseries_correlation():

   model_time = []
   model_time.append( ['ACCESS-CM2','piControl_r1i1p1f1_gn', '095001-144912', 6000] ) 
   model_time.append( ['CanESM5','piControl_r1i1p1f1_gn', '520101-620012', 12000] )
   model_time.append( ['GFDL-CM4','piControl_r1i1p1f1_gn', '015101-065012', 6000 ] )
#   model_time.append( ['IPSL-CM6A-LR','piControl_r1i1p1f1_gr', '185001-384912', 18000 ] )    # tau lat, lon is missing ; hfds timeseries is shorter
#   model_time.append( ['MIROC6','piControl_r1i1p1f1_gn', '320001-369912', 6000 ] )           # hfds timeseries is shorter
   model_time.append( ['MPI-ESM1-2-LR','piControl_r1i1p1f1_gn', '185001-284912', 12000 ] )
   model_time.append( ['NorESM2-LM','piControl_r1i1p1f1_gn', '160001-210012', 6000] )     # actually has 6012 
   model_time.append( ['HadGEM3-GC31-LL','piControl_r1i1p1f1_gn', '225001-384912'] )
   model_time.append( ['HadGEM3-GC31-MM','piControl_r1i1p1f1_gn', '185001-234912'] )

   for clist in model_time : 

      c_model   = clist[0] 
      c_descrip = clist[1]
      c_period  = clist[2]

      for cvar in [ 'tos'] : # , 'hfds', 'tos' ] : # 'tauu', 'hfds', 'tos' ] : 

         c_mth_cal_strt_list = ['06']

         for c_mth_cal_strt in c_mth_cal_strt_list : 
            field_timeseries_correlation(c_model, c_descrip, c_period, cvar, c_mth_cal_strt)

run_field_timeseries_correlation()    


