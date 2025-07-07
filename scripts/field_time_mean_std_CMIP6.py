
# -*- coding: iso-8859-1 -*-
'''
field_time_mean_std  

Date: 22 June 2020
Purpose: Calculate the time mean and standard deviation from a time-series of annual mean fields      
Author: Mike Bell   

'''

def field_time_mean_std_CMIP6(c_model, c_descrip, c_period, cvar, c_mth_start):

   import sys
   import netCDF4
   import numpy as np
   from utilities import rd_ncdf_var_check_one, create_2D_file, create_3D_file 

   dirin = '/data/users/frmk/CMIP6/'   # CMIP6AMIP or CMIP6HighResMIP

   ll_convert_Kelvin_to_Centigrade = False    #  Used to convert mean tos from Kelvin to Centigrade (important for calculation of standard devn)
   
   cyr_st  = '00'     # usually 1960 for era5  ;   1950 or 1980 for ocean_model ; 1985 for DEEPC;  1850 for ar766
   iyr_st = int(cyr_st) 
   iyr_end = 499   # the last year used IS iyr_end ; usually 2019 for era5; 2049 for ocean_model ; 2016 for DEEPC; 2348 for ar766

   ll_running_mean = True
   len_time_window = 5

  
   ll_omit_coords = False 

#----------------------------------------------------------------------------------------------

#  accomodate file naming anomalies
   if cvar == 'tauu' and c_model == 'GFDL-CM4'   : c_descrip = 'piControl_r1i1p1f1_gr1'    
   
   if not ll_running_mean :
      len_time_window = 1   # to ensure this is the case 
      str_tw = ''
   else : 
      str_tw = '_rf_'+ str(len_time_window) 
   
   file_in = cvar+'_'+ c_model +'_' + c_descrip +'_' + c_period + str_tw + '_00'+c_mth_start+'01-'+str(iyr_end)+c_mth_start+'01.nc'

   file_var = dirin+c_model+'/'+ file_in 

   g_var = netCDF4.Dataset(file_var,'r')
   nav_lon =   rd_ncdf_var_check_one (g_var, 'nav_lon')         
   nav_lat =   rd_ncdf_var_check_one (g_var, 'nav_lat')         
   fields = g_var.variables[cvar][:]
   g_var.close()

   print ( ' np.shape(fields)  = ', np.shape(fields) ) 
   npt, npj, npi = np.shape(fields)
   
   depth = 0.0 
   
   sum_var    = np.zeros_like( nav_lat ) 
   sum_var_sq = np.zeros_like( nav_lat ) 

   for iyr in range( npt ) :             

      var = fields[iyr]
      if ll_convert_Kelvin_to_Centigrade :
         if varname == 'tos'  : 
            var = var - 273.15   # original data are in Kelvin 
    
      sum_var = sum_var + var  
      sum_var_sq = sum_var_sq + var*var 

   r_len = 1.0 / float(npt) 
   mean_var = sum_var * r_len 
   mean_var_sq = sum_var_sq * r_len 
  
   variance = mean_var_sq - mean_var * mean_var 
   std_devn = np.sqrt(variance) 

   file_mean = dirin+c_model+'/'+ 'mean_'+ file_in 
   print ( 'outputing file_mean = ', file_mean)
   create_2D_file ( nav_lat, nav_lon, depth, mean_var, 'mean_'+cvar, file_mean, ll_omit_coords )  

   file_stdvn = dirin+c_model+'/'+ 'std_'+ file_in 
   create_2D_file ( nav_lat, nav_lon, depth, std_devn, 'std_'+cvar, file_stdvn, ll_omit_coords )  
         
def run_field_time_mean_std():

   model_time = []
#   model_time.append( ['ACCESS-CM2','piControl_r1i1p1f1_gn', '095001-144912', 6000] ) 
#   model_time.append( ['CanESM5','piControl_r1i1p1f1_gn', '520101-620012', 12000] )
#   model_time.append( ['GFDL-CM4','piControl_r1i1p1f1_gn', '015101-065012', 6000 ] )
#   model_time.append( ['IPSL-CM6A-LR','piControl_r1i1p1f1_gr', '185001-384912', 18000 ] )    # tau lat, lon is missing 
   model_time.append( ['MIROC6','piControl_r1i1p1f1_gn', '320001-399912', 6000 ] )
#   model_time.append( ['MPI-ESM1-2-LR','piControl_r1i1p1f1_gn', '185001-284912', 12000 ] )
#   model_time.append( ['NorESM2-LM','piControl_r1i1p1f1_gn', '160001-210012', 6000] )     # actually has 6012 
#   model_time.append( ['HadGEM3-GC31-LL','piControl_r1i1p1f1_gn', '225001-384912'] )
#   model_time.append( ['HadGEM3-GC31-MM','piControl_r1i1p1f1_gn', '185001-234912'] )

   for clist in model_time : 

      c_model   = clist[0] 
      c_descrip = clist[1]
      c_period  = clist[2]

      for cvar in [ 'tos' ] : # 'tauu', 'hfds', 'tos' ] : 

         c_mth_cal_strt_list = ['06']

         for c_mth_cal_strt in c_mth_cal_strt_list : 
            field_time_mean_std_CMIP6(c_model, c_descrip, c_period, cvar, c_mth_cal_strt)

run_field_time_mean_std()    

