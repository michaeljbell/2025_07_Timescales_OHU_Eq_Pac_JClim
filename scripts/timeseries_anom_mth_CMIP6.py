# -*- coding: iso-8859-1 -*-
'''
timeseries_anom_mth_CMIP6  

Date: 7 Sept 2020
Purpose: Calculate the timeseries of monthly anomalies from a timeseries of monthly values      
Author: Mike Bell   

'''

def timeseries_anom_mth_CMIP6(c_model, c_expt_type, n_times):

   import sys
   import math         
   import numpy as np
   from timeseries_read_write import timeseries_read, timeseries_write
   from utilities_CMIP6 import cmip6_file_names
   import numpy.linalg as linalg

   ll_subtract_linear_trend = True   #  subtract linear trend? 

   cc_mth_seas = 'mth'   # mth seas

# Variables 'hfds', 'tauuo', 'tos', 'depsq', 'taux_dt_dz', 't20d_sq'. 
   c_var = 'hfds'
   area_name  = '-90.0-lat-90.0_-180.0-lon-180.0'
   c_type      = 'integral'   # mean integral
   
   c_dir, c_descrip, c_period, c_end = cmip6_file_names(c_model, c_expt_type, c_var)

   infile_stub   = '../Time_series/'+c_model+'/' + c_model + '_' + c_descrip + '_' + c_period + '_'

   file_var_in  = infile_stub + c_var +'_'+ area_name +'_'+ c_type + '.txt'  
   
   file_var_out = infile_stub + c_var +'_'+ area_name +'_'+ c_type      
   if ll_subtract_linear_trend :    
      file_var_out = file_var_out + '_anom_minus_trend.txt'     
   else : 
      file_var_out = file_var_out + '_anom.txt'     

   file_std = open( '../std_out/anom_mth'+'.txt','w')

###----------------------------------------------------------------------------------
   
# 1. Read the time-series into numpy arrays of length itim_end - itim_st + 1 

   itim_st = 0 
   itim_end = n_times - 1     #  uses full length of time-series - other options perfectly possible 

   var_ts_in  = timeseries_read(file_var_in,  itim_st, itim_end, file_std ) 
   ntimes = len(var_ts_in) 

   if cc_mth_seas == 'mth' :
      n_per_year = 12  
      out_title = 'monthly anomalies'
   elif cc_mth_seas == 'seas' :
      n_per_year = 4 
      out_title = 'seasonal anomalies'
       
   nyears =  int( ntimes / n_per_year )  
   
   var_ts_out = np.zeros_like ( var_ts_in ) 
   
   for iper in range (n_per_year) : 

      mean_value = 0.0 
      for iyr in range (nyears) : 
         mean_value = mean_value + var_ts_in[iper + iyr*n_per_year]  

      mean_value = mean_value / float(nyears) 
      for iyr in range (nyears) : 
         var_ts_out[iper+iyr*n_per_year] = var_ts_in[iper + iyr*n_per_year] - mean_value    

   if ll_subtract_linear_trend : 
      matrix = np.zeros ( (ntimes,2) )     
      matrix[:,0] =  1.0 
      matrix[:,1] =  np.linspace(0.0,1.0,ntimes)         

      x_soln = linalg.lstsq( matrix, var_ts_out, rcond = None )[0] 

      var_soln = np.zeros(ntimes) 
      var_soln[:] =  matrix[:,0] * x_soln[0] + matrix[:,1] * x_soln[1] 
      
      print ( c_model ) 
      print ( 'constant subtracted; this should be almost identically zero:',  x_soln[0] )
      print ( 'linear trend subtracted; difference between start and end: x_soln[1] = ',  x_soln[1])

      var_ts_out = var_ts_out - var_soln

   
   outfile=open(file_var_out,'w')
   timeseries_write (outfile, out_title, var_ts_out, itim_st, itim_end)

def run_timeseries_anom_mth() :

   c_model_list = ['ACCESS-CM2', 'CanESM5', 'GFDL-CM4', 'IPSL-CM6A-LR', 'MIROC6', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM' ]  # 'MIROC6'
   c_expt_type = 'piControl' # '4xCO2' or 'piControl'

   n_times = 6000   #   500 or 150 years of monthly means 

   for c_model in c_model_list : 

      timeseries_anom_mth_CMIP6(c_model, c_expt_type, n_times)
 
run_timeseries_anom_mth()
