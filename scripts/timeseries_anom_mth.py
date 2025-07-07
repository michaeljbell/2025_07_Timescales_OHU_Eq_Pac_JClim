# -*- coding: iso-8859-1 -*-
'''
timeseries_anom_mth  

Date: 7 Sept 2020
Purpose: Calculate the timeseries of monthly anomalies from a timeseries of monthly values      
Author: Mike Bell   

'''

def timeseries_anom_mth():

 import sys
 import math         
 import numpy as np
 from timeseries_read_write import timeseries_read, timeseries_write
 import numpy.linalg as linalg

 for expt_id in ['DEEPC_v05'] : # ['era5', 'u-co779'] : # , 'u-co766'] :
 
# [ 'u-aj354', 'u-ak141', 'u-ak144', 'u-ar599']
# [ 'u-bh325', 'u-bh326', 'u-bh327', 'u-bh328' ] :    # 'u-ak306', 'u-am120',  'u-ar766'  'DEEPC_v05' 'en4' 'era5' 
# [ 'era5' ] 
  for c_var in ['hfds'] : #    [ 'tauuo', 'pmsl_e_m_w' ] : # 'nshf', 'tos', 'slhf', 'ssr' ] : 


# Variables 'hfds', 'tauuo', 'tos', 'depsq', 'taux_dt_dz', 't20d_sq'. 
#   c_var = 'tos'

   ll_subtract_linear_trend = True   #  subtract linear trend? 
   cc_mth_seas = 'mth'   # mth seas
  
   infilestub = expt_id 
   
   if cc_mth_seas == 'mth' :
     if expt_id in [ 'DEEPC_v05', 'en4',  'ORAS5' ] :  #  'era5', has been in this list in the past
     
        itim_file_st = 1985 ; itim_file_end = 2016        #   1985 or 1993 ; 2013 or 2016 ; 1960 and 2018
        itim_st = 1         ; itim_end = 384              #    384  or 348 or 252; 708
     else :
        itim_file_st = 1982 ; itim_file_end = 2008  # 2048 or 2049  1960 2014 
        itim_st      = 0 ; itim_end      = 323      # 1200 or 1188    659   
     infilestub = infilestub+'_1m_' + ''    # + '01_' or + '' is non-standard 

   elif cc_mth_seas == 'seas' :
     if expt_id == 'u-ar766' :
       itim_file_st = 1850 ; itim_file_end = 2348 
       itim_st = 1         ; itim_end = 1996    
     else: 
       itim_file_st = 1950 ; itim_file_end = 2049 
       itim_st = 1         ; itim_end = 396 
     infilestub = infilestub+'_1s_03_'

   infolder   = '../Time_series/'+expt_id+'/'

   area_name  = '-90.0-lat-90.0_-180.0-lon-180.0'   # '-5.0-lat-5.0_35.0-lon-105.0' or '-90.0-lat-90.0_-180.0-lon-180.0'
#   area_name  = '236-jj-248_80-ii-150'           # '236-jj-261_80-ii-150'  

   c_type      = 'integral'   # mean integral

   file_var_in  = infolder+infilestub+str(itim_file_st)+'_'+str(itim_file_end)+'_'+c_var+'_'+area_name+'_'+c_type+'.txt'     
   file_var_out = infolder+infilestub+str(itim_file_st)+'_'+str(itim_file_end)+'_'+c_var+'_'+area_name+'_'+c_type     

   if ll_subtract_linear_trend :    
      file_var_out = file_var_out + '_anom_minus_trend.txt'     
   else : 
      file_var_out = file_var_out + '_anom.txt'     

   file_std = open( '../std_out/anom_mth'+'.txt','w')

###----------------------------------------------------------------------------------
   
# 1. Read the time-series into numpy arrays of length itim_end - itim_st + 1 

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
      
      print ( 'constant subtracted; this should be almost identically zero:',  x_soln[0] )
      print ( 'linear trend subtracted; difference between start and end: x_soln[1] = ',  x_soln[1])

      var_ts_out = var_ts_out - var_soln
   
   outfile=open(file_var_out,'w')
   timeseries_write (outfile, out_title, var_ts_out, itim_st, itim_end)

timeseries_anom_mth()
