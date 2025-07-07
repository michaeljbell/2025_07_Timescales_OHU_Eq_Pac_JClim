# -*- coding: iso-8859-1 -*-
'''
field_running_mean  

Date: 22 June 2020
Purpose: Calculate running means of annual mean fields      
Author: Mike Bell   

'''

def field_running_mean_CMIP6(c_model, c_descrip, c_period):

   import netCDF4
   import numpy as np
   from utilities import period_string, rd_ncdf_var_check_one, create_2D_timeseries_file
   
   cvar='tauu'  # tos, sos, hfds, empmr, tauuo, G 

   len_run = 5    # length (in years) of running mean

   data_dir = '/data/users/frmk/CMIP6/' 
   infolder=data_dir+ c_model+'/'
                   
   if cvar == 'tauu' and c_model == 'GFDL-CM4'   : 
      c_descrip_fld = 'piControl_r1i1p1f1_gr1'    
   else :
      c_descrip_fld = c_descrip    

   filestub= cvar+'_'+ c_model +'_' + c_descrip_fld +'_' + c_period

   cyr_st  = '00'     # usually 1960 for era5  ;   1950 or 1980 for ocean_model ; 1985 for DEEPC;  1850 for ar766
   iyr_st = int(cyr_st) 
   iyr_end = 499   # the last year used (there are iyr_end fields in the input file) 
   nyears_in = iyr_end 
   
   infiledates = cyr_st + '0601-' + str(iyr_end) + '0601'

   ll_omit_coords = False 

#----------------------------------------------------------------------------------------------

   file_var = infolder + filestub + '_' + infiledates + '.nc' 
   file_out = infolder + filestub + '_rf_' + str(len_run) + '_' + infiledates + '.nc'   # it is simpler to keep the same dates

# Mesh-mask file used 
   g = netCDF4.Dataset(file_var,'r')

   nav_lat = rd_ncdf_var_check_one(g,'nav_lat')
   nav_lon = rd_ncdf_var_check_one(g,'nav_lon')
   (npj, npi) = np.shape(nav_lat)

   print (' np.shape(nav_lat) = ', np.shape(nav_lat) )

   nyears_out = nyears_in - len_run + 1 #  length of series unchanged if len_run = 1

   fields_out = np.zeros( (nyears_out, npj, npi) )
   times_out = np.zeros( nyears_out )
   
   var_win = np.zeros( (len_run, npj, npi) ) 

   ioffset = int(len_run/2)

   r_len_run = 1.0 / float(len_run) 

   var = rd_ncdf_var_check_one(g,cvar) 

   for iyr in range( iyr_st, iyr_end ) :                	     
      var_win[iyr%len_run] = 0.0 
      for iwin in range (len_run) :  
         var_win[iwin] = var_win[iwin] + var[iyr]

      if iyr >= iyr_st + len_run - 1 : 
         iout = (iyr+1)%len_run
         out_field = var_win[iout] * r_len_run    # divide by the number of running values   
         iyr_out = iyr - ( len_run - 1)           # first ouput to iyr if len_run = 1 

         fields_out[iyr_out] = out_field
         times_out[iyr_out]  = iyr + ioffset - ( len_run - 1)   # times are not actually used in field_timeseries_correlation_CMIP6.py
         print ( ' c_model, iyr = ', c_model, iyr )  
            
   create_2D_timeseries_file ( nav_lat, nav_lon, times_out, fields_out, cvar, file_out )  
         
def run_field_running_mean_CMIP6() :

   model_time = []
   model_time.append( ['ACCESS-CM2','piControl_r1i1p1f1_gn', '095001-144912', 6000] ) 
   model_time.append( ['CanESM5','piControl_r1i1p1f1_gn', '520101-620012', 12000] )
   model_time.append( ['GFDL-CM4','piControl_r1i1p1f1_gn', '015101-065012', 6000 ] )
   model_time.append( ['IPSL-CM6A-LR','piControl_r1i1p1f1_gr', '185001-384912', 18000 ] )    # tau lat, lon is missing 
   model_time.append( ['MIROC6','piControl_r1i1p1f1_gn', '320001-399912', 6000 ] )
   model_time.append( ['MPI-ESM1-2-LR','piControl_r1i1p1f1_gn', '185001-284912', 12000 ] )
   model_time.append( ['NorESM2-LM','piControl_r1i1p1f1_gn', '160001-210012', 6000] )     # actually has 6012 
   model_time.append( ['HadGEM3-GC31-LL','piControl_r1i1p1f1_gn', '225001-384912'] )
   model_time.append( ['HadGEM3-GC31-MM','piControl_r1i1p1f1_gn', '185001-234912'] )

   for clist in model_time : 

      c_model   = clist[0] 
      c_descrip = clist[1]
      c_period  = clist[2]

      field_running_mean_CMIP6(c_model, c_descrip, c_period)

run_field_running_mean_CMIP6()
