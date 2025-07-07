# -*- coding: iso-8859-1 -*-
'''
field_annual_mean_timeseries  

Date: 13 April 2020
Purpose: Calculate annual means of fields from monthly means (as a check)     
Author: Mike Bell   

'''

def annual_mean_timeseries(grid_type, varname_list):

   import netCDF4
   import numpy as np
   from utilities import period_string, days_in_month, rd_ncdf_var_check_one, create_list_2D_file
   
   expt_id = 'bg555'   #   am120 or ak306
   config_res = 'eORCA1'   
   expt_series = 'CMIP6HighResMIP' #  CMIP6HighResMIP or GC5 

   c_mth_start = '06'                     # calendar month of start of the year 
   i_mth_cal_st = int(c_mth_start)        # January has index 1 in this function 
   
   data_dir = '/data/users/frmk/'+expt_series+'/'
   infolder=data_dir+'u-'+expt_id                  
   filestub='surface_'+grid_type+'/nemo_'+expt_id 
   infilestub = infolder+'/onm/'+filestub+'o_1m_'
   outfilestub= infolder+'/ony_'+c_mth_start+'/'+filestub+'o_1y_'    
   
   iyr_st  = 1931 
   iyr_end = 1949

#----------------------------------------------------------------------------------------------
   varname_out_list = []

   for varname in varname_list : 
      if varname == 'tos_con' : 
         varname_out_list.append('tos') 
      if varname == 'sozotaux' : 
         varname_out_list.append('tauuo') 
      else : 
         varname_out_list.append(varname)    

   print( 'varname_list = ', varname_list)
   print( 'varname_out_list = ', varname_out_list)

#----------------------------------------------------------------------------------------------
# read latitudes and longitudes
#   infilemesh = data_dir + 'MESH_MASK/mesh_mask_'+ config_res + '_'+ expt_series + '.nc' 
#   g = netCDF4.Dataset(infilemesh,'r')
#   nav_lat = rd_ncdf_var_check_one(g, 'nav_lat')
#   nav_lon = rd_ncdf_var_check_one(g, 'nav_lon')
#   g.close()
#   print ('shape(nav_lat) = ', np.shape(nav_lat) )


#----------------------------------------------------------------------------------------------

   deptht = np.zeros (1) 

   ll_ann_means = False  # working with monthly means 

   for iyr in range( iyr_st, iyr_end+1 ) :   
   
   
      sum_fields_list = [] 
         
      for imth_cal in range ( i_mth_cal_st, i_mth_cal_st+12 ) :  

         date_str = period_string( iyr, imth_cal, 1 )
       
         infile=infilestub+date_str+'_grid-'+grid_type+'.nc'

         g = netCDF4.Dataset(infile,'r')
         var_list = []
         for varname in varname_list : 
            var = rd_ncdf_var_check_one(g, varname)
            var_list.append(var)  
         g.close()
#         print ('shape(var) = ', np.shape(var) )

         ndays = days_in_month( iyr, imth_cal-1, calendar_type ='Gregorian')   # 
    
         if imth_cal == i_mth_cal_st : 
            print ('reading infile = ', infile)
            sum_ndays  =  ndays
            sum_fields_list = []
            for item in range ( len (varname_list) ) :
               sum_fields_list.append ( var_list[item] * float(ndays) )
         else: 
            sum_ndays  = sum_ndays  + ndays
            for item in range ( len (varname_list) ) :
               sum_fields_list[item] = sum_fields_list[item] + var_list[item] * float(ndays) 
	  
      r_ndays = 1.0 / float(sum_ndays) 
      for item in range ( len (varname_list) ) :
         sum_fields_list[item] = sum_fields_list[item] * r_ndays    # divide by the number of days in year   

      file_out = outfilestub+str(iyr)+c_mth_start+'01-'+str(iyr+1)+c_mth_start+'01_grid-'+grid_type+'.nc'
      print ('writing file_out = ', file_out)

# generate dimensions 
      ny, nx = np.shape( sum_fields_list[0] ) 
      x = np.arange ( nx ) 
      y = np.arange ( ny ) 
      
      create_list_2D_file ( x, 'x', y, 'y', sum_fields_list, varname_out_list, file_out, depth=deptht, ll_time=False )        

         
#----------------------------------------------------------------------------------------------

def run_annual_mean_timeseries(): 

   annual_mean_timeseries('T', ['tos', 'hfds', 't20d'] )  # tos_con 
   annual_mean_timeseries('U', ['tauuo'] )     # tauuo, sozotaux
   
run_annual_mean_timeseries()   
   
