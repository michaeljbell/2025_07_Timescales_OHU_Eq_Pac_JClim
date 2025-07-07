# -*- coding: iso-8859-1 -*-
'''
field_annual_means_ERA5_DEEPC_EN4  

Date: 11 October 2022
Purpose: Calculate the annual means or seasonal means of ERA5 or DEEPC surface fields or EN4 isotherm depths or 3D temperature fields  ; 
         EN4 t20d fields are calculated by field_EN4_isotherm_depths.py   
Author: Mike Bell   

'''

def field_annual_seasonal_means_ERA5_DEEPC_EN4(c_mth_strt, cc_timemean_type):

#   c_mth_strt             # start month of period output;  between 1 and 12 ('01' for Jan; '12' for Dec) 
#   cc_timemean_type       #  Annual or Monthly or Single_seasons

   import sys        
   import netCDF4
   import math
   import numpy as np
   from constants import ra  
   from utilities import period_string, days_in_month, rd_ncdf_var_check_one, create_list_2D_file, create_list_2D_3D_file

###----------------------------------------------------------------------------------
# Start of user inputs 


   dataset = 'ORAS5' #    'ERA5' or 'DEEPC' or 'EN4' or 'GloSea5' or 'ORAS5'
   
   ll_3D_fields = False   #  only relevant for EN4 

   ll_verbose = True

#----------------------------------------------------------------------------------------------

   i_mth_strt = int(c_mth_strt)  
   if cc_timemean_type == 'Annual' :  
      y_m = 'y_' 
      y_m_out = y_m + c_mth_strt
      i_mth_end = i_mth_strt + 11 
      n_mths = 12 
   elif cc_timemean_type == 'Single_seasons' : 
      y_m = 's_' 
      y_m_out = 's'
      i_mth_end =  i_mth_strt+2
      n_mths = 3 
   if cc_timemean_type == 'Monthly' :  
      y_m = 'm_' 
      y_m_out = 'm' 
      i_mth_end =  i_mth_strt
      n_mths = 1 
      
   if dataset == 'ERA5' :    

      list_varnamein  =  ['sst'] #  [   'ewss',   'nsss',  'sst',  'slhf',   'sshf',    'ssr',    'str' ]
      list_varnameout =  ['tos'] #  [  'tauuo',  'tauvo',  'tos',  'slhf',   'sshf',    'ssr',    'str' ]
      list_filenamein =  ['tos_'] # [ 'tauuo_', 'tauvo_', 'tos_', 'nshfs_', 'nshfs_', 'nshfs_', 'nshfs_' ]

      iyr_st  = 1960    # 1960 is earliest
      iyr_end = 2018    # 2019 is latest   
      infolder='/data/users/frmk/ERA5_pm90/'    # this is also the root of the outfolder 

      infilestub  = 'era5_'                     # +filein+'_'
      infileend  = '_pm90.nc'                   # this is also the outfile end 
      outfilestub = 'on'+y_m_out+'/'+'era5_1'+ y_m 

   elif dataset == 'DEEPC' :

      list_varnamein  =  [ 'Fmass' ]
      list_varnameout =  [ 'hfds' ]
      list_filenamein =  [ '' ]

      iyr_st  = 1985    # 1985 is earliest
      iyr_end = 2016    # 2016; 2017 only has data to November  
      infolder='/data/users/frmk/DEEPC_v05/'     

      infilestub  = 'DEEP-C_Fmass_v5_198501-201711'
      infileend  = '.nc'
      outfilestub = 'on'+y_m_out+'/DEEPC_v05_1'+ y_m 

   elif dataset == 'EN4' :

     if ll_3D_fields : 
      list_varnamein  =  [ 'temperature' ]
      list_varnameout =  [ 'thetao' ]
      list_filenamein =  [ '' ]
      dirD = '/3D/'
     else  : 
      list_varnamein  =  [ 't20d' ]
      list_varnameout =  [ 't20d' ]
      list_filenamein =  [ 't20d.' ]
      dirD = '/2D/'

      iyr_st  = 2016    # 1960 is earliest
      iyr_end = 2016    # 2019 is last  
      infolder='/data/users/frmk/EN4/'     

      infilestub  = 'EN.4.2.2.f.analysis.g10.'   
      infileend   = '.nc'
      outfilestub = 'on'+y_m_out+dirD+ infilestub 

   elif dataset == 'GloSea5' :

      list_varnamein  =  [ 't20d' ]
      list_varnameout =  [ 't20d' ]
      list_filenamein =  [ '' ]
      dirD = '/2D/'

      iyr_st  = 1993    # 1993 is earliest
      iyr_end = 2015    # 2016 is last  
      infolder='/data/users/frmk/GloSea5/'     

      infilestub  = 'onm/2D/monthmean_t20d.'   
      infileend   = '.nc'
      outfilestub = 'on' + y_m_out + '/2D/GloSea5_1' + y_m 
   elif dataset == 'ORAS5' :

      list_varnamein  =  [ 'so20chgt' ]
      list_varnameout =  [ 't20d' ]
      list_filenamein =  [ '' ]
      dirD = '/surface_T/'

      iyr_st  = 1985    # 1985 is earliest
      iyr_end = 1993    # 2014 is last  
      infolder='/data/users/frmk/ORAS5/'     

      infilestub  = 'onm/surface_T/so20chgt_control_monthly_highres_2D_'   
      infileend   = '_CONS_v0.1.nc'
      outfilestub = 'on' + y_m_out + '/surface_T/ORAS5_1' + y_m 

   len_listnames = len( list_varnamein ) 
   
#---------------------------------------------------------------------------------------------
# End of user inputs 
#---------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------
# Read data from the mesh mask file  
#---------------------------------------------------------------------------------------------

   meshfilename =  infolder + 'mesh_mask.nc'
   g = netCDF4.Dataset(meshfilename,'r')
   nav_lon =   rd_ncdf_var_check_one (g, 'nav_lon')         
   nav_lat =   rd_ncdf_var_check_one (g, 'nav_lat')         
   tmask   =   rd_ncdf_var_check_one (g, 'tmask')         
   g.close()

   if not ll_3D_fields and dataset == 'EN4' :
      tmask = tmask[0]
  
   dep = 0.0 

#----------------------------------------------------------------------------------------------
# loop through the time-varying fields calculating and outputting seasonal or annual means 
#---------------------------------------------------------------------------------------------

   for iyr in range ( iyr_st, iyr_end + 1  ) :   

    if ll_3D_fields : 
      list_varout = []
    else : 
      list_varout = [nav_lat, nav_lon]
      list_var_names_out = ['nav_lat','nav_lon'] 

    for item in range ( len_listnames ) : 
 
     ll_first  = True 
     ndays_sum = 0

     varnamein  = list_varnamein[item] 
     varnameout = list_varnameout[item]
     filein     = list_filenamein[item]
      
     for imth in range ( i_mth_strt, i_mth_end + 1 ) :  

       if imth < 13 : 
         iyr_input = iyr
         imth_input = imth-1      # imth_input runs from 0 to 11
       else : 
         iyr_input = iyr+1
         imth_input = imth-13	 

       if dataset == 'ERA5' :    
         idec_st = ( iyr_input // 10 ) * 10 # integer division 
         itime = imth_input + 12 * (iyr_input - idec_st)  
         date_str = str(idec_st) + '_' + str(idec_st+9)
       elif dataset == 'DEEPC' : 
         iyr_off = iyr_input - 1985   # 1985 is first year in data set
         itime = imth_input + 12 * iyr_off
         date_str = ''
       elif dataset == 'EN4' : 
         imth_cal = imth_input + 1 
         if imth_cal > 9 : 
           imth_str = str(imth_cal)
         else : 
            imth_str = '0'+str(imth_cal) 
         date_str = str(iyr_input)+imth_str 
         itime = 0 
       elif dataset == 'GloSea5' : 
         date_str = period_string( iyr, imth, 1 )
         itime = 0 
       elif dataset == 'ORAS5' : 
         date_str = period_string( iyr, imth, 1 )
         date_str = date_str[0:6]
         itime = 0 
   
       infile=infolder + infilestub + filein + date_str + infileend
       if ll_verbose: print ('iyr, imth, imth_input, iyr_input, infile = ', iyr, imth, imth_input, iyr_input, infile)
 	
       g = netCDF4.Dataset(infile,'r')
       var = g.variables[varnamein][itime] 
       g.close()

       var = var * tmask   #   multiply by land-sea mask 

# The ERA5 tauuo data contain accumulations over 1 day  
# I've divided them by the number of seconds in a day (following the documentation on the web site)  
       if dataset == 'ERA5' :
          if varnamein == 'ewss' or varnamein == 'nsss' or filein == 'nshfs_' :     
             var = var / 86400.    # divide by the number of seconds in a day (original units are N m^-2 s or J m^-2 )

       if dataset == 'ERA5' and varnamein == 'sst' :     
          var = var - 273.15    # convert Kelvin to Celsius
   
       if dataset == 'EN4' :
          if varnameout == 't20d' :     
             var = var.data
             var = np.where(var>0.0, var, 0.0 )    # convert "missing data values" to zeros 
          if varnameout == 'thetao' :     
             var = var - 273.15    # convert Kelvin to Celsius

       ndays_this_mth = days_in_month(iyr_input, imth_input, calendar_type = 'Gregorian' )  
       var       = var * ndays_this_mth 

#       print ( ' varnamein,  np.shape (var) = ', varnamein, np.shape (var) )  

       if ll_first :
          ll_first = False 
          var_sum   =  var
       else : 
          var_sum   = var_sum   + var

       ndays_sum = ndays_sum + ndays_this_mth

     varout = var_sum / ndays_sum 

#     print ( ' np.shape (varout) = ', np.shape (varout) )  

     list_varout.append(varout) 
    
    period_str = period_string( iyr, i_mth_strt, n_mths )
    file_out = infolder+outfilestub+period_str+infileend

    if ll_3D_fields : 
      create_list_2D_3D_file ( nav_lat, nav_lon, [], [], list_var_out, list_varnameout, False, 0, file_out)
    else: 
      ny, nx = np.shape( list_varout[0] ) 
      x = np.arange ( nx ) 
      y = np.arange ( ny ) 
      print ( 'ny, nx = ', ny, nx )  
      list_var_names_out = list_var_names_out +  list_varnameout
      create_list_2D_file ( x, 'x', y, 'y', list_varout, list_var_names_out, file_out, depth=0.0, ll_time=False )        
      print (' written ', file_out )
      
#----------------------------------------------------------------------------------------------------------------
def run_field_annual_seasonal_means_ERA5_DEEPC_EN4():

#  cc_timemean_type = 'Monthly' 
#  cmth_list = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12' ]   # use to generate all months in year

#  cc_timemean_type = 'Single_seasons' 
#  cmth_list = ['03', '06', '09', '12' ]   # use to generate standard seasons 
#  cmth_list = ['12']

  cc_timemean_type = 'Annual' 
  cmth_list = ['06'] 

  for c_mth_strt in cmth_list : 
     field_annual_seasonal_means_ERA5_DEEPC_EN4(c_mth_strt,cc_timemean_type)

run_field_annual_seasonal_means_ERA5_DEEPC_EN4()
