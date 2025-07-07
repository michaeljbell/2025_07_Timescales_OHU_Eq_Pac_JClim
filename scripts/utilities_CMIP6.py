# -*- coding: iso-8859-1 -*-
'''
CMIP6 utilities: contains cmip6_file_names, lat_lon_area_mask_cmip6, cmip6_convert_int_to_str4, cmip6_int_to_str_mth

Date: 31 March 2021
Author: Mike Bell   

'''

def cmip6_file_names(c_model, c_expt_type, c_var): 

   import sys 
   
   if c_model in ['ACCESS-CM2', 'CanESM5', 'GFDL-CM4', 'IPSL-CM6A-LR', 'MIROC6', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM' ] : 
     print( 'c_model = ', c_model ) 
   else : 
     print ( ' cmip6_file_names: stopping because of unrecognised c_model = ', c_model ) 
     sys.exit()

   if c_var in [ 'hfds', 'tauu', 'tos' ] :
     print( 'c_var = ', c_var ) 
   else : 
     print ( ' cmip6_file_names: stopping because of unrecognised c_var = ', c_var ) 
     sys.exit()

   if c_expt_type == 'piControl'  :  
     if c_model == 'ACCESS-CM2' : 
        c_period = '095001-144912'
     elif c_model == 'CanESM5' : 
        c_period = '520101-620012'
     elif c_model == 'GFDL-CM4' : 
        c_period = '015101-065012'   
     elif c_model == 'IPSL-CM6A-LR' : 
        if c_var == 'hfds':              # the period covered differs from this !!
           c_period = '185001-209912'
        else :
           c_period = '185001-384912'
     elif c_model == 'MIROC6' : 
        if c_var == 'hfds':
           c_period = '320001-369912'
        else :
           c_period = '320001-399912'
     elif c_model == 'MPI-ESM1-2-LR' : 
        c_period = '185001-284912'
     elif c_model == 'NorESM2-LM' : 
        c_period = '160001-210012'
     elif c_model == 'HadGEM3-GC31-LL' : 
        c_period = '225001-384912'
     elif c_model == 'HadGEM3-GC31-MM' : 
        c_period = '185001-234912'

     c_descrip = 'piControl_r1i1p1f1_gn'  
     if c_model == 'IPSL-CM6A-LR' :  c_descrip = 'piControl_r1i1p1f1_gr'
     if c_var == 'tauu' and c_model == 'GFDL-CM4'   : c_descrip = 'piControl_r1i1p1f1_gr1'    

     c_end = '_combined'
     if c_var == 'tauu' and c_model == 'ACCESS-CM2' : c_end = ''

     c_dir = '/scratch/jonathab/data_for_mike_cmip6/'
       
   elif c_expt_type == '4xCO2' :
     c_period = '' 
     c_descrip = 'abrupt-4xCO2_r1i1p1f1' 
     if c_model in [ 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM' ] :
        c_descrip = 'abrupt-4xCO2_r1i1p1f3' 
     c_end = 'full' 
     c_dir = '/scratch/jonathab/data_for_mike_cmip6/4xCO2/'


   else : 
     print ( ' cmip6_file_names: stopping because of unrecognised c_expt_type = ', c_expt_type ) 
     sys.exit()
     

   if c_var == 'tos' : 
      c_dir = c_dir + 'sst/tos_Omon_'
   elif c_var == 'hfds' :
      c_dir = c_dir + 'nshf/hfds_Omon_'
   elif c_var == 'tauu' :
      c_dir = c_dir + 'wind_stress/tauu_Amon_'

   return c_dir, c_descrip, c_period, c_end  
   
#------------------------------------------------------------------------------------------

def lat_lon_area_mask_CMIP6(cvar, c_model, c_descrip, c_period): 

   import sys 
   import numpy as np 
   import netCDF4
   from utilities import rd_ncdf_var_check_one
   from constants import ra 
   
   dirin  = '/scratch/jonathab/data_for_mike_cmip6/'

   if cvar in [  'hfds', 'tos' ] :
   
      c_descrip_area = c_descrip
      if c_model == 'IPSL-CM6A-LR' :
         c_descrip_area = 'piControl_r1i1p1f1_gn'   #  accomodates a file naming anomaly  

      file_area = dirin + 'area/areacello_Ofx_' + c_model +'_' + c_descrip_area + '.nc' 
      g_area = netCDF4.Dataset(file_area,'r')
      
      if c_model == 'GFDL-CM4' : 
        clat = 'lat'      ;  clon = 'lon'  
      elif c_model == 'IPSL-CM6A-LR' : 
        clat = 'nav_lat'      ;  clon = 'nav_lon'  
      else : 
        clat = 'latitude' ;  clon = 'longitude'  
   
      nav_lat  = rd_ncdf_var_check_one(g_area, clat)
      nav_lon  = rd_ncdf_var_check_one(g_area, clon)
      nav_area = rd_ncdf_var_check_one(g_area, 'areacello')

      comb = '_combined'

      if c_model == 'GFDL-CM4' : 
        clat = 'lat'      ;  clon = 'lon'  
      
      file_tos = dirin + 'sst/tos_Omon_' + c_model +'_' + c_descrip +'_' + c_period + comb + '.nc' 
      g_tos = netCDF4.Dataset(file_tos,'r')
      tos = g_tos.variables['tos'][0]
      land_mask = np.where( tos > 1E10, 0.0, 1.0 )     
      
   elif cvar in [  'tauu' ] :
      if c_model != 'ACCESS-CM2' :
         comb = '_combined' 
      else : 
         comb = ''

      if c_model == 'GFDL-CM4' : 
         c_descrip = 'piControl_r1i1p1f1_gr1'   #  accomodates a file naming anomaly 

      file_area = dirin + 'wind_stress/tauu_Amon_' + c_model +'_' + c_descrip +'_' + c_period  + comb + '.nc'
      g_area = netCDF4.Dataset(file_area,'r')

      if c_model == 'IPSL-CM6A-LR' :
        lat_cen = rd_ncdf_var_check_one(g_area, 'lat') 
        lon_cen = rd_ncdf_var_check_one(g_area, 'lon') 
        print (' np.shape(lat_cen) = ', np.shape(lat_cen) ) 
        print (' np.shape(lon_cen) = ', np.shape(lon_cen) ) 
        dlat    = np.zeros_like (lat_cen) 
        dlat[1:-1] = 0.5 * np.radians ( lat_cen[2:] - lat_cen[:-2] ) ; dlat[0] = dlat[1] ; dlat[-1] = dlat[-2]
        dlon    = np.zeros_like (lon_cen) 
        dlon[1:-1] = 0.5 * np.radians ( lon_cen[2:] - lon_cen[:-2] )  ; dlon[0] = dlon[1] ; dlon[-1] = dlon[-2]
      
      else : 
        lats = rd_ncdf_var_check_one(g_area, 'lat_bnds') 
        lons = rd_ncdf_var_check_one(g_area, 'lon_bnds') 
        print (' np.shape(lats) = ', np.shape(lats) ) 
        print (' np.shape(lons) = ', np.shape(lons) ) 

        dlat    = np.radians ( lats[:,1] - lats[:,0] ) 
        lat_cen =       0.5* ( lats[:,1] + lats[:,0] ) 

        dlon    = np.radians ( lons[:,1] - lons[:,0] )
        lon_cen =       0.5* ( lons[:,1] + lons[:,0] ) 

      nav_lon, nav_lat = np.meshgrid( lon_cen, lat_cen )
      dlon_msh, dlat_msh = np.meshgrid( dlon, dlat )

      nav_lat_rad = np.radians ( nav_lat)        
      cos_nav_lat = np.cos ( nav_lat_rad ) 
      
      nav_area =  ra * ra * dlon_msh * dlat_msh * cos_nav_lat
      
      land_mask = np.ones_like(nav_lon)
   
   else :
      print ('stopping because variable choice not recognised, cvar = ', cvar ) 
      sys.exit()

   g_area.close()

  
   return nav_lat, nav_lon, nav_area, land_mask   

#------------------------------------------------------------------------------------------------

def cmip6_convert_int_to_str4(Iyr) : 

# converts Iyr (integer) to a string of length 4 (e.g. 13 is converted to '0013') 

   import sys 
   Cyr = str(Iyr)
   LenCyr = len(Cyr)
   LenZeros = 4 - LenCyr
   if LenZeros > 0 : Cyr = LenZeros*'0'+Cyr
   if LenZeros < 0 : 
     print ( ' convert_int_to_str4: stopping because the year number in the filename is greater than 9999: Iyr = ', Iyr ) 
     sys.exit()  
   return Cyr    

#------------------------------------------------------------------------------------------------

def cmip6_int_to_str_mth(month) :

# input:   month is integer month (between 1 and 12) 
# returns: c_mth month as string of length 2 

   if month < 10 : 
      c_mth = '0'+str(month)
   else : 
      c_mth = str(month)
   return c_mth
   
#------------------------------------------------------------------------------------------------
