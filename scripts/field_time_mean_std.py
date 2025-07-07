
# -*- coding: iso-8859-1 -*-
'''
field_time_mean_std  

Date: 22 June 2020
Purpose: Calculate the time mean and standard deviation from a time-series of annual mean fields      
Author: Mike Bell   

'''

def field_time_mean_std(grid_type, varname, c_mth_cal_strt):

#   grid_type        # 'T' or 'U' or 'F'
#   varname          # tos, sos, hfds, empmr, tauuo, tauvo, t20d, uo, thetao, ssr, str, slhf, sshf, nshf, tauuo, tauvo
#   c_mth_cal_strt   # '3', '6', '9' or '12'  

   import netCDF4
   import numpy as np
   from utilities import date_string, period_string, rd_ncdf_var_check_one, create_2D_file, create_3D_file, rd_ncdf_var_check_one

   data_dir = '/data/users/frmk/GC345/'   # CMIP6AMIP or CMIP6HighResMIP or GC5
   data_type = 'ocean_model' # 'ocean_model' or 'era5' or 'DEEPC' or 'EN4' or 'atmos_model'
   ll_3D_field = False

   ll_running_mean = False
   len_run = 1 
   
   iyr_st  = 1979   # usually 1960 for era5  ;   1950 or 1980 for ocean_model ; 1985 for DEEPC;  1850 for ar766 
   iyr_end = 2077   # the last year used IS iyr_end ; usually 2018 for era5; 2049 for ocean_model ; 2016 for DEEPC; 2348 for ar766

#   i_mth_cal_strt = 12 #  12 or 6  # start month of annual mean  between 1 and 12 (1 for Jan; 12 for Dec) (set to 1 for Monthly means)
   nmonths_per_period = 1   # 12 for annual means, 3 for seasonal means
   str_imth = '' #  '_' + c_mth_cal_strt #  '' or '_' + c_mth_cal_strt
   file_msy = 'm'   # 's' or 'm' or 'y' 
   dir_msy = file_msy + str_imth
     
   i_mth_cal_strt = int(c_mth_cal_strt) 

   if data_type == 'atmos_model' : 
      ll_convert_Kelvin_to_Centigrade = True    #  Used to convert mean tos from Kelvin to Centigrade (important for calculation of standard devn)
   else :                                       # applies for ocean_model, era5 and EN4 (era5 and EN4 were converted to Centigrade earlier) 
      ll_convert_Kelvin_to_Centigrade = False   

   if ll_3D_field :
      grid_surface = 'grid_'
   else :
      grid_surface = 'surface_'
   
   if data_type == 'ocean_model' : 
      expt_id = 'ai599'   #   'am120' or 'ak306'  aj368   ar766  co779
      config_res = 'eORCA025' # 'eORCA1', eORCA025, eORCA12
      inmesh = data_dir+'MESH_MASK/mesh_mask_'+config_res+'_GC5.nc'  #  -GO6 or _GC5
      infolder=data_dir+'u-'+expt_id + '/on'         
      filestub='/'+ grid_surface+grid_type+'/nemo_'+expt_id+'o_1'+file_msy+'_'   # + varname + '_' is optional 
      fileend='.nc'
      co_depth = 'gdept_1d'     # depends on grid type - but lower case 
   elif data_type == 'atmos_model' : 
      expt_id = 'bh325'   #   'am120' or 'ak306'  aj368 ar766
      config_res = 'N216' # 'eORCA1', eORCA025, eORCA12
      inmesh = data_dir+'MESH_MASK/mesh_mask_'+config_res+'_'+grid_type+'.nc'  
      infolder=data_dir+'u-'+expt_id + '/an'                 
      filestub='/'+ grid_surface+grid_type+'/'+expt_id+'_1'+file_msy+'_'  # + varname + '_' ;  + grid_surface+grid_type+'/'+ is optional ; 1y or 1s 
      fileend='.nc' # '' or '_grid-T.nc'
      co_depth = 'gdept_1d'     # depends on grid type - but lower case 
   elif data_type == 'era5' :
      infolder='/data/users/frmk/ERA5_pm90/on'      # _pm30 is optional
      inmesh = '/data/users/frmk/ERA5_pm90/mesh_mask.nc'  
      filestub='/era5_1'+file_msy +'_' #  + varname + '_' # + varname + '_' is optional 
      fileend='_pm90.nc'                         # _pm30 is optional
   elif data_type == 'DEEPC' :
      infolder='/data/users/frmk/DEEPC_v05/on'    
      inmesh = '/data/users/frmk/DEEPC_v05/mesh_mask.nc'  
      filestub='/DEEPC_v05_1'+file_msy+'_'
      fileend='.nc'                       
   elif data_type == 'EN4' :
      infolder='/data/users/frmk/EN4/on'    
      inmesh = '/data/users/frmk/EN4/mesh_mask.nc'  
      filestub='/2D/EN.4.2.2.f.analysis.g10.' #   + varname + '.'   #  + varname + '.' is optional 
      fileend='.nc'                       
      co_depth = 'depth'      

    
   str_period_out = period_string( iyr_st, i_mth_cal_strt, nmonths_per_period, iyr_end ) 
    
#   if data_type == 'ocean_model' or data_type == 'atmos_model' : # or data_type == 'atmos_model' : 
#      period_string = str(iyr_st)+'1201-'+str(iyr_end+1)+'1201'
#   else :
#      period_string = str(iyr_st)+'_'+str(iyr_end)

   ll_omit_coords = False 

#----------------------------------------------------------------------------------------------

   if varname == 'sozotaux' : 
      varname_out = 'tauuo'
   elif varname == 'tos_con' : 
      varname_out = 'tos'
   else :
      varname_out = varname

   if not ll_running_mean :  
      infilestub = infolder+dir_msy+filestub    #   +varname+'_'
   else : 
      infilestub = infolder+dir_msy+'_rf'+filestub+varname+'_' + str(len_run) + '_'

   imth = 0 

   g = netCDF4.Dataset(inmesh,'r')

   if ll_3D_field : 
      depth = rd_ncdf_var_check_one (g, co_depth) 
   else : 
      depth = 0.0 

   nav_lon =   rd_ncdf_var_check_one (g, 'nav_lon')         
   nav_lat =   rd_ncdf_var_check_one (g, 'nav_lat')         
   g.close()

   (npj, npi) = np.shape(nav_lat)
   print ('shape(nav_lat) = ', np.shape(nav_lat) )
   
   sum_var    = np.zeros_like( nav_lat ) 
   sum_var_sq = np.zeros_like( nav_lat ) 

   for iyr in range( iyr_st, iyr_end+1 ) :             

      date_str = period_string( iyr, i_mth_cal_strt, nmonths_per_period ) 


      if data_type == 'ocean_model' :  
         infile=infilestub+date_str+'_grid-'+grid_type+fileend
      elif data_type == 'atmos_model' :   
         infile=infilestub+date_str+fileend #  '_grid-'+grid_type+fileend
      else :
#         date_str = str(iyr)
         infile=infilestub+date_str+fileend
	    
      g = netCDF4.Dataset(infile,'r')
      var = rd_ncdf_var_check_one(g, varname) 
      g.close()

      print ('data_type, date_str, infile, varname, np.shape(var) = ', data_type, date_str, infile, varname, np.shape(var))

      if ll_convert_Kelvin_to_Centigrade :
         if varname == 'tos'  : 
            var = var - 273.15   # original data are in Kelvin 
    
      sum_var = sum_var + var  
      sum_var_sq = sum_var_sq + var*var 

   r_len = 1.0 / (iyr_end + 1 - iyr_st) 
   mean_var = sum_var * r_len 
   mean_var_sq = sum_var_sq * r_len 
  
   variance = mean_var_sq - mean_var * mean_var 
   std_devn = np.sqrt(variance) 

   file_mean = infilestub+'mean_'+varname_out+'_'+str_period_out+'_grid-'+grid_type+fileend
   print ( 'outputing file_mean = ', file_mean)
   
   if ll_3D_field :
      create_3D_file ( nav_lat, nav_lon, depth, mean_var, 'mean_'+varname_out, file_mean )  

   else :
      create_2D_file ( nav_lat, nav_lon, depth, mean_var, 'mean_'+varname_out, file_mean, ll_omit_coords )  

   file_stdvn = infilestub+'std_'+varname_out+'_'+str_period_out+'_grid-'+grid_type+fileend
   if ll_3D_field :
      create_3D_file ( nav_lat, nav_lon, depth, std_devn, 'std_'+varname_out, file_stdvn )  

   else :
      create_2D_file ( nav_lat, nav_lon, depth, std_devn, 'std_'+varname_out, file_stdvn, ll_omit_coords )  
         
def run_field_time_mean_std():

  grid_type = 'T' 
  varname_list = ['tos', 'hfds'] # 't20d', 'slhf', 'sshf', 'ssr', 'str', 'nshf', 'tos' ] # 'tos', 'tos_con' , 'nshf'
#  varname_list = ['nshf']   # hfds or nshf
#  varname_list = ['tos',  'hfds', 't20d' ] 
#  varname_list = ['t20d']

  c_mth_cal_strt_list = [  '01',  '02',  '03',  '04',  '05',  '06',  '07',  '08',  '09',  '10',  '11',  '12' ] # ['06'] # [ '03', '06', '09', '12' ] # , '12' ]

#  c_mth_cal_strt_list = ['06']

  for c_mth_cal_strt in c_mth_cal_strt_list : 
#     for varname in varname_list : 
#        field_time_mean_std(grid_type, varname, c_mth_cal_strt)

     field_time_mean_std( 'U', 'sozotaux', c_mth_cal_strt)      #  tauuo sozotaux   
#     field_time_mean_std( 'V', 'tauvo', c_mth_cal_strt)      # 

run_field_time_mean_std()    

