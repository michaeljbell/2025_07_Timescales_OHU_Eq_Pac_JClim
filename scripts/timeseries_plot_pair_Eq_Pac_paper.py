# -*- coding: iso-8859-1 -*-
'''
timeseries_sum_diff

Date: 19 June 2020
Purpose: Output two time-series on a single plot     
Author: Mike Bell   

'''

def mean_std ( var ) : 
   import numpy as np
   var_mean = np.mean(var)
   var_std  = np.std(var, ddof=1)
   str_mean_std = '; mean = '+"{:.3f}".format(var_mean) +'; stdvn = ' + "{:.4f}".format(var_std)
   return str_mean_std

def timeseries_sum_diff():

   import sys
   import math         
   import numpy as np
   import matplotlib.pyplot as plt
   from timeseries_read_write import timeseries_read, timeseries_write

   ll_CMIP6 = True 

   if ll_CMIP6 : 
      expt_id = 'HadGEM3-GC31-LL_'
      expt_id_list      = [ expt_id+'piControl', expt_id+'ssp245',  expt_id+'1pctCO2', expt_id+'abrupt-4xCO2' ] 
      itim_st_list      = [   165,                165,                   0,                0 ]
      itim_end_list     = [   249,                249,                  80,               80 ]
      itim_off_set_list = [  1850,               1850,                2015,             2015 ]
   else : 
      expt_id_list      = [ 'bg555', 'bm758', 'ar766' ]
      itim_st_list      = 3*[ 1950 ] 
      itim_end_list     = 3*[ 2013 ]
      itim_off_set_list = 3*[    0 ]

   list_color  = ['tab:green', 'tab:purple', 'tab:cyan', 'tab:red' ]    
   list_style  = ['solid',     'dashed',     'dashdot',  'dotted'  ]

   if ll_CMIP6 : 
      outplotfolder = '../Time_series/CMIP6/'
   else : 
      expt_id_out = expt_id_list[0]
      outplotfolder = '../Time_series/u-'+expt_id_out+'/'

   l_plot_ts_legend = False
   ll_diff = False
   ll_ann_mean = True 

   if ll_ann_mean: 
      y_m = 'y_06' #   'y_06'
      ch_x_axis = 'Years'
   else : 
      y_m = 'm' 
      ch_x_axis = 'Months'

   outfilestub = ''
   for item in range( len(expt_id_list) ) : 
      outfilestub = outfilestub + expt_id_list[item] + '_' + str(itim_st_list[item]) + '_'

# equatorial Pacific areas

   lon_min_A =   130.0  
   lon_max_A =   -70.0 
   if ll_diff : 
     lon_min_A =   160.0 #  130.0,  160.0  
     lon_max_A =  -150.0 #  -70.0, -150.0 
     lon_min_B =   -150.0
     lon_max_B =   -90.0

   lat_min=     -3.0
   lat_max =     3.0

   c_type = 'mean'    # 'mean', 'integral'

   c_head_lat   = '_'+str(lat_min)+'-lat-'+str(lat_max) 
   c_head_lon_A = '_'+str(lon_min_A)+'-lon-'+str(lon_max_A)+'_'+c_type
   c_head_out =  c_head_lat + c_head_lon_A
   if ll_diff : 
     c_head_lon_B = '_'+str(lon_min_B)+'-lon-'+str(lon_max_B)+'_'+c_type
     c_head_out =  c_head_out + c_head_lon_B

   file_stem_out =outplotfolder+outfilestub+c_head_out     
   file_std = open('../std_out/timeseries_plot.txt','w') 
   
###----------------------------------------------------------------------------------
   
# 1. Read the time-series into numpy arrays of length itim_end - itim_st + 1 
   
   cvar_list                  =  [ 'tauuo' ]        # 'hfds',  'tos', 'tauuo'
   multiplication_factor_list =  [ 1.0   ]     # 1.0E-15, 1.0,   1.0

   nitems = len(cvar_list)

   for item in range(nitems) :   
      cvar = cvar_list[item]
      multiplication_factor = multiplication_factor_list[item]

# 2. Plot out the two time-series
      label_size = 30
      plt.figure(1,figsize=(18,4.5),dpi=600)
#      plt.subplot(111) 

      for item in range ( len(expt_id_list) ) :
         expt_id      = expt_id_list[item]
         itim_st      = itim_st_list[item] 
         itim_end     = itim_end_list[item]
         itim_off_set = itim_off_set_list[item]

         timestub = '_1'+ y_m +'_'+ str(itim_st) + '_'+str(itim_end)+'_'

         if ll_CMIP6 : 
            infolder = '../Time_series/CMIP6/'
         else : 
            infolder='../Time_series/u-'+expt_id+'/u-'
         infilestub = expt_id+ timestub
         fileA  =infolder+infilestub+cvar+c_head_lat+c_head_lon_A+'.txt'     
         ts     = timeseries_read(fileA,  itim_st, itim_end, file_std ) 

         if ll_diff : 
            fileB  =infolder+infilestub+cvar+c_head_lat+c_head_lon_B+'.txt'     
            tsB    = timeseries_read(fileB,  itim_st, itim_end, file_std ) 
            ts     = ts - tsB

         ts = ts * multiplication_factor

         print ( ' np.size(ts) =', np.size(ts) ) 

         times = np.arange ( itim_st + itim_off_set, itim_end + itim_off_set + 1 )
      
#         col = list_color[ item % 4 ]
#         sty = list_style[ int(item/4) ]

         col = list_color[ item ]
         sty = list_style[ item ]

         str_mean_std_ts = mean_std(ts) 
          
         plt.plot(times,ts, color = col, linestyle = sty, linewidth=4, label=expt_id+str_mean_std_ts)

      plt.tick_params(axis='x', labelsize=label_size )      
      plt.tick_params(axis='y', labelsize=label_size )

      title = cvar + ' for latitudes ' + c_head_lat + '; longitudes ' + c_head_lon_A + '; multiplication_factor = ' + str(multiplication_factor)
      if ll_diff : 
         print ( title) 
         title = title +' minus '+c_head_lon_B

      plt.grid(linewidth=3)
      if l_plot_ts_legend :
         plt.title( title ) 
         plt.legend()
         plt.xlabel(ch_x_axis, fontsize=label_size)

      file_out = file_stem_out+'_'+cvar +'_'+ str(multiplication_factor) + '_legend'+str(l_plot_ts_legend)+'.png'
      plt.savefig(file_out) 
      plt.close('all') 
      print ( 'generated file ', file_out ) 
         
timeseries_sum_diff()
