# -*- coding: iso-8859-1 -*-
'''
run_lsq_fit_plots  uses  lsq_fit_read and lsq_fit_plots

Date: 8 March 2023
Purpose: Plot out dependence of least squares fits on time-averaging    
Author: Mike Bell   

'''

#------------------------------------------------------------------------------------------
# To build the file this script processes do the following steps    
# copy lsq_fit_prep.sh into the ../lsq_fit_plots/expt_id2 directory 
# F="filename_for_input"
# ./lsq_fit_prep.sh $F 
# 
# or 
# 
# 1. cd ../lsq_fit_plots/expt_id2   
# 2. F="filename to use as input"  remove the .txt part of the name 
# 3. grep rsq_subset "$F".txt > Summaries/"$F"_rsq.txt   *** see further notes below 
# 4. Edit the new file and include 3 lines at the start taken from the bottom of timeseries_lsq_fit_multiple.py 
#    n_itim_st, n_len_average
#    itim_st     
#    len_average 
#  
#    for example for files built from seasonal data
#       4, 8 
#       1, 2, 3, 4
#       1, 2, 4, 8, 12, 20, 28, 40
#
#  In step 3 one can use alternative inputs to rsq_subset
#    rsq_subset              _rsq.txt
#    x_soln_subset.ifit.0    _a0.txt
#    x_soln_subset.ifit.1    _a1.txt
#    std_subset              _std.txt


def lsq_fit_read(file_in, i_fit_var_plot) :  

   import sys 
   import numpy as np       

   print ( 'opening file_in = ', file_in ) 
   f1  = open(file_in,'r')

   line1 = f1.readline()      
   print ( line1 )
   
   [ c_itim_st, c_len_tim_avg ] = line1.split(', ')
   n_itim_st     = int(c_itim_st)
   n_len_tim_avg = int(c_len_tim_avg) 
   print ( ' n_itim_st, n_len_tim_avg = ', n_itim_st, n_len_tim_avg ) 
   
   line2 = f1.readline()      
   print ( line2 )
   c_tim_st = line2.split(', ')
   list_tim_st = []
   for item in range ( n_itim_st ) :
      list_tim_st.append( c_tim_st[item].strip() ) 
   print ( ' c_tim_st = ', c_tim_st ) 
   print ( ' list_tim_st = ', list_tim_st ) 

   line3 = f1.readline()      
   print ( line3 )
   c_len_tim_avg = line3.split(', ')
   list_len_tim_avg = []
   for item in range ( n_len_tim_avg ) :
      list_len_tim_avg.append( c_len_tim_avg[item].strip() ) 
   print ( ' c_len_tim_avg = ', c_len_tim_avg ) 
   print ( ' list_len_tim_avg = ', list_len_tim_avg ) 

   values = np.zeros( (n_itim_st, n_len_tim_avg) )
    
   for it_st in range ( n_itim_st ) :
      for it_avg in range ( n_len_tim_avg ) :
          for ifit in range (3) :
             linen = f1.readline()      
             if ifit == i_fit_var_plot :
                c_len_typ, c_rsq, c_it_in, c_value = linen.split(', ')  
                print ( 'c_len_typ, c_rsq, c_it_in, c_value = ', c_len_typ, c_rsq, c_it_in, c_value ) 

                if c_it_in.strip() != list_len_tim_avg[it_avg] : 
                   print( ' stopping because of inconsistency; c_it_in, list_len_tim_avg[it_avg] = ', c_it_in, list_len_tim_avg[it_avg] ) 
                   sys.exit()

                values[it_st, it_avg] = float(c_value)
                print ( 'it_st, it_avg, values[it_st, it_avg] = ', it_st, it_avg, values[it_st, it_avg] ) 


   return list_tim_st, list_len_tim_avg, values 

def lsq_fit_plots(c_fit_typ, i_fit_var_plot):

   import sys
   import math         
   import numpy as np
   import matplotlib.pyplot as plt

#------------------------------------------------------------------------------------------
# 1. build file names - most of this is copied from timeseries_lsq_fit_multiple.py 

   expt_id1='DEEPC_v05'   # 'DEEPC_v05'    #'u-ak306' # 'u-ak306' # 'u-am120' (4xCO2),  'u-ar766' (control run)  
   expt_id2='era5'        # 'era5'         #'u-ak306' #
   expt_id3='en4'         # 'en4'          #'u-ak306' #

#   expt_id1 = 'u-ak306' ; expt_id2 = expt_id1 ; expt_id3 = expt_id1

   expt_id1_2 = expt_id1+'_'+expt_id2
   expt_id_123 = expt_id1+'_'+expt_id2+'_'+expt_id3

# Variables 'hfds', 'tauuo', 'tos', 'depsq', 'taux_dt_dz', 't20d_sq'. 
   c_var1 = 'hfds'
   c_var2 = 'tauuo'  
   c_var3 = 't20d' 
   
   c_var = [ 'var1', 'var2', 'var3' ] 
   c_fit_plot = ['var2_only', 'var3_only', 'full'] 

   c_var_name = [ c_var1, c_var2, c_var3 ] 

   if c_fit_typ == 'std' :
      c_fit_var      = c_var[i_fit_var_plot] 
      c_fit_var_name = c_var_name[i_fit_var_plot] 

   else :
      c_fit_var = c_fit_plot[i_fit_var_plot]

 
   
   c_lat1 = '-5.0-lat-5.0'   #   '_-5   o '-5
   c_lat2 = '-3.0-lat-3.0'
   c_lat3 = '-5.0-lat-5.0'

   c_lat1 = '5.0'
   c_lat2 = '3.0'
   c_lat3 = c_lat1

   c_lon1='_130.0-lon--70.0_'  #   '_35.0-lon-105.0_' 
   c_lon2='_130.0-lon--70.0_'  #  '_40.0-lon-45.0_'   
   c_lon3='_135.0-lon-137.0_'  #   '_35.0-lon-105.0_'   

   ts_type1 = 'integral_anom'      #    anom, mean, integral integral_anom
   ts_type2 = 'mean_anom'      #    anom, mean, integral mean_anom
   ts_type3 = 'mean_anom'      #    anom, mean, integral mean_anom

# name of file output
   area_var_name = c_var1+'_'+c_lat1+c_lon1+ts_type1+'_'+c_var2+'_'+c_lat2+c_lon2+c_var3+'_'+c_lat3+c_lon3  

   cc_ann_seas_mth_mean = 'mth'
   ll_convert_season_labels_to_months = False #  True makes plot consistent with CMIP6 plots 

   if cc_ann_seas_mth_mean == 'ann' : 
      y_m = 'y_06'
      ch_x_axis = '' # 'Years'
      itim_file_st = 1950         ; itim_file_end = 2049            # 1951, 2050    ; 1851, 2348

   elif cc_ann_seas_mth_mean == 'seas' : 
      y_m = 's_03' 
      ch_x_axis = '' # 'Season'
      if expt_id1 == 'u-ar766' :
         itim_file_st = 1850 ; itim_file_end = 2348 
      else :
         itim_file_st = 1950 ; itim_file_end = 2049 

   elif cc_ann_seas_mth_mean == 'mth'  : 
      y_m = 'm_01'                            # was 'm_01'  or 'm' 
      ch_x_axis = '' # 'Months'
      if expt_id1 in [ 'DEEPC_v05', 'en4', 'era5' ] :
         itim_file_st  = 1985      #  1950     1985
         itim_file_end = 2016      #  2009     2016
      else :
         itim_file_st  = 1950      #  1950     1985
         itim_file_end = 2049      #  2009     2016

   else : 
      print ( ' stopping: user error specifying cc_ann_seas_mth_mean = ', cc_ann_seas_mth_mean ) 
      sys.exit()

   str_times= str(itim_file_st)+'_'+str(itim_file_end)

# do not choose both ll_extended_mean and ll_running_mean = True 
   ll_running_mean = False
   ll_extended_mean = True   
   ll_dh_dt         = False 
#------------------------------------------------------------------------------------------
# end of user inputs 
#------------------------------------------------------------------------------------------

# labels for plots generated
   label_var1 = c_var1 +'_'+ c_lat1 + '_' + c_lon1   
   label_var2 = c_var2 +'_'+ c_lat2 + '_' + c_lon2   
   label_var3 = c_var3 +'_'+ c_lat3 + '_' + c_lon3

#------------------------------------------------------------------------------------------
# 2. Read the summary file as input 

   infolder='../lsq_fit_plots/'+expt_id1_2+'/Summaries/'

   file_nam_in = expt_id1_2+'_1'+y_m+'_'+str_times+'_'+area_var_name 
   if ll_running_mean: 
      file_nam_in = file_nam_in + 'run_means_'
   elif ll_extended_mean :    
      if ll_dh_dt : 
         file_nam_in = file_nam_in + 'ext_ddt_' 
      else :
         file_nam_in = file_nam_in + 'ext_means_' 
      
   file_nam_in = file_nam_in + 'lsq_fit_'+c_fit_typ 

   file_in = infolder + file_nam_in +'.txt'

   list_tim_st, list_len_tim_avg, values   = lsq_fit_read(file_in, i_fit_var_plot) 
         
   if ll_convert_season_labels_to_months : 
      for item in range ( len(list_len_tim_avg) ) : 
          iseasons = int (list_len_tim_avg[item]) 
          cmonths  = str( 3 * iseasons)
          list_len_tim_avg[item] = cmonths
	  
#------------------------------------------------------------------------------------------
# 3. Plot outputs as required (user can control what to output here) 

   fig = plt.figure(figsize=(15,8))
   ax1 = plt.subplot2grid((1,1), (0,0))

   plotfilename = infolder + file_nam_in +'_'+ c_fit_var + '.png'

   col_type = [ 'r', 'g', 'cyan', 'darkviolet' ] 

#    plt.plot(time_axis,v_seq*1.E-4,color='k',linewidth=2, linestyle='-')   factor 1.E-4 removed from following 3 lines (24/04/2019)

   for item in [0, 1, 2, 3] : 
       plt_lab =  str(list_tim_st[item])
       col     =  col_type[item] 
       vals = values[item]        
       plt.plot(list_len_tim_avg,vals, color=col, linewidth=8, label = plt_lab)

#   plt.title (file_nam_in+ ' ' + c_fit_var+'\n', fontsize=10) # + c_fit_var

   if c_fit_typ == 'rsq' : 
      ax1.set_ylim([0,1.0])
   elif c_fit_typ == 'std' : 
      if c_fit_var_name == 'hfds' : 
         ymax = 2.5e14  ; ymin = 5.0e13 
      elif c_fit_var_name == 'tauuo' : 
         ymax = 0.008   ; ymin = 0.003
      elif c_fit_var_name == 't20d' : 
         ymax = 20.     ; ymin = 5.0
      ax1.set_ylim([ymin,ymax])
   
   plt.xlabel(ch_x_axis, fontsize=30)
#    plt.ylabel(c_fit_typ, fontsize=30)
   plt.tick_params(axis='both', size = 30, length = 5, width = 3, labelsize=45, grid_linewidth = 20 )      

   ax1.spines['left'].set_linewidth(5)
   ax1.spines['bottom'].set_linewidth(5)
   
   plt.locator_params(axis='y', nbins=4)
   
   plt.grid(linewidth=3)

#   plt.legend()     # there are multiple legends for each line (because it is plotted for each segment so labels & legends are commented out 

   fig.savefig(plotfilename, transparent=False, bbox_inches='tight', pad_inches=0.1)
   plt.close('all')

###----------------------------------------------------------------------------------
   
def run_lsq_fit_plots():

# choose which std or fit to plot (there are three options in each case 
   for c_fit_typ in [ 'rsq', 'std' ] :  #  'rsq', 'a0', 'a1', 'std'   see top of this script   
     for i_fit_var_plot in [ 0, 1, 2 ] :  
        lsq_fit_plots(c_fit_typ, i_fit_var_plot) 

run_lsq_fit_plots()
