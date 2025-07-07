# -*- coding: iso-8859-1 -*-
'''
run_lsq_fit_plot_revised_CMIP6_with_reanalysis  uses  lsq_fit_read and lsq_fit_plots

Date: 8 March 2023
Purpose: Plot out dependence of least squares fits on time-averaging. All CMIP6 models and reanalysis on one plot  
Author: Mike Bell   

'''

#------------------------------------------------------------------------------------------
# To build the file this script processes do the following steps    
#
# edit list_files_CMIP6.py then run it   
# edit lsq_fit_prep_rsq_CMIP6.sh or lsq_fit_prep_reg_coeff_CMIP6.sh and run it  
#
#

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
#                print ( 'c_len_typ, c_rsq, c_it_in, c_value = ', c_len_typ, c_rsq, c_it_in, c_value ) 

                if c_it_in.strip() != list_len_tim_avg[it_avg] : 
                   print( ' stopping because of inconsistency; c_it_in, list_len_tim_avg[it_avg] = ', c_it_in, list_len_tim_avg[it_avg] ) 
                   sys.exit()

                values[it_st, it_avg] = float(c_value)
#                print ( 'it_st, it_avg, values[it_st, it_avg] = ', it_st, it_avg, values[it_st, it_avg] ) 


   return list_tim_st, list_len_tim_avg, values 

def lsq_fit_plots(c_model, c_descrip, c_period, area_var_name, time_mean_type, c_fit_typ, i_fit_var_plot):

   import sys 
   import numpy as np       

#------------------------------------------------------------------------------------------
# 1. Read the summary file as input 

   infolder='../lsq_fit_plots/'+c_model+'/' # '/' or + '/Summaries/'
   file_nam_in = c_model + '_' + c_descrip + '_' + c_period + '_' + area_var_name + time_mean_type
   file_nam_in = file_nam_in + 'lsq_fit_'+c_fit_typ 
   file_in = infolder + file_nam_in +'.txt'

   list_tim_st, list_len_tim_avg, values   = lsq_fit_read(file_in, i_fit_var_plot) 

# find the means of the values at the different start times 
   n_tim_starts = len(list_tim_st)
   n_time_averages = len(list_len_tim_avg)
   vals = np.zeros(n_time_averages) 
   for ist in range( n_tim_starts ) : 
      for iav in range (n_time_averages) : 
         vals[iav] = vals[iav] + values[ist,iav]        
   vals = vals / float(n_tim_starts)

   return list_len_tim_avg, vals 

###----------------------------------------------------------------------------------
   
def run_lsq_fit_plot_revised_CMIP6_with_reanalysis():

   import sys 
   import numpy as np       
   import matplotlib.pyplot as plt
   from utilities_CMIP6 import cmip6_file_names

# 1. build file names - most of this is copied from timeseries_lsq_fit_multiple.py 

# Variables 'hfds', 'tauuo', 'tos', 'depsq', 'taux_dt_dz', 't20d_sq'. 
   c_var1 = 'hfds'   # 'hfds' or 'tos'
   c_var2 = 'hfds'   #   'tauu' or 'tos' 
   c_var2_reanal = 'hfds'   # 'tos' or 'tauuo'
   c_var3 = 'hfds' 
   
   l_std_plot_min_zero = True  #  set this to True to set minimum value on y-axis for std plot to zero 
   
   c_var = [ 'var1', 'var2', 'var3' ] 
   c_fit_plot = ['var2_only', 'var3_only', 'full'] 
   
   c_lat1 = '90.0' # '5.0' or '90.0'
   c_lat2 = '5.0'   # '3.0' or '5.0'
   c_lat3 = '5.0'
   c_lat2_reanal = c_lat2  # c_lat2 or c_lat3

   c_lon1= '_-180.0-lon-180.0_' # '_-180.0-lon-180.0_'  or '_130.0-lon--70.0_'
   c_lon2= '_130.0-lon--70.0_' # '_135.0-lon--75.0_'
   c_lon3= '_130.0-lon--70.0_'   
   c_lon2_reanal = c_lon2 
#   c_lon2 = c_lon3 

   n_times_in_series = 6000  #   500 years of monthly data 


   ts_type1 = 'integral_anom_minus_trend'      #     integral_anom or  integral_anom_minus_trend ; anom, mean, integral
# seem not to be used 
#   ts_type2 = 'mean_anom_minus_trend'      #    anom, mean, integral mean_anom
#   ts_type3 = 'mean_anom_minus_trend'      #    anom, mean, integral mean_anom

# name of file output
   area_var_name = c_var1+'_'+c_lat1+c_lon1+ts_type1+'_'+c_var2+'_'+c_lat2+c_lon2+c_var3+'_'+c_lat3+c_lon3  
   area_var_name_reanal = c_var1+'_'+c_lat1+c_lon1+ts_type1+'_'+c_var2_reanal+'_'+c_lat2_reanal+c_lon2_reanal+c_var3+'_'+c_lat3+c_lon3  

# do not choose both ll_extended_mean and ll_running_mean = True 
   ll_running_mean = False
   ll_extended_mean = True   
   
   ll_legend = False
   
# build list of data sets to use as input 

   c_model_list = ['ACCESS-CM2',  'CanESM5', 'GFDL-CM4', 'IPSL-CM6A-LR', 'MIROC6', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM' ]  
   c_expt_type = 'piControl' # '4xCO2' or 'piControl'

   model_time = []   
   for c_model in c_model_list : 
      c_dir, c_descrip, c_period, c_end = cmip6_file_names(c_model, c_expt_type, c_var1)
      model_time.append( [c_model,c_descrip, c_period, n_times_in_series] ) 
   
#   model_time.append( ['ACCESS-CM2','piControl_r1i1p1f1_gn', '095001-144912', 6000] ) 
#   model_time.append( ['CanESM5','piControl_r1i1p1f1_gn', '520101-620012', 12000] )
#   model_time.append( ['GFDL-CM4','piControl_r1i1p1f1_gn', '015101-065012', 6000 ] )
#   model_time.append( ['HadGEM3-GC31-LL','piControl_r1i1p1f1_gn', '225001-384912', 24000] )
#   model_time.append( ['HadGEM3-GC31-MM','piControl_r1i1p1f1_gn', '185001-234912', 6000] )
#   model_time.append( ['IPSL-CM6A-LR','piControl_r1i1p1f1_gr', '185001-384912', 18000] )    # tau lat, lon is missing 
#   model_time.append( ['MIROC6','piControl_r1i1p1f1_gn', '320001-399912', 6000 ] )
#   model_time.append( ['MPI-ESM1-2-LR','piControl_r1i1p1f1_gn', '185001-284912', 12000 ] )
#   model_time.append( ['NorESM2-LM','piControl_r1i1p1f1_gn', '160001-210012', 6000] )     # actually has 6012 

   file_out_stem = '../lsq_fit_plots/AllCMIP6/piControl_' + area_var_name

#------------------------------------------------------------------------------------------
# end of user inputs 
#------------------------------------------------------------------------------------------

   if ll_running_mean: 
      time_mean_type = 'run_means_'
   elif ll_extended_mean :    
      time_mean_type = 'ext_means_' 

   file_out_stem = file_out_stem + time_mean_type

   if ll_legend :    
      file_out_stem = file_out_stem + 'legend_'
      
# choose which std or fit to plot (there are three options in each case 
   for c_fit_typ in [ 'rsq'] :  #  'rsq',  'std', 'reg_coeff_0'   see top of this script   
      for i_fit_var_plot in [ 0, 1 ] :    # [0, 1] for 'rsq' or 'a0' ; [0, 1, 2] for 'std' or [0] for 'std'

         if c_fit_typ == 'std' :
            c_fit_var = c_var[i_fit_var_plot] 
         else :
            c_fit_var = c_fit_plot[i_fit_var_plot]

         file_out = file_out_stem + c_fit_var
  
#------------------------------------------------------------------------------------------
# 3. Plot outputs as required (user can control what to output here) 

         fig = plt.figure(figsize=(15,8))
         ax1 = plt.subplot2grid((1,1), (0,0))
	 
#         col_type = [ 'r', 'g', 'cyan', 'darkviolet' ]   # not used yet

         ll_CMIP6 = True 

# calculate and plot the CMIP6 data; one line per model 
         for clist in model_time : 
            c_model   = clist[0] 
            c_descrip = clist[1]
            c_period  = clist[2]
            n_times = clist[3]

            c_label = c_model
            if c_label == 'HadGEM3-GC31-LL' : c_label = 'HadGEM3-LL'
            if c_label == 'HadGEM3-GC31-MM' : c_label = 'HadGEM3-MM'

            list_len_tim_avg,vals = lsq_fit_plots(c_model, c_descrip, c_period, area_var_name, time_mean_type, c_fit_typ, i_fit_var_plot) 
            plt.plot(list_len_tim_avg, vals, linewidth=8, label = c_label)

# now calculate and plot the reanalysis line
         
         ll_CMIP6 = False 
         if c_var1 == 'hfds' : 
            c_period  = '1985_2016' 
            c_descrip = '1m'        # '1m' or '1m_01'
            if c_var2 == 'hfds' : 
               c_model   = 'DEEPC_v05_DEEPC_v05'
            else : 
               c_model   = 'DEEPC_v05_era5'
         else : 
            c_model   = 'era5_era5'
            c_period  = '1985_2016'#  '1960_2018' or '1985_2016'
            c_descrip = '1m_01'    # '1m' or '1m_01' 

         list_len_tim_avg,vals = lsq_fit_plots(c_model, c_descrip, c_period, area_var_name_reanal, time_mean_type, c_fit_typ, i_fit_var_plot) 
         plt.plot(list_len_tim_avg, vals, linewidth=8, linestyle = ':', label = 'obs product')
   
         if c_fit_typ == 'rsq' : 
            ax1.set_ylim([0,1.0])
         elif c_fit_typ == 'std' and l_std_plot_min_zero :
            ax1.set_ylim(ymin=0.0)
#      ax1.set_ylim(ymax=0.9)    # used just for tos max value (figure 9) 
 
         if ll_legend :
            plt.legend()
 
         ch_x_axis = 'months'
         plt.xlabel(ch_x_axis, fontsize=30)
#   plt.ylabel(c_fit_typ, fontsize=20)

         plt.tick_params(axis='both', size = 30, length = 5, width = 3, labelsize=45, grid_linewidth = 20 )      

         ax1.spines['left'].set_linewidth(5)
         ax1.spines['bottom'].set_linewidth(5)

         plt.locator_params(axis='y', nbins=4)
   
         if not ll_legend : 
            plt.grid(linewidth=3)

#   plt.legend()     # there are multiple legends for each line (because it is plotted for each segment so labels & legends are commented out 

         plotfilename = file_out + c_fit_typ 
         if ll_legend : 
            plotfilename = plotfilename + '_legend'
         plotfilename = plotfilename + '.png'
         print( 'plotfilename = ', plotfilename)
         fig.savefig(plotfilename, transparent=False, bbox_inches='tight', pad_inches=0.1)
         plt.close('all')

run_lsq_fit_plot_revised_CMIP6_with_reanalysis()
