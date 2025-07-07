# -*- coding: iso-8859-1 -*-
'''
timeseries_lsq_fit_multiple  

Date: 22 April 2016
Purpose: Least squares fit to time-series: Regression of var1 using var2 and var3. Splits time-series into multiple chunks 
Author: Mike Bell   

07/03/2023:   revised to allow subsets; monthly and seasonal data; loop over starting and averaging periods  

'''

def prt_mean_std ( var, cvar, file_out ) : 
   import numpy as np
   var_mean = np.mean(var)  
   var_std  = np.std(var, ddof=1)
   print ( ' variable = ', cvar, '; mean = ', var_mean, '; standard devn = ', var_std, file = file_out) 
   return

def timeseries_lsq_fit_multiple(file_std, c_model, c_expt_type, n_times, itim_st, ll_running_mean, len_time_window, ll_extended_mean, len_average, cc_ann_seas_mth_chk):

   import sys
   import math         
   import numpy as np
   import matplotlib.pyplot as plt
   import numpy.linalg as linalg
   from timeseries_read_write import timeseries_read, timeseries_write
   from timeseries_running_mean import timeseries_running_mean
   from timeseries_extended_mean import timeseries_extended_mean
   from utilities_CMIP6 import cmip6_file_names

   expt_id=c_model 

   l_subset = False    # do not use subsetting as well as running means or time windows    
   if l_subset :
      len_time_subset = 30   #   typically 40 or or 50 or 100   ; 
   
   l_plot = False    #   True to output plots - that would generate a lot of output - file names have not been set up for that

# Variables 'hfds', 'tauuo', 'tos', 'depsq', 'taux_dt_dz', 't20d_sq'. 
   c_var1 = 'hfds' 
   c_var2 = 'hfds'  
   c_var3 = 'hfds'
   
   c_lat1 = '90.0'
   c_lat2 = '5.0'
   c_lat3 = '5.0'

   c_lon1='_-180.0-lon-180.0_'
#   c_lon1='_130.0-lon--70.0_'   
#   c_lon1='_-150.0-lon--90.0_'   
#   c_lon2='_135.0-lon--75.0_'  
   c_lon3='_130.0-lon--70.0_'   
#   c_lon3='_-150.0-lon--90.0_'   
   c_lon2 = c_lon3 

   ts_type1 = 'integral_anom_minus_trend'      #    anom, mean, integral
   ts_type2 = 'integral_anom_minus_trend'  #    anom, mean, integral integral_anom
   ts_type3 = 'integral_anom_minus_trend'      #    anom, mean, integral integral_anom_minus_trend

   itim_end = n_times - 1     #  uses full length of time-series - other options perfectly possible 

# name of file output
   area_var_name = c_var1+'_'+c_lat1+c_lon1+ts_type1+'_'+c_var2+'_'+c_lat2+c_lon2+c_var3+'_'+c_lat3+c_lon3  

   cc_ann_seas_mth_mean = 'mth'
   
   if cc_ann_seas_mth_mean != cc_ann_seas_mth_chk : 
      print ( 'Stopping because of inconsistent specification of cc_ann_seas_mth_mean at top and bottom of script.  ')
      print ( ' cc_ann_seas_mth_mean, cc_ann_seas_mth_chk = ', cc_ann_seas_mth_mean, cc_ann_seas_mth_chk ) 
      sys.exit()


   if cc_ann_seas_mth_mean == 'mth'  : 
      ch_x_axis = 'Months'
   else : 
      print ( ' stopping: user error specifying cc_ann_seas_mth_mean = ', cc_ann_seas_mth_mean ) 
      sys.exit()

   infolder='../Time_series/'
   outfolder='../lsq_fit_plots/'+expt_id+'/'

   c_dir, c_descrip1, c_period1, c_end = cmip6_file_names(c_model, c_expt_type, c_var1)
   c_dir, c_descrip2, c_period2, c_end = cmip6_file_names(c_model, c_expt_type, c_var2)
   c_dir, c_descrip3, c_period3, c_end = cmip6_file_names(c_model, c_expt_type, c_var3)

   c_in_file_stub1 = c_model + '/' + c_model + '_' + c_descrip1 + '_' + c_period1 + '_'
   c_in_file_stub2 = c_model + '/' + c_model + '_' + c_descrip2 + '_' + c_period2 + '_'
   c_in_file_stub3 = c_model + '/' + c_model + '_' + c_descrip3 + '_' + c_period3 + '_'


   ll_legend = False   # true implies legend is written on plots (tends to overwrite the timeseries ...)  

# following are now set in an outer wrapper function 
# do not choose both ll_extended_mean and ll_running_mean = True 

#   ll_running_mean = False
#   len_time_window = 0    # integer number 0, 5, 10, 20 

#   ll_extended_mean = False   
#   len_average = 0

   ll_write_files = False # write out extended mean files to aid debugging 
   
#------------------------------------------------------------------------------------------
# end of user inputs 
#------------------------------------------------------------------------------------------

   file_var1=infolder+c_in_file_stub1+c_var1+'_-'+c_lat1+'-lat-'+c_lat1+c_lon1+ts_type1    
   file_var2=infolder+c_in_file_stub2+c_var2+'_-'+c_lat2+'-lat-'+c_lat2+c_lon2+ts_type2   
   file_var3=infolder+c_in_file_stub3+c_var3+'_-'+c_lat3+'-lat-'+c_lat3+c_lon3+ts_type3   

# labels for plots generated
   if l_plot : 
      label_var1 = c_var1 +'_'+ c_lat1 + '_' + c_lon1   
      label_var2 = c_var2 +'_'+ c_var2 + '_' + c_lon2   
      label_var3 = c_var3 +'_'+ c_var3 + '_' + c_lon3

#-----------------------------------------------------------------------------------------------------

   file_nam_plots = c_model + '_' + c_descrip1 + '_' + c_period1 + '_' + area_var_name

   if l_subset : 
      file_nam_plots = file_nam_plots + 'subsets_' + str(len_time_subset) + '_'

   if ll_running_mean: 
      file_nam_plots = file_nam_plots + 'run_means_'
      if l_plot :
         file_nam_plots = file_nam_plots + str(len_time_window) + '_'
      
   elif ll_extended_mean :    
      file_nam_plots = file_nam_plots + 'ext_means_' 
      if l_plot :
         file_nam_plots = file_nam_plots + str(len_average) + '_'

   if file_std == 'Null' : 
      file_std = open( '../lsq_fit_plots/'+expt_id+'/' + file_nam_plots + 'lsq_fit.txt','w') 
      ll_print_inputs = True

      print ( ' c_var1 = ', c_var1, '; c_var2 = ', c_var2, '; c_var3 = ', c_var3 , file = file_std )
      if l_plot : print ( ' label_var1 = ', label_var1, '; label_var2 = ', label_var2, '; label_var3 = ', label_var3, file = file_std )
      print ( ' area_var_name = ', area_var_name, '; cc_ann_seas_mth_mean = ', cc_ann_seas_mth_mean, file = file_std ) 
      print ( ' n_times, itim_st = ', n_times, itim_st, file = file_std ) 
   else :
      ll_print_inputs = False

###----------------------------------------------------------------------------------
   
# 1. Read the time-series into numpy arrays of length itim_end - itim_st + 1 

   var1_ts_full  = timeseries_read(file_var1+'.txt',  itim_st, itim_end, file_std, ll_print = ll_print_inputs ) 
   var2_ts_full  = timeseries_read(file_var2+'.txt',  itim_st, itim_end, file_std, ll_print = ll_print_inputs ) 
   var3_ts_full  = timeseries_read(file_var3+'.txt',  itim_st, itim_end, file_std, ll_print = ll_print_inputs ) 

   len_ts_full = itim_end - itim_st + 1    #   revised 07/03/2023  

   if l_subset :
      n_subsets = int(len_ts_full/len_time_subset) 
   else : 
      n_subsets = 1 
      len_time_subset = len_ts_full

# print each time 
   print ( ' itim_st, itim_end = ', itim_st, itim_end, file = file_std ) 
   print ( ' ll_running_mean= ', ll_running_mean, '; len_time_window = ', len_time_window,file = file_std ) 
   print ( ' ll_extended_mean= ', ll_extended_mean, '; len_average = ', len_average,file = file_std ) 
      
   print ( ' n_subsets = ',  n_subsets, file = file_std ) 
   
   rsq_subset    = np.zeros( (3,n_subsets) ) 
   x_soln_subset = np.zeros( (3,2,n_subsets) ) 
   std_subset    = np.zeros( (3,n_subsets) )

   for isubset in range ( n_subsets ) : 
    itim_st_sub  = isubset * len_time_subset 
    itim_end_sub = itim_st_sub + len_time_subset 
    var1_ts = var1_ts_full[itim_st_sub:itim_end_sub+1]
    var2_ts = var2_ts_full[itim_st_sub:itim_end_sub+1]
    var3_ts = var3_ts_full[itim_st_sub:itim_end_sub+1]
   
    mean_var1 = np.mean ( var1_ts) 
    mean_var2 = np.mean ( var2_ts) 
    mean_var3 = np.mean ( var3_ts)

    if isubset == 0 : print ( ' np.size(var1_ts) =', np.size(var1_ts) , file = file_std) 

###----------------------------------------------------------------------------------
# 1.B do any time-filtering 

    if ll_running_mean : 
      var1_ts = timeseries_running_mean(var1_ts, len_time_window)
      var2_ts = timeseries_running_mean(var2_ts, len_time_window)
      var3_ts = timeseries_running_mean(var3_ts, len_time_window)
      ioffset_st = int(len_time_window/2)      
      itim_st  = itim_st + ioffset_st
      itim_end = itim_end - (len_time_window - 1 - ioffset_st)    # this ensures that itim_end - itim_st is reduced by len_time_window - 1
      times = np.arange ( itim_st, itim_end+1 )  

    elif ll_extended_mean :    
      var1_ts = timeseries_extended_mean(var1_ts, len_average)  
      var2_ts = timeseries_extended_mean(var2_ts, len_average)  
      var3_ts = timeseries_extended_mean(var3_ts, len_average)  
      ioffset_st = int(len_average/2)      
      itim_st  = itim_st + ioffset_st
      n_times  = len (var1_ts) 
      itim_end = itim_st + (n_times - 1) * len_average 
      times = np.arange( itim_st, itim_end+1, len_average ) 

      if ll_write_files : 
         outfile=open(file_var1+'_em_'+str(len_average),'w')
         timeseries_write (outfile, 'extended means', var1_ts, itim_file_st, itim_file_end-1)   # -1 temporarily 

    else : 
      times = np.arange ( itim_st, itim_end+1 )  

###----------------------------------------------------------------------------------

# 2. Store the standard deviations and plot out the timeseries    

    std_subset[0,isubset] = np.std(var1_ts, ddof=1)
    std_subset[1,isubset] = np.std(var2_ts, ddof=1)
    std_subset[2,isubset] = np.std(var3_ts, ddof=1)
      
    if l_plot :  

      plt.figure(1)
      plt.subplot(311) 
      plt.plot(times,var1_ts, linewidth=4, color='r', label=label_var1)
      plt.legend()
      plt.subplot(312) 
      plt.plot(times,var2_ts, linewidth=4, color='b', label=label_var2)
      plt.legend()
      plt.subplot(313) 
      plt.plot(times,var3_ts, linewidth=4, color='g', label=label_var3)
      plt.legend()
      plt.xlabel(ch_x_axis, fontsize=12)
      plt.savefig(infolder+file_nam_plots+'timeseries.png') 
      plt.close('all') 

###----------------------------------------------------------------------------------
   
# 3. Do least square fits 
   
# 3.1 Build the matrix   var1_ts = var2_ts . x_0 + var3_ts . x_1 + const . x_2 + epsilon  

    ntimes = np.size( times )
    
#    print ('ntimes, times[0], times[-1] =', ntimes, times[0], times[-1]) 
  
    matrix = np.zeros ( (ntimes,2) )     

    ch_fit_array = ['var2_only', 'var3_only', 'full']

    for ifit in range(3) : 
      ch_fit = ch_fit_array[ifit]
      if ch_fit == 'var2_only' :
           matrix[:,0] =  var2_ts
           matrix[:,1] =  mean_var2         # using 1.0 here gave incorrect results  
      elif ch_fit == 'var3_only' :
           matrix[:,0] =  var3_ts 
           matrix[:,1] =  mean_var3          # using 1.0E17 here gave incorrect results    
      elif ch_fit == 'full' : 
           matrix = np.zeros ( (ntimes,3) ) 
           mean_prod = math.sqrt( abs ( mean_var2 * mean_var3 ) ) 
           matrix[:,0] =  var2_ts 
           matrix[:,1] =  var3_ts               
           matrix[:,2] =  mean_prod            

# 3.2 Find the least squares solution 
      x_soln = linalg.lstsq( matrix, var1_ts, rcond = None )[0] 
      
      var1_soln = np.zeros(ntimes) 
      var1_soln[:] =  matrix[:,0] * x_soln[0] + matrix[:,1] * x_soln[1] 
      if ch_fit == 'full' : var1_soln[:] = var1_soln[:] + matrix[:,2] * x_soln[2]
      var1_err = var1_ts - var1_soln
      mean_var1_err_sq = np.mean( var1_err * var1_err) 
      mean_var1 = np.mean( var1_ts ) 
      var1_var  =  var1_ts - mean_var1
      mean_var1_var_sq = np.mean( var1_var * var1_var ) 

      ratio =  mean_var1_err_sq / mean_var1_var_sq
      if ratio < 1.0 : 
         rsq = 1.0 - ratio
      else : 
         rsq = 0.0 
	 
      rsq_subset[ifit,isubset]        = rsq 
      x_soln_subset[ifit,0:2,isubset] = x_soln[0:2] 

###----------------------------------------------------------------------------------
              
# 4.3 Produce summaries of solution  as comma separated strings

   if ll_running_mean : 
      len_tim_avg = len_time_window
   elif ll_extended_mean :    
      len_tim_avg = len_average
   else : 
      len_tim_avg = 0 

   print ( 'outfolder, file_nam_plots = ', outfolder, file_nam_plots ) 

   for ifit in range (3) :     
      print ( ' len_tim_avg, ifit, ch_fit_array[ifit]   ', len_tim_avg, ',', ifit,ch_fit_array[ifit], file = file_std ) 
      print ( '   len_tim_avg, x_soln_subset[ifit,0],   ', len_tim_avg, ',', str(x_soln_subset[ifit,0]).strip('[]'), file = file_std ) 
      print ( '   len_tim_avg, x_soln_subset[ifit,1],   ', len_tim_avg, ',', str(x_soln_subset[ifit,1]).strip('[]'), file = file_std ) 
      print ( '   len_tim_avg, rsq_subset[ifit],        ', len_tim_avg, ',', str(rsq_subset[ifit]).strip('[]'), file = file_std ) 
      if l_subset : prt_mean_std ( rsq_subset[ifit], 'rsq', file_std ) 

   file_names = [ file_var1, file_var2, file_var3 ] 
   for ivar in range (3) :     
      print ( ' variable : ', file_names[ivar], file = file_std )
      print ( 'len_tim_avg,  std_subset[ivar],   ', len_tim_avg, ',', str(std_subset[ivar]).strip('[]'), file = file_std ) 
      if l_subset : prt_mean_std ( std_subset[ivar], 'std[ivar]', file_std ) 

   return file_std

###----------------------------------------------------------------------------------

def run_timeseries_lsq_fit_multiple():

# do not choose both ll_extended_mean and ll_running_mean = True 

# for monthly and seasonal data, itim_st determines the month or season when time-averaging starts 

# if cc_ann_seas_mth_mean == 'seas' and 'm_03' set itim_st  = 1 for MAM; 2 for JJA; 3 for SON; 4 for DJF)    
# len_average = 4 * n_years  to give n_years average 


# default settings 
   ll_extended_mean = True   
   len_average = 0

   ll_running_mean = False
   len_time_window = 0    # integer number 0, 5, 10, 20 

   itim_st = 1

   c_model_list = ['ACCESS-CM2', 'CanESM5', 'GFDL-CM4', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'IPSL-CM6A-LR', 'MIROC6', 'MPI-ESM1-2-LR', 'NorESM2-LM' ]  
#    c_model_list = ['ACCESS-CM2', 'CanESM5' ] #  'GFDL-CM4', 'IPSL-CM6A-LR', 'MIROC6', 'MPI-ESM1-2-LR' ] 

   c_expt_type = 'piControl' # '4xCO2' or 'piControl'
       
   for c_model in c_model_list :
      
# section to loop over multiple choices of time averaging (CMIP6 data are monthly means)  

      cc_ann_seas_mth_chk = 'mth'   

      file_std = 'Null'        # this is only used as input when the function is called when ll_open_out_file = False

      n_times = 6000  #   first 500 years of data

      for itim_st in [ 3, 6, 9, 12 ] : #   [ 3, 6, 9, 12 ]  :
         for len_average in [ 1, 2, 3, 4, 6, 12, 24, 48, 60, 84, 120 ] : #  [ 1, 2, 3, 4, 6, 12, 24, 48, 60, 84, 120 ] : 

            file_std = timeseries_lsq_fit_multiple(file_std, c_model, c_expt_type, n_times, itim_st, ll_running_mean, len_time_window, ll_extended_mean, len_average, cc_ann_seas_mth_chk)


run_timeseries_lsq_fit_multiple()

