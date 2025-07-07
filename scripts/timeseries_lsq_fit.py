# -*- coding: iso-8859-1 -*-
'''
timeseries_lsq_fit  

Date: 22 April 2016
Purpose: Least squares fit to time-series: Regression of var1 using var2 and var3    
Author: Mike Bell   

'''

def prt_mean_std ( var, cvar, file_out ) : 
   import numpy as np
   var_mean = np.mean(var)  
   var_std  = np.std(var)
   print ( ' variable = ', cvar, '; mean = ', var_mean, '; standard devn = ', var_std, file = file_out) 
   return

def timeseries_lsq_fit():

   import sys
   import math         
   import numpy as np
   import matplotlib.pyplot as plt
   import numpy.linalg as linalg
   from timeseries_read_write import timeseries_read, timeseries_write
   from timeseries_running_mean import timeseries_running_mean
   from timeseries_extended_mean import timeseries_extended_mean

   ll_CMIP6 = True 
   
   if ll_CMIP6 : 
      expt_id1 = 'HadGEM3-GC31-LL_1pctCO2'   #  _1pctCO2 or _piControl
      expt_id2 = expt_id1 
      expt_id3 = expt_id1 
      infolder='../Time_series/CMIP6/'
      outfolder='../lsq_fit_plots/CMIP6/'
   else : 
# decided to change to use u-ak306 as experiment name; renaming files from nemo_ak306 to u-ak306  ; era5 ; DEEPC_v05 ; en4
      expt_id1='DEEPC_v05' #'DEEPC_v05'    # 'u-ak306' # 'u-am120' (4xCO2),  'u-ar766' (control run)  
      expt_id2='era5' #'era5'
      expt_id3='ORAS5' #'en4' 
      infolder='../Time_series/'
      outfolder='../lsq_fit_plots/'+expt_id1_2+'/'

   expt_id1_2  = expt_id1+'_'+expt_id2
   expt_id_123 = expt_id1+'_'+expt_id2+'_'+expt_id3

   l_plot_ts  = False    #   True to output timeseries plots 
   l_plot_lsq_fit = False    #   True to output least square fit plots  
   
   l_plot = l_plot_ts or l_plot_lsq_fit

# Variables 'hfds', 'tauuo', 'tos', 'depsq', 'taux_dt_dz', 't20d_sq'. 
   c_var1 = 'hfds'
   c_var2 = 'tauuo'  
   c_var3 = 'tos'    # 't20d' or 'tos'
   
   c_lat1='_-5.0-lat-5.0'  
   c_lat2='_-3.0-lat-3.0'  
   c_lat3='_-5.0-lat-5.0'  

   c_lon1='_130.0-lon--70.0_'  
   c_lon2='_130.0-lon--70.0_'
   c_lon3='_135.0-lon-137.0_'   
   c_lon3 = c_lon1 


   ts_type1 = 'integral'   #    anom, mean, integral integral_anom
   ts_type2 = 'mean'       #    anom, mean, integral mean_anom
   ts_type3 = 'mean'       #    anom, mean, integral mean_anom

   l_sq_var2 = False
   l_sq_var3 = False
   

# name of file output
   area_var_name = c_var1+'_'+c_lat1+c_lon1+ts_type1+'_'+c_var2+'_'+c_lat2+c_lon2+c_var3+'_'+c_lat3+c_lon3  

   ll_ann_mean = True    
   if ll_ann_mean: 
      y_m = 'y_06'    # 'y_06' '1s_06'
      ch_x_axis = 'Years'
      itim_st = 1993         ; itim_end = 2013            # 1950, 1985 ;  2048, 2016    ; 1851, 2348
      itim_file_st = itim_st ; itim_file_end = itim_end   

      if ll_CMIP6 :
         itim_st =   0         ; itim_end = 40
         itim_file_st = 0    ; itim_file_end = 139

   else : 
      y_m = 'm' 
      ch_x_axis = 'Months'
#      itim_st = 11 ; itim_end = 1198    #  479 or 60       # this aligns the months with the annual means I had been using (to check the results)  
#      itim_st = 0 ; itim_end = 1199    #  479 or 60    
      itim_file_st = 1 ; itim_file_end = 720   
      itim_st = 6      ; itim_end = 720      # 1950 - 2009 (50 years)     # itim_st is very important 3, 6, 9 or 12 

   str_times_in  = str(itim_file_st)+'_'+str(itim_file_end)
   str_times_out = str(itim_st)+'_'+str(itim_end)
   infilestub = '_1'+ y_m +'_'+str_times_in

   ll_legend = False   # true implies legend is written on plots (tends to overwrite the timeseries ...)  

# do not choose both ll_extended_mean and ll_running_mean = True 

   ll_running_mean = False
   len_time_window = 0    # integer number 0, 5, 10, 20 

   ll_extended_mean = False   
   len_average = 0

   ll_write_files = True # write out extended mean files to aid debugging 
   
#------------------------------------------------------------------------------------------
# end of user inputs 
#------------------------------------------------------------------------------------------
   u_opt = ''  # 'u-' o''
   if ll_CMIP6 : 
      list_infolder = 3*[ infolder ] 
   else : 
      infu =  infolder+u_opt
      list_infolder = [ infu+expt_id1, infu+expt_id2, infu+expt_id3 ]  

   file_var1=list_infolder[0]+'/'+expt_id1+infilestub+'_'+c_var1+c_lat1+c_lon1+ts_type1    # +'u-' optional 
   file_var2=list_infolder[1]+'/'+expt_id2+infilestub+'_'+c_var2+c_lat2+c_lon2+ts_type2    # +'u-' optional 
   file_var3=list_infolder[2]+'/'+expt_id3+infilestub+'_'+c_var3+c_lat3+c_lon3+ts_type3    # +'u-' optional 

# labels for plots generated
   if l_plot : 
      label_var1 = c_var1 +'_'+ c_lat1 + '_' + c_lon1   
      label_var2 = c_var2 +'_'+ c_lat2 + '_' + c_lon2   
      label_var3 = c_var3 +'_'+ c_lat3 + '_' + c_lon3

#------------------------------------------------------------------------------------------

   file_nam_plots = expt_id_123+'_1'+y_m+'_'+str_times_out+'_'+area_var_name 
   file_nam_123  = expt_id_123+'_1'+y_m+'_'+str_times_out+'_'+area_var_name 
   
   if ll_running_mean: 
      file_nam_plots = file_nam_plots + 'tw_'+str(len_time_window) + '_'
   elif ll_extended_mean :    
      file_nam_plots = file_nam_plots + 'em_'+str(len_average) + '_' 


   if ll_CMIP6 : 
      file_std = open( '../lsq_fit_plots/CMIP6/' + file_nam_plots + 'lsq_fit.txt','w') 
   else : 
      file_std = open( '../lsq_fit_plots/'+expt_id1_2+'/' + file_nam_plots + 'lsq_fit.txt','w') 
   
   print ( ' c_var1 = ', c_var1, '; c_var2 = ', c_var2, '; c_var3 = ', c_var3 , file = file_std )
   print ( ' l_sq_var2 = ', l_sq_var2, ' l_sq_var3 = ', l_sq_var3, file = file_std )
   if l_plot : print ( ' label_var1 = ', label_var1, '; label_var2 = ', label_var2, '; label_var3 = ', label_var3, file = file_std )
   print ( ' area_var_name = ', area_var_name, '; ll_ann_mean = ', ll_ann_mean, file = file_std ) 
   print ( ' itim_st, itim_end = ', itim_st, itim_end, file = file_std ) 
   print ( ' ll_running_mean= ', ll_running_mean, '; len_time_window = ', len_time_window,file = file_std ) 
   print ( ' ll_extended_mean= ', ll_extended_mean, '; len_average = ', len_average,file = file_std ) 
   
###----------------------------------------------------------------------------------
   
# 1. Read the time-series into numpy arrays of length itim_end - itim_st + 1 

   var1_ts  = timeseries_read(file_var1+'.txt',  itim_st, itim_end, file_std ) 
   var2_ts  = timeseries_read(file_var2+'.txt',  itim_st, itim_end, file_std ) 
   var3_ts  = timeseries_read(file_var3+'.txt',  itim_st, itim_end, file_std ) 

   if l_sq_var2 : 
      print ( ' squaring var2 ' ) 
      var2_ts = var2_ts * var2_ts 

   if l_sq_var3 : 
      print ( ' squaring var3 ' ) 
      var3_ts = var3_ts * var3_ts 
   
   mean_var1 = np.mean ( var1_ts) 
   mean_var2 = np.mean ( var2_ts) 
   mean_var3 = np.mean ( var3_ts)
   
   print ( ' np.size(var1_ts) =', np.size(var1_ts) , file = file_std) 


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

# 2. Plot out the timeseries    

   
   if l_plot_ts :  

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
      plt.xlabel(ch_x_axis, fontsize=15)
      plt.savefig(infolder+file_nam_plots+'timeseries.png') 
      plt.close('all') 

   prt_mean_std ( var1_ts, 'var1', file_std ) 
   prt_mean_std ( var2_ts, 'var2', file_std ) 
   prt_mean_std ( var3_ts, 'var3', file_std ) 
   
# 3. Do least square fits 
   
# 3.1 Build the matrix   var1_ts = var2_ts . x_0 + var3_ts . x_1 + const . x_2 + epsilon  

   ntimes = np.size( times )
   print ( ' ntimes =', ntimes , file = file_std) 
  
   matrix = np.zeros ( (ntimes,2) )     

   for ich_fit in ['var2_var3', 'var2_only', 'var3_only', 'full'] : 
    
      if ich_fit == 'var2_var3' : 
           matrix[:,0] =  var2_ts 
           matrix[:,1] =  var3_ts   
      elif ich_fit == 'var2_only' :
           matrix[:,0] =  var2_ts
           matrix[:,1] =  mean_var2         # using 1.0 here gave incorrect results  
      elif ich_fit == 'var3_only' :
           matrix[:,0] =  var3_ts 
           matrix[:,1] =  mean_var3          # using 1.0E17 here gave incorrect results    
      elif ich_fit == 'full' : 
           matrix = np.zeros ( (ntimes,3) ) 
           mean_prod = math.sqrt( abs ( mean_var2 * mean_var3 ) ) 
           matrix[:,0] =  var2_ts 
           matrix[:,1] =  var3_ts               
           matrix[:,2] =  mean_prod            

# 3.2 Find the least squares solution 
      x_soln = linalg.lstsq( matrix, var1_ts, rcond = None )[0] 
      
      var1_soln = np.zeros(ntimes) 
      var1_soln[:] =  matrix[:,0] * x_soln[0] + matrix[:,1] * x_soln[1] 
      if ich_fit == 'full' : var1_soln[:] = var1_soln[:] + matrix[:,2] * x_soln[2]
      var1_err = var1_ts - var1_soln
      mean_var1_err_sq = np.mean( var1_err * var1_err) 
      mean_var1 = np.mean( var1_ts ) 
      var1_var  =  var1_ts - mean_var1
      mean_var1_var_sq = np.mean( var1_var * var1_var ) 

      if ll_write_files : 
         if ll_CMIP6 : 
            file_out = infolder + file_nam_123 + ich_fit +'.txt' 
         else : 
            file_out = infolder +u_opt+expt_id1+'/'+ file_nam_123 + ich_fit +'.txt' 
         outfile=open(file_out,'w')     
         timeseries_write (outfile, 'least squares fit', var1_ts, itim_st, itim_end)    #  was itim_file_st, itim_file_end
      
      if ich_fit == 'var2_only' : 
        offset = mean_var2 * x_soln[1]
        if l_plot : c_out_fit = label_var2+'_only'
      elif ich_fit == 'var3_only' :
        offset = mean_var3 * x_soln[1]
        if l_plot : c_out_fit = label_var3+'_only'
      elif ich_fit == 'full' :
        offset = mean_prod * x_soln[2]
        if l_plot : c_out_fit = 'full'
      elif ich_fit == 'var2_var3' :
        if l_plot : c_out_fit = c_var2+'_'+c_var3
   
      ratio =  mean_var1_err_sq / mean_var1_var_sq
      print ( ' ich_fit =  ', ich_fit, file = file_std ) 
      print ( '    x_soln[0]  = ', x_soln[0], file = file_std ) 
      print ( '    x_soln[1]  = ', x_soln[1], file = file_std ) 
      if  ich_fit == 'full' :  print ( '    x_soln[2]  = ', x_soln[2], file = file_std ) 
      print ( '    mean_var1  = ', mean_var1, file = file_std )   
      print ( '    mean square error, variance, ratio = ', mean_var1_err_sq, mean_var1_var_sq, ratio, file = file_std ) 
      if ich_fit != 'var2_var3' : 
         print ( '    offset = ', offset, file = file_std )
      if ratio < 1.0 : 
         rsq = 1.0 - ratio
         r = math.sqrt(rsq) 
         print ( '    rsq, r = ', rsq , r,  file = file_std ) 
        
# 3.3 Produce time-series and scatter plot of solution  

      tick_size = 20 

      if l_plot_lsq_fit :    
         plt.plot(times,var1_ts, linewidth=4, color='r', label=label_var1+' actual')
         plt.plot(times,var1_soln, linewidth=4, color='b', label=c_out_fit+' lsq fit') # linestyle=':', 

         plt.title("Actual and lsq fit time-series " + c_out_fit, fontsize=16) 
         plt.xlabel(ch_x_axis, fontsize=16)
         plt.xticks(fontsize=tick_size)
         plt.yticks(fontsize=tick_size)
         if ll_legend : 
            plt.legend()
         plt.savefig(outfolder + file_nam_plots + ich_fit +'_timeseries.png') 
         plt.close('all') 

         plt.scatter(var1_soln,var1_ts, linewidth=4, color='k' )
         plt.title("Scatter plot " + c_out_fit, fontsize=16) 
         plt.xlabel(label_var1+' fit', fontsize=12)
         plt.ylabel(label_var1+' actual', fontsize=12)
         plt.xticks(fontsize=tick_size)
         plt.yticks(fontsize=tick_size)
         plt.savefig(outfolder + file_nam_plots + ich_fit + '_scatter_plot.png') 
         plt.close('all') 

   print ( 'outfolder, file_nam_plots = ', outfolder, file_nam_plots ) 

# 4. Produce coherence plot of var1 and var2    

#   plt.cohere( var1_ts, var2_ts, NFFT=60,  Fs=12 ) 
#   plt.title('Coherence; ' + expt_id + ' ' + str(itim_end) , fontsize=16) 
#   plt.xlabel('frequency (cycles per year)', fontsize=12)
#   plt.ylabel('coherence', fontsize=12)
#   plt.savefig(infolder + file_nam_plots+ '_coherence_plot.png') 
#   plt.close('all') 
   
timeseries_lsq_fit()
