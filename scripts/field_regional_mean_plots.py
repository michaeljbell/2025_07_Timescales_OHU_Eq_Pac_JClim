# -*- coding: iso-8859-1 -*-
'''
field_regional_mean_plots.py 

Date: 6 Sept 2021
Purpose: Calculate and plot regional means of DEEPC or exptid net surface heat fluxes    
Author: Mike Bell   
'''

#----------------------------------------------------------------------------------------------------------

def field_regional_mean_plots() :

 import sys        
 import math as math
 import numpy as np
 import netCDF4
 import matplotlib.pyplot as plt
 from user import data_dir 
 
 config_res = 'eORCA1'

 inmesh = data_dir+'/MESH_MASK/mesh_mask_'+config_res+'-GO6.nc'  
 in_region_masks = data_dir+'/MESH_MASK/region_masks_'+config_res+'-GO6.nc'

 expt_id = 'DEEPC_v5'    # 'DEEPC' (this is v2) or 'DEEPC_v5' or 'ak306' or 'aj368' (i.e. any model expt_id) ; aj368 ; ce853 ;  ay490 

 if expt_id == 'DEEPC' : 
   hf_var = 'mean_Fmass'
   date_str =  '2000-2009_Annual'      #     '1985-1994_Annual'  
   file_in = '/data/users/frmk/DEEPC/DEEPC_SFC_NET_v02_'   + date_str + '_' + config_res+'-GO6.nc'
 elif expt_id == 'DEEPC_v5' : 
   hf_var = 'mean_Fmass'
   date_str =  '2000-2009_Annual'      #     '1985-1994_Annual'  
   file_in = '/data/users/frmk/DEEPC_v05/DEEP-C_Fmass_v5_' + date_str + '_' + config_res+'-GO6.nc'

 else :  
   hf_var = 'hfds'
   if expt_id != 'ce853' : 
     date_str = '2000-2009_Annual'
     file_in = data_dir+'u-'+expt_id +'/dec/surface_T/nemo_'+ expt_id +'o_hfds_'+ date_str + '.nc'
   else : 
     date_str = '20000101-20100101'
     file_in = data_dir+'u-'+expt_id +'/dec/surface_T/nemo_'+ expt_id+'o_'+hf_var+'_' + date_str + '_grid-T.nc'

 ll_all_labels = True  # controls labels on final plot
 label_size = 20
 
 outdir = '../region_plots/'+expt_id+'/'

# latbdys copied from classes_sections_plot.py section 1 

 latbdys = [-90.0, -60.0, -45.0, -30.0, -15.0, 15.0, 30.0, 45.0, 60.0, 90.0]   # boundaries of the latitude bands (section 11)  
 n_latbdys = len(latbdys) 
 n_latbands = n_latbdys - 1
 
 print ('latbdys, n_latbands = ', latbdys, n_latbands ) 
 
 str_latbdys = []
 for j_lat in range (n_latbands) : 
    str_latbdys.append( str(latbdys[j_lat]) )     #  strings of latitudes of inter-region boundaries

 c_basins   = ['tmaskatl', 'tmaskpac','tmaskind'] 
 n_basins = len ( c_basins ) 
 
#-------------------------------------------------------------------------------------------------------------
# 2.1 Read the mesh mask fields  
#------------------------------------------------------------------------------------------------------------- 

 g = netCDF4.Dataset(inmesh,'r')

 deptht    = g.variables["gdept_1d"][0] 
    
 nav_lat_t = g.variables["nav_lat"][:]        #    longitude is not needed
 maskt_srf = g.variables["tmask"][0,0]     
 print ( ' maskt_srf.dtype = ', maskt_srf.dtype ) 

 e1t   = g.variables["e1t"][0]
 e2t   = g.variables["e2t"][0]

 if isinstance(nav_lat_t, np.ma.MaskedArray) :  nav_lat_t = nav_lat_t.data
 if isinstance(maskt_srf, np.ma.MaskedArray) :  maskt_srf = maskt_srf.data
 if isinstance(e1t, np.ma.MaskedArray) : 
    e1t = e1t.data
    e2t = e2t.data 

 g.close()
 
#-------------------------------------------------------------------------------------------------------------
# 2.2 Read the basin mask fields 
#------------------------------------------------------------------------------------------------------------- 
 g = netCDF4.Dataset(in_region_masks, 'r')
 basin_masks = []
 for cbasin in c_basins : 
    basin = g.variables[cbasin][0]    
    if isinstance(basin, np.ma.MaskedArray) : 
       basin = basin.data
    basin_masks.append( basin )   
 print ( ' cbasin, np.shape(basin)   = ', cbasin, np.shape(basin)  ) 
 
#-------------------------------------------------------------------------------------------------------------
# 3. Read the DEEPC or model fields interpolated to the selected grid and 
#    calculate fields with land/sea mask imposed multiplied by the horizontal area  
#-------------------------------------------------------------------------------------------------------------

 g = netCDF4.Dataset(file_in,'r')
 field_in = g.variables[hf_var][0] 
 g.close()

 if isinstance(field_in, np.ma.MaskedArray) : field_in = field_in.data

 print ( ' np.shape(field_in) = ', np.shape(field_in) )   

 field_in = e1t * e2t * field_in * maskt_srf 

#-------------------------------------------------------------------------------------------------------------
# 4. Loop over the basins and latitude regions calculating sums for each region
#-------------------------------------------------------------------------------------------------------------

 area_sums = np.zeros( (n_basins, n_latbands) ) 

 for i_basin in range ( n_basins ) : 
   cbasin = c_basins[i_basin]
   field_out = field_in * basin_masks[i_basin] 

   for j_lat in range ( n_latbands ) : 
      lat_min = latbdys[j_lat]         
      lat_max = latbdys[j_lat+1]
      field_latband = np.where( nav_lat_t >= lat_min, field_out, 0.0) 
      field_latband = np.where( nav_lat_t < lat_max, field_latband, 0.0) 
      area_sums[i_basin, j_lat] = np.sum(field_latband) 	
      print (' i_basin, j_lat, area_sums[i_basin, j_lat]= ', i_basin, j_lat, area_sums[i_basin, j_lat] ) 

#-------------------------------------------------------------------------------------------------------------
# 5. Output the sums on a bar chart (based on section 15 of classes_sections_plot.py  
#-------------------------------------------------------------------------------------------------------------

# The number of basins is hard-wired in this plot 

 x = np.arange( n_latbands )
 width = 0.2
 fig, ax = plt.subplots()
 rects1 = ax.bar(x - width, area_sums[0], width, color='r', label = c_basins[0])
 rects2 = ax.bar(x        , area_sums[1], width, color='b', label = c_basins[1])
 rects3 = ax.bar(x + width, area_sums[2], width, color='g', label = c_basins[2])
 ax.set_xticks(x-0.5)
 ax.set_xticklabels(str_latbdys,rotation=30) 
 ax.tick_params(axis='both', labelsize=label_size )
 y_lims = [-0.5E15,2.1E15] 
 plt.ylim(y_lims)  
  
 if ll_all_labels :
   ax.set_ylabel('heat flux')
   ax.legend()
   ax.set_title( expt_id+' latband sum: '+ hf_var + ' ' + ' '  + date_str  )
 outfile = outdir+ expt_id+'_bar_chart_latband_' + hf_var + '_' + date_str
 if ll_all_labels : 
    outfile = outfile + '_labelled.png'
 else : 
    outfile = outfile + '_no_labels.png'
 
 plt.tight_layout()    # added to ensure axis labels fit within the plot
 plt.savefig(outfile)
 plt.close()   

#----------------------------------------------------------------------------------------------------------
field_regional_mean_plots()
