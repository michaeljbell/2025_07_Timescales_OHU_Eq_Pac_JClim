# -*- coding: iso-8859-1 -*-
'''
field_plot_difference_CMIP6 

Date: 10 April 2020
Purpose: Plot out geographical maps of difference between two fields that can be defined on different grids   
Author: Mike Bell   

'''
def field_plot_differences_CMIP6(c_model, c_descrip, c_period):

### originally based on a script provided by Daley Calvert
###
### options include: difference plots (from a constant field or one with the same time string) 
###                  annual mean or decadal plots (monthly plots could also be produced) 

   import sys        
   import numpy as np
   import netCDF4 as nc4
   from matplotlib import cm
   from matplotlib.colors import ListedColormap, to_rgba
   import matplotlib.pyplot as plt
   import matplotlib.ticker as mticker
   import cartopy
   import cartopy.img_transform
   import cartopy.crs as ccrs
   from utilities import regrid_field, set_19_distinct_colors, rd_ncdf_var_check_one
   

   var = 'hfds'    #  hfds mean_Fmass tos, sos, hfds, empmr, G, tauuo, sos, zos, corrln_hfds_tauuo, mean_t20d rhop_tnd, mean_Fmass
   varfilename=var # 'std_dec_means_'+var 
   varname = var + '' # + '_hfds' # + '' #  + '_tauu'  #  +'_std'      #   +'_std'      #  +'_pot'  #  +'o_'
   outvar= var        # net_empmr_flux hfds sst  tauuo  G, sos, zos rho_tnd

   vmax   =  190.    # contour max & min values    -190, 190, 20 for hfds; 0, 30, 16 for SST ; -0.29, 0.29, 30 for taux ; -0.038, 0.038, 20 for G 
                      #                             -5.0, 5.0, 21 for SST differences ; -95, 95, 20 for hfds differences ; -0.15, 0.15, 30 for taux diff
   vmin   = - vmax
   nsteps =   20    # number of contour levels;  contour interval = (vmax - vmin) / (nsteps - 1) 

   ll_color_sequential = False  # False implies 19 color map with white in middle is used 
   if ll_color_sequential: 
      tab_sequential = 'coolwarm'   # 'jet'
   else : 
      white_range = 0          # number of central values that are whited out 


   figsize=(8, 5)
   label_size = 15   # 10 is original # 14 used for figure 2 ; 12 for figure 1 

   ll_include_title_lats_lons = False

   multiplication_factor = 1.0    # can be used to avoid horrible labels on colour bars  

   c_mdp = c_model +'_' + c_descrip +'_' + c_period
   
   c_rf = 'rf_5_'   #  '' or 'rf_5_' 
   
   c_mth_start = '06'
   
   cyr_st  = '00'     # usually 1960 for era5  ;   1950 or 1980 for ocean_model ; 1985 for DEEPC;  1850 for ar766
   iyr_st = int(cyr_st) 
   iyr_end = 499   # the last year used IS iyr_end ; usually 2019 for era5; 2049 for ocean_model ; 2016 for DEEPC; 2348 for ar766

   c_per_mean = cyr_st +c_mth_start+'01-'+str(iyr_end)+c_mth_start+'01'

   c_ts = '_hfds_-5.0-lat-5.0_130.0-lon--70.0' #   '_tauu_-3.0-lat-3.0_135.0-lon--75.0_'  # '' # '_hfds_-5.0-lat-5.0_130.0-lon--70.0'
   varname = varname #    + '_tauu'

   infolder='/data/users/frmk/CMIP6/'+ c_model +'/'  #         # dec/ 
   infilestub =  varfilename+'_' + c_mdp + '_' + c_rf +  c_per_mean   + c_ts    
   
#   infilestub = 'hfds_HadGEM3-GC31-LL_1pctCO2_008001-008912'  


   l_pcolormesh = False    # if fields are very noisy you need to use pcolormesh to see the negative values

   central_longitude = 180.0  # 0.0 or 180. 

   ( lon_min, lon_max, lat_min, lat_max ) = (-180, 180, -90, 90)      # extents of the plots 
#   ( lon_min, lon_max, lat_min, lat_max ) = (-180, 180, -15, 15)

# see https://scitools.org.uk/cartopy/docs/latest/crs/projections.html for options
# class cartopy.crs.Mercator(central_longitude=0.0, min_latitude=-80.0, max_latitude=84.0)                  
   projection_out=ccrs.PlateCarree(central_longitude=central_longitude)    # PlateCarree - for output plot

   difference_plot = False
   c_sum_diff      = 'diff'        # 'sum' or 'dif' 
      
   if difference_plot : 
      infolder2 ='/data/users/frmk/DEEPC_v05/ony_06/'
      infilestub2='DEEPC_v05_1y_20000601-20010601'     
      varname2 = 'hfds' 

   outdir='../Surface_plots/'+c_model+'/'     

### end of user inputs 
###----------------------------------------------------------------------------------

   projection_in=ccrs.PlateCarree()    # PlateCarree - for input fields (lat-lon grid) 
   
   for cseason in [ 'Annual'] :     # a dummy loop (to avoid changing existing spacing)  

      infilename = infilestub 
      infile=infolder+infilename+'.nc'
      
      print ('infile = ', infile)
      
      data = nc4.Dataset(infile,'r')
      lat1   = rd_ncdf_var_check_one(data, 'nav_lat')
      lon1   = rd_ncdf_var_check_one(data, 'nav_lon')
      field1 = rd_ncdf_var_check_one(data,   varname)   
      data.close()

      print ( ' np.shape(field1) = ', np.shape(field1), field1.dtype ) 
      print ( 'field1_min = ', np.min(field1), 'field1_max = ', np.max(field1) )
            
      if difference_plot : 

         infilename2 = infilestub2 
         infile2=infolder2+infilename2+'.nc'  
         print ( 'infile2 = ', infile2 )  

         g = nc4.Dataset(infile2,'r')
         lat2   = rd_ncdf_var_check_one(g, 'nav_lat')
         lon2   = rd_ncdf_var_check_one(g, 'nav_lon')
         field2 = rd_ncdf_var_check_one(g,  varname2)
         g.close()
         print ( ' np.shape(field2) = ', np.shape(field2) ) 

# Set the figure size
      fig= plt.figure(figsize=figsize)

# set the axes and calculate the coordinates and regridded field on the output projection   
      ax = plt.axes(projection=projection_out)
      ax.set_extent( [lon_min, lon_max, lat_min, lat_max], projection_out )    
      
      x, y, field1_rgd = regrid_field(ax, fig, lon1, lat1, field1, projection_in, lon_min, lon_max, lat_min, lat_max, projection_out )
      if difference_plot : 
         x, y, field2_rgd = regrid_field(ax, fig, lon2, lat2, field2, projection_in, lon_min, lon_max, lat_min, lat_max, projection_out )

      if difference_plot: 
         if c_sum_diff == 'sum' : 
            field_rgd  = field1_rgd + field2_rgd 
         else : 
            field_rgd  = field1_rgd - field2_rgd 

      else : 
         field_rgd = field1_rgd

      print ( ' field_rgd_max, field_rgd_min = ', np.amax(field_rgd), np.amin(field_rgd) ) 

      field_rgd = field_rgd * multiplication_factor

      if ll_color_sequential :  
         tab20b = cm.get_cmap(tab_sequential, nsteps)
         newcolors = tab20b(np.linspace(0,1,nsteps))
         white = to_rgba('w')
         newcolors[0,:]=white
         newcmp=ListedColormap(newcolors) 
      else : 
          newcmp=set_19_distinct_colors(white_range)
      
      levels = np.linspace( vmin, vmax, nsteps) 
      if l_pcolormesh : 
         c = plt.pcolormesh(x, y, field_rgd, vmin=vmin, vmax=vmax,  cmap=newcmp) # 'tab20b'
      else : 
         c = plt.contourf(x, y, field_rgd, vmin=vmin, vmax=vmax, levels = levels , cmap=newcmp)  # 'tab20b'
      
# Add coastlines, contour labels and a colour bar
      plt.gca().coastlines()
 #     plt.clabel(c, inline=False, colors='k')
      cb=plt.colorbar(c, orientation='horizontal', extend='both')
      cb.ax.tick_params(labelsize = label_size) 
 
# Add latitude and longitude labels and then matching gridlines
      if ll_include_title_lats_lons :
         plt.xlabel("latitude", fontsize=label_size)
         plt.ylabel("longitude", fontsize=label_size)
      ax.set_xticks(np.linspace(lon_min, lon_max, 7), crs=projection_out)
      ax.set_yticks(np.linspace(lat_min, lat_max, 7), crs=projection_out) 
      ax.tick_params(axis='both', labelsize=label_size )

      xlabels = np.linspace(lon_min, lon_max, 7)
      if central_longitude == 180.0 :
         for ilab in range (7) :   
            xlabels[ilab] = xlabels[ilab] + 180
            if xlabels[ilab] > 180 : xlabels[ilab] = xlabels[ilab] - 360
      ylabels = np.linspace(lat_min, lat_max, 7)

      lat_lon_labelsize = 15 
      ax.set_xticklabels(xlabels,fontsize=lat_lon_labelsize)
      ax.set_yticklabels(ylabels,fontsize=lat_lon_labelsize)

      g1 = ax.gridlines(color='k')
      g1.xlocator = mticker.FixedLocator( np.linspace(lon_min, lon_max, 7) )
      g1.ylocator = mticker.FixedLocator( np.linspace(lat_min, lat_max, 7) )

# set the filename 

      if difference_plot : 
         if c_sum_diff == 'sum' : 
            outfile=outdir+'_'+infilename+'_sum_'+varname+'_'+infilename2
         else : 
            outfile=outdir+'_'+infilename+'_min_'+varname+'_'+infilename2
      else : 
         outfile=outdir+'_'+infilename
	 
      outfile = outfile +'_titles_'+str(central_longitude)+'_'+str(ll_include_title_lats_lons)+'.png'  

      title = infilename  

      if difference_plot : 
         title2 = infilename2  

         if c_sum_diff == 'sum' : 
            title  = title + ' plus ' + title2 
         else  : 
            title  = title + ' minus ' + title2 

      if multiplication_factor != 1.0 : 
         title = title + ' multiplied by ' + str(multiplication_factor)

      if ll_include_title_lats_lons :
         plt.title(title, fontsize=8) 

# can just hack in the output filename here 
#      outfile = outdir+'nemo_ak306o_hfds_1950_1960_'+simth+'_'+simthp1+'.png'	 
      print ( 'outfile = ', outfile ) 
      plt.savefig(outfile)
      plt.close()   

def run_field_plot_differences_CMIP6():   

   model_time = []
#   model_time.append( ['ACCESS-CM2','piControl_r1i1p1f1_gn', '095001-144912', 6000] ) 
#   model_time.append( ['CanESM5','piControl_r1i1p1f1_gn', '520101-620012', 12000] )
#   model_time.append( ['GFDL-CM4','piControl_r1i1p1f1_gn', '015101-065012', 6000 ] )
#   model_time.append( ['IPSL-CM6A-LR','piControl_r1i1p1f1_gr', '185001-384912', 18000 ] )    # tau lat, lon is missing 
#   model_time.append( ['MIROC6','piControl_r1i1p1f1_gn', '320001-399912', 6000 ] )
#   model_time.append( ['MPI-ESM1-2-LR','piControl_r1i1p1f1_gn', '185001-284912', 12000 ] )
#   model_time.append( ['NorESM2-LM','piControl_r1i1p1f1_gn', '160001-210012', 6000] )     # actually has 6012 
   model_time.append( ['HadGEM3-GC31-LL','piControl_r1i1p1f1_gn', '225001-384912'] )
#   model_time.append( ['HadGEM3-GC31-MM','piControl_r1i1p1f1_gn', '185001-234912'] )

   for clist in model_time : 

      c_model   = clist[0] 
      c_descrip = clist[1]
      c_period  = clist[2]

      field_plot_differences_CMIP6(c_model, c_descrip, c_period)


run_field_plot_differences_CMIP6()
