# -*- coding: iso-8859-1 -*-
'''
field_plot_difference  

Date: 10 April 2020
Purpose: Plot out geographical maps of difference between two fields that can be defined on different grids   
Author: Mike Bell   

'''
def field_plot_differences():

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
   
   data_type = 'DEEPC' # 'ocean_expt', 'era5', 'DEEPC' , 'en4', cmems, 'atmos_expt'
   grid_type = 'T'         # 'T' , 'U' , 'V', 'F' 

#   varfilename='o_'          #  or varfilename = 'o_'+ varname + '_' 
   var = 'corrln_hfds'    #  hfds mean_Fmass tos, sos, hfds, empmr, G, tauuo, sos, zos, corrln_hfds_tauuo, mean_t20d rhop_tnd, mean_Fmass
   varfilename= var      # 'std_dec_means_'+var  'corrln_tos_tauuo_' var+'_10_std_'+var ; 'std_'+var 

   var_ts = '_tauuo'  # ''  or '_tauuo'
   c_ts   = '_-3.0-lat-3.0_130.0-lon--70.0' # '_-3.0-lat-3.0_130.0-lon--70.0'   # '' 

   varname =  var + var_ts  #  var  ;   +'_std'  # 'anom_'+var # 'std_' + var      #  +'_pot'  #  +'o_'  used to pick out the field 
   outvar  =  var + var_ts  # var or ''  
   
   c_rf =     ''        #   '_rf'   #  ''
   c_rf_len = ''        #   '_tw_5' #  '' 

   vmax   =  0.95  # contour max & min values    -190, 190, 20 for hfds; 0, 30, 16 for SST ; -0.29, 0.29, 30 for taux ; -0.038, 0.038, 20 for G 
#                                                   -5.0, 5.0, 21 for SST differences ; -95, 95, 20 for hfds differences ; -0.15, 0.15, 30 for taux diff
   vmin   = - vmax
#   vmin   =  0.0 
   nsteps =   20      # number of contour levels;  contour interval = (vmax - vmin) / (nsteps - 1) 

   ll_title = False

   ll_convert_Kelvin_to_Centigrade = False #  this is a hack 

   multiplication_factor = 1.0    # used to avoid horrible labels on colour bars  

# Contour the data
# for diverging colormaps, cmap=, 'BrBG' and 'RdBu' are good options. 'coolwarm' is not bad; 'bwr'
# cmaps['Sequential'] = [
#            'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
#            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
#            'GnBu', 'PuBu', 'YlGnBu', 'PuGn', 'BuGn', 'YlGn']

   ll_color_sequential = False  # False implies 19 color map with white in middle is used 
   if ll_color_sequential: 
      tab_sequential = 'jet'   # 'jet'
   else : 
      white_range = 0          # number of central values that are whited out 
      
   pos_neg = 'pos' 
   str_interval = '' #  var_ts + c_ts + '_mean_'+pos_neg+'_grid-' + grid_type # or '' 

   if data_type == 'ocean_expt' : 
      expt_id = 'ak306'   # ak306 am120 ay490 aj368 ap409  era5_pm30
      infolder='/data/users/mike.bell/CMIP6HighResMIP/u-'+expt_id+'/ony_06'+c_rf+'/surface_'+grid_type+'/'  
      infilestub = 'nemo_'+ expt_id  + 'o_1y_' + varfilename + c_rf_len + c_ts+ '_19500601-20490601' + str_interval + '.nc'  #  _grid-T     
      l_filestub  = True   # True => filename  is set from filestub (no dates appended to it)    var = 'corrln_hfds_hfds'    #  hfds mean_Fmass tos, sos, hfds, empmr, G, tauuo, sos, zos, corrln_hfds_tauuo, mean_t20d rhop_tnd, mean_Fmass

   elif data_type == 'atmos_expt' : 
      expt_id = 'bhENS'   # ak306 am120 ay490 aj368 ap409  era5_pm30
      infolder='/data/users/mike.bell/CMIP6AMIP/u-'+expt_id+'/anm/surface_'+grid_type+'/'    # dec or ony # 'surface_'+grid_type+'/'
      infilestub = expt_id +'_1m_'  + varfilename +'_19601101-21791201_grid-'+grid_type+'.nc'  #     + outvar+'_'  ; _grid-'+grid_type+'
      l_filestub  = True   # True => filename  is set from filestub (no dates appended to it) 

   elif data_type == 'era5' : 
      expt_id = 'era5' 
      infolder='/data/users/mike.bell/ERA5_pm90/ony_06/'     # ony or onm  _pm30 is optional 
      infilestub = 'era5_1y_' +varfilename + var_ts + c_ts +'_19600601-20190601_pm90.nc' # +outvar+'_'       #   -T_pm30 
      l_filestub  = True   # True => filename  is set from filestub (no dates appended to it) 

   elif data_type == 'DEEPC' : 
      expt_id = 'DEEPC' 
      infolder='/data/users/mike.bell/DEEPC_v05/ony_06/'     # ony or onm
      infilestub = 'DEEPC_v05_1y_'+varfilename + var_ts + c_ts +'_19850601-20170601.nc'    # _grid-T
      l_filestub  = True   # True => filename  is set from filestub (no dates appended to it) 

   elif data_type == 'en4' : 
      expt_id = 'en4' 
      infolder='/data/users/mike.bell/EN4/ony_06/2D/'     # ony or onm
      infilestub = 'EN.4.2.2.f.analysis.g10.'+varfilename+ c_ts +'_19850601-20170601.nc'   # temperature or t20d
      l_filestub  = True   # True => filename  is set from filestub (no dates appended to it) 

   elif data_type == 'cmems' : 
      expt_id = 'cmems'
      var = 'eastward_stress_bias' 
      infolder='/data/users/mike.bell/'     # ony or onm
      infilestub = 'cmems_obs-wind_glo_phy_nrt_l4_0.125deg_PT1H_1669307499197.nc'
      l_filestub  = True   # True => filename  is set from filestub (no dates appended to it) 


#   expt_id = 'am916' 
#   infolder='/data/users/mike.bell/CMIP6HighResMIP/u-'+expt_id+'/instant/rhop_tnd/'    # dec or cen
#   infilestub = 'u-'+expt_id+'o_1m_'+varfilename + '_reg_623_745_19991201_19991231_grid-T'  
#   infilestub = 'u-'+expt_id+'o_1m_'+varfilename + '_reg_0_1205_19991201_19991231_instant_grid-T'  
 
#   expt_id = 'DEEPC'
#   infolder ='/data/users/mike.bell/DEEPC/'
#   infilestub='DEEPC_SFC_NET_v02_'     #   2000-2010'    # _MERRA_v02_

#   expt_id = 'DEEPC'
#   infolder ='/data/users/mike.bell/DEEPC/'
#   infilestub='DEEPC_SFC_NET_v02_'     #   2000-2010'    # _MERRA_v02_

   one_plot=True 
   annual_plots = True
   iyr_dec_step = 5 
   if one_plot : 
      iyr_dec_first  = 1960      # 1960 or 1950 
      iyr_dec_end    = 2019      # 2019 or 2049 or 2050
      iyr_dec_step   = iyr_dec_end + 1 - iyr_dec_first

   elif annual_plots : 
      iyr_dec_first  = 1960
      iyr_dec_end    = 1990    # actually the year after the last one 

   else :    # decadal plots 
      iyr_dec_first  = 195
      iyr_dec_end    = 2010

   l_pcolormesh = False    # if fields are very noisy you need to use pcolormesh to see the negative values

   field1_3d = False
   k_lev = 0   # only relevant if field1_3d or field2_3d = True     
   var_depth = 'depth'    # to read depths of levels 

   central_longitude = 180.0  # 0.0 or 180. 

   ( lon_min, lon_max, lat_min, lat_max ) = (-180, 180, -90, 90)      # extents of the plots 
#   ( lon_min, lon_max, lat_min, lat_max ) = (-180, 180, -30, 30)
#   ( lon_min, lon_max, lat_min, lat_max ) = (130, -70, -30, 30)        # these are longitudes relative to greenwich 
#   ( lon_min, lon_max, lat_min, lat_max ) = (-60, 120, -30, 30)     # these are longitudes relative to 180.  


# see https://scitools.org.uk/cartopy/docs/latest/crs/projections.html for options
# class cartopy.crs.Mercator(central_longitude=0.0, min_latitude=-80.0, max_latitude=84.0)                  
   projection_out=ccrs.PlateCarree(central_longitude=central_longitude)    # PlateCarree - for output plot

#  lsm_overlay_plot is intended to read in a land-sea mask and multiply field1 by it. Do not set difference_plot and lsm_overlay_plot = True !
   difference_plot = False
   lsm_overlay_plot = True
   c_sum_diff = 'diff'                # 'sum' or 'dif' 

   if difference_plot and lsm_overlay_plot :
      print ('  stopping because both difference_plot and lsm_overlay_plot are set to True; ') 
      sys.exit()
      
   if difference_plot or lsm_overlay_plot : 
      
#      l_filestub2 = True   # True => filename2 is set from filestub2 (for the sum/difference case: difference_plot = True) 
#      expt_id2 = 'am916'
#      infolder2= '/data/users/mike.bell/CMIP6HighResMIP/u-'+expt_id2+'/instant/rhop_tnd/'
#      infilestub2 = 'u-'+expt_id2+'o_1m_'+varfilename+'_reg_623_745_19991201_19991231_grid-T' 
#      infilestub2 = 'u-'+expt_id2+'o_1m_'+varfilename + '_reg_0_1205_19991201_19991231_instant_grid-T'  
#      varname2 = 'rhop_tnd_sal'

      l_filestub2  = True
      data_type2   = 'nav_lat'    # 'cmems' reads in lat and lon differently 
      field2_3d    = False     #  this is the default 

      expt_id2 = 'DEEPC_MERRA'                                 # 'DEEPC' or 'DEEPC_MERRA'
      infolder2 ='/data/users/mike.bell/DEEPC/'
      infilestub2='DEEPC_SFC_NET_MERRA_v02_'    #  'DEEPC_SFC_NET_v02_' or 'DEEPC_SFC_NET_MERRA_v02_'
      varname2 = 'mean_Fmass' 

      expt_id2 = 'ak306'   # ak306 am120 ay490 aj368
      infolder2='/data/users/mike.bell/CMIP6HighResMIP/u-'+expt_id2+'/ony/surface_'+grid_type+'/'     # dec or ony
      infilestub2 = 'nemo_'+expt_id2 +'o_1y_'+varfilename+'_'     #   1y_'  'o_hfds_' or '_hfds_'    
      varname2 = varname 

      expt_id2 = 'cmems'
      infolder2='/data/users/mike.bell/'     # ony or onm
      infilestub2 = 'cmems_obs-wind_glo_phy_nrt_l4_0.125deg_PT1H_1669307499197'
      varname2 = 'eastward_stress_bias' 

      expt_id2 = 'ERA5'                                 # 'DEEPC' or 'DEEPC_MERRA'
      infolder2 ='/data/users/mike.bell/ERA5_pm90/'
      infilestub2='mesh_mask'    #  'DEEPC_SFC_NET_v02_' or 'DEEPC_SFC_NET_MERRA_v02_'
      varname2 = 'tmask' 

      expt_id2 = 'CMIP6HighResLSM'   # ak306 am120 ay490 aj368
      infolder2='/data/users/mike.bell/CMIP6HighResMIP/MESH_MASK/'     # dec or ony
      infilestub2 = 'mesh_mask_eORCA025-GO6'    #   1y_'  'o_hfds_' or '_hfds_'    
      varname2 = 'tmask' 
      field2_3d = True

   if difference_plot : 
      outdir='../Surface_plots/'+outvar+'/'+expt_id+'_'+expt_id2+'/'     
   else : 
      outdir='../Surface_plots/'+outvar+'/'+expt_id +'/'    

### end of user inputs 
###----------------------------------------------------------------------------------

   projection_in=ccrs.PlateCarree()    # PlateCarree - for input fields (lat-lon grid) 
   
   for iyr_dec in range (iyr_dec_first, iyr_dec_end, iyr_dec_step ) : 

#    for imth in range (1):    
    for cseason in [ 'Annual'] :     # , 'Spring', 'Summer', 'Fall', 'Winter' ] : 
      if one_plot :
         iyr  = iyr_dec_first
         iyp1 = iyr_dec_end 
      elif annual_plots: 
         iyr  = iyr_dec
         iyp1 = iyr_dec + 1
      else : 
         iyr  = iyr_dec*10
         iyp1 = (iyr_dec+1)*10

      siyr = str(iyr )
      siyp1 = str(iyp1)
      if l_filestub :
        infilename = infilestub 
      else : 
        infilename = infilestub +siyr+'1001-'+siyp1+'1101'+str_interval +'_pm30.nc' # + '.nc' # 'grid-'+grid_type+'.nc'
	
#        infilename = infilestub +siyr+'-'+siyp1+'_'+cseason
      
# can just hack in the infile  and infilename (used to name the output file) here 
#      infilename = 'DEEPC_v05_1y_hfds_anom_hfds_1985_2016_tauuo_-3.0-lat-3.0_130.0-lon--70.0_mean_'+pos_neg+'_grid-T.nc'

      infile=infolder+infilename
      print ('infile = ', infile)
      
      data = nc4.Dataset(infile,'r')
      if data_type == 'cmems' :
         lat = data.variables['lat'][:] 
         lon = data.variables['lon'][:] 
         lon1, lat1 = np.meshgrid( lon, lat )    
      else : 
         lon1 =   rd_ncdf_var_check_one (data, 'nav_lon')         
         lat1 =   rd_ncdf_var_check_one (data, 'nav_lat')         

      field1 =   rd_ncdf_var_check_one (data, varname)         
      if field1_3d : 
         field1 = field1[k_lev]   
         depth1 = data.variables[var_depth][0][k_lev]
      data.close()
      print ( ' np.shape(lat1), np.shape(field1) = ', np.shape(lat1), np.shape(field1), field1.dtype ) 
      if field1_3d : print ( 'k_lev, depth1 = ', k_lev, depth1 ) 
      print ( 'field1_min = ', np.min(field1), 'field1_max = ', np.max(field1) )

      if ll_convert_Kelvin_to_Centigrade :
        field1 = field1 - 273.15   # original data are in Kelvin 
            
      if difference_plot or lsm_overlay_plot : 

         if l_filestub2 :
            infilename2 = infilestub2 
         else : 
            infilename2 = infilestub2 +siyr+'1201-'+siyp1+'1201_grid-'+grid_type+''   
#            infilename2 = infilestub2 +siyr+'-'+siyp1+'_'+cseason
         infile2=infolder2+infilename2+'.nc'  
         print ('infile2 = ', infile2)

         g = nc4.Dataset(infile2,'r')
         if data_type2 == 'latitude' :
            lat = g.variables['latitude'][:] 
            lon = g.variables['longitude'][:] 
            lon2, lat2 = np.meshgrid( lon, lat )    
         elif data_type2 == 'lat' :
            lat = g.variables['lat'][:] 
            lon = g.variables['lon'][:] 
            lon2, lat2 = np.meshgrid( lon, lat )    
         else : 
            lat2 = g.variables['nav_lat'][:]
            lon2 = g.variables['nav_lon'][:]

         if field2_3d : 
            field2_all_levels = rd_ncdf_var_check_one(g,varname2) 
            field2 = field2_all_levels[k_lev]
         else : 
            field2 = rd_ncdf_var_check_one(g,varname2)
         g.close()
         print ( ' np.shape(field2) = ', np.shape(field2) ) 

         if ll_convert_Kelvin_to_Centigrade :
            field2 = field2 - 273.15   # original data are in Kelvin 


# Set the figure size
      fig= plt.figure(figsize=(8, 5))
      label_size = 15   # 10 is original # 14 used for figure 2 ; 12 for figure 1 
      ll_include_title_lats_lons = False

# set the axes and calculate the coordinates and regridded field on the output projection   

      ax = plt.axes(projection=projection_out)
      ax.set_extent( [lon_min, lon_max, lat_min, lat_max], crs=projection_out )    
      
      x, y, field1_rgd = regrid_field(ax, fig, lon1, lat1, field1, projection_in, lon_min, lon_max, lat_min, lat_max, projection_out )
      if difference_plot or lsm_overlay_plot : 
         x, y, field2_rgd = regrid_field(ax, fig, lon2, lat2, field2, projection_in, lon_min, lon_max, lat_min, lat_max, projection_out )

      if difference_plot: 
         if c_sum_diff == 'sum' : 
            field_rgd  = field1_rgd + field2_rgd 
         else : 
            field_rgd  = field1_rgd - field2_rgd 

      elif lsm_overlay_plot : 
            field_rgd  = field1_rgd * field2_rgd
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
      label_size_colorbar = 15 
      cb = plt.colorbar(c, orientation='horizontal', extend='both')
      cb.ax.tick_params(labelsize = label_size_colorbar) 
 
# Add latitude and longitude labels and then matching gridlines
      if ll_include_title_lats_lons :
         plt.xlabel("latitude", fontsize=label_size)
         plt.ylabel("longitude", fontsize=label_size)

      xlabels = np.linspace(lon_min, lon_max, 7)
      if central_longitude == 180.0 :
         for ilab in range (7) :   
            xlabels[ilab] = xlabels[ilab] + 180
            if xlabels[ilab] > 180 : xlabels[ilab] = xlabels[ilab] - 360
      ylabels = np.linspace(lat_min, lat_max, 7)
      
      print ('xlabels = ', xlabels) 

      ax.set_xticks(np.linspace(lon_min, lon_max, 7), crs=projection_out)
      ax.set_yticks(np.linspace(lat_min, lat_max, 7), crs=projection_out)
      ax.tick_params(axis='both', labelsize=label_size )

      lat_lon_labelsize = 15 
      ax.set_xticklabels(xlabels,fontsize=lat_lon_labelsize)
      ax.set_yticklabels(ylabels,fontsize=lat_lon_labelsize)
  
      g1 = ax.gridlines(color='k', crs=projection_out)
      g1.xlocator = mticker.FixedLocator( np.linspace(lon_min, lon_max, 7) )
      g1.ylocator = mticker.FixedLocator( np.linspace(lat_min, lat_max, 7) )

# set the filename 
      if field1_3d : 
         st_dep = '{:>6.0f}'.format(depth1)
         st_dep = st_dep.replace(" ","")
         st_lev = '_lev_'+str(k_lev) + '_dep_'+ st_dep     
         print ( 'st_lev = ', st_lev ) 
      else : 
         st_lev = ''

      if difference_plot : 
         if c_sum_diff == 'sum' : 
            outfile=outdir+varname+'_'+infilename+'_sum_'+varname+'_'+infilename2+st_lev
         else : 
            outfile=outdir+varname+'_'+infilename+'_min_'+varname+'_'+infilename2+st_lev
      else : 

         outfile=outdir+varname+'_'+infilename+st_lev  

      title = cseason + ' ' + infilename 
      if difference_plot : 
         if c_sum_diff == 'sum' : 
            title  = title + ' plus ' + infilename2 
         else  : 
            title  = title + ' minus ' + infilename2  
      if field1_3d :  
         title = title + '; ' + st_lev   

      if multiplication_factor != 1.0 : 
         title = title + ' multiplied by ' + str(multiplication_factor)

      if ll_title : 
         plt.title(title, fontsize=8) 
      else :    
         outfile = outfile+'_notitle'

      outfile = outfile+'.png'

# can just hack in the output filename here 
#      outfile = outdir+'nemo_ak306o_hfds_1950_1960_'+simth+'_'+simthp1+'.png'	 
      print ( 'outfile = ', outfile ) 
      plt.savefig(outfile)
      plt.close()   
   
field_plot_differences()
