# -*- coding: iso-8859-1 -*-
'''
utilities: contains test_opts, date_string, period_string, days_in_month, create_2D_file, create_list_2D_file, create_list_2D_3D_file, create_3D_file, 
                    check_mask, rd_ncdf_var_check_one, sum_values_tri, sum_in_class, set_ORCA_halos, regrid_field, set_19_distinct_colors 

Date: 31 March 2021
Author: Mike Bell   
Purpose: enable user to set options in a similar manner to a Fortran namelist (this script has not been exploited much)

'''

def test_opts(opts, keys_req, keys_req_typ, keys_opt, keys_opt_typ, keys_opt_def): 

   import sys 
   
# check that user has set dictionary opts correctly 

   l_exit = False

# check the keys the user must set first 
   for key, ktype in zip(keys_req, keys_req_typ)  : 
      if not ( key in opts ) : 
         print ( ' key '+ key + ' is required in opts input ')  
         l_exit = True     
      else : 
         if ktype != type( opts[key] ) : 
            print ( 'key ' + key + ' has the wrong data type ') 
            l_exit = True 

# check or set the optional keys next  	    
   for key, ktype, key_def in zip ( keys_opt, keys_opt_typ, keys_opt_def ) : 
      if not ( key in opts ) : 
         opts[key] = key_def 
      else : 
         if ktype != type( opts[key] ) : 
            print ( 'key ' + key + ' has the wrong data type ' ) 
            l_exit = True 

   if l_exit : 
      sys.exit() 
   
   return opts 

#---------------------------------------------------------------------------------------------- 

def date_string( ll_ann_mean, iyr, imth, ll_out_all = False, iday=0 ) :

# Purpose: Construct date string for CMIPHiRes file names     

# if iday is not None the data are daily means 
# otherwise if ll_ann_mean = True the data are annual means; else they are monthly data

# note that imth is assumed to run from 0 to 11 (0 for Jan; 11 for Dec)  
# iday is the day (between 1 and 31)  
   
  if iday != 0 : 
       siyr = str(iyr) 
       simth = str(imth+1) 
       if imth+1< 10: simth = '0'+simth 
       siday = str(iday) 
       if iday < 10 : siday = '0'+siday 
       date_str = siyr+simth+siday+'-'+siyr+simth+siday   
  else : 

    if ll_ann_mean : 
       iyp1 = iyr + 1
       siyr = str(iyr )
       siyp1 = str(iyp1)
       simth = str(12)
       simthp1 = str(12) 
    else : 
       simth = str(imth+1) 
       simthp1 = str(imth+2)
       siyr   = str(iyr) 
       siyp1 = str(iyr)   # only if imth = 11 is it iyr + 1  
       if imth < 9 : simth   = '0'+str(imth+1)
       if imth < 8 : simthp1 = '0'+str(imth+2)
       if imth == 11 : 
           simthp1 = '01' 
           siyp1   = str(iyr + 1)  

    date_str = siyr+simth+'01-'+siyp1+simthp1+'01'

  if ll_out_all : 
    return date_str, siyr, siyp1, simth, simthp1
  else :   
    return  date_str      

#---------------------------------------------------------------------------------------------- 

def period_string( iyr_st, imth_cal_st, nmths, iyr_end = 'None' ) :

# Purpose: Construct date string similar to CMIPHiRes file names     

# iyr_st       year of start of period   
# imth_cal_st  calendar month (1-12) of start of period 
# nmths        number of months in the period (e.g. 12 for an annual mean) 


# 6 March 2023 extended script so that imth_cal_st may be > 12. This uses fact that input parameters are private to the function
#               

  if imth_cal_st > 12:
     imth_cal_st = imth_cal_st - 12 
     iyr_st = iyr_st+1

  str_iyr_st = str(iyr_st) 

  str_imth_cal_st = str(imth_cal_st) 
  if imth_cal_st < 10 : 
    str_imth_cal_st = '0' + str_imth_cal_st 

  if iyr_end == 'None' : 
     iyr_end = iyr_st
  
  if imth_cal_st+nmths< 13:
       str_iyr_end  = str(iyr_end) 
       imth_end = imth_cal_st+nmths
  else : 
       str_iyr_end  = str(iyr_end+1) 
       imth_end = imth_cal_st+nmths-12

  str_imth_end = str(imth_end) 
  if imth_end < 10 : 
     str_imth_end = '0'+str_imth_end

  period_str = str_iyr_st+str_imth_cal_st+'01-'+str_iyr_end+str_imth_end+'01'

  return  period_str      

#----------------------------------------------------------------------------------------------

def days_in_month(iyr, imth, calendar_type) : 

# iyr is the year
# imth is the month (0 is Jan, 11 is Dec)
# calendar_type:  'Gregorian' or '360'  

  import sys 

  if calendar_type == '360' :
    ll_360 = 'True' 
  elif calendar_type == 'Gregorian' :
    ll_360 = 'Fale'
  else: 
    print ('Stopping beause of user error in determining days_in_month: iyr, imth, calendar_type = ', iyr, imth, calendar_type ) 
    sys.exit()
  
  if ll_360 : 
    ndays = 30 
  else : 
    if imth == 3 or imth == 5 or imth == 8 or imth == 10 :  # April, Jun, Sep, Nov 
       ndays = 30
    elif imth == 1 :  # Feb 
       ndays = 28     # leap year not taken into account 
    else :
       ndays = 31. 
    if imth == 1 and iyr % 4 == 0 : 
       ndays = 29

  return ndays
#----------------------------------------------------------------------------------------------

def create_2D_file ( lat, lon, dep, field, varname, file_out, ll_omit_coords=False ) :

      import netCDF4
      import numpy as np

      npj, npi = np.shape(lat) 
      g = netCDF4.Dataset(file_out,'w')
      time = g.createDimension("t", 1)
      x = g.createDimension("x", npi)
      y = g.createDimension("y", npj)
      z = g.createDimension("z", 1) 
      if not ll_omit_coords :
         depths     = g.createVariable("gdept_1d","f4",("t","z"))
         latitudes  = g.createVariable("nav_lat","f4",("t","y","x"))
         longitudes = g.createVariable("nav_lon","f4",("t","y","x"))
      field_out  = g.createVariable(varname,"f4",("t","y","x"))
      if not ll_omit_coords :
         depths[:] = dep 
         latitudes[:]  = lat
         longitudes[:] = lon
      field_out[:]  = field 
      g.close()
#----------------------------------------------------------------------------------------------
def create_list_2D_file ( xdim, cxdim, ydim, cydim,  field2d_list, varname_list, file_out, depth = None, ll_time = True ) :

      import netCDF4
      import numpy as np

      npi = len(xdim) 
      npj = len(ydim) 

      g = netCDF4.Dataset(file_out,'w')
      if ll_time :
         time = g.createDimension("time", 1)

      x = g.createDimension("x", npi)
      y = g.createDimension("y", npj)
      if depth is not None :
         z = g.createDimension("z", 1)     
      
      if ll_time :
         xdims  = g.createVariable(cxdim,"f4",("time","x"))
         ydims  = g.createVariable(cydim,"f4",("time","y"))
         if depth is not None :
            zdims  = g.createVariable("gdept_1d","f4",("time","z"))
         xdims[:] = xdim
         ydims[:] = ydim 
         if depth is not None :
            zdims[:] = depth 
      else :
         xdims  = g.createVariable(cxdim,"f4",("x"))
         ydims  = g.createVariable(cydim,"f4",("y"))
         if depth is not None :
            zdims  = g.createVariable("gdept_1d","f4",("z"))
         xdims[:] = xdim
         ydims[:] = ydim 
         if depth is not None :
            zdims[:] = depth 

      list_field_out_2D = [] 
      for varname in varname_list :
         if ll_time :
            list_field_out_2D.append ( g.createVariable(varname,"f4",("time","y","x")) )         
         else :
            list_field_out_2D.append ( g.createVariable(varname,"f4",("y","x")) )         

#      print ( ' len (field2d_list) = ', len (field2d_list) ) 

      for i in range ( len (field2d_list) )  :       
         if ll_time :
            list_field_out_2D[i][:]  = field2d_list[i] 
         else :
            list_field_out_2D[i][:]  = field2d_list[i] 

#         print ( 'i, min, max = ', i, np.min(field2d_list[i]), np.max(field2d_list[i]) ) 
	 
      g.close()
#----------------------------------------------------------------------------------------------

def create_list_2D_3D_file ( lat, lon, dep, list_fields_2D, list_varname_2D, 
                                            list_fields_3D, list_varname_3D, 
                                            l_time_included, time, file_out ) :
      import netCDF4
      import numpy as np

      npj,npi = np.shape(lat) 
      npk = len(dep) 
      g = netCDF4.Dataset(file_out,'w')
      t = g.createDimension("time", 1)
      x = g.createDimension("x", npi)
      y = g.createDimension("y", npj)
      z = g.createDimension("z", npk) 
# changed to work with cdfmoc tool 
#       depths     = g.createVariable("gdept_1d","f4",("t","z"))
      depths     = g.createVariable("depth","f4",("time","z"))
      if l_time_included: 
         times = g.createVariable("time_centered","f4",("time"))
      latitudes  = g.createVariable("nav_lat","f4",("time","y","x"))
      longitudes = g.createVariable("nav_lon","f4",("time","y","x"))

      list_field_out_2D = []      
      for varname in list_varname_2D :  
         list_field_out_2D.append ( g.createVariable(varname,"f4",("time","y","x")) ) 

      list_field_out_3D = []      
      for varname in list_varname_3D :  
         list_field_out_3D.append ( g.createVariable(varname,"f4",("time","z","y","x")) ) 

      depths[:] = dep 
      if l_time_included:
         times[:] = time 
      latitudes[:]  = lat
      longitudes[:] = lon

      for i in range ( len (list_fields_2D) )  :       
         list_field_out_2D[i][:]  = list_fields_2D[i] 

      for i in range ( len (list_fields_3D) )  :       
         list_field_out_3D[i][:]  = list_fields_3D[i] 

      g.close()
      
#----------------------------------------------------------------------------------------------

def create_3D_file ( lat, lon, dep, field, varname, file_out ) :

      import netCDF4
      import numpy as np

      npk, npj, npi = np.shape(field) 
      g3d = netCDF4.Dataset(file_out,'w')
      time = g3d.createDimension("t", 1)
      z = g3d.createDimension("z", npk)
      y = g3d.createDimension("y", npj)
      x = g3d.createDimension("x", npi)
      depths = g3d.createVariable("gdept_1d","f4",("t","z"))
      latitudes  = g3d.createVariable("nav_lat","f4",("t","y","x"))
      longitudes = g3d.createVariable("nav_lon","f4",("t","y","x"))
      field_out  = g3d.createVariable(varname,"f4",("t","z","y","x"))
      depths[:] = dep 
      latitudes[:]  = lat
      longitudes[:] = lon
      field_out[:]  = field 
      g3d.close()

#----------------------------------------------------------------------------------------------

def create_2D_timeseries_file ( lat, lon, time, field, varname, file_out ) :

      import netCDF4
      import numpy as np

      npt, npj, npi = np.shape(field) 
      print ('npt, npj, npi = ', npt, npj, npi ) 
      g3d = netCDF4.Dataset(file_out,'w')
      t = g3d.createDimension("t", npt)
      y = g3d.createDimension("y", npj)
      x = g3d.createDimension("x", npi)
      times      = g3d.createVariable("times","f4", ("t"))
      latitudes  = g3d.createVariable("nav_lat","f4",("y","x"))
      longitudes = g3d.createVariable("nav_lon","f4",("y","x"))
      field_out  = g3d.createVariable(varname,"f4",("t","y","x"))
      for j in range (npj) :
        latitudes[j]  = lat[j]
        longitudes[j] = lon[j]
      for jt in range (npt) :
        times[jt] = time[jt]
        field_out[jt]  = field[jt]
      g3d.close()
      return
      
#----------------------------------------------------------------------------------------------
def check_mask ( field, mask, min_mdi, mask_name, file_std ) : 
      import numpy as np  
      field_mask = np.ma.getmask(field)
      inverted_mask = np.logical_not ( field_mask )
      agree = np.equal( inverted_mask, mask )  # all should agree ! 
      all_agree = np.all( agree )

      print ( 'mask_name, all agree = ', mask_name, all_agree, file = file_std ) 
      print ( 'np.shape(field) = ', np.shape(field), file = file_std ) 

      if not all_agree : 
          bad_pts = np.where ( agree, 0, 1)          
          n_bad_pts = np.sum ( bad_pts ) 
          print ( 'n_bad_pts = ', n_bad_pts, file = file_std )
      
      
      return

#----------------------------------------------------------------------------------------------

def rd_ncdf_var_check_one(g, varin) :

# reads a variable from a netcdf file; if the length of the first dimension is 1 it returns the variable with that dimension omitted 

# g        is handle for a file already opened by netCDF4
# varin    is the name of the variable that is to be returned


   import netCDF4
   variable = g.variables[varin][:]   # [:] or [0]
   if len (variable) == 1 :      
      return variable[0]
   else :
      return variable

#------------------------------------------------------------------------------------------------
import numba
@numba.jit(nopython=True, parallel=True)
def sum_values_tri( nvals, nclasses, int_1d, fld_frac, val_1d, sum_values_class) :  

   for i in range (nvals) : 
      sum_values_class[ int_1d[i] ]  = sum_values_class[ int_1d[i] ] + val_1d[i] * ( 1.0 - fld_frac[i] )   
      if int_1d[i] < nclasses - 1 :
         sum_values_class[int_1d[i] + 1] = sum_values_class[int_1d[i] + 1] + val_1d[i] * fld_frac[i]    

   return sum_values_class

#------------------------------------------------------------------------------------------------
def sum_in_class( field, values, class_min, class_width, nclasses ): 

# inputs
# ------
# field                 numpy array of floats: field values that are compared with those defined by the class ; 0 at masked points  
# values                numpy array of floats: the values (volumes or fluxes) to be assigned to the classes. 
#                       Must be 0.0 for field values that are out of range of the classes 
# class_min             central value of the lowest class
# class_width           width of each class (assumed to be uniform) 
# nclasses              number of classes

# Uses triangular algorithm and @numba.jit

# returns 
# -------
# sum_class     1D numpy array of length nclasses containing the sum of the values assigned to each class  

   import numpy as np

   recip_class_width = 1.0 / class_width

   fld_norm = (field - class_min ) * recip_class_width    # normalised values

   val_1d   = np.ndarray.flatten(values)
   f_norm_1d = np.ndarray.flatten(fld_norm)        #   normalised field values in a flattened 1D array
     
   int_1d = f_norm_1d.astype(int)                  # the interval in which the field sits
   int_1d = np.minimum ( int_1d , nclasses - 1 ) 
   int_1d = np.maximum ( int_1d, 0 )  
   
   fld_frac = f_norm_1d - ( class_min + int_1d * class_width)    
   fld_frac = np.minimum ( fld_frac, 1.0 )          
   fld_frac = np.maximum ( fld_frac, 0.0 ) 

   sum_in_class = np.zeros(nclasses)     # 
   nvals = len(val_1d) 
   sum_in_class = sum_values_tri(nvals, nclasses, int_1d, fld_frac, val_1d, sum_in_class )

   return sum_in_class   
#------------------------------------------------------------------------------------------------
def set_ORCA_halos(field, tuvf, isign, pole_type ): 

# inputs
# ------
# field                 2D numpy array of floats on an ORCA grid including halos (halo values might be incorrect on input)  
# tuvf                  string: 't', 'u', 'v' or 'f' depending on the grid on which field is valid 
# isign                 +1. for true scalars (e.g. tracers, divergence); -1 for true vectors (velocities) or pseudo scalars (vorticity) 
#                       to set halos for a tracer field that has been averaged onto the U grid, use tuvf = 'u' and isign = +1.                     
# pole_type             'f' or 't' depending on the model configuration: ORCA2, ORCA025 and ORCA12 use a 't' point pivot; ORCA1 an 'f' type

# returns 
# -------
# field                 2D numpy array of floats on an ORCA grid with the halos set correctly   

# The derivation of this code from the Fortran is described in ../@Notes/lbc_nfd_generic_convert_to_python.txt 


   import sys
   import numpy as np

# check that the inputs for pole_type and tuvf are recognised
   if pole_type != 't' and pole_type != 'f' : 
      print ( ' set_ORCA_halos: pole_type is set incorrectly: it should be \'t\' or \'f\'; it is ', pole_type )   
      sys.exit()

   if tuvf != 't' and tuvf != 'u' and tuvf != 'v' and tuvf != 'f' and tuvf != 'w' : 
      print ( ' set_ORCA_halos: tuvf is set incorrectly: it should be \'t\', \'u\', \'v\' or \'f\'; it is ', tuvf )   
      sys.exit()

# set the cyclic boundary conditions 
   field[:, 0] = field[:,-2]
   field[:,-1] = field[:, 0]

   jpj, jpi = np.shape(field)     # jpi is required below 
   
   if pole_type == 't' :       

      if tuvf == 't' or tuvf == 'w' : 
         field[-1, 1:] = isign * field [ -3, -1:0:-1 ]   
         field[-1, 0]  = isign * field [ -3, 2 ]           
         field[-2, int(jpi/2):] = field[-2, int((jpi+1)/2): 0: -1]   
	 
      elif tuvf == 'u' :      # isign = +1. for vectors like velocity
         field[-1, :-1] = isign * field[-3, :0:-1]  	 
         field[-1, 0 ]  = isign * field[-3, 1]
         field[-1, -1]  = isign * field[-3, -2]
         field[-2, int(jpi/2)-1:-1] = - field[-2, int((jpi+1)/2):0:-1]   
	 
      elif tuvf == 'v' : 
         field[-2, 1:] = isign * field [ -3, -1:0:-1 ]     # rhs corrected from , -1:1:-1] to -1:0:-1] on 31/03/21 
         field[-1, 1:] = isign * field [ -4, -1:0:-1 ]     # same correction as above
         field[-1, 0 ] = isign * field [ -4, 2 ] 

      elif tuvf == 'f' : 
         field[-2, :-1] = isign * field[-3, -1:0:-1]
         field[-1, :-1] = isign * field[-4, -1:0:-1]
         field[-1,0]  = isign * field[-4,1]
         field[-1,-1] = isign * field[-4,-2]

   if pole_type == 'f' :  

      if tuvf == 't' or tuvf == 'w' : 
         field[-1, : ] = isign * field[ -2, ::-1 ]

      elif tuvf == 'u' : 
         field[-1, :-1] = isign * field[-2, -2::-1]  	 
         field[-1,-1]   = isign *  field[-2,-3]	 

      elif tuvf == 'v' : 
         field[-1, :] = isign * field[-3, ::-1] 
         field[-2, int(jpi/2):] = isign * field[-2, int((jpi-1)/2)::-1] 

      elif tuvf == 'f' : 
         field[-1, :-1] = isign * field[-3, -2::-1]
         field[-1,-1] = isign * field[-3, -3]
         field[-2, int(jpi/2):-1] = isign * field[-2, int((jpi-1)/2)-1::-1] 

   return field 

#------------------------------------------------------------------------------------------------
def regrid_field(ax, fig, lon, lat, field, projection_in, lon_min, lon_max, lat_min, lat_max, projection_out ): 

  import cartopy
  import matplotlib.pyplot as plt

# Determine number of pixels in the subplot
  bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
  nx = bbox.width * fig.get_dpi()
  ny = bbox.height * fig.get_dpi()
  nx = round( nx)     # added 2 March 2022
  ny = round (ny)   
  print ( 'nx, ny = ', nx, ny ) 

# Reproject the data onto a regular grid (with dimensions set by the number of pixels in the subplot, as above)

  x, y = cartopy.img_transform.mesh_projection(projection_in, nx, ny, x_extents=[lon_min, lon_max], y_extents=[lat_min, lat_max] )[:2]

  field_rgd = cartopy.img_transform.regrid(field, lon, lat, projection_in, projection_out, x, y)

  return x, y, field_rgd

#------------------------------------------------------------------------------------------------

def set_19_distinct_colors(white_range=0):


# no inputs
# ---------

# returns 
# -------
# newcmp                 a colour map with 19 colours intended to be distinct with light colours near the middle and white in the middle, 
#                        transitioning from blue to green through white to red then purple   

    import numpy as np
    from matplotlib import cm
    from matplotlib.colors import ListedColormap, to_rgba
    import copy 

    tab20b = cm.get_cmap('tab20b', 19)
    newcolors = tab20b(np.linspace(0,1,19))
    white = to_rgba('w')
    newcolors[9,:]=white
# move 8th colour to 4th and shift colours 4 to 7 up into 5 to 8  
    temp_colour = copy.deepcopy( newcolors[8,:] ) 
    newcolors[8,:] = newcolors[7,:]
    newcolors[7,:] = newcolors[6,:]
    newcolors[6,:] = newcolors[5,:]
    newcolors[5,:] = newcolors[4,:]
    newcolors[4,:] = temp_colour

# move 
    new_scheme = True
    if new_scheme :
       temp_colour = copy.deepcopy( newcolors[0,:] ) 
       newcolors[0,:] = newcolors[3,:]
       newcolors[3,:] = temp_colour
       temp_colour = copy.deepcopy( newcolors[2,:] ) 
       newcolors[2,:] = newcolors[1,:]
       newcolors[1,:] = temp_colour
       temp_colour = copy.deepcopy( newcolors[11,:] ) 
       newcolors[11,:] = newcolors[14,:]
       newcolors[14,:] = temp_colour
       temp_colour = copy.deepcopy( newcolors[12,:] ) 
       newcolors[12,:] = newcolors[13,:]
       newcolors[13,:] = temp_colour


    for icolor in range ( 1, white_range+1 ) : 
       newcolors[9+icolor,:]=white 
       newcolors[9-icolor,:]=white 
    
    newcmp=ListedColormap(newcolors) 
   
    return newcmp  

#------------------------------------------------------------------------------------------------
