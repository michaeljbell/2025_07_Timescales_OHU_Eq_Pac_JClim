# -*- coding: iso-8859-1 -*-
'''
field_build_decad_means_CMIP6_v2.py 

contains cmip6_champ_file_name, field_build_decad_means_CMIP6 

Date: 15 Feb 2025
Author: Mike Bell   

10 Apr 2025 extended to allow annual means with start_month different from January 

'''
   
def cmip6_champ_file_name(c_model, c_expt_type, c_var, start_month, start_year, end_year, file_std, c_indices='missing'): 

# Purpose: convert CMIP6 datasets into decadal means [or shorter period means] with a common time-axis 
#          1850 in historical and ssp files is mapped to year 0 on the common time-axis; 
#          The start of 'piControl', '1pctCO2' and 'abrupt-4xCO2' integrations are also mapped to year 0     

# Inputs 
# c_model       string: the name of the model data to be read: e.g.  'HadGEM3-GC31-LL'
# c_expt_type   string: the eperiment type. One of: 'abrupt-4xCO2', '1pctCO2', 'piControl', 'historical', 'ssp245', ssp585'
# c_var         string: the short name of the variable: e.g. 'tos', 'tauuo', 'hfds'

# start_month   integer: the start month for the period (usually 1 but for annual means often use 6) 
# start_year    integer: the start year in the common time-axis, i.e. relative to the start of the experiment (start of first year = 0 )    
# end_year      integer: the end year     "              "                "                           "

#Optional inputs
# c_indices     string: when provided overrides the "default" settings specified in this function (to enable future loops over c_indices)  

# Outputs: (matching lists of the same length)
# list_filenames            list of strings: the names of the files in the champ directory to be read 
# list_start_month_offsets  list of integers: the number of months after the start of the file of the first month in the file to be read
# list_months_to_read       list of integers: the number of months in the file to be read


# Errors detected
#  c_model, c_expt_type, c_var not in the list of variables expected
# the data to the end of end_year is not within the file returned 

# Mike Bell 6 December 2024

   import sys 
   from utilities_CMIP6 import cmip6_convert_int_to_str4

# 0. Initialise outputs
   list_filenames = []
   list_start_month_offsets = []
   list_nmonths_to_read=[]
   
# 1. Set default values solely from c_model    
#    period = -1 means all the data is in one dataset; 
#                 institute,        version,    period(yrs), DateStartPIC, yearsPIC   FileStartOther  YearsOther  
   if c_model == 'ACCESS-CM2':
     itemlist = ['CSIRO-ARCCSS',  'gn', 'v20191109',   -1,      '095001',     500,      '095001',      150 ] 
   elif c_model == 'CanESM5' : 
     itemlist = ['CCCma',         'gn', 'v20190429',  200,      '520101',     800,      '185001',      151 ] 
   elif c_model == 'GFDL-CM4' : 
     itemlist = ['NOAA-GFDL',     'gn', 'v20180701',   20,      '015101',     500,      '185001',      160 ] 
   elif c_model == 'HadGEM3-GC31-LL' :  
     itemlist = ['MOHC',          'gn', 'v20190620',  100,      '185001',     500,      '185001',      150 ]
   elif c_model == 'HadGEM3-GC31-MM' :  
     itemlist = ['MOHC',          'gn', 'v20200302',   20,      '185001',     500,      '185001',      150 ]
   elif c_model == 'IPSL-CM6A-LR' : 
     itemlist = [ 'IPSL',         'gn', 'v20200326',   -1,      '185001',     500,      '185001',      300 ] 
   elif c_model == 'MIROC6' : 
     itemlist = ['MIROC',         'gn', 'v20210129',  100,      '320001',     800,      '185001',      150 ]
   elif c_model == 'MPI-ESM1-2-LR' : 
     itemlist = ['MPI-M',         'gn', 'v20190710',   20,      '185001',    1000,      '185001',      140 ]  
   elif c_model == 'NorESM2-LM' :  
     itemlist = ['NCC',           'gn', 'v20230616',   10,      '160001',     300,      '000101',      140 ]   # startPIC corrected to 160001 from 170001 18 Feb 2025
   else : 
     print ( ' cmip6_file_names: stopping because of unrecognised c_model = ', c_model ) 
     sys.exit()

   c_institute        = itemlist[0]
   c_grid             = itemlist[1]
   c_version          = itemlist[2]
   PeriodYears        = itemlist[3]
   c_date_start_PIC   = itemlist[4]
   YearsPIC           = itemlist[5]
   c_file_start_Other = itemlist[6]
   YearsOther         = itemlist[7]

   print (' c_model = ', c_model, file = file_std)
   print (' itemlist = ', itemlist, file = file_std)

# 2. Set default values from  c_expt_type and c_var
#

   if c_expt_type == 'piControl' : 
      c_date_start     = c_date_start_PIC
      c_file_start     = c_date_start
      YearsInFile      = YearsPIC
   elif c_expt_type in ['1pctCO2', 'abrupt-4xCO2' ] :  
      c_file_start     = c_file_start_Other
      if c_model == 'MIROC6' : c_file_start = c_date_start_PIC
      c_date_start     = c_file_start
      YearsInFile      = YearsOther
   elif c_expt_type =='historical' : 
      c_file_start     = '185001'
      c_date_start     = '185001'
      YearsInFile      = 165
   elif c_expt_type in ['ssp245', 'ssp585' ] : 
      c_file_start     = '201501'
      c_date_start     = '185001'
      YearsInFile      = 2101 - 1850
      YearsInFile2     = 2101 - 1850

   else : 
     print ( ' cmip6_champ_file_name: stopping because of unrecognised c_expt_type = ', c_expt_type, file = file_std ) 
     sys.exit()

   if c_var in [ 'hfds', 'tauuo', 'tos', 'areacello' ] :
     print( 'c_var = ', c_var, file = file_std ) 
   else : 
     print ( ' cmip6_champ_file_name: stopping because of unrecognised c_var = ', c_var, file = file_std ) 
     sys.exit()

# 3. Revise default values to take account of dataset specific variations in the file naming  

   if c_model == 'ACCESS-CM2' : 
      if c_expt_type in ['abrupt-4xCO2', 'historical', 'ssp245' ]  : c_version = 'v20191108' 
      if c_expt_type == 'piControl': c_version = 'v20191112'   
      if c_expt_type == 'ssp245' : c_version = 'v20191108'
      if c_expt_type == 'ssp585' : c_version = 'v20210317'

   elif c_model == 'CanESM5' : 
      if c_expt_type in [ 'ssp245', 'ssp585' ]  : 
           if  start_year < 2101 : 
              PeriodYears = -1  
              YearsInFile = 2101 - 1850 
              c_file_start2 = '210101'
           else : 
              c_file_start = '210101'

   elif c_model == 'GFDL-CM4' : 
      if c_expt_type in [ 'abrupt-4xCO2', '1pctCO2' ] : 
         c_date_start = '000101'
         c_file_start = '000101'
      if c_expt_type == '1pctCO2' : YearsInFile = 150
      if c_expt_type == 'piControl'    : 
         if c_var == 'tos' : c_version = 'v20190201'     # this is a bit weird ?? perhaps other downloads failed? check the list_full

   elif c_model in [ 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM'] : 
      if  c_indices == 'missing' : 
        if c_expt_type == 'piControl' :    
           c_indices = 'r1i1p1f1'
        else :
           c_indices = 'r1i1p1f3'

      if c_model == 'HadGEM3-GC31-LL' : 
        if c_expt_type == 'piControl'  : c_version = 'v20211103'
        if c_expt_type == 'historical' : c_version = 'v20190624'
        if c_expt_type == 'ssp245'     : c_version = 'v20190908'
        if c_expt_type == 'ssp585'     : c_version = 'v20200114'
        if c_expt_type in [ 'ssp245', 'ssp585' ] :
           PeriodYears = -1  
           if  start_year < 200 : 
              YearsInFile = 2050 - 1850 
              c_file_start2 = '205001'
              YearsInFile2 = 2101 - 1850
           else : 
              c_file_start = '205001'
              YearsInFile = 2101 - 1850
        if c_var == 'areacello'        :  c_version = 'v20190709'
      else : 
        if c_expt_type == 'piControl'    : c_version = 'v20191204'
        if c_expt_type == 'historical'   : c_version = 'v20191207'
        if c_expt_type == 'abrupt-4xCO2' : c_version = 'v20200302'
        if c_expt_type == '1pctCO2'      : c_version = 'v20200228'
        if c_expt_type == 'ssp585'       : c_version = 'v20200515'
        if c_expt_type in [ 'ssp245', 'ssp585' ] :
           if  start_year < 180 : 
              PeriodYears = -1 
              YearsInFile = 2030 - 1850
           else : 
              c_file_start = '203001'
              YearsInFile = 2101 - 1850
        if c_var == 'areacello'          : c_version = 'v20200108'

   elif c_model == 'IPSL-CM6A-LR' : 
      if c_expt_type == 'piControl':   
         if c_var == 'tos' : c_version = 'v20200326'
         if c_var == 'hfds': c_version = 'v20181022'
         if c_var == 'tauuo': c_version = 'v20181123'
      elif c_expt_type == 'historical':  
         c_version = 'v20180803'
         YearsInFile = 165
      elif c_expt_type == '1pctCO2' : 
         c_version = 'v20180727' 
         YearsInFile = 150
      elif c_expt_type == 'abrupt-4xCO2' : 
         c_version = 'v20190118'
      elif c_var=='tauuo' : 
         if c_expt_type =='abrupt-4xCO2' : c_version = 'v20190118'
         if c_expt_type =='piControl'    : c_version = 'v20181123'
      if c_expt_type == 'ssp245'     : c_version = 'v20190119'
      if c_expt_type == 'ssp585'     : c_version = 'v20190903'

   elif c_model == 'MIROC6' :
      if c_var == 'areacello'      : c_version = 'v20190311'    
      if c_expt_type in [ 'ssp245', 'ssp585' ] :
         PeriodYears = -1   
         if c_var == 'tos': 
            c_version = 'v20190627'
         else :
            c_version = 'v20210129'
      if c_expt_type == 'piControl' and c_var == 'tos' :
            c_version = 'v20181212'
      if c_expt_type == 'abrupt-4xCO2' and c_var == 'tos' :
            c_version = 'v20190705'
      if c_expt_type =='historical' and c_var == 'tos':
            c_version = 'v20181212'

   elif c_model == 'MPI-ESM1-2-LR' : 
      if c_var == 'tauuo': c_version = 'v20200909'
      if c_expt_type =='historical' :  YearsInFile = 165
      if c_expt_type in [ 'ssp245', 'ssp585' ] :
         c_file_start = '201501'
         if c_var == 'tauo':  c_version = 'v20200909' 
   
   elif c_model == 'NorESM2-LM' : 
      if   c_expt_type =='abrupt-4xCO2' :  c_version = 'v20210118'
      elif c_expt_type in [ 'ssp245', 'ssp585' ] :       
         c_file_start = '202101'                     # not easy to build the file for 2020-2029 
         c_version = 'v20191108'
      else : 
         if c_indices == 'missing'  : c_indices = 'r1i1p4f1'

      if c_expt_type =='historical' : 
         c_file_start   = '185001'
         c_date_start   = '185001'
         YearsInFile = 165
	 
   print (' start_year, end_year, c_date_start, c_file_start, PeriodYears, YearsInFile = ', start_year, end_year, c_date_start, c_file_start, PeriodYears, YearsInFile, file = file_std)

# 4. Set c_indices if not set previously (this must be done last!)
   if c_indices == 'missing' : c_indices = 'r1i1p1f1'

# 5. Construct the list of dates in the filenames, the list_start_month_offsets and the lists_months_to_read

#  find list_CDateOut             list of strings (len4) each is the date-string in the filename of the file to be read 
#       list_start_month_offsets  list of integers: the number of months after the start of the file of the first month in the file to be read
#       list_months_to_read       list of integers: the number of months in each file to be read

   list_CDateOut = []    # list_start_month_offsets and list_months_to_read were initialised at the top of the function (as they are function outputs) 

   if c_var in [ 'areacello' ] : 
      list_start_month_offsets.append( 0 ) 
      list_nmonths_to_read.append( 1 )   
      list_CDateOut.append( '' )
      Otype = 'Ofx'   
   else : 
   
      Otype = 'Omon' 
      if PeriodYears == -1 :    # all the data are in one file 
        c_file_start_out = c_file_start[:4]

        IyrFirstFile = int(c_file_start_out) 
        IYrSeriesStart = int(c_date_start[:4])

        list_start_month_offsets.append( (start_year + IYrSeriesStart - IyrFirstFile) * 12 + start_month - 1 )

        IyrEnd = IYrSeriesStart + YearsInFile - 1 
        CyrEnd  = cmip6_convert_int_to_str4(IyrEnd)
        list_CDateOut.append( '_'+ c_file_start_out +'01-'+CyrEnd+'12' )

        if start_month != 1: 
           end_year_out = end_year + 1
        else : 
           end_year_out = end_year

        if IyrEnd < IYrSeriesStart + end_year_out :   # in this case we need to read a second file
           IyrStart = IYrSeriesStart + start_year      
           list_nmonths_to_read.append( (IyrEnd - IyrStart + 1)*12 - (start_month - 1) ) # this reads to the end of the file
        else :                                        # in this case we only need to read one file
           list_nmonths_to_read.append( (end_year - start_year + 1)*12 ) 

        if IyrEnd < IYrSeriesStart + end_year_out :   # in this case we need to read a second file 
        #   (we assume it is enough to read two files; that requires end_year - start_year < PeriodYears )
        # for this case we need c_file_start2 and YearsInFile2 to be set ; this has only been done for HadGEM3-LL

          c_file_start_out2 = c_file_start2[:4]
          IyrEnd2 = IYrSeriesStart + YearsInFile2 - 1 
          CyrEnd2  = cmip6_convert_int_to_str4(IyrEnd2)
          list_CDateOut.append( '_'+ c_file_start_out2 +'01-'+CyrEnd2+'12' )

          list_start_month_offsets.append(0)

          IyrStartFile2 = int(c_file_start_out2) 
          IyrEndRead = IYrSeriesStart + end_year
          n_months_to_read = (IyrEndRead - IyrStartFile2 + 1)*12 + (start_month - 1)
          list_nmonths_to_read.append( n_months_to_read )

      else : # the data are in chunks. 
        IyrFirstFile = int(c_file_start[:4])
        IYrSeriesStart = int(c_date_start[:4])
        n_periods_offset = int ( (start_year + IYrSeriesStart - IyrFirstFile) / PeriodYears )        # integer arithmetic ! 
        print ( ' start_year, IYrSeriesStart, IyrFirstFile, n_periods_offset = ', start_year, IYrSeriesStart, IyrFirstFile, n_periods_offset, file = file_std )

        IyrStartFile = IyrFirstFile + n_periods_offset*PeriodYears  # year of the start file we will read 
        c_file_start_out = cmip6_convert_int_to_str4(IyrStartFile)
        n_start_month_offset =  (start_year + IYrSeriesStart - IyrStartFile) * 12 +  start_month - 1
        if n_start_month_offset < 0 : 
           print (' n_start_month_offset is negative. Something has gone wrong so stopping. n_start_month_offset = ', n_start_month_offset ) 
           print (' n_start_month_offset is negative. Something has gone wrong so stopping. n_start_month_offset = ', n_start_month_offset, file = file_std ) 
           sys.exit()
        list_start_month_offsets.append( n_start_month_offset )

        IyrEnd =  IyrStartFile + PeriodYears - 1 

        if IyrEnd > IYrSeriesStart + YearsInFile - 1 :     # the last file sometimes ends before the end of the period used for the other files 
           IyrEnd = IYrSeriesStart + YearsInFile - 1
        CyrEnd  = cmip6_convert_int_to_str4(IyrEnd)
        list_CDateOut.append( '_'+ c_file_start_out +'01-'+CyrEnd+'12' )

        if start_month != 1: 
           end_year_out = end_year + 1
        else : 
           end_year_out = end_year

        IyrStart = IYrSeriesStart + start_year      
        if IyrEnd < IYrSeriesStart + end_year_out :   # in this case we need to read a second file
           list_nmonths_to_read.append( (IyrEnd - IyrStart + 1)*12 - (start_month - 1) ) # this reads to the end of the file
        else :                                        # in this case we only need to read one file
           list_nmonths_to_read.append( (end_year - start_year + 1)*12 ) 
 
        if IyrEnd < IYrSeriesStart + end_year_out :   # in this case we need to read a second file (we assume it is enough to read two files; that requires end_year - start_year < PeriodYears )
          n_periods_offset = n_periods_offset + 1
          IyrStartFile = IyrStartFile + PeriodYears
          c_file_start_out = cmip6_convert_int_to_str4(IyrStartFile)
          list_start_month_offsets.append( 0 )    # the first field to be read is the first field in the file

          IyrEndRead = IYrSeriesStart + end_year
          n_months_to_read = (IyrEndRead - IyrStartFile + 1)*12 + (start_month - 1)
          if n_months_to_read > PeriodYears*12 : 
             print (' n_months_to_read is greater than the number of months in the file. Something has gone wrong so stopping. n_months_to_read = ', n_months_to_read ) 
             print (' n_months_to_read is greater than the number of months in the file. Something has gone wrong so stopping. n_months_to_read = ', n_months_to_read, file = file_std ) 
             sys.exit() 
          list_nmonths_to_read.append( n_months_to_read ) 

          IyrEnd = IyrStartFile + PeriodYears - 1 
          if IyrEnd > IYrSeriesStart + YearsInFile - 1 : 
             IyrEnd = IYrSeriesStart + YearsInFile - 1
          CyrEnd = cmip6_convert_int_to_str4(IyrEnd)
          list_CDateOut.append( '_'+ c_file_start_out +'01-'+CyrEnd+'12' )

        if end_year_out > YearsInFile : 
           print ( ' cmip6_champ_file_name: stopping because you would read beyond the end of the file: end_year_out, YearsInFile = ', end_year_out, YearsInFile, file = file_std ) 
           sys.exit()
      
      print (  'list_start_month_offsets, list_nmonths_to_read, list_CDateOut = ', list_start_month_offsets, list_nmonths_to_read, list_CDateOut, file = file_std ) 
   
# 6. Construct the filename
   if c_expt_type in [ 'ssp245', 'ssp585' ]  : 
      cfolder = 'ScenarioMIP'
   else :
      cfolder = 'CMIP'
   
   for CDateOut in list_CDateOut :
      c_dir = '/data/users/managecmip/champ/CMIP6/'+cfolder+'/'+c_institute+'/'+c_model+'/'+c_expt_type+'/'+c_indices+'/'+Otype+'/'+c_var+'/'+c_grid+'/'+c_version
      c_name=  c_var+'_'+Otype+'_'+c_model+'_'+c_expt_type+'_'+c_indices+'_'+c_grid+CDateOut+'.nc'
      list_filenames.append( c_dir + '/'+ c_name )
     
   return list_filenames, list_start_month_offsets, list_nmonths_to_read  
   
#------------------------------------------------------------------------------------------

def field_build_decad_means_CMIP6() : 


   import sys 
   import netCDF4
   import numpy as np
   from utilities import create_2D_file, rd_ncdf_var_check_one
   from utilities_CMIP6 import cmip6_convert_int_to_str4, cmip6_int_to_str_mth

# User choices start 

   datadir = '/data/users/mike.bell/'
   dirout = datadir + 'CMIP6/'

   nyears_per_output =  1      #     either 1 to output annual means or 10 to output decadal means 
   start_month       =  6      #     integer: start month of output periods ; usually 1 (for decadal means) or 6 for annual means focused on eq Pacific

   for c_expt_type in [ 'abrupt-4xCO2', 'piControl' ] : # 'abrupt-4xCO2', '1pctCO2', 'ssp245', 'ssp585', 'historical' ] : # [ 'abrupt-4xCO2', '1pctCO2', 'ssp245', 'ssp585', 'piControl', 'historical', 'areacello' ] :

      filename_std = '../std_out/build_CMIP6_decade_means/'+c_expt_type+'.txt'
      file_std = open(filename_std, 'w') 
      print ( 'writing debug output to filename_std = ', filename_std ) 

      for c_model in [ 'ACCESS-CM2', 'CanESM5', 'GFDL-CM4', 'HadGEM3-GC31-LL', 'HadGEM3-GC31-MM', 'IPSL-CM6A-LR', 'MIROC6', 'MPI-ESM1-2-LR', 'NorESM2-LM'   ] : 
         for c_var in [ 'tos'] : # [ 'areacello'] :

# User choices end 

           ll_calculate = True
	 
           if c_expt_type in ['ssp245', 'ssp585' ] :
             nyear_start = 165     #     Note that calculations for NorESM are incorrect for nyear_start = 170  ;  was 180 
             nyear_end   = 180     #     user choice for length of time-series to use as input                  ;  was 250
             if c_model ==  'HadGEM3-GC31-MM' and c_expt_type == 'ssp245' :              ll_calculate = False 
             if c_model ==  'MIROC6' and c_expt_type == 'ssp585' :                       ll_calculate = False 
             if c_model ==  'ACCESS-CM2' and c_var == 'tauuo' :                          ll_calculate = False
             if c_model ==  'CanESM5' and c_var == 'tauuo' and c_expt_type == 'ssp585' : ll_calculate = False

           if c_expt_type == 'piControl' : 
             nyear_start = 150     # normally 0 !
             if c_model ==  'NorESM2-LM' : #   and nyear_start < 10 : 
                print ( 'not generating data for NorESM2-LM for nyear_start < 10 = ' )      
                nyear_start = 10  # first file only covers 9 years; the time series is actually offset ...   
             nyear_end   = 160 #   250 is available

           if c_expt_type in [ 'abrupt-4xCO2' ] :
             nyear_start =   0
             nyear_end   = 120
             if c_model ==  'MIROC6' and c_var != 'tos' :  ll_calculate = False

           if c_expt_type in [ '1pctCO2' ] :
             nyear_start = 0
             nyear_end   = 140 #  for temporary testing   140
             if c_model ==  'MIROC6'                    :  ll_calculate = False

           if c_expt_type == 'historical' :
             nyear_start = 150     #  0 or 30     #     Note that calculations for NorESM are incorrect for nyear_start = 170 
             nyear_end   = 164     #     user choice for length of time-series to use as input  (was 160) 
             if  c_model ==  'HadGEM3-GC31-MM' : nyear_end = 120     #   periods 120-130 and 130-140 are difficult because the files are irregular 
	        
           if ll_calculate == False : 
              print ( 'not generating data for c_model, c_var, c_expt_type = ', c_model, c_var, c_expt_type )
              print ( 'not generating data for c_model, c_var, c_expt_type = ', c_model, c_var, c_expt_type, file = file_std )

           if c_var == 'areacello' :  
               nyears_per_output = nyear_end    # so the areacello is only build once
               if c_expt_type != 'piControl' : 
                  print ( 'not generating areacello data for c_expt_type = ', c_expt_type ) 
                  ll_calculate = False  
	 
           if ll_calculate : 
             for start_year in range(nyear_start, nyear_end, nyears_per_output) : 
	       
               end_year = start_year + nyears_per_output - 1    	    # does not include addition of 1 for start_month != 1 
                  
               list_filenames, list_start_month_offsets, list_nmonths_to_read = cmip6_champ_file_name(c_model, c_expt_type, c_var, start_month, start_year, end_year, file_std)	       
               print ( ' list_filenames, list_start_month_offsets, list_nmonths_to_read = ', list_filenames, list_start_month_offsets, list_nmonths_to_read ) 
               print ( ' list_filenames, list_start_month_offsets, list_nmonths_to_read = ', list_filenames, list_start_month_offsets, list_nmonths_to_read, file = file_std ) 

               nmonths_total = 0 
               ll_first_file = True  #   included on 02 March 2025; decadal means built from two files were previously wrong 

               for (c_filename, start_month_offset, nmonths_to_read) in zip( list_filenames, list_start_month_offsets, list_nmonths_to_read ) :  
                 nmonths_total = nmonths_total + nmonths_to_read 
                 g_file = netCDF4.Dataset(c_filename,'r')         

                 if c_model == 'GFDL-CM4' : 
                   clat = 'lat'      ;  clon = 'lon'  
                 elif c_model == 'IPSL-CM6A-LR' : 
                   clat = 'nav_lat'      ;  clon = 'nav_lon'  
                 else : 
                   clat = 'latitude' ;  clon = 'longitude'  
   
                 nav_lat  = rd_ncdf_var_check_one(g_file, clat)
                 nav_lon  = rd_ncdf_var_check_one(g_file, clon)

                 print ( ' nmonths_to_read = ', nmonths_to_read, file = file_std ) 
                 for imth in range ( nmonths_to_read ) :    #  typically a decade of months
	       
                   if c_var == 'areacello' :
                     field = g_file.variables[c_var][:]
                   else : 
                     field = g_file.variables[c_var][imth+start_month_offset]
                   if ll_first_file : 
                     print ( ' c_model, c_expt_type, c_var, np.shape(field) = ', c_model, c_expt_type, c_var, np.shape(field), file = file_std )  
                     field_mean = field	       
                     ll_first_file = False 
                   else : 
                     field_mean = field_mean + field

               recip = 1.0 /  float( nmonths_total )
               print ('recip = ', recip, file = file_std ) 
               field_mean = recip * field_mean 
	       
               if c_var == 'areacello' :
                  TimeStringOut= ''
               else : 
                  cStYr  = cmip6_convert_int_to_str4(start_year)

                  if start_month != 1 : 
                     end_year_out = end_year+1
                  else : 
                     end_year_out = end_year

# strings for the start and end months in the output file
                  cEndYr = cmip6_convert_int_to_str4(end_year_out)
                  c_mth_start = cmip6_int_to_str_mth( start_month)
                  end_month = (start_month + 10) % 12 + 1       # start_month =1 has end_month = 12 
                  c_mth_end   = cmip6_int_to_str_mth( end_month)         
                  TimeStringOut = '_' + cStYr + c_mth_start +'-'+cEndYr + c_mth_end
                  if nyears_per_output == 10 : 
                     dir_period = '/dec/'
                  elif nyears_per_output == 1 : 
                     dir_period = '/1y_'+c_mth_start+'/'
               dep = 0.0 
	       
               file_out = dirout+c_model+dir_period+ c_var+'_'+ c_model +'_' + c_expt_type + TimeStringOut + '.nc'  
               print ( ' writing file_out = ', file_out, file = file_std ) 
               create_2D_file ( nav_lat, nav_lon, dep, field_mean, c_var, file_out ) 
     
      file_std.close()  

field_build_decad_means_CMIP6()
