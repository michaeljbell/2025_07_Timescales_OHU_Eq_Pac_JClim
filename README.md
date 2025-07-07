# README file for git@github.com:michaeljbell/2025\_07\_Timescales\_OHU\_Eq\_Pac\_JClim

## Introduction 

This github repository contains the scripts used to produce the figures in: <br> 

Bell, M. J. and J. A. Baker, 2025 On what Timescales do zonal surface wind stresses drive heat uptake in the equatorial Pacific ocean? J Climatology, accepted. 
 
This README file contains the following sections: 

- Data set inputs. This describes how to access the main data sets used as inputs to the scripts
- Short descriptions of the scripts used to generate each figure 
- Description of the organisation and naming of the scripts and data sets 

## 1) Data set inputs

### a) CMIP6HighResMIP

These datasets are described on https://code.metoffice.gov.uk/trac/ukcmip6/wiki

### b) ERA5 

These datasets were downloaded as monthly mean fields (in netcdf format) from the Copernicus web site: 
 https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means?tab=download

### c) DEEP-C

These data were downloaded from https://researchdata.reading.ac.uk/347/

### d) CMIP6

There are two sources of input files. Both were downloaded from ESGF (Earth System Grid Federation) described in https://doi.org/10.1016/j.future.2013.07.002 using the MetOffice's managecmip tool 

1. 500-year timeseries of piControl. These files contain all months for the first 500 years of the piControl simulation for a specific variable. The (slightly irregular) filenames for these files are constructed by cmip6\_file\_names within utilities\_CMIP6.py.

2. Shorter timeseries for piControl, abrupt-4xCO2, 1pctCO2 and ssp245 integrations. These files are read by field\_build\_annual\_decade\_means\_CMIP6.py

## 2) Short descriptions of the scripts used to generate each figure 

Note that all of these scripts include user input sections. These sections need to be modified by the user to generate specific plots. 
This library does not include details of the specific user choices made to generate all the plots in the paper. All the scripts are in the scripts directory. 

### Figure 1 Time-mean DEEP-C v5 1985 - 2016 NSHF

DEEPC-v5 - can refer to web-site 

- field\_annual\_seasonal\_means\_ERA5\_DEEPC\_EN4.py
- field\_time\_mean\_std.py br>

a) by latitude band

field\_regional\_mean\_plots.py  from 2023\_03\_NSHF\_Zenodo      Figure 3a) of Bell, Storkey and Nurser

b) in equatorial Pacific between 30S and 30N. 

field\_plot\_difference.py

### Figure 2 - conceptual schematic 

### Figure 3 Time-series of correlations between either NSHFs or SSTs and zonal wind stresses for HadGEM3 and DEEP-C

for a) and c) 

- field\_annual\_mean\_timeseries.py to form annual means from monthly means 
- field\_area\_mean\_timeseries\_orig.py to calculate timeseries for the chosen areas
- timeseries\_anom\_mth.py to calculate timeseries of monthly anomalies from the time-mean
- timeseries\_lsq\_fit.py to calculate and plot the least squares fits 

for b) and d):

- field\_annual\_seasonal\_means\_ERA5\_DEEPC\_EN4.py to form annual means from ERA5 and DEEP-C data
- field\_area\_mean\_timeseries\_orig.py to calculate timeseries for the chosen areas
- timeseries\_anom\_mth.py to calculate timeseries of monthly anomalies from the time-mean 
- timeseries\_lsq\_fit.py to calculate and plot the least squares fits

### Figure 4: (a) - (c) Correlations between HadGEM3 100-year NSHF, taux and t20dw time-series as a function of the length of the time-mean period   
 
- field\_area\_mean\_timeseries\_orig.py to calculate monthly-mean timeseries for the chosen areas
- timeseries\_anom\_mth.py to calculate timeseries of monthly anomalies from the time-mean 
- timeseries\_lsq\_fit\_multiple.py to calculate the correlation coefficients as a function of time-mean period and start month 
- lsq\_fit\_prep.sh to prepare the file for input to: 
- lsq\_fit\_plots.py to plot out the correlation coefficients

### Figure 4: (d) - (f) Correlations between DEEP-C NSHF, ERA5 taux and EN4 t20dw as in (a) - (c) 

- field\_annual\_seasonal\_means\_ERA5\_DEEPC\_EN4.py to form monthly means from ERA5 and DEEP-C data inputs
- field\_area\_mean\_timeseries\_orig.py to calculate monthly-mean timeseries for the chosen areas
- timeseries\_anom\_mth.py to calculate timeseries of monthly anomalies from the time-mean 
- timeseries\_lsq\_fit\_multiple.py to calculate the correlation coefficients as a function of-  time-mean period and start month 
- lsq\_fit\_prep.sh to prepare the file for input to: 
- lsq\_fit\_plots.py to plot out the correlation coefficients

### Figure 4: (g) - (i) standard deviations of HadGEM3 100-year timeseries and reanalyses timeseries of NSHF, taux and t20dw as a  function of the length of the time-mean period 

Scripts as for Figure 4 (a) - (f) then
- lsq\_fit\_prep.sh to prepare the file for input to: 
- lsq\_fit\_plots\_std\_LL\_and\_reanal.py

### Figure 5: as for Figure 4 except for SST rather than NSHF. 

Same scripts as Figure 4

### Figure 6: Correlations between pairs of timeseries as a function of time-mean period for CMIP6 models 

- field\_area\_mean\_timeseries\_CMIP6.py  build 500-year long timeseries of area-means for chosen variables and areas for each CMIP6 model
- timeseries\_anom\_mth\_CMIP6.py to calculate timeseries of monthly anomalies from the time-mean
- timeseries\_lsq\_fit\_multiple\_CMIP6.py to calculate the correlation coefficients
- lsq\_fit\_prep\_CMIP6.sh to prepare the files for input to the next script. Uses header\_mth.txt as input copied either from header\_mth\_4\_6.txt or header\_mth\_4\_11.txt.
- lsq\_fit\_plot\_revised\_CMIP6\_with\_reanalysis.py to plot out the panels of figure 6. 

### Figure 7: Geographical plots of time-mean NSHF for each CMIP6 model 

- field\_time\_mean\_std\_CMIP6.py - calculates time-mean of 500 years of monthly data for each CMIP6 model
- field\_plot\_difference\_CMIP6.py - plots out the time-mean field for each model 

### Figure 8: Correlations between taux or tos annual-mean timeseries and monthly fields of other variables for CMIP6HiRes and reanalysis data

- field\_area\_mean\_timeseries\_orig.py to calculate the annual mean timeseries for the selected field and area
- field\_timeseries\_correlation\_v2.py to calculate the fields of correlations 
- field\_plot\_difference.py to plot the fields 
 
### Figure 9: Correlations between taux monthly-mean anomaly timeseries and monthly anomaly fields of tos for CMIP6 models 
 
- field\_area\_mean\_timeseries\_CMIP6.py to calculate the annual mean timeseries for the selected field and area 
- field\_time\_mean\_std\_CMIP6.py to calculate fields of time-means and standard deviations for each model 
- field\_timeseries\_correlation\_CMIP6.py to calculate the fields of correlations for each model
- field\_plot\_difference\_CMIP6.py to plot the fields 
 
### Figure 10: Correlations as in Figure 9 when fields and timeseries are passed through a 5-year running mean filter
 
- field\_area\_mean\_timeseries\_CMIP6.py to calculate the annual mean timeseries for the selected field and area 
- field\_running\_mean\_CMIP6.py to calculate 5-year running mean fields
- field\_time\_mean\_std\_CMIP6.py
- field\_timeseries\_correlation\_CMIP6.py to calculate the fields of correlations for each model
- field\_plot\_difference\_CMIP6.py to plot the fields 

### Figure 11: Correlation and Standard deviations of global-mean tos timeseries and correlation with tos equatorial Pacific timeseries as functions of time-mean period for CMIP6 models
 
- field\_area\_mean\_timeseries\_CMIP6.py to calculate the global or equatorial Pacific mean, annual mean timeseries for tos for each CMIP6 model
- timeseries\_anom\_mth\_CMIP6.py to calculate timeseries of monthly anomalies from the time-mean and remove linear trend 
- timeseries\_lsq\_fit\_multiple\_CMIP6.py to calculate the correlation coefficients
- lsq\_fit\_prep\_CMIP6.sh to prepare the files for input to: 
- lsq\_fit\_plot\_revised\_CMIP6\_with\_reanalysis.py to plot out the panels of figure 6. 
 
### Figure 12: Timeseries of tos, hfds, taux and tos gradient in equatorial Pacific in HadGEM3 CMIP6 piControl, abrupt-4xCO2, 1pctCO2 and ssp245 simulations 
 
- field\_build\_annual\_decade\_means\_CMIP6\_v2.py to generate annual mean fields for HadGEM3 CMIP6 data set 
- field\_area\_mean\_timeseries\_orig.py to generate area means of the annual means
- timeseries\_plot\_pair\_Eq\_Pac\_paper.py to plot out the timeseries
 
### Figure 13: Scatter plots of selected area-mean, time-mean differences from piControl of various fields from CMIP6 integrations
 
- field\_build\_annual\_decade\_means\_CMIP6\_v2.py to generate decadal mean fields for each CMIP6 data set 
- scatter\_plot\_CMIP6\_decadal\_means\_v3.py calculate area-means, time-means and differences from piControl and generate scatter plots for a list of CMIP6 integrations
 
### Figure B1: Timeseries of HadGEM3 HighResMIP correlations as a function of time-mean period for alternative area-mean 
 
 Same scripts as Figure 4 (a) - (c).
 
 ## 3) Description of the organisation and naming of the scripts and data sets 

The scripts assume that the following top-level directories have been set up: 

scripts 	- the main python processing scripts
lsq\\_fit\_plots  - outputs from least square fits (used to generate Figures 3 to 6) 
Surface\_plots   - geographical plots (as in Figures 1 and 7 to 11 
Time\_series 	- time-series of area-means. But also plots of outputs
std\_out   	- output of diagnostics (from some scripts)  

The file naming convention for the scripts is broadly: INPUT\_ACTION\_OUTPUT.py
Example: field\_area\_mean\_timeseries\_orig.py   ; fields are the inputs; area\_meaning is the action; timeseries is the output

There is a utilities.py script that is widely used. It includes functions to: 

generate dates for filenames: date\_string, period\_string, days\_in\_month
create netcdf files: create\_2D\_file, create\_list\_2D\_file, create\_list\_2D\_3D\_file
read fields from netcdf files discarding the first dimension if it is of unit length: field rd\_ncdf\_var\_check\_one

utilities\_CMIP6.py contrains utility scripts that are widely used for processing of the CMIP6 data.

In the scripts, the lowest level functions are defined before the higher level ones. The sections are: 

1) a header describing the name and purpose ; history of important changes 
2) section where user choices are made
3) the "meat" of the script - this is usually split into sections each with a header


The names of the output files (png or timeseries) are intended to be self-describing. Unlike pp or netcdf files there is no header inside the file that contains loads of metadata. The metadata is in the filename itself. This has pros and cons - but is not a bad choice for research code.

