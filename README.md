# README file for git@github.com:michaeljbell/2025_07_Timescales_OHU_Eq_Pac_JClim

## Introduction 

This github repository contains the scripts used to produce the figures in: <br> 

Bell, M. J. and J. A. Baker, 2025 On what Timescales do zonal surface wind stresses drive heat uptake in the equatorial Pacific ocean?
 
This README file contains the following sections: 

- Data set inputs. This describes how to access the main data sets used as inputs to the scripts
- Short descriptions of the scripts used to generate each figure 
- The organisation of the scripts and data sets 

## 1) Data set inputs

### a) CMIP6HighResMIP - 

### b) ERA5 

### c) DEEP-C

### d) CMIP6

There are two sources of input files. Both were downloaded from the CMIP6 facility: 

1. 500-year timeseries of piControl. These files contain all months for the first 500 years of the piControl simulation for a specific variable. The (slightly irregular) filenames for these files are constructed by cmip6_file_names within utilities_CMIP6.py.

2. Shorter timeseries for piControl, abrupt-4xCO2, 1pctCO2 and ssp245 integrations. These files are read by field_build_annual_decade_means_CMIP6.py

## 2) Utility scripts 

utilities.py and utilities_CMIP6.py are utility scripts that are widely used by the other scripts 

## 2) Short descriptions of the scripts used to generate each figure 

Note that all of these scripts include user input sections. These sections need to be modified by the user to generate specific plots. 
This library does not include details of the specific user choices made to generate all the plots in the paper.

### Figure 1 Time-mean DEEP-C v5 1985 - 2016 NSHF

DEEPC-v5 - can refer to web-site 

- field_annual_seasonal_means_ERA5_DEEPC_EN4.py
- field_time_mean_std.py br>

a) by latitude band

field_regional_mean_plots.py  from 2023_03_NSHF_Zenodo      Figure 3a) of Bell, Storkey and Nurser

b) in equatorial Pacific between 30S and 30N. 

field_plot_difference.py

### Figure 2 - conceptual schematic 

### Figure 3 Time-series of correlations between either NSHFs or SSTs and zonal wind stresses for HadGEM3 and DEEP-C

for a) and c) 

- field_annual_mean_timeseries.py to form annual means from monthly means 
- field_area_mean_timeseries_orig.py to calculate timeseries for the chosen areas
- timeseries_anom_mth.py to calculate timeseries of monthly anomalies from the time-mean
- timeseries_lsq_fit.py to calculate and plot the least squares fits 

for b) and d):

- field_annual_seasonal_means_ERA5_DEEPC_EN4.py to form annual means from ERA5 and DEEP-C data
- field_area_mean_timeseries_orig.py to calculate timeseries for the chosen areas
- timeseries_anom_mth.py to calculate timeseries of monthly anomalies from the time-mean 
- timeseries_lsq_fit.py to calculate and plot the least squares fits

### Figure 4: (a) - (c) Correlations between HadGEM3 100-year NSHF, taux and t20dw time-series as a function of the length of the time-mean period   
 
- field_area_mean_timeseries_orig.py to calculate monthly-mean timeseries for the chosen areas
- timeseries_anom_mth.py to calculate timeseries of monthly anomalies from the time-mean 
- timeseries_lsq_fit_multiple.py to calculate the correlation coefficients as a function of time-mean period and start month 
- lsq_fit_prep.sh to prepare the file for input to: 
- lsq_fit_plots.py to plot out the correlation coefficients

### Figure 4: (d) - (f) Correlations between DEEP-C NSHF, ERA5 taux and EN4 t20dw as in (a) - (c) 

- field_annual_seasonal_means_ERA5_DEEPC_EN4.py to form monthly means from ERA5 and DEEP-C data inputs
- field_area_mean_timeseries_orig.py to calculate monthly-mean timeseries for the chosen areas
- timeseries_anom_mth.py to calculate timeseries of monthly anomalies from the time-mean 
- timeseries_lsq_fit_multiple.py to calculate the correlation coefficients as a function of-  time-mean period and start month 
- lsq_fit_prep.sh to prepare the file for input to: 
- lsq_fit_plots.py to plot out the correlation coefficients

### Figure 4: (g) - (i) standard deviations of HadGEM3 100-year timeseries and reanalyses timeseries of NSHF, taux and t20dw as a  function of the length of the time-mean period 

Scripts as for Figure 4 (a) - (f) then
- lsq_fit_prep.sh to prepare the file for input to: 
- lsq_fit_plots_std_LL_and_reanal.py

### Figure 5: as for Figure 4 except for SST rather than NSHF. 

Same scripts as Figure 4

### Figure 6: Correlations between pairs of timeseries as a function of time-mean period for CMIP6 models 

- field_area_mean_timeseries_CMIP6.py  build 500-year long timeseries of area-means for chosen variables and areas for each CMIP6 model
- timeseries_anom_mth_CMIP6.py to calculate timeseries of monthly anomalies from the time-mean
- timeseries_lsq_fit_multiple_CMIP6.py to calculate the correlation coefficients
- lsq_fit_prep_CMIP6.sh to prepare the files for input to the next script. Uses header_mth.txt as input copied either from header_mth_4_6.txt or header_mth_4_11.txt.
- lsq_fit_plot_revised_CMIP6_with_reanalysis.py to plot out the panels of figure 6. 

### Figure 7: Geographical plots of time-mean NSHF for each CMIP6 model 

- field_time_mean_std_CMIP6.py - calculates time-mean of 500 years of monthly data for each CMIP6 model
- field_plot_difference_CMIP6.py - plots out the time-mean field for each model 

### Figure 8: Correlations between taux or tos annual-mean timeseries and monthly fields of other variables for CMIP6HiRes and reanalysis data

- field_area_mean_timeseries_orig.py to calculate the annual mean timeseries for the selected field and area
- field_timeseries_correlation_v2.py to calculate the fields of correlations 
- field_plot_difference.py to plot the fields 
 
### Figure 9: Correlations between taux monthly-mean anomaly timeseries and monthly anomaly fields of tos for CMIP6 models 
 
- field_area_mean_timeseries_CMIP6.py to calculate the annual mean timeseries for the selected field and area 
- field_time_mean_std_CMIP6.py to calculate fields of time-means and standard deviations for each model 
- field_timeseries_correlation_CMIP6.py to calculate the fields of correlations for each model
- field_plot_difference_CMIP6.py to plot the fields 
 
### Figure 10: Correlations as in Figure 9 when fields and timeseries are passed through a 5-year running mean filter
 
- field_area_mean_timeseries_CMIP6.py to calculate the annual mean timeseries for the selected field and area 
- field_running_mean_CMIP6.py to calculate 5-year running mean fields
- field_time_mean_std_CMIP6.py
- field_timeseries_correlation_CMIP6.py to calculate the fields of correlations for each model
- field_plot_difference_CMIP6.py to plot the fields 

### Figure 11: Correlation and Standard deviations of global-mean tos timeseries and correlation with tos equatorial Pacific timeseries as functions of time-mean period for CMIP6 models
 
- field_area_mean_timeseries_CMIP6.py to calculate the global or equatorial Pacific mean, annual mean timeseries for tos for each CMIP6 model
- timeseries_anom_mth_CMIP6.py to calculate timeseries of monthly anomalies from the time-mean and remove linear trend 
- timeseries_lsq_fit_multiple_CMIP6.py to calculate the correlation coefficients
- lsq_fit_prep_CMIP6.sh to prepare the files for input to: 
- lsq_fit_plot_revised_CMIP6_with_reanalysis.py to plot out the panels of figure 6. 
 
### Figure 12: Timeseries of tos, hfds, taux and tos gradient in equatorial Pacific in HadGEM3 CMIP6 piControl, abrupt-4xCO2, 1pctCO2 and ssp245 simulations 
 
- field_build_annual_decade_means_CMIP6_v2.py to generate annual mean fields for HadGEM3 CMIP6 data set 
- field_area_mean_timeseries_orig.py to generate area means of the annual means
- timeseries_plot_pair_Eq_Pac_paper.py to plot out the timeseries
 
### Figure 13: Scatter plots of selected area-mean, time-mean differences from piControl of various fields from CMIP6 integrations
 
- field_build_annual_decade_means_CMIP6_v2.py to generate decadal mean fields for each CMIP6 data set 
- scatter_plot_CMIP6_decadal_means_v3.py calculate area-means, time-means and differences from piControl and generate scatter plots for a list of CMIP6 integrations
 
### Figure B1: Timeseries of HadGEM3 HighResMIP correlations as a function of time-mean period for alternative area-mean 
 
 Same scripts as Figure 4 (a) - (c).

