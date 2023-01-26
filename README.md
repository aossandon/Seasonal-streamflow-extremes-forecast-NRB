# Seasonal-streamflow-extremes-forecast-NRB
This repository provides the dataset for the Narmada River Basin (NRB), R, and STAN scripts for implementing space-time Bayesian hierarchical modeling (BHM) framework for two flood risk attributes  - seasonal daily maximum flow and the number of events that exceed a threshold during a season (NEETM) - at five gauge locations on the Narmada River Basin proposed in _Ossandón et al. (2022)_. 

## Dataset for the NRB
This dataset contains the files with time series of potential covariates, seasonal daily maximum streamflow and NEETM across the Narmada River basin network (five gauges), India, for the calibration period (2003-2018). It also contains basic information (longitude, latitude, and area) for the gauges considered here. The potential covariates comprise PWPR, Niño 3.4, Niño 1+2, IOD large-scale climate indices, and the the Eastern Pacific Cold Tongue index (EPCT) defined here from March, April, May, or June depending on the lead time (0-, 1-, 2- or 3-month lead time). The daily observed streamflow data were obtained from the [India Water Resources Information System](https://indiawris.gov.in/wris/#/) (India-WRIS).
## Scripts
### R
- Library.R: Install and load all the packages required. Also, it contains all the personal functions created for this implementation. 
- Calibration_Nst_GEV_1cov_sgst_best_dif_month_lead.R: Fit the best candidate BHMC for different lead times (1 to 10-day lead time). As output, it generates the files _Par_NameCovariates_SgSt.rds_ and _post_Q_cali_NameCovariates_SgSt.rds_ with the posterior distribution of the model parameters and posterior seasonal maximum streamflow forecasts, respectively. These files are saved on the subdirectory _\results\GEV_kmth_ Where k indicates the lead time.
- Get_forecast_2021.R: Generate the ensemble streamflow forecast for the peak monsoon season 2021. Before running this script, the script Par_bayesian_model_gamma_best_mod_kLT.rds must be run for 1 to 5-days lead time. 
### STAN
- bayesian_model_gamma_2step_diff_st_2_best_kdlt.STAN: Contain stan code for the model structure of the best candidate BHMC for each lead time. Where k indicates the lead time.
## References
Ossandón, Á., Rajagopalan, B., & Kleiber, W. (2022). Forecasting magnitude and frequency of seasonal streamflow extremes using a Bayesian hierarchical framework. Under review in Water Resources Research
