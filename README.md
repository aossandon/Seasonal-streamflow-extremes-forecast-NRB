# Seasonal-streamflow-extremes-forecast-NRB
This repository provides the dataset for the Narmada River Basin (NRB), R, and STAN scripts for implementing space-time Bayesian hierarchical modeling (BHM) framework for two flood risk attributes  - seasonal daily maximum flow and the number of events that exceed a threshold during a season (NEETM) - at five gauge locations on the Narmada River Basin proposed in _Ossandón et al. (2022)_. 

## Dataset for the NRB
This dataset contains the files with time series of potential covariates, seasonal daily maximum streamflow and NEETM across the Narmada River basin network (five gauges), India, for the calibration period (2003-2018). It also contains basic information (longitude, latitude, and area) for the gauges considered here. The potential covariates comprise PWPR, Niño 3.4, Niño 1+2, IOD large-scale climate indices, and the the Eastern Pacific Cold Tongue index (EPCT) defined here from March, April, May, or June depending on the lead time (0-, 1-, 2- or 3-month lead time). The daily observed streamflow data were obtained from the [India Water Resources Information System](https://indiawris.gov.in/wris/#/) (India-WRIS).
## Scripts
### R
- **Library.R**: Install and load all the packages required. Also, it contains all the personal functions created for this implementation. 
- **Calibration_Nst_GEV_1cov_sgst_best_dif_month_lead.R**: Fit the best candidate BHMC for different lead times (1 to 10-day lead time). As output, it generates the files **Par_NameCovariates_SgSt.rds** and **post_Q_cali_NameCovariates_SgSt.rds** with the posterior distribution of the model parameters and posterior seasonal maximum streamflow forecasts, respectively. These files are saved on the subdirectory **\results\GEV_kmth** Where k indicates the lead time.
- **Calibration_Nst_Poisson_1cov_sgst_best_dif_month_lead.R**: Fit the best candidate BHMC for different lead times (1 to 10-day lead time). As output, it generates the files **Par_NameCovariates_SgSt.rds** and **post_Q_cali_NameCovariates_SgSt.rds** with the posterior distribution of the model parameters and posterior seasonal maximum streamflow forecasts, respectively. These files are saved on the subdirectory **\results\POI_kmth** Where k indicates the lead time.
### STAN
- **gev_multi_cop_St.STAN**: Contains stan code for the model structure of a stationary BHM considering GEV margins.
- **gev_multi_cop_1cov_Sgst.STAN**: Contains stan code for the model structure of a nonstationary BHM considering GEV margins with one covariate for the location paramater.
- **gev_multi_cop_2cov_Sgst.STAN**: Contains stan code for the model structure of a nonstationary BHM considering GEV margins with two covariates for the location paramater.
- **Poisson_multi_st.STAN**: Contains stan code for the model structure of a stationary BHM considering Poisson margins.
- **Poisson_multi_1cov.STAN**: Contains stan code for the model structure of a stationary BHM considering Poisson margins  with one covariate for the distribution paramater.
- **Poisson_multi_2cov.STAN**: Contains stan code for the model structure of a stationary BHM considering Poisson margins  with two covariates for the distribution paramater.
## References
Ossandón, Á., Rajagopalan, B., & Kleiber, W. (2022). Forecasting magnitude and frequency of seasonal streamflow extremes using a Bayesian hierarchical framework. Under review in Water Resources Research
