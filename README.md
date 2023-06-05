# Modelling pesticide transport during rainfall runoff events with OpenLISEM-pesticide
This repository contains the code and secondary data used to initialize and calibrate 
OpenLISEM-pesticide v1. This model is designed to simulate the runoff, erosion and 
pesticide transport in particulate and dissolved phase during rainfall-runoff events. 
For model simulations and calibration an observational dataset in a small agricultural 
catchment in Limburg, The Netherlands, is used. 

**Authors:**  
 <a href="https://orcid.org/0000-0001-7460-1915">M.C. Commelin (1) <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, 
 <a href="https://orcid.org/0000-0003-1499-1047">J.G. Wesseling (1) <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>,
 <a href="https://orcid.org/0000-0001-6051-8619">J.E.M. Baartman (1) <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>,
 <a href="https://orcid.org/0000-0002-3523-4830">V.G. Jetten (2) <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>,  
1: Soil Physics and Land Management group, Wageningen University, 6708 WG Wageningen, The Netherlands  
2: Faculty of Geo-Information Science and Earth Observation (ITC), University of Twente, Enschede, The Netherlands.     
Contact: meindert.commelin@wur.nl  

**Related publication:**  
Simulating event based pesticide transport with runoff and erosion: OpenLISEM-pesticide v.1  
[add doi]  

**Related datasets and software:**  
The code and executable for OpenLISEM-pesticide on 4TU repository:  
DOI: [add link to OpenLISEM-pesticide on repository]  

The code on github: https://github.com/vjetten/openlisem/tree/lisem_pest  

The observational dataset:  
Runoff, erosion and pesticide transport in a small catchment in Limburg - The Netherlands  
DOI: https://doi.org/10.4121/19690684  


The TCRP v1.1 model is used as additional software:  
Title: Tillage-controlled Runoff Pattern model (TCRP) v1.1; a GIS model to adjust runoff 
patterns for tillage induced roughness in agricultural catchments  
DOI: https://doi.org/10.4121/22153622   
Original TCRP version: https://ees.kuleuven.be/eng/geography/modelling/tcrp/   


## Performing all calculations and analysis
To run the code and analysis for this research the following software and data is required:
 - R, download from: https://cran.r-project.org/
 - RStudio, download from: https://www.rstudio.com/
 - The OpenLISEM-pesticide model
 - The TCRP model
 - the observational dataset  

The observational dataset should be downloaded, and added in a subfolder 'data/'. Then open the R project file: 'openlisem_pesticide_model_development.Rproj'. 
The R script with the full workflow for this study is: 'paper_all_steps.R' In this document further guidance on setup and choices is given.
Initial spatial analysis and transformations have been done separatly and the resulting dataset is provided in 'spatial_data_SL'.
The secondary data used for calculations in this research are stored in 'ext_data/'. Functions in R and intermediate data files are provided in 'sources/'

Required packages can be installed with: 'install.packages("package name")'.

## Session info  
The calculations for the original manuscript where done using the following software and package versions in R.

**R version 4.2.3 (2023-03-15)**

**Versions of loaded packages:** 
_hydroGOF(v.0.4-0)_, _extrafont(v.0.19)_, _gridExtra(v.2.3)_, _cowplot(v.1.1.1)_, _zoo(v.1.8-11)_, _pander(v.0.6.5)_, _forcats(v.1.0.0)_, _stringr(v.1.5.0)_, _dplyr(v.1.1.0)_, _purrr(v.1.0.1)_, _readr(v.2.1.4)_, _tidyr(v.1.3.0)_, _tibble(v.3.2.0)_, _ggplot2(v.3.4.1)_, _tidyverse(v.2.0.0)_, _raster(v.3.6-20)_, _sp(v.1.6-0)_, _sf(v.1.0-11)_, _hms(v.1.1.3)_ and _lubridate(v.1.9.2)_
