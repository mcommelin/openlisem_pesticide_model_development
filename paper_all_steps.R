#' This R code contains the data processing, model runs and figures for: 
#' 
#' Simulating event based pesticide transport with runoff and erosion 
#' OpenLISEM-pesticide v.1
#' 
#' For more information see the readme file in this repository.

# General options --------------------------------------------------------------
#'Executing code in RStudio works best by selecting 1 or multiple lines and than
#'press 'Ctrl + Enter' this will execute the specific lines. This script is 
#'designed to work your way through with this method.
#'
#'To run this code smoothly make choices on the following settings:
#'1. Do you want to get verbose feedback on reading data files?
#'default = FALSE, change to TRUE if you want this data in the console
show_data_info <- FALSE
#'2. set the location where PCRaster is installed.
#'For installation instructions see: 
#'https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_project/install.html
#' e.g. Windows = "C:/MyPrograms/Miniconda3", Linux = "~/miniconda3"
miniconda_dir <- "C:/MyPrograms/Miniconda3"
#' and give the name of the conda environment where PCRaster is installed
#' e.g. 'lisem'
pcr_env <- "lisem"
#'3. give the full path to the OpenLISEM executable
#'See readme files for download information.
# e.g.: "C:/Programs/OpenLISEM/lisem_bin/Lisem.exe"
lisem_dir <- "C:/Users/MC/Werk/OpenLISEM/lisem_bin/Lisem.exe"
#'4. set the directory where the TCRP model scripts are stored
#'See readme files for download information.
# e.g. "C:/Programs/TCRP_v1.1/models"
tcrp_dir <- "C:/Users/MC/Werk/Programs/TCRP_v1.1/models"

# Initialization ---------------------------------------------------------------
# load packages
library(lubridate)
library(hms)
library(sf)
library(raster)
library(tidyverse)
library(pander)

#general settings
#1.
options(readr.show_col_types = show_data_info)
#2.
source("sources/pcrasteR.R")
set_pcraster(miniconda = miniconda_dir, env = pcr_env)

# load functions within this study
# these source files contain much of the main calculations and choices, for a
# better understanding check these functions in detail!
source("sources/functions_dataset_lisem.R") # functions to prepare OpenLISEM input
source("sources/functions_lisem_results.R") # functions for figures and results

#save version of packages use to run the code - added to readme.md
pander(sessionInfo())

# for the correct font on figures we use the package 'extrafont' if you install
# this for the first time, you also need to import the font database with:
#     font_import()
# after this reload extrafont with: 
#     library(extrafont)

# OpenLISEM for runoff and erosion ---------------------------------------------

## TCRP model ------------------------------------------------------------------
#' Adjust the flow paths in the catchment for tillage directions with the 
#' TCRP model. As input we use maps with fields, tillage direction and the 
#' height of this direction (e.g. potato ridges are 25 cm high). When the data
#' belonging to this study is downloaded, the model should work directly.
#' As input for the TCRP model the height of the tillage patterns per crop is
#' needed. This file is stored in the directory '~/sources/tcs_files'
#' - crop_roughness.csv : for each crop with tillage provide RR and RO
#'      RR = random roughness (cm)
#'      RO = oriented roughness (cm)

# only execute once!
# this function also creates a folder structure and converts the input GeoTIFF
# raster maps to the PCRaster .map format.
run_tcrp(tcrp_dir)

## Fixed input -----------------------------------------------------------------
#' The next step is to produce all input maps from variables that we will fix.
#' These will not be used for calibration. Several input files are called. These
#' are stored in the directory '~/sources/tcs_files'
#' The files are:
#' - vars_landuse.csv : containing LAI and RR for all landuse types. For
#'        arable crops the LAI is time dependend and will be calculated based on
#'        a simple crop model from SWAP. The nodata value in this files is 
#'        always '-77'
#' - crop_images_yyyy.txt : for each event observation data on crop variables:
#'        CH (crop height, m), cover (% surface covered with vegetation), 
#'        nresi (the mannings n added by residues), nveg (the mannings n added 
#'        by vegetation)
#' - crops.csv : a file with start and end date of each crop in the catchment
#'        for each growing season.
#' - lu_yyyy_croptype.csv : for each season a crop_type - lu_nr table, provided 
#'        with the original dataset
#' - homogeneous_vars.txt : a text file providing the values for all variables
#'        that are modelled homogeneous over the whole catchment.

fixed_input_lisem()

## Rainfall data ---------------------------------------------------------------
#' OpenLISEM also requires rainfall input data. In the dataset 2 sources of 
#' rainfall are available. These sources are measured at several locations.
#'  
#' Choose which source of rainfall data should be used for each event:
#'  1 = KNMI radar : available for all events, the mean of the 6 clossest 
#'                  points to the catchment, risk of wrong timing and under
#'                  estimation of intensity
#'  2 = Tipping bucket : tipping bucket data collected at the outlet and 2 point 
#'                  in the catchment, this timeseries is not complete and can not 
#'                  always be used.The tips are linearly interpolated and if 
#'                  available might provide better estimates of rainfall than KNMI radar.
#' And choose the method of spatial representation for each event:
#'  1 = mean of points covering the catchment
#'  2 = Link every cell to the closest measurement point with Voronoi polygons.

rainfall_input(rain_source = c(2,1), method = c(1,1))

## Initial run OpenLISEM -------------------------------------------------------
#' The variables that will be used for calibration of runoff and erosion are 
#' produced with the next function. A part of the input is still fixed, but
#' other settings can be updated to calibrate the model. The input files are:
#' - vars_landuse.csv : containing values for mannings n, initial soil
#'        moisture, hydraulic conductivity of the upper soil layer and root 
#'        cohesion of the crops. This file contains the values that are fixed, 
#'        all values that are calibrated should have '-77' in this file.
#' - homogeneous_vars.txt : a text file providing the values for cohesion 
#'        variables that are modeled homogeneous over the whole catchment.
#' - calibration_file : the name of this file is input for the following
#'        function, so multiple files can be made to store different settings.
#'        Three examples are given in the 'sources/tcs_files' folder with:
#'        - tcs_initial_calibration.txt : this are the NOT calibrated settings
#'            of the first run.
#'        - manual_calibration.txt : this file contains the settings after 
#'            manual calibration, this was the starting point for pesticide 
#'            simulations.

cal_file <- "tcs_initial_calibration.txt"
# this function makes the input raster maps
initial_set <- sed_water_cal_input_lisem(calibration_file = 
                                           paste0(tcs_dir, cal_file))

# run settings file (here you can set the timestep)
run_set <- "manual_run_settings.txt"
# this function will make the runfile for OpenLISEM
make_runfile_lisem(run_in = paste0(tcs_dir, run_set))

# run OpenLISEM with initial settings (~10 min runtime)
# the wait command defines if R wait for OpenLISEM to finish running.
# if you have a strong PC using 'c(FALSE, TRUE)' decrease the computation time
# if you want to see the user interface during the model run, set GUI = TRUE
wait <- c(TRUE, TRUE)
for (i in seq_along(event_sel)) {
  run_file <- dir(main_dirs[i], pattern = ".run$")
  run_lisem(main_dirs[i], run_file, lisem_dir, wait[i], GUI = FALSE)
}

# export initial run performance to table
events <- ymd(c("2019-05-28", "2020-08-16"))
base_dir <- paste0("lisem_runs/", year(events), "/", 
                   str_remove_all(as.character(date(events)), "-"))
initial_perf <- figures_results_paper(figure = 1, events = events,
                                      base_dir = base_dir)
# this data is stored after we also have the calibrated results.

## Calibration for Runoff and erosion ------------------------------------------

# create input for the manual calibrated runs
# this will overwrite the input and result directories of the initial runs!!!
cal_file <- "tcs_manual_calibration.txt"
manual_set <- sed_water_cal_input_lisem(calibration_file = 
                                          paste0(tcs_dir, cal_file))

# make runfile
run_set <- "manual_run_settings.txt"
make_runfile_lisem(run_in = paste0(tcs_dir, run_set))

# run OpenLISEM with calibrated settings
wait <- c(TRUE, TRUE) 
for (i in seq_along(event_sel)) {
  run_file <- dir(main_dirs[i], pattern = ".run$")
  run_lisem(main_dirs[i], run_file, lisem_dir, wait[i], GUI = FALSE)
}

# combine input settings table for initial and calibrated model runs.
initial_set <- initial_set %>%
  mutate(version = "initial")
manual_set <- manual_set %>%
  mutate(version = "manual")
comb_set <- bind_rows(initial_set, manual_set) %>%
  pivot_longer(cols = ksat1:d50, values_to = "value", names_to = "var") %>%
  mutate(date = year(date),
         value = round(value, digits = 2)) %>%
  pivot_wider(names_from = c(date, version), values_from = value) %>%
  select(var, crop_type, `2019_initial`, `2019_manual`, 
         `2020_initial`, `2020_manual`) %>%
  arrange(var)
# save input settings to the results directory.
write_csv(comb_set, "results/settings_runoff_erosion.csv")


### make figure 2 --------------------------------------------------------------
events <- ymd(c("2019-05-28", "2020-08-16"))
base_dir <- paste0("lisem_runs/", year(events), "/", 
                   str_remove_all(as.character(date(events)), "-"))
manual_perf <- figures_results_paper(figure = 2, fig_name = "manual",
                                     events = events, base_dir = base_dir)

# show performance table
p1 <- initial_perf %>%
  mutate(date = event_sel,
         version = "initial")
p2 <- manual_perf %>%
  mutate(date = event_sel,
         version = "manual")

comb_perf <- bind_rows(p1, p2) %>%
  select(version, date, everything())
# save performance statistics of initial and calibrated simulations
write_csv(comb_perf, "results/performance_runoff_erosion.csv")

# OpenLISEM for pesticides -----------------------------------------------------
#' at this stage we have simulations for runoff and erosion, which perform 
#' adequatly. The next step is to add the pesticide simulations.

## Initial simulations ---------------------------------------------------------

#' This function produces the pesticide maps, and adjusts the runfile to model
#' the correct rainfall event. The following file and settings are required:
#' - pesticide_calibration_file : the name of this file is input for the 
#'        following function, so multiple files can be made to store different 
#'        settings. Three examples are given in the 'sources/tcs_files' folder 
#'        with:
#'        - pesticide_initial_calibration.txt : this are the NOT calibrated 
#'            settings of the first run with pesticide values, based on the
#'            manual calibration of runoff and erosion.
#'        - pesticide_manual_calibration.txt : this file contains the settings 
#'            after manual calibration, based on the CRS calibration of runoff
#'            and erosion.
#' - run_settings_file : a file with the timestep for the pesticide model
 

# to evaluate the model performance first prepare the observations
# this only needs to be done once
pesticides <- c("Glyphosate", "Metobromuron")
load_observations(pesticides = pesticides)

#'For this study we use simulate 2 different pesticides, in the loop below both
#'pesticides are simulated for both events, which is 4 model runs in total.
#'The function 'compound_pest_input()' makes the dataset and runfile for
#'pesticide simulations. the first option is 'pest_in' which three options:
# No pesticides = 1
# Glyphosate = 2
# Metobromuron = 3

# we run the loop, for Glyphosate and Metobromuron.
pesticides <- c(2,3)

# Should map series for pesticides be saved? This takes more memory and
# computation time, but gives the oppurtunity to analyse dynamics in more detail
# afterwards.
report_maps <- FALSE # set to TRUE to save maps
# input files for OpenLISEM runs
cal_file <- "pesticide_initial_calibration.txt"
run_set <- "manual_run_settings.txt"

# loop OpenLISEM with initial settings (result = setting & performance)
ini_set_pest <- vector("list", length = length(pesticides))
pest_ini_perf <- vector("list", length = length(pesticides))
events <- ymd(c("2019-05-28", "2020-08-16"))
base_dir <- paste0("lisem_runs/", year(events), "/", 
                   str_remove_all(as.character(date(events)), "-"))

for (i in seq_along(pesticides)) {
  
  # create input and runfile
  ini_set_pest[[i]] <- compound_pest_input(pesticides[i], 
                        pest_cal_file = paste0(tcs_dir, cal_file),
                        report = report_maps, run_in = paste0(tcs_dir, run_set))
  
# run lisem
  for (j in seq_along(event_sel)) {
    run_file <- dir(main_dirs[j], pattern = ".run$")
    run_lisem(main_dirs[j], run_file, lisem_dir, GUI = FALSE)
  }
  
  # store performance
  pest_ini_perf[[i]] <- figures_results_paper(figure = 3, fig_name = "manual",
                                          events = events, base_dir = base_dir)
  
}
initial_settings_pest <- bind_rows(ini_set_pest)
initial_perf_pest <- bind_rows(pest_ini_perf)

## Calibrated results ----------------------------------------------------------

#' we repeat the loop over two events and two pesticides, but now with 
#' calibrated settings. If you want to try yourself different setting. This
#' can be change in 'sources/tcs_files/'
#' find the file 'pesticide_manual_calibration.txt' save it with a different
#' name and adjust the name below
cal_file <- "pesticide_manual_calibration.txt"

# Should map series for pesticides be saved?
report_maps <- FALSE # set to TRUE to save maps
run_set <- "manual_run_settings.txt"
pesticides <- c(2,3)
pest_names <- c("Glyphosate", "Metobromuron")
#loop over pesticides and events
man_set_pest <- vector("list", length = length(pesticides))
pest_man_perf <- vector("list", length = length(pesticides))
res_error <- vector("list", length = length(pesticides) * length(event_sel))
res_hydrograph <- vector("list", length = length(pesticides) *length(event_sel))
res_totals <- vector("list", length = length(pesticides) * length(event_sel))
for (i in seq_along(pesticides)) {
  man_set_pest[[i]] <- compound_pest_input(pesticides[i], 
                                  pest_cal_file = paste0(tcs_dir, cal_file),
                                  report = report_maps, 
                                  run_in = paste0(tcs_dir, run_set))
  
  events <- ymd(c("2019-05-28", "2020-08-16"))
  base_dir <- paste0("lisem_runs/", year(events), "/", 
                     str_remove_all(as.character(date(events)), "-"))
  

  # run lisem
  for (j in seq_along(event_sel)) {
    run_file <- dir(main_dirs[j], pattern = ".run$")
    run_lisem(main_dirs[j], run_file, lisem_dir, GUI = FALSE)
    #store results for later
    x <- matrix(seq(1,4,1), nrow = 2)
    k <- x[j,i]
    res_totals[[k]] <- store_result_total(main_dirs[j])
    res_hydrograph[[k]] <- store_result_hydrograph(main_dirs[j]) %>%
      mutate(compound = pest_names[i])
    res_error[[k]] <- store_result_error(main_dirs[j])
    # save list to .Rdata file, can be used for further analysis later
    save(res_error, res_hydrograph, res_totals, file = "results/pest_man.Rdata")
  }
  
  # store performance
  pest_man_perf[[i]] <- figures_results_paper(figure = 3, fig_name = "manual",
                                          events = events, base_dir = base_dir)
}
manual_settings_pest <- bind_rows(man_set_pest)
manual_perf_pest <- bind_rows(pest_man_perf)


#combine settings and performance of initial and calibrated simulations

# combine setting table
s1 <- initial_settings_pest %>%
  mutate(version = "initial")
s2 <- manual_settings_pest %>%
  mutate(version = "manual")
comb_set <- bind_rows(s1, s2) %>%
  pivot_longer(cols = conc_ads:kd, values_to = "value", names_to = "var") %>%
  mutate(date = year(date),
         value = round(value, digits = 2)) %>%
  pivot_wider(names_from = c(date, version), values_from = value) %>%
  select(compound, var, large, `2019_initial`, `2019_manual`, `2020_initial`, `2020_manual`) %>%
  arrange(var)
# save the pesticide settings
write_csv(comb_set, "results/settings_pesticides.csv")

# performance
dates <- rep(event_sel, 2)
# show performance table
p1 <- initial_perf_pest %>%
  mutate(date = dates,
         version = "initial")
p2 <- manual_perf_pest %>%
  mutate(date = dates,
         version = "manual")

comb_perf <- bind_rows(p1, p2) %>%
  select(version, date, everything())
# save the performance statistics of the pesticide simulations
write_csv(comb_perf, "results/performance_pesticides.csv")

### make figure 6 --------------------------------------------------------------
events_all <- rep(ymd(c("2019-05-28", "2020-08-16")), 2)
base_dir_all <- paste0("lisem_runs/", year(events_all), "/", 
                   str_remove_all(as.character(date(events_all)), "-"))
load("results/pest_man.Rdata")
pest_plots_man <- vector("list", length = length(res_error))
add_legend <- c(F,F,T,F)
for (i in seq_along(res_error)) {
  pest_plots_man[[i]] <- figures_results_paper(figure = 5,
                                               events = events_all[i],
                                               base_dir = base_dir_all[i],
                                               list_pos = i, 
                                               legend = add_legend[i])
}

plot <- plot_grid(pest_plots_man[[3]], pest_plots_man[[4]], 
                  pest_plots_man[[1]], pest_plots_man[[2]],
                  nrow = length(res_error),
                  labels = c("A: MET 28-05", "B: MET 16-08", 
                             "C: GLY 28-05", "D: GLY - 16-08"),
                  align = "hv", axis = "rl",
                  label_fontfamily = "Times New Roman",
                  label_x = 0.2, label_y = 1, label_size = 11) 

ggsave(plot = plot, filename = "images/figure6.tiff", 
       width = 170, height = 240, units = "mm",
       dpi = 600, device = "tiff")

## Sensitivity analysis --------------------------------------------------------

# the content of the event folder should be the calibrated pesticide run!
# We create the correct here!
report_maps <- FALSE # set to TRUE to save maps
cal_file <- "pesticide_manual_calibration.txt"
run_set <- "manual_run_settings.txt"
pesticides <- c(3)
compound_pest_input(pesticides[1], pest_cal_file = paste0(tcs_dir, cal_file),
                      report = report_maps, run_in = paste0(tcs_dir, run_set))

# create a new directory to execute the sensitivty analysis
sens_dir <- "lisem_runs/sensitivity_oat/"
# link to the file with the input for all model simulations
run_vars <- "sources/tcs_files/sensitivity_analysis_input.csv"
# this is a list with standard map names of variables that will be adjusted
maps_names <- c("pcmixsoil.map", "pestmixdep.map", "pestsoildep1.map", 
                "ksat1.map", "coh.map", "n.map", "psi1.map")
# select the event as base for OAT analysis : 2019-05-28
event <- event_sel[1]

#do sensitivity analysis
#...... WARNING ...................................................
# This function only works if the runfile from the selected event
# runs without warning or error message in OpenLISEM with GUI!
# Executing this function takes about 2 - 4 hours, the results are
# already given in 'results/sensitivity.Rdata'
#..................................................................
oat_sensitivity_OLP(run_vars, map_names, sens_dir, event)

# analyse and present results of the sensitivity analysis.
sens_vars <- read_delim(run_vars, delim = "\t")

# collect the values of the base run to compare with the sensitivity results
events <- ymd(c("2019-05-28"))
base_dir <- paste0("lisem_runs/", year(events), "/", 
                   str_remove_all(as.character(date(events)), "-"))
run_file <- dir(base_dir, pattern = ".run$")
run_lisem(main_dirs[1], run_file, lisem_dir, GUI = FALSE)

pest_perf_base <- figures_results_paper(figure = 3, fig_name = "oat_sens",
                                        events = events, base_dir = base_dir) %>%
  mutate(run = 0)

# loop over sensitivity results to calculate performance statistics
sens_perf <- vector("list", length = nrow(sens_vars))
base_dir <- "lisem_runs/sensitivity_oat/20190528"
for(i in 1:nrow(sens_vars)) {
  sens_perf[[i]] <- figures_results_paper(figure = 4,
                                          events = events, base_dir = base_dir,
                                          list_pos = i)
}

# combine everything into one table
perf_sens_tab <- bind_rows(sens_perf) %>%
  mutate(run = 1:nrow(.)) %>%
  left_join(sens_vars, by = "run") %>%
  bind_rows(pest_perf_base) %>%
  arrange(run) %>%
  mutate(across((PMerr:PMwerr), ~signif(., digits = 1))) %>%
  mutate(across((DP_tot_sim:PP_mod_sim), ~round(., digits = 0)),
         across((KGEs_pest:RMSEw_pest), ~round(., digits = 2)),
         rel_DP = round(DP_tot_sim / DP_tot_sim[1] * 100),
         rel_PP = round(PP_mod_sim / PP_mod_sim[1] * 100)) %>%
  select(var, val, everything(), -compound, - run)
# save the results of the sensitivity analysis.
write_csv(perf_sens_tab, "results/oat_sensitivity_analysis.csv")
