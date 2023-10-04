#prepare a input dataset for lisem based on temporal constraints

# Initialization ----------------------------------------------------------

library(lubridate)
library(hms)
library(sf)
library(raster)
library(tidyverse)

# functions
# discharge functions
source("sources/discharge_calculations.R")

# function to run OpenLLISEM-pesticide from the command line
# dir = the directory where the runfile is located
# run_file = runfile name (with .run extension!)
# lisem_dir = full path name to Lisem executable
# GUI = show the user interface (default = TRUE)
# wait = should R wait until Lisem is finished? (default = TRUE)

run_lisem <- function(dir, run_file, lisem_dir, GUI = TRUE, wait = TRUE) {
  projwd <- getwd()
  run <- paste0(projwd, "/", dir, run_file)
  if (!GUI) {gui <-  "-ni"} else {gui <- ""}
  command <- paste0(gui, " -r ", run)
  print("Running OpenLISEM...")
  t_s <- Sys.time()
  system2(lisem_dir, command, wait = wait)
  t_e <- Sys.time()
  t_r <- t_e - t_s
  print(paste0("Run finished!"))
  return(t_r)
}

#' Give all words in a character column CAPS at start of word.
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

nl_tz <- locale(tz = "Etc/GMT-1")

## load data - make folders --------------

#load data for the two events for this study
ev_dates <- read_csv("sources/ev_dates.csv")
event_sel <- ymd(c("2019-05-28", "2020-08-16"))
seasons <- c("2019", "2020")
tcs_dir <- "sources/tcs_files/"
main_dirs <- c("lisem_runs/2019/20190528/", "lisem_runs/2020/20200816/")

event_pos <- c(3, 17)
# check if lisem directory exists, otherwise make
if (!dir.exists("lisem_runs/")) {
  dir.create("lisem_runs/")
}
for (j in seq_along(seasons)) {
  # check if lisem directories exist, otherwise make
  if (!dir.exists(paste0("lisem_runs/", seasons[j]))) {
    dir.create(paste0("lisem_runs/", seasons[j]))
  }
  # check if tcrp directory exists, otherwise make
  if (!dir.exists(paste0("lisem_runs/", seasons[j], "/tcrp/"))) {
    dir.create(paste0("lisem_runs/", seasons[j], "/tcrp/"))
  }
  # check if event directories exist, otherwise make
  event <- event_sel[j]
  if (!dir.exists(str_c("lisem_runs/", year(event), "/", 
                        str_remove_all(as.character(date(event)), "-")))) {
    dir.create(str_c("lisem_runs/", year(event), "/", 
                     str_remove_all(as.character(date(event)), "-")))
  }
}

# load rainfall
knmi_p <- read_csv("ext_data/KNMI_rain.csv") %>%
  group_by(timestamp) %>%
  summarize(mean(P))
names(knmi_p) <- c("timestamp", "P_m_knmi")
knmi_p <- mutate(knmi_p, P_i_knmi = P_m_knmi * 12)

#cumulative rainfall
knmi_day <- knmi_p %>%
  mutate(date = date(timestamp)) %>%
  group_by(date) %>%
  summarise(P = sum(P_m_knmi, na.rm = T)) %>%
  mutate(Pcum = cumsum(P))

# Run TCRP ------------------------------------------------------

run_tcrp <- function(tcrp_dir) {
#To run the TCRP model and prepare the input for OpenLISEM the dataset
# 'spatial_data_SL' must be placed in the main project directory. (default)

# for each growing season, convert the tif maps to PCRaster format
# loop over years
for (j in seq_along(seasons)) {
  # map names we want to load and convert
  map_dir <- paste0("spatial_data_SL/rasters", seasons[j])
  files <- unique(dir(map_dir, pattern = ".tif$"))
  maps <- str_remove(files, ".tif")
  season_dir <- paste0("lisem_runs/", seasons[j], "/tcrp/")
#change .tif to .asc and make a mask map.
for (i in seq_along(maps)) {
  ras <- raster(paste0(map_dir, "/", files[i]))
  store <- paste0(season_dir, maps[i], ".asc")
  writeRaster(ras, store, format = "ascii", NAflag = -9999, overwrite = TRUE)
  # save location attributes of main map to make a mask map
  if (maps[i] == "dem") {
    nc <- ras@ncols
    nr <- ras@nrows
    cs <- xres(ras)
    xulc <- ras@extent@xmin
    yulc <- ras@extent@ymax
    #make new pcr map
    newmap(mapname = "mask.map", nrows = nr, ncols = nc, datatype = "S",
           xulc = xulc, yulc = yulc, cellsize = cs, dir = season_dir)
  }
}
for (i in seq_along(maps)) {  
# convert to .map and resample
  map_in <- paste0(maps[i], ".asc")
  map_tmp <- "tmp.map"
  map_out <- paste0(maps[i], ".map")
  subdir <- season_dir
  asc2map(map_in = map_in, map_out = map_tmp,
          sub_dir = subdir)
  resample_location(map_in = "tmp.map", map_out = map_out, dir = subdir)
  file.remove(paste0(subdir, map_tmp))
  file.remove(paste0(subdir, map_in))
}
}

# make additional maps and run TCRP model
layers <- paste0("landuse_", seasons)
for (i in seq_along(layers)) {
  lu <- st_read(dsn = "spatial_data_SL/fields.gpkg", layer = layers[i])
  rough <- read_csv(paste0(tcs_dir, "crop_roughness.csv"))
  lu_tbl <- lu %>%
    select(lu_nr, ro) %>%
    st_drop_geometry() %>%
    mutate(ro = if_else(is.na(ro), 0, ro)) %>%
    distinct() %>%
    rename_with(~as.character(seq(0,1,1)))
  rough_name <- paste0("lisem_runs/", seasons[i], "/tcrp/roughness.tbl")
  write_delim(lu_tbl, rough_name, col_names = T)
}

# prepare tcrp input maps
script <- "prepare_tcrp.mod"
script_dir <- "sources/pcr_scripts"
for (i in seq_along(seasons)) {
work_dir <- paste0("lisem_runs/", seasons[i], "/tcrp")
pcr_script(script = script, work_dir = work_dir, script_dir = script_dir)
}
#execute tcrp scripts
tcrp_order <- c("LddTill", "PitRem", "CirkRem", "KruisRem", "Logit", "LddRes")
script_dir <- tcrp_dir
for (j in seq_along(seasons)) {
  for (i in seq_along(tcrp_order)) {
    script <- paste0(tcrp_order[i], ".mod")
    work_dir <- paste0("lisem_runs/", seasons[j], "/tcrp")
    pcr_script(script = script, work_dir = work_dir, script_dir = script_dir,
               script_path_rel = FALSE)
  }
}

} # end function run_tcrp()

# OpenLISEM input ---------------------------------------------------------

fixed_input_lisem <- function() {
# make lisem input from tcrp output
script <- "tcrp_to_lisem.mod"
script_dir <- "sources/pcr_scripts"
for (i in seq_along(seasons)) {
  work_dir <- paste0("lisem_runs/", seasons[i], "/tcrp")
  pcr_script(script = script, work_dir = work_dir, script_dir = script_dir)
}

# copy files to lisem_runs
fcopy <- readLines("sources/tcs_files/copy_maps.txt")[3:14]
for (j in seq_along(seasons)) {
  event <- event_sel[j]
  # check if lisem directories exist, otherwise make
  new_dir <- paste0("lisem_runs/", seasons[j], "/", 
                    str_remove_all(as.character(date(event)), "-"))
  cop <- paste0("lisem_runs/", seasons[j], "/tcrp", "/", fcopy)
  file.copy(cop, new_dir, overwrite = TRUE)
}

## Input fixed variables -------------

# data on start and end of each crop per season
crops <- read_csv(paste0(tcs_dir, "crops.csv")) %>%
  mutate(cropstart = dmy(cropstart),
         cropend = dmy(cropend))

# LAI based on simple crop model from SWAP
lai <- read_csv(paste0(tcs_dir, "lai_crops.txt"))

# fixed crops variables
lu_input <- read_csv(paste0(tcs_dir, "vars_landuse.csv")) %>%
  select(crop_type, LAI, rr) %>%
  filter(crop_type != "wheeltracks")
#main_dirs <- c("lisem_runs/2019/20190528/", "lisem_runs/2020/20200816/")

## RR data
rough <- read_csv(paste0(tcs_dir, "crop_roughness.csv"))
# calculate vars for RR decay function by Potter 1990
clay <- 3.5    # clay content of the soil
oc <- 3.6      # organic matter content
b <- 63 + 62.7 * log(oc) + 15.7 * clay - 0.25 * clay^2

cropsel <- unique(crops$crop_type)
aggr_stab = 14

# loop over events to make a unique input table for each event
# and make all pcraster maps
script <- "fixed_input_maps.mod"
script_dir <- "sources/pcr_scripts"

for (j in seq_along(event_sel)) {
  event <- event_sel[j]
  rr_tbl <- c()
  
  # data on cover and crop height
  crop_images <- read_csv(paste0(tcs_dir, "crop_images_", seasons[j], ".txt"))
  
  ### LAI interpolate ---------------------------------------------------------
  lai_event <- c()
  for (i in seq_along(cropsel)) {
    lai_crop <- lai %>%
      filter(crop_type == cropsel[i])
    crops_event <- crops %>%
      filter(year(cropend) == year(event)) %>%
      filter(crop_type == cropsel[i])
    
    lai_date <- tibble(date = seq(from = crops_event$cropstart[1], 
                                  to = crops_event$cropend[1], by = 1),
                       crop_stage = round(seq(from = 0, to = 2, 
                                              length.out = length(date)), digits = 2)) %>%
      bind_rows(lai_crop)%>%
      arrange(crop_stage) %>%
      mutate(LAI = round(na.approx(LAI, rule = 2), digits = 2),
             crop_type = cropsel[i]) %>%
      filter(date == event)
    
    lai_event <- bind_rows(lai_event, lai_date)
  }
  # from images it is clear that the potatoes crop on 20200816 is already dead, so LAI = 0.85
  if (event == ymd("2020-08-16")) {
    lai_event <- lai_event %>%
      mutate(LAI = if_else(crop_type == "potato", 0.85, LAI))
  }
  
  # add cover, crop height and calculate canopy storage
  # TODO: soft-code the LAI for grass and fraction of trees in apple var calculations.
  lai_event <- lai_event %>%
    full_join(crop_images, by = "crop_type") %>%
    mutate(LAIn = if_else(crop_type == "apple", 0.33 * LAI + 0.67 * 3, LAI),
           LAIn = if_else(is.na(LAIn), 1, LAIn),
           Smax = 0.935 + 0.498 * LAIn - 0.00575 * LAIn^2) %>%
    mutate(ch = if_else(crop_type == "apple", 0.33 * ch + 0.67 * 0.1, ch),
           cover = if_else(crop_type == "apple", 1, cover),
           crop_type = if_else(crop_type == "wheatw", "cereals", crop_type)) %>%
    select(crop_type, LAIn, ch, cover, Smax)
  
  lu_tmp1 <- lu_input %>%
    left_join(lai_event, by = "crop_type") %>%
    mutate(LAI = if_else(LAI == -77, LAIn, LAI)) %>%
    select(-LAIn)
  
  ## RR calculation -------------
  # only for arable crops (wheat, potato, sugarbeet)
  rr_adj <- c()
  rr_crops <- c("potato", "sugarbeet", "cereals")
  for (k in seq_along(rr_crops)) {
    crops_event <- crops %>%
      filter(year(cropend) == year(event)) %>%
      mutate(crop_type = if_else(crop_type == "wheatw", "cereals", crop_type)) %>%
      filter(crop_type == rr_crops[k]) %>%
      mutate(pcum_til = filter(knmi_day, date == cropstart)$Pcum,
             pcum_ev = filter(knmi_day, date == event)$Pcum,
             pcum = pcum_ev - pcum_til,
             rri = filter(rough, crop_type == rr_crops[k])$rr,
             rrt = round(rri * exp(-(pcum/b)^0.6), digits = 1),
             rrt = if_else(rrt < 0.5, 0.5, rrt)) %>%
      select(crop_type, rrt)
    rr_adj <- bind_rows(rr_adj, crops_event)
    
  }
  lu_tmp2 <- lu_tmp1 %>%
    left_join(rr_adj, by = "crop_type") %>%
    mutate(rr = if_else(!is.na(rrt), rrt, rr)) %>%
    select(-rrt)
  
  # combine the crop parameters with the landuse classes to make a PCR table
  layers <- paste0("landuse_", year(event_sel))
  lu <- st_read(dsn = "spatial_data_SL/fields.gpkg", layer = layers[j])
  lu_crop_type <- read_csv(paste0(tcs_dir, "lu_", seasons[j], "_croptype.csv"))
  
  lu_out <- lu %>%
    select(-rr, -crop_type) %>%
    left_join(lu_crop_type, by = "lu_nr") %>%
    left_join(lu_tmp2, by = "crop_type") %>%
    st_drop_geometry() %>%
    filter(!is.na(crop_type)) %>%
    select(lu_nr, LAI, Smax, ch, cover, rr) %>%
    mutate(LAI = if_else(is.na(LAI), 0, LAI),
           ch = if_else(is.na(ch), 0, ch),
           cover = if_else(is.na(cover), 0, cover),
           Smax = if_else(is.na(Smax), 0, Smax),
           aggr = aggr_stab) %>%
    rename_with(~as.character(seq(0,length(.)-1,1))) %>%
    mutate(across(everything(), ~round(., digits = 1))) %>%
    distinct()
  
  write_delim(lu_out, paste0(main_dirs[j],"fixed.tbl"), col_names = T)
  
  #make homogeneous.tbl
  hom_tbl <- read_csv(paste0(tcs_dir, "homogeneous_vars.txt"), skip = 2) %>%
    pivot_wider(names_from = "var", values_from = "val") %>%
    mutate(area = 1) %>%
    select(area, everything()) %>%
    rename_with(~as.character(seq(0,length(.)-1,1)))
  
  write_delim(hom_tbl, paste0(main_dirs[j],"homogeneous.tbl"), col_names = T)
  
  # execute pcr script
  work_dir <- main_dirs[j]
  pcr_script(script = script, work_dir = work_dir, script_dir = script_dir)
}

} # end function fixed_input_lisem()


# Rainfall input ----------------------------------------------------------

rainfall_input <- function(rain_source = c(1,1), method = c(1,1)) {
## Rainfall files -------------
# combine different rainfall sources (KNMI_radar and CR6 - tipping bucket)
wb_data <- read_csv("data/WB_data.csv") # data from waterboard (Ha)
cr6 <- read_csv("data/CR6.csv", col_types = "TddddddddddddddD") # from CR6 (Hb)
crest_nap <- 79.331 # m + NAP as measured by the waterboard

cr1000_rain <- read_csv("data/cr1000_rain.csv") # 2 additional tipping buckets

# use Q_calc to calculate discharge for each event and store in list
# the dicharge calculations can be found in the 'sources' folder for more details.
Q_tcor <- Q_calc(ev_dates, wb_data, cr6, crest_nap)
df_dat <- Q_tcor[[1]] # data with discharge for each event
df_tcor <- Q_tcor[[2]] # based on synchronization between different time series a time correction is performed

# load rainfall
knmi_radar <- read_csv("ext_data/KNMI_rain.csv") %>%
  mutate(P = round(P * 12, digits = 1)) %>% # calculate intensity
  filter(loc_id > 84) # point 83 and 84 are not used in the catchment

# script for id.map
script_dir <- "sources/pcr_scripts"

for (i in seq_along(event_sel)) {
  event <- event_sel[i]
  event_time <- filter(ev_dates, date == event)
  spatial <- method[i]
  
  rainstart <- knmi_p %>%
    filter(timestamp >= event_time$start_t & timestamp <= event_time$end_t) %>%
    filter(!is.na(P_i_knmi)) %>%
    filter(P_i_knmi > 0) %>%
    summarise(rainstart = min(timestamp))
  
  # rain knmi
  if (rain_source[i] == 1) { 
    source_rain <- "KNMI radar"
    if (spatial == 2) {
      rain <- knmi_radar %>%
        filter(timestamp >= event_time$start_t & timestamp <= event_time$end_t) %>%
        pivot_wider(names_from = loc_id, values_from = P) %>%
        mutate(mins = (as.numeric(timestamp) - as.numeric(rainstart$rainstart))/60) %>%
        arrange(mins) %>%
        filter(mins >= 0) %>%
        select(mins, everything(), -timestamp)
    }
    if(spatial == 1) {
      rain <- knmi_radar %>%
        filter(loc_id == 85 | loc_id == 86) %>%
        filter(timestamp >= event_time$start_t & timestamp <= event_time$end_t) %>%
        mutate(mins = (as.numeric(timestamp) - as.numeric(rainstart$rainstart))/60) %>%
        group_by(mins) %>%
        summarise(P = mean(P)) %>%
        arrange(mins) %>%
        filter(mins >= 0)
    }
    
  }
  if (rain_source[i] == 2) {
    source_rain <- "Tipping buckets"
    
    # join tipping bucket rain data
    rain_cr6 <- df_dat[[event_pos[i]]] %>%
      select(timestamp, rain_tot) %>%
      mutate(loc_id = 89)
    rain_cr1000 <- cr1000_rain %>%
      filter(timestamp >= event_time$start_t & timestamp <= event_time$end_t) %>%
      select(timestamp, loc_id, rain_tot)
    
    rain_raw <- bind_rows(rain_cr6, rain_cr1000)
    stations <- unique(rain_raw$loc_id)
    #add time corrections for the 2 cr1000 tipping buckets. There is drift on the timer
    # however we don't know how much. So we synchronize the first 12 mm/h minute! not sure if that is correct
    tcor_tip <- c(0, -7, -4)
    
    # make loop for tipping interpolation
    df_rain <- vector("list", length = length(stations))
    for (j in seq_along(stations)) {
      
      
      rain_station <- rain_raw %>%
        filter(loc_id == stations[j]) %>%
        mutate(timestamp = timestamp + minutes(tcor_tip[j]))
      
      rainstart_int <- rain_station %>%
        filter(!is.na(rain_tot)) %>%
        filter(rain_tot > 0) %>%
        summarise(rainstart = min(timestamp))
      dat <- rain_station %>%
        mutate(mins = (as.numeric(timestamp) - as.numeric(rainstart_int$rainstart))/60,
               rain_int = na.approx(rain_tot, rule = 2),
               tip_min = if_else(rain_int > 0.0, mins, NaN)) %>%
        filter(!is.na(tip_min)) %>%
        mutate(P_i = (rain_int * 60) / (tip_min - lag(tip_min)),
               P_i = if_else(is.na(P_i), rain_int * 60, P_i)) %>%
        select(timestamp, P_i)
      
      df_rain[[j]] <- rain_station %>%
        full_join(dat, by = "timestamp") %>%
        mutate(Pi_tip = if_else(timestamp > rainstart_int$rainstart, 
                                P_i, 0),
               mins = (as.numeric(timestamp) - as.numeric(rainstart$rainstart))/60) %>%
        arrange(timestamp) %>%
        mutate(P_i = na.locf0(Pi_tip, fromLast = TRUE),
               P_i = if_else(is.na(P_i), 0, P_i),
               P_i = round(P_i, digits = 2)) %>%
        arrange(mins) %>%
        filter(mins >= 0) %>%
        select(mins, P_i) %>%
        mutate(loc_id = stations[j])
    }
    if (spatial == 2) {
      rain <- bind_rows(df_rain) %>%
        pivot_wider(names_from = loc_id, values_from = P_i) %>%
        filter(if_all(everything(), ~ !is.na(.)))
    }
    if(spatial == 1) {
      rain <- bind_rows(df_rain) %>%
        filter(!is.na(P_i)) %>%
        group_by(mins) %>%
        summarise(P_i = round(mean(P_i), digits = 2)) %>%
        filter(cumsum(P_i) != 0) %>%
        mutate(mins = mins - (min(mins) - 1))
    }
  }
  
  if (spatial == 2) {
    #make table for id maps
    loc_station <- tibble(loc_id = names(rain)[-1]) %>%
      mutate(station_nr = seq(1, nrow(.), 1)) %>%
      rename_with(~as.character(seq(0,length(.)-1,1)))
    
    write_delim(loc_station, paste0(main_dirs[i],"rain_id.tbl"), col_names = T)
  }
  # check if event directories exist, otherwise make
  if (!dir.exists(str_c("lisem_runs/", year(event), "/", 
                        str_remove_all(as.character(date(event)), "-")))) {
    dir.create(str_c("lisem_runs/", year(event), "/", 
                     str_remove_all(as.character(date(event)), "-")))
  }
  # location to store rain.txt
  store <- str_c("lisem_runs/", year(event), "/", 
                 str_remove_all(as.character(date(event)), "-") , "/rain.txt")
  write_delim(rain, store, col_names = F)
  #adjust header for correct input
  colnr <- as.character(ncol(rain))
  station_names <- paste(paste("station ", seq(1, ncol(rain) - 1, 1)), collapse = " \n")
  
  header <- readLines("sources/tcs_files/rain_header.txt") %>%
    str_replace("<<date>>", paste0(as.character(event), ", Rainfall source = ",
                                   source_rain)) %>%
    str_replace("<<cols>>", colnr) %>%
    str_replace("<<station_names>>", station_names)
  writeLines(header, store)
  write.table(rain, store, append = T, col.names = F, row.names = F)
  
  ## make the correct id.map based on the rainfall files
  work_dir <- main_dirs[i]
  if (spatial == 1) {
    script <- "rain_id_one.mod"
    pcr_script(script = script, work_dir = work_dir, script_dir = script_dir)
  }
  if (spatial == 2) {
    rain_maps_tif <- c("KNMI_rain_area", "Tipping_rain_area")
    ras <- raster(paste0("spatial_data_SL/", rain_maps_tif[rain_source[i]], ".tif"))
    store <- paste0(main_dirs[i], rain_maps_tif[rain_source[i]], ".asc")
    writeRaster(ras, store, format = "ascii", NAflag = -9999, overwrite = TRUE)
    # convert to pcraster map
    map_in <- paste0(rain_maps_tif[rain_source[i]], ".asc")
    map_tmp <- "tmp.map"
    map_out <- paste0("id_raw.map")
    subdir <- main_dirs[i]
    asc2map(map_in = map_in, map_out = map_tmp,
            sub_dir = subdir)
    resample_location(map_in = "tmp.map", map_out = map_out, dir = subdir)
    file.remove(paste0(subdir, map_tmp))
    file.remove(paste0(subdir, map_in))
    # execute pcr script
    script <- "rain_id_maps.mod"
    pcr_script(script = script, work_dir = work_dir, script_dir = script_dir)
    file.remove(paste0(subdir, map_out))
  }
}

} # end function rainfall input


## Calibration input -------------------------------------------------------
sed_water_cal_input_lisem <- function(calibration_file) {
  
### Depletion model theta I --------------------------------------------------

# daily temp
knmi_t <- read_csv("ext_data/knmi_daily.txt", skip = 15) %>%
  select(YYYYMMDD, TG) %>%
  rename(date = YYYYMMDD, temp = TG) %>%
  mutate(date = ymd(date),
         temp = temp / 10) %>%
  full_join(knmi_day, by = "date") %>%
  mutate(APImax = 30 * 10 * 0.44,
         APImin = 30 * 10 * 0.025,
         rn = row_number(),
         API = if_else(rn == 1, APImax, NaN),
         k = 0.87 + 0.012*(20 - temp),
         k = if_else(k > 1, 1, k))
API_vals <- knmi_t

for (i in 2:nrow(knmi_t)) {
  API_vals <- API_vals %>%
    mutate(API = if_else(!is.na(API), API, lag(API) * k + P),
           API = if_else(API > APImax, APImax, API),
           API = if_else(API < APImin, APImin, API))
}

API_vals <- API_vals %>%
  mutate(theta_i = API / 300)

# fixed crops variables
lu_input <- read_csv(paste0(tcs_dir, "vars_landuse.csv")) %>%
  select(-rr, -LAI)
main_dirs <- c("lisem_runs/2019/20190528/", "lisem_runs/2020/20200816/")

# values to load calibration files per event
skip_events <- which(grepl("^\\[", readLines(calibration_file)))

hom_tbl <- read_csv(paste0(tcs_dir, "homogeneous_vars.txt"), skip = 2)

# and make all pcraster maps
script <- "dynamic_input_maps.mod"
script_dir <- "sources/pcr_scripts"

lu_settings <- vector("list", length = length(event_sel))

for (i in seq_along(event_sel)) {
  # load data
  crop_images <- read_csv(paste0(tcs_dir, "crop_images_", seasons[i], ".txt")) %>%
    mutate(crop_type = if_else(crop_type == "wheatw", "cereals", crop_type))
  fixed <- read_delim(paste0(main_dirs[i],"fixed.tbl"), col_names = T) %>%
    rename_with(~c("lu_nr", "LAI", "Smax", "ch", "cover", "rr", "aggr"))
  
  #calibration file
  cal_file <- read_csv(calibration_file, skip = skip_events[i], n_max = 5)
  
  #cal single values
  cal_vars <- read_csv(calibration_file, skip = skip_events[i]+7, n_max = 3)
  
  coh_bare = filter(hom_tbl, var == "coh_bare")$val * 
    filter(cal_vars, var == "coh_bare_cal")$val # cohesion of bare soil
  coh_wheels = filter(hom_tbl, var == "coh_wheels")$val * 
    filter(cal_vars, var == "coh_wheels_cal")$val # cohesion compacted wheel tracks
  d50_val = filter(hom_tbl, var == "d50")$val * 
    filter(cal_vars, var == "d50_cal")$val
  
  lu_tmp <- lu_input %>%
     left_join(crop_images, by = "crop_type") %>%
     select(-cover, -ch)
  
  # combine the crop parameters with the landuse classes to make a PCR table
  layers <- paste0("landuse_", year(event_sel))
  lu <- st_read(dsn = "spatial_data_SL/fields.gpkg", layer = layers[i])
  lu_crop_type <- read_csv(paste0(tcs_dir, "lu_", seasons[i], "_croptype.csv"))
  
  lu_out <- lu %>%
    st_drop_geometry() %>%
    select(lu_nr) %>%
    left_join(lu_crop_type, by = "lu_nr") %>%
    left_join(lu_tmp, by = "crop_type") %>%
    left_join(fixed, by = "lu_nr") %>%
    filter(!is.na(crop_type)) %>%
    # initial values
    mutate(n = if_else(n == -77, rr/100 + nresi + nveg*cover, n),
           thetai1 = if_else(thetai1 == -77, 
                     filter(API_vals, date == event_sel[i])$theta_i, thetai1)) %>%
    # calibration
    left_join(cal_file, by = "crop_type") %>%
    mutate(n = if_else(!is.na(n_cal), n * n_cal, n),
           thetai1 = if_else(!is.na(thetai1_cal), thetai1 * thetai1_cal, thetai1),
           ksat1 = if_else(!is.na(ksat1_cal), ksat1 * ksat1_cal, ksat1),
           psi1 = if_else(!is.na(psi1_cal), psi1 * psi1_cal, psi1),
           cohadd = if_else(!is.na(cohadd_cal), cohadd * cohadd_cal, cohadd),
           coh = coh_bare + cohadd * cover,
           coh = if_else(cohadd == -1, -1, coh),
           d50 = d50_val)
  # table for paper
  lu_settings[[i]] <- lu_out %>%
    select(crop_type, ksat1, n, thetai1, psi1, coh, d50, rr) %>%
    distinct() %>%
    filter(crop_type %in% cal_file$crop_type)
  lu_wheel <- filter(lu_settings[[i]], crop_type == "potato")
  # table for OpenLISEM
  lu_out <- lu_out %>%
    select(lu_nr, ksat1, n, thetai1, psi1, coh, d50) %>%
    rename_with(~as.character(seq(0,length(.)-1,1))) %>%
    mutate(across(everything(), ~round(., digits = 5))) %>%
    distinct()
  
  wheeltrack_out <- lu_input %>%
    filter(crop_type == "wheeltracks") %>%
    select(crop_type, ksat1, psi1) %>%
    mutate(coh = coh_wheels,
           ksat1 = ksat1 * filter(cal_file, crop_type == "wheeltracks")$ksat1_cal,
           psi1 = psi1 * filter(cal_file, crop_type == "wheeltracks")$psi1_cal,
           n = filter(cal_file, crop_type == "wheeltracks")$n_cal) 
  # table for paper
  lu_settings[[i]] <- bind_rows(lu_settings[[i]], wheeltrack_out) %>%
    mutate(n = if_else(is.na(rr), n * (lu_wheel$rr/100), n),
           d50 = if_else(is.na(d50), lu_wheel$d50, d50),
           thetai1 = if_else(is.na(thetai1), lu_wheel$thetai1, thetai1),
           date = date(event_sel[i])) %>%
    select(-rr) %>%
    arrange(crop_type)
  # table for OpenLISEM
  wheeltrack_out <- wheeltrack_out %>%
    mutate(crop_type = 1) %>%
    rename_with(~as.character(seq(0,length(.)-1,1)))
  
  write_delim(lu_out, paste0(main_dirs[i],"dynamic.tbl"), col_names = T)
  write_delim(wheeltrack_out, paste0(main_dirs[i],"wheels.tbl"), col_names = T)
  
  #execute pcr script
  work_dir <- main_dirs[i]
  pcr_script(script = script, work_dir = work_dir, script_dir = script_dir)
  
}
settings <- bind_rows(lu_settings)
return(settings)

} # end function sed_water_cal_input_lisem()

# Pesticide input lisem ---------------------------------------------------

compound_pest_input <- function(comp_in, pest_cal_file, report = FALSE, 
                                run_in) {

  # which compound?
  compound_sel <- c("Glyphosate", "Metobromuron")
  if (comp_in == 1) {
    compound_in = "No_pesticide"
  } else {
    compound_in <- compound_sel[comp_in - 1]
  }
  #load data
  main_dirs <- c("lisem_runs/2019/20190528/", "lisem_runs/2020/20200816/")
  

if (comp_in > 1) {  
  #' field numbers for three main fields.
fields_large <- tibble(large = c(rep("C", 9), rep("A", 5), rep("B", 7), rep("D", 2), "E", "G"),
                       field_nr = c(1009, 1020, 1013, 2024, 2019, 2023, 3024, 3012, 3022, 
                                    1021, 2017, 2009, 3014, 3015,
                                    1007, 1011, 2003, 2013, 3019, 3016, 3023,
                                    2020, 2021, 2002, 2016))

# soil sample data
soil_samp <- read_csv("data/Soil_samples.csv", lazy = F)

#' Load all pesticide concentrations and application data 
#' Pesticide concentration in W and S of all samples, calculations from raw data
#' are performed in pesticide_calculations.R
lc_all_data <- read_csv("sources/lc_all_data.csv", na = c("-999", "-888", "-777", "NA"))
# Pesticide application data
pest_appl <- read_csv("data/Pest_application.csv") %>%
  left_join(fields_large, by = "field_nr")
# description of applied pesticides - link to active ingredients
pest_desc <- read_csv("data/Pest_description.csv") %>%
  mutate(ai = capwords(ai),
         ai_dose = if_else(ai_dose_unit == "%", ai_dose * 10, ai_dose),
         ai = if_else(ai == "Glyphosate-isopropylammonium", "Glyphosate", ai)) %>%
  rename(compound = ai)

#' difference in molecular weight adjust for that in application rate!!!!!

# chemical characteristics of active ingredients
compound_char <- read_csv("ext_data/AI_characteristics.csv") %>%
  mutate(koc_class = if_else(koc < 15, "very mobile", "mobile"),
         koc_class = if_else(koc > 75 & koc < 500, "moderately mobile", koc_class),
         koc_class = if_else(koc > 500 & koc < 4000, "slightly mobile", koc_class),
         koc_class = if_else(koc > 4000, "non-mobile", koc_class),
         dt50_class = if_else(dt50 < 30, "non-persistent", "moderateley persistent"),
         dt50_class = if_else(dt50 > 100 & dt50 < 365, "persistent", dt50_class),
         dt50_class = if_else(dt50 > 365, "very persistent", dt50_class),
         Sw_class = if_else(Sw < 50, "Low", "Moderate"),
         Sw_class = if_else(Sw > 500, "High", Sw_class))

# analysed AS
analysed <- str_remove(str_subset(names(lc_all_data), "conc"), "conc_")

# remove terbuthylazine from analysis, is not applied or detected, we added it to analysis
# because it is often detected on fields, however does not have a function in this research
analysed_comp <- tibble(compound = analysed) %>%
  filter(compound != "Terbuthylazine")

#' AI analysed
ai_analysed <- analysed_comp %>%
  mutate(analysed = rep(1, nrow(analysed_comp))) %>%
  mutate(compound = str_remove(compound, "conc_")) %>%
  left_join(compound_char, by = "compound") %>%
  select(compound, analysed, type)

## Catchment data --------------------------------
nms <- str_subset(names(lc_all_data), "conc_")
# dates of field campaigns related to events
catch_dates <- ymd(c("2019-05-28", "2020-08-20"))
pest_catch <- lc_all_data %>%
  filter(samp_type == "C")%>%
  left_join(soil_samp, by = "pest_ID") %>%
  left_join(fields_large, by = "field_nr") %>%
  select(pest_ID, field_nr, sub_loc, everything()) %>%
  pivot_longer(starts_with("conc_"),
               names_to = "compound", values_to = "conc") %>%
  mutate(compound = str_remove(compound, "conc_"),
         conc = if_else(conc == -999, 0, conc)) %>%
  semi_join(ai_analysed, by = "compound") %>%
  left_join(compound_char, by = "compound")

# loop to make pest.tbl

#' chemical characteristics for compound based on focus modelling study choices 
#' in the EFSA report (page 72).

# values to load calibration files per event
skip_events <- which(grepl("^\\[", readLines(pest_cal_file)))

# bulkdensity of soil
comp_rho <- 1550

# compound specific settings
comp_name <- compound_in

# and make all pcraster maps
script <- "pesticide_maps.mod"
script_dir <- "sources/pcr_scripts"
}

# loop over events
  pest_settings <- vector("list", length = length(event_sel))

for(i in seq_along(event_sel)) {
  if(comp_in > 1) {
  #load calibration file
  comp_vals <- read_csv(pest_cal_file, skip = skip_events[i], n_max = 2) %>%
    filter(name == comp_name)
  
  conc_cal_rat <- comp_vals$conc_cal_rat
  comp_kfilm <- comp_vals$kfilm
  comp_beta <- comp_vals$ERbeta
  comp_kd <- comp_vals$kd
  comp_ERM <- comp_vals$ERmax
  
  #cal single values
  cal_vars <- read_csv(pest_cal_file, skip = skip_events[i]+4, n_max = 2)
  
  zs_val <- cal_vars$val[1] # pesticide depth in soil layer 1
  zm_val <- cal_vars$val[2] # mixing layer 
  
  # specific choices with samples based on understanding of the study
  # 1 - add sample 'C_B_LC' from field B to field A, this is within the stream
  #   line of field B and the measured values correspond with that instead of
  #   other measurements of field B
  
  # 1. filter for date and compound 
  comp_catch <- pest_catch %>%
    filter(compound == comp_name) %>%
    filter(date == catch_dates[i]) %>%
    distinct(an_name, .keep_all = TRUE) %>%
    group_by(large) %>%
    mutate(large = if_else(str_detect(an_name, "20200820_C_B-LC"), "A", large)) %>%
    # 2. calculate value for every large field
    #' taking the mean is not the only possibility other options are
    #' max or min but also a manual decision for each field based on the
    #' description for the samples
    summarise(conc_field = mean(conc, na.rm = T),
              conc_field = if_else(is.na(conc_field), 0, conc_field),
              conc_field = if_else(large == "C", conc_field / 2, conc_field))
  # the concentration of pesticides on the apple orchard is divided by 2 because
  # the pesticide only are sprayed under the trees and not on the grass sections
  
  #'the field concentration are a combination of dissolved and sorbed pesticide
  #'in the soil, to calculate the separate concentration we use a sample soil mass 
  #'of 2 grams with 30% soil moisture.
  
  # 3. link lu_nr and make pest.tbl
  layers <- paste0("landuse_", year(event_sel))
  pest_lu <- st_read(dsn = "spatial_data_SL/fields.gpkg", layer = layers[i]) %>%
    st_drop_geometry() %>%
    left_join(comp_catch, by = "large") %>%
    select(lu_nr, conc_field, large) %>%
    mutate(conc_dis = if_else(is.na(conc_field), 0, conc_field / (0.3 + 0.7*comp_kd)),
           conc_ads = round(conc_dis * comp_kd / 1000 * conc_cal_rat, digits = 3), # to mg/kg
           conc_dis = round(conc_dis / 1000 * conc_cal_rat, digits = 6), # to mg/L
           zm = zm_val,
           zs = zs_val,
           conc_zs = round(conc_ads / 10, digits = 3)) 
  # table for paper
  pest_settings[[i]] <- pest_lu %>%
    select(large, conc_ads, zm, zs) %>%
    distinct() %>%
    filter(!is.na(large)) %>%
    arrange(large) %>%
    mutate(kfilm = comp_kfilm,
           beta = comp_beta,
           date = date(event_sel[i]),
           compound = comp_name)
  
  pest_lu <- pest_lu %>%
    select(lu_nr, conc_ads, conc_dis, zm, zs, conc_zs) %>%
    rename_with(~as.character(seq(0,length(.)-1,1))) %>%
    distinct()
  
  write_delim(pest_lu, paste0(main_dirs[i], "pest.tbl"), col_names = T)
  
  #make pcr maps
  work_dir <- main_dirs[i]
  pcr_script(script = script, work_dir = work_dir, script_dir = script_dir)
  
  
  }
  #make a table with sample field and total n per event
  comp_tab <- pest_catch %>%
    filter(compound == comp_name) %>%
    filter(date %in% catch_dates) %>%
    distinct(an_name, .keep_all = TRUE) %>%
    group_by(large) %>%
    mutate(large = if_else(str_detect(an_name, "20200820_C_B-LC"), "A", large)) %>%
    filter(!str_detect(name, "^slib potatoe$")) %>%
    ungroup() %>%
    group_by(date) %>%
    mutate(n_tot = n()) %>%
    select(date, n_tot, large) %>%
    distinct(large, .keep_all = TRUE) %>%
    group_by(date) %>%
    arrange(date, large) %>%
    mutate(large = paste0(large, collapse = " "))%>%
    distinct()
  
  if (!dir.exists("results/")) {
    dir.create("results/")
  }
      
  write_csv(comp_tab, "results/table_samples.csv")

# Adjust runfile lisem ----------------------------------------------------

  
  # numerical parameters
  dt_set <- read_csv(run_in, skip = 2)
  dt <- dt_set$val[1]
  
  # check if lisem result directory exists, otherwise make
  if (!dir.exists(paste0(main_dirs[i], "res"))) {
    dir.create(paste0(main_dirs[i], "res"))
  }
  # if it exists, remove and make new
  if (dir.exists(paste0(main_dirs[i], "res"))) {
    unlink(paste0(main_dirs[i], "res"), recursive = TRUE)
    dir.create(paste0(main_dirs[i], "res"))
  }
  
#load template
run_temp <- readLines(paste0(tcs_dir, "manual.run"))
proj_wd <- getwd()

#adjust paths based on system
# home directory
run_temp <- str_replace_all(run_temp, "^Map Directory=<path>", 
                            paste0("Map Directory=", proj_wd, "/", main_dirs[i]))
# result directory
run_temp <- str_replace_all(run_temp, "^Result Directory=<path>", 
                            paste0("Result Directory=", proj_wd, "/", main_dirs[i], "res/"))
# rain files
run_temp <- str_replace_all(run_temp, "^Rainfall Directory=<path>", 
                            paste0("Rainfall Directory=", proj_wd, "/", main_dirs[i]))
run_temp <- str_replace_all(run_temp, "^Rainfall Map Directory=<path>", 
                            paste0("Rainfall Map Directory=", proj_wd, "/", main_dirs[i]))
# set timestep
run_temp <- str_replace_all(run_temp, "Timestep=0", paste0("Timestep=", dt)) # Timestep model 

# adjust data for each run
pest_inc <- 0
if (comp_in > 1) {pest_inc <- 1}
run_temp <- str_replace_all(run_temp, "Include Pesticides=0", paste0("Include Pesticides=", pest_inc))

if (report) {report_pest <- 1}

if (comp_in > 1) {
  report_pest <- 0
  if (report) {report_pest <- 1}  
  run_temp <- str_replace_all(run_temp, "Report Pesticides=0", paste0("Report Pesticides=", report_pest))

#adjust pesticide values
run_temp <- str_replace_all(run_temp, "Pesticide name=foobicide", paste0("Pesticide name=", comp_name)) # compound name
run_temp <- str_replace_all(run_temp, "Kfilm pesticide=0.0", paste0("Kfilm pesticide=", comp_kfilm)) # initial estimate
run_temp <- str_replace_all(run_temp, "Kd pesticide=0.0", paste0("Kd pesticide=", comp_kd)) # 
run_temp <- str_replace_all(run_temp, "ERbeta pesticide=-0.2", paste0("ERbeta pesticide=", comp_beta)) # 
run_temp <- str_replace_all(run_temp, "Rho mixing layer=0.0", paste0("Rho mixing layer=", comp_rho)) # mean bulkdensity
run_temp <- str_replace_all(run_temp, "ERmax pesticide=7.4", paste0("ERmax pesticide=", comp_ERM))

}

# check if a runfile exists, remove
if (length(dir(main_dirs[i], ".run$"))>0) {
  file.remove(dir(main_dirs[i], ".run$", full.names = TRUE))
}

writeLines(run_temp, paste0(main_dirs[i], compound_in, ".run"))
}
  settings <- bind_rows(pest_settings)
  return(settings)  
  
} # end function compound_pest_input()


#..............................................
#................................................

make_runfile_lisem <- function(run_in, name_run = "No_pesticide") {
  # Adjust runfile lisem ----------------------------------------------------
for (i in seq_along(event_sel)) {  
  # numerical parameters
  dt_set <- read_csv(run_in, skip = 2)
  dt <- dt_set$val[1]
  
  # check if lisem result directory exists, otherwise make
  if (!dir.exists(paste0(main_dirs[i], "res"))) {
    dir.create(paste0(main_dirs[i], "res"))
  }
  # if it exists, remove and make new
  if (dir.exists(paste0(main_dirs[i], "res"))) {
    unlink(paste0(main_dirs[i], "res"), recursive = TRUE)
    dir.create(paste0(main_dirs[i], "res"))
  }
  
  #load template
  run_temp <- readLines(paste0(tcs_dir, "manual.run"))
  proj_wd <- getwd()
  
  #adjust paths based on system
  # home directory
  run_temp <- str_replace_all(run_temp, "^Map Directory=<path>", 
                              paste0("Map Directory=", proj_wd, "/", main_dirs[i]))
  # result directory
  run_temp <- str_replace_all(run_temp, "^Result Directory=<path>", 
                              paste0("Result Directory=", proj_wd, "/", main_dirs[i], "res/"))
  # rain files
  run_temp <- str_replace_all(run_temp, "^Rainfall Directory=<path>", 
                              paste0("Rainfall Directory=", proj_wd, "/", main_dirs[i]))
  run_temp <- str_replace_all(run_temp, "^Rainfall Map Directory=<path>", 
                              paste0("Rainfall Map Directory=", proj_wd, "/", main_dirs[i]))
  # set timestep
  run_temp <- str_replace_all(run_temp, "Timestep=0", paste0("Timestep=", dt)) # Timestep model 
  
  # adjust data for each run

  # check if a runfile exists, remove
  if (length(dir(main_dirs[i], ".run$"))>0) {
    file.remove(dir(main_dirs[i], ".run$", full.names = TRUE))
  }
  
  writeLines(run_temp, paste0(main_dirs[i], name_run, ".run"))
}
  
} # end function make_runfile_lisem()


# Sensitivity analysis ----------------------------------------------------


oat_sensitivity_OLP <- function(run_vars, map_names, sens_dir, event) {

# check if sensitivity analysis directory exists, remove and make new.
if (dir.exists(sens_dir)) {
  unlink(sens_dir, recursive = TRUE)
}
dir.create(sens_dir)

year <- year(event)
sens_wd <- paste0(sens_dir, str_remove_all(as.character(date(event)), "-"), "/")
#load input for runs
run_input <- read_delim(run_vars,
                      delim = "\t")

# copy base run maps and runfile to location for sensitivity analysis
if (!dir.exists(sens_wd)) {
  dir.create(sens_wd)
  file.copy(str_c("lisem_runs/", year, "/", 
                  str_remove_all(as.character(date(event)), "-")), 
            str_c("lisem_runs/sensitivity_oat/"),
            recursive = TRUE)
}
# remove results folder and make new
if (dir.exists(paste0(sens_wd, "res"))) {
  unlink(paste0(sens_wd, "res"), recursive = TRUE)
  dir.create(paste0(sens_wd, "res"))
}


# adjust runfile to new folder path names
run_name <- dir(sens_wd, ".run$")[1]
#load template
run_temp <- readLines(paste0(sens_wd, run_name))
proj_wd <- getwd()

#adjust paths based on system
# home directory
run_temp <- str_replace_all(run_temp, "^Map Directory=.*", 
                            paste0("Map Directory=", proj_wd, "/", sens_wd))
# result directory
run_temp <- str_replace_all(run_temp, "^Result Directory=.*", 
                            paste0("Result Directory=", proj_wd, "/", sens_wd, "res/"))
# rain files
run_temp <- str_replace_all(run_temp, "^Rainfall Directory=.*", 
                            paste0("Rainfall Directory=", proj_wd, "/", sens_wd))
run_temp <- str_replace_all(run_temp, "^Rainfall Map Directory=.*", 
                            paste0("Rainfall Map Directory=", proj_wd, "/", sens_wd))
writeLines(run_temp, paste0(sens_wd, "oat_", run_name))

# make base version of all files that will be changed!!
run_base <- run_temp # store the base runfile for the loop
# update runfile name
run_name <- paste0("oat_", run_name)

#make base maps
maps_sens_base <- paste0("base_", maps_names)

file.copy(paste0(sens_wd, maps_names), paste0(sens_wd, maps_sens_base))
#make base rainfile
file.copy(paste0(sens_wd, "rain.txt"), paste0(sens_wd, "rain_base.txt"))


# adjust % or value and select datatype in OLP
sdat <- run_input %>%
  mutate(v_type = if_else(str_detect(val, "%"), "p", "v"),
         val = as.numeric(str_remove_all(val, "\\+|%|v")),
         d_type = if_else(var %in% c("kfilm", "ERbeta", "ERmax", "kd"), "r", "m"),
         d_type = if_else(var == "p", "p", d_type))

var_to_map <- tibble(var = c("conc_ads", "zm", "ksat", "coh", "n", "psi"),
                     map_base = maps_sens_base,
                     map_out = maps_names)

# loop over all runs
# 1 prepare input
# 2 run lisem
# 3 store output and clean res folder

# START LOOP .........................................
res_oat_totals <- vector("list", length = nrow(sdat))
res_oat_hydrograph <- vector("list", length = nrow(sdat))
res_oat_error <- vector("list", length = nrow(sdat))
for (i in seq_along(sdat$run)) {
#1 prepare input
# set everything to default!!
  #rainfall
  file.copy(paste0(sens_wd, "rain_base.txt"), paste0(sens_wd, "rain.txt"), overwrite = TRUE)
  #maps
  file.copy(paste0(sens_wd, maps_sens_base), paste0(sens_wd, maps_names), overwrite = TRUE)
  #runfile
  writeLines(run_base, paste0(sens_wd, run_name))
  
# check if lisem result directory exists, otherwise make
if (!dir.exists(paste0(sens_wd, "res"))) {
  dir.create(paste0(sens_wd, "res"))
}

# adjust maps
if (sdat$d_type[i] == "m") {
  # find correct map
  map_in <- filter(var_to_map, var == sdat$var[i])$map_base
  map_out <- filter(var_to_map, var == sdat$var[i])$map_out
  # pcrcalc for %
  if (sdat$v_type[i] == "p") {
    com <- paste0(map_out, "=", map_in, "+(", map_in, "*", sdat$val[i]/100, ")")
  } else # pcrcalc for v
    { 
    com <- paste0(map_out, "=", map_in, "*0+", sdat$val[i])
  }
  pcrcalc(com, sens_wd)
  
}

  # always adjust DP to PP depending on Kd and PP_conc
  kd_val = 3.5
  if (sdat$var[i] == "kd") {
    kd_val = if_else(sdat$v_type[i] == "p", 3.5 + (sdat$val[i]/100 * 3.5),
                     sdat$val[i])
  }
  com <- paste0("pcmixwat.map=pcmixsoil.map/",kd_val)
  pcrcalc(com, sens_wd)
  
# adjust runfile
if (sdat$d_type[i] == "r") {
  var_in <- capwords(sdat$var[i])
  val_base <- as.numeric(str_remove(na.omit(
    str_extract(run_base, paste0("^", var_in, " pesticide=.*"))),
                          "^.*pesticide="))
  val_in <- if_else(sdat$v_type[i] == "p", val_base + 
                      ((sdat$val[i]/100) * val_base), sdat$val[i])
  run_temp <- run_base
  run_temp <- str_replace_all(run_temp, paste0("^", var_in, " pesticide=.*"), 
                              paste0(var_in, " pesticide=", as.character(val_in)))
  writeLines(run_temp, paste0(sens_wd, run_name))
}
# adjust rainfall file
if (sdat$d_type[i] == "p") {
  header <- readLines(paste0(sens_wd, "rain_base.txt"))[1:5]
  rain_in <- read_delim(paste0(sens_wd, "rain_base.txt"), 
                        delim = " ", skip = 5, col_names = FALSE)
  rain_out <- rain_in %>%
    mutate(X2 = X2 + (X2 * (sdat$val[i] / 100)))
  store <- paste0(sens_wd, "rain.txt")
  writeLines(header, store)
  write.table(rain_out, store, append = T, col.names = F, row.names = F)
}

#2 run lisem
  print(paste0("Run number ", i, " out of ", nrow(sdat)))
# run lisem initial
lisem_dir <- "C:/Users/MC/Werk/OpenLISEM/lisem_bin/Lisem.exe"
run_time <- run_lisem(sens_wd, run_name, lisem_dir, GUI = FALSE)
runs_to_go <- nrow(sdat) - i
pred_end <- Sys.time() + runs_to_go * run_time
print(paste0(runs_to_go, " runs to go. Predicted end time = ", pred_end))

#3 store output

#totals
filetot <- paste0(sens_wd, "res/totals.csv")
res_oat_totals[[i]] <- read_csv(filetot) %>%
  rename_with( ~c("var", "val"))
#hydrograph

file <- paste0(sens_wd, "res/hydrographs_1.csv")
hy_names <- readLines(file)[2] %>%
  str_split(",", simplify = TRUE) %>%
  str_remove_all(" |#")
res_oat_hydrograph[[i]] <- read_csv(file, skip = 2) %>%
  rename_with(~hy_names) %>%
  mutate(mins = round(Time * 24 * 60, digits = 2)) %>%
  distinct()

#errors
file_err_pest <- paste0(sens_wd, "res/error_pest.txt")
err_var <- as.numeric(readLines(file_err_pest)[2])
err_nms <- readLines(file_err_pest)[3:(err_var+2)]
res_oat_error[[i]]  <- read_delim(file_err_pest, skip = err_var+2, col_names = F, delim = " ") %>%
  select(-X1) %>%
  rename_with(~err_nms) %>%
  rename(runstep = 'run step') %>%
  mutate(time = runstep / 60)

# remove results folder and make new
if (dir.exists(paste0(sens_wd, "res"))) {
   unlink(paste0(sens_wd, "res"), recursive = TRUE)
  dir.create(paste0(sens_wd, "res"))
  
# save list to .RDate file, can be used for further analysis later
  save(res_oat_error, res_oat_hydrograph, res_oat_totals, file = "sensitivity.Rdata")
}
} #end sensitivity loop
} # end sensitivity function

# functions to save results in .Rdata list files

store_result_total <- function(work_dir) {
  #totals
  filetot <- paste0(work_dir, "res/totals.csv")
  res <- read_csv(filetot) %>%
    rename_with( ~c("var", "val"))
  return(res)
}
store_result_hydrograph <- function(work_dir) {
  file <- paste0(work_dir, "res/hydrographs_1.csv")
  hy_names <- readLines(file)[2] %>%
    str_split(",", simplify = TRUE) %>%
    str_remove_all(" |#")
  res <- read_csv(file, skip = 2) %>%
    rename_with(~hy_names) %>%
    mutate(mins = round(Time * 24 * 60, digits = 2)) %>%
    distinct()
  return(res)
}
store_result_error <- function(work_dir) {  
  #errors
  file_err_pest <- paste0(work_dir, "res/error_pest.txt")
  err_var <- as.numeric(readLines(file_err_pest)[2])
  err_nms <- readLines(file_err_pest)[3:(err_var+2)]
  res  <- read_delim(file_err_pest, skip = err_var+2, col_names = F, delim = " ") %>%
    select(-X1) %>%
    rename_with(~err_nms) %>%
    rename(runstep = 'run step') %>%
    mutate(time = runstep / 60)
  return(res)
}

# end of dataset functions