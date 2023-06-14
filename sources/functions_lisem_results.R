# analyse the results of lisem runs and compare these with observed data

# Initialization ----------------------------------------------------------

library(lubridate)
library(hms)
library(sf)
library(raster)
library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(extrafont)
library(hydroGOF)

#functions
# discharge functions
source("sources/discharge_calculations.R")

nl_tz <- locale(tz = "Etc/GMT-1")

#function to evaluate model performance and make figure

# choose which results you want to evaluate
# Only runoff = 1
# Runoff and erosion = 2
# Runoff, erosion and pesticides = 3

load_observations <- function(pesticides = NULL) {


#load generic data
ev_dates <- read_csv("sources/ev_dates.csv")

# selected events
event_sel <- ymd(c("2019-05-28", "2020-08-16"))
event_pos <- c(3, 17)
seasons <- c("2019", "2020")
tcs_dir <- "sources/tcs_files/"

# Observations ------------------------------------------------------------

## load data ------------------------------
wb_data <- read_csv("data/WB_data.csv") # data from waterboard (Ha)
cr6 <- read_csv("data/CR6.csv", col_types = "TddddddddddddddD") # from CR6 (Hb)
crest_nap <- 79.331 # m + NAP as measured by the waterboard

## Calculate discharge -----------------------------------------------------

# use Q_calc to calculate discharge for each event and store in list
# the dicharge calculations can be found in the 'sources' folder for more details.
Q_tcor <- Q_calc(ev_dates, wb_data, cr6, crest_nap)
df_dat <- Q_tcor[[1]] # data with discharge for each event
df_tcor <- Q_tcor[[2]] # based on synchronization between different time series a time correction is performed

## Precipitation data ------------------------------------------------------
#' this data is provided as secondary source, it is clipped from the 5 minute radar
#' of KNMU for 2018 - 2020.
knmi_p <- read_csv("ext_data/KNMI_rain.csv") %>%
  group_by(timestamp) %>%
  summarize(mean(P))
names(knmi_p) <- c("timestamp", "P_m_knmi")
knmi_p <- mutate(knmi_p, P_i_knmi = P_m_knmi * 12)
#add rainfall from knmi to event data, also keep rainfall by cr6
name <- str_remove_all(ev_dates$date, "-")
for (j in seq_along(event_sel)) {
  i <- event_pos[j]
  raintmp <- knmi_p %>%
    filter(timestamp >= ev_dates$start_t[i] & timestamp <= ev_dates$end_t[i]) 
  df_dat[[i]] <- left_join(df_dat[[i]], raintmp, by = "timestamp") %>%
    mutate(P_i_cr6 = rain_tot * 60)
           #mins = (as.numeric(timestamp) - as.numeric(min(timestamp)))/60)
  rainstart <- df_dat[[i]] %>%
    filter(!is.na(P_i_knmi)) %>%
    filter(P_i_knmi > 0) %>%
    summarise(rainstart = min(timestamp))
  
  df_dat[[i]]<- df_dat[[i]] %>%
    mutate(mins = (as.numeric(timestamp) - as.numeric(rainstart$rainstart))/60,
           Q_peak_min = mins)
}

Q_data <- bind_rows(df_dat[event_pos])
store <- "lisem_runs/Q_data.csv"
write_csv(Q_data, store)


## Erosion data ------------------------------------------------------------
# add tss estimate for 20200812 (made in pesticide calculations.R)
# which was not available due to loss of samples.
ev_samples <- read_csv("data/Event_samples.csv")
tss_est <- read_csv("sources/tss_20200812.csv", lazy = FALSE) %>%
  left_join(ev_samples, by = "pest_ID") %>%
  select(-sed_conc.y) %>%
  rename(sed_conc = sed_conc.x)

# add tss estimate to the event samples data
ev_samples <- ev_samples %>%
  anti_join(tss_est, by = "pest_ID") %>%
  bind_rows(tss_est)

# make sediment runoff data for each event
df_er <- vector("list", length = length(name))
for (j in seq_along(event_sel)) {
  i <- event_pos[j]
  ev_start <- ev_dates$start_t[i]
  ev_end <- ev_dates$end_t[i]
  df_er[[i]] <- ev_samples %>%
    filter(timestamp > ev_start & timestamp < ev_end) %>%
    select(timestamp, sed_conc) %>%
    mutate(timestamp = timestamp + df_tcor[[i]]) %>%
    group_by(timestamp) %>%
    summarise(sed_conc = mean(sed_conc)) %>%
    left_join(df_dat[[i]], by = "timestamp") %>%
    mutate(load = sed_conc * Q_int) %>%
    select(timestamp, sed_conc, load)
}

data <- bind_rows(df_er[event_pos])
store <- "lisem_runs/erosion_data.csv"
write_csv(data, store)

## Pesticide data ----------------------------------------------------------

#' Load all pesticide concentrations
lc_all_data <- read_csv("sources/lc_all_data.csv", na = c("-999", "-888", "-777", "NA"))

# analysed AS
analysed <- str_remove(str_subset(names(lc_all_data), "conc"), "conc_")

# remove terbuthylazine from analysis, is not applied or detected, we added it to analysis
# because it is often detected on fields, however does not have a function in this research
analysed_comp <- tibble(compound = analysed) %>%
  filter(compound != "Terbuthylazine")



#outlet pesticide data
#' pesticide data for each event
pest_ts <- ev_samples %>%
  mutate(date = date(timestamp),
         sed_gr = if_else(is.na(sed_gr), w_d_cup - w_e_cup, sed_gr)) %>%
  select(date, pest_ID, sed_gr, timestamp) %>%
  group_by(timestamp) %>%
  summarise_all(mean, na.rm = T) %>%
  group_by(pest_ID) %>%
  summarise_all(mean, na.rm = T) %>%
  distinct() %>%
  filter(!is.na(pest_ID))

lc_all <- lc_all_data %>%
  left_join(pest_ts, by = "pest_ID") %>%
  filter(!is.na(date)) %>%
  select(-an_name, -week, -conc_Terbuthylazine) %>%
  replace(is.na(.), 0) %>%
  select(date, timestamp, sed_gr, samp_type, pest_ID, everything())

df_data <- vector("list", length = length(pesticides))
for (k in seq_along(pesticides)) {
  # select compounds 
  comp <- pesticides[k] 

#' mass transported for one specific compound in water and sed.
loq = if_else(comp == "Glyphosate", 0.01, 0.0025)

df_pest <- vector("list", length = length(name))
for (j in seq_along(event_sel)) {
  i <- event_pos[j]
  ev_start <- ev_dates$start_t[i]
  ev_end <- ev_dates$end_t[i]
  df_pest[[i]] <- lc_all %>%
    filter(timestamp > ev_start & timestamp < ev_end) %>%
    mutate(timestamp = timestamp + df_tcor[[i]])
if (nrow(df_pest[[i]]) > 0) {  
  nms <- names(df_pest[[i]])
  #oob_lim <-  5 * df_totals[[i]]$undisc_w[1] #limit to end barplots, so the lower values are still visible
 # data with timestamp to compare with model output
  Q_tot <- df_dat[[i]] %>%
    arrange(timestamp) %>%
    mutate(secs = as.numeric(lead(timestamp) - timestamp)*60,
           Q_m3 = Q_int * secs) %>%
    left_join(df_er[[i]], by = "timestamp") %>%
    arrange(timestamp) %>%
    mutate(sed_int = na.approx(sed_conc, maxgap = 10),
           Sed_kg = sed_int * Q_m3) %>%
    mutate(Q_m3 = if_else(Q_int > 0.012, Q_m3, 0),
           Q_duration = if_else(Q_int > 0.012, 1, 0)) 
 
  df_pest[[i]] <- pivot_longer(df_pest[[i]], nms[6:length(nms)],
                               names_to = "compound", values_to = "conc") %>%
     mutate(compound = str_replace(compound, "conc_", "")) %>%
     pivot_wider(values_from = "conc", names_from = "samp_type", names_prefix = "conc_") %>%
     select(-date) %>%
     filter(compound == comp) %>%
     full_join(Q_tot, by = "timestamp") %>%
     arrange(timestamp) %>%
     mutate(conc_S = conc_S / 1000, # change to mg/kg
            conc_W = conc_W / 1000, # mg/L
            PP_mg = conc_S * Sed_kg, # mg
            DP_mg = conc_W * Q_m3 * 1000,   # mg
            DP_loq = if_else(!is.na(conc_W), loq * Q_m3 * 1000, NaN), #mg 
            compound = comp) 
}
}
df_data[[k]] <- bind_rows(df_pest[event_pos])
}
data <- bind_rows(df_data)

store <- "lisem_runs/pesticide_data.csv"
write_csv(data, store)

# store table with samples at outlet and combine with fields for methods
tab <- read_csv("results/table_samples.csv") %>%
  mutate(date = if_else(date == "2020-08-20", ymd("2020-08-16"), date)) 
dat <- data %>%
  filter(!is.na(sed_conc)) %>%
  mutate(date = date(timestamp)) %>%
  group_by(date) %>%
  summarise(ss_n = n(),
            p_n = sum(!is.na(PP_mg))) %>%
  full_join(tab, by = "date")

write_csv(dat, "results/table_samples.csv")

} # end of function load_observations()

#.....................................................................
#.....................................................................
#' This function handels all possible output from the OpenLISEM results folder
#' (MC) maybe it is not very clear any more with all the if statement, sorry...
#' The function has the following input options:
#' figure = (1-5)
#'    1 = initial sed + water; only table, 
#'    2 = sed + water with figure draft version
#'    3 = pesticides + figure 6 in manuscript
#'    4 = sensitivity only table
#'    5 = erosion + runoff, figure 3 in manuscript
#' fig_name = ""; fill anything you like to add to the name. besides that as
#'                default the event date and pesticide type will be added.
#' events: expects a vector with dates for evaluated events.
#' basedir: the location where the '/res' folder of the OpenLISEM runs is 
#'          located
#' list_pos = NULL: if instead from OpenLISEM output the results are loaded
#'                  from a .Rdata file (e.g. sensitivity analysis), this value
#'                  indicates the position in the list file to read.
#' legend = TRUE: default a legend will be added to figures, for some combined
#'                figure it is better to not repeat the legend every time.


figures_results_paper <- function(figure = 1, fig_name = "",
                                  events, base_dir, list_pos = NULL, 
                                  legend = TRUE) {
  # evaluate
  # 1 = water, 2 = 1 + sed, 3 = 2 + pesticides, 4 = 3 with specific output
  
  figure_eval <- c(2, 2, 3, 4, 5)
  evaluate <- figure_eval[figure]
  
  if (figure == 4) {figure <- 1}
  
  if (fig_name != "")  {
    fig_name <- paste0(fig_name, "_")
  }
  #load generic data
  ev_dates <- read_csv("sources/ev_dates.csv")
  
  # selected events
  event_sel <- events 
  tcs_dir <- "sources/tcs_files/"
  
  # select compounds in lisem run
  if (figure == 5) {
    load("results/pest_man.Rdata")
    hy_list <- res_hydrograph[[list_pos]]
    totals_list <- res_totals[[list_pos]]
    error_list <- res_error[[list_pos]]
    comp <- hy_list$compound[1]
  } else {
  comp <- dir(base_dir[1], "^[A-Z].*run$") %>%
    str_remove(".run")
  }
  leg <- legend
  #if evaluate is 4 - sensitivity analysis. this function is part of a loop
  # select the correct event based on the 'list_pos'
  if (evaluate == 4) {
    load("sensitivity.Rdata")
    hy_list <- res_oat_hydrograph[[list_pos]]
    totals_list <- res_oat_totals[[list_pos]]
    error_list <- res_oat_error[[list_pos]]
  }
  
  if (comp == "No_pesticide" & evaluate == 3) {evaluate == 2}
  # load data ------------------
  Q_dat <-  read_csv("lisem_runs/Q_data.csv")
  if (evaluate > 1) {
    er_dat <- read_csv("lisem_runs/erosion_data.csv")
  }
  if (evaluate > 2) {
    pest_dat <- read_csv("lisem_runs/pesticide_data.csv") %>%
      filter(compound == comp)
  }
  if (figure > 1) {
  # Define different color pallets for figures
  # 2 colors for water and sediment OR DP and PP
  colors1 <- c("#E6BB00", "#0072B2")
  
  colors2 <- c("#E6BB00", "#0072B2", "#CC0000", "#8df55b")
  
  # 12 color version
  colors3 <- c("#7465AC", "#E6BB00", "#56B4E9",
               "#50ECCC", "#009E73", "#F0E442", "#885522",
               "#0072B2", "#CC0000", "#8df55b", "#BFC816", "#999999")
  
  my_theme = theme(
    text = element_text(size = 10, family = "Times New Roman"),
    axis.title.x = element_text(size = 10, family = "Times New Roman"),
    axis.text.x = element_text(size = 10, family = "Times New Roman"),
    axis.title.y = element_text(size = 10, family = "Times New Roman"),
    plot.margin = margin(6, 6, 6, 0)
  )
  }
  # Model results -----------------------------------------------------------
  
  margin <- c(0,0)
  xlimits <- c(100, 170) # x axis limits per event
  mod_totals <- vector("list", length = length(event_sel))
  mb_error <- vector("list", length = length(event_sel))
  plot_list <- vector("list", length = length(event_sel))
  plot_rain <- vector("list", length = length(event_sel))
  comb_plot <- vector("list", length = length(event_sel))

  for (j in seq_along(event_sel)) {
    rain_file <- paste0(base_dir[j], "/rain.txt")
    rain_type <- readLines(rain_file)[1] %>%
      str_remove("^.*, ")
    skipnr <- as.numeric(readLines(rain_file)[2]) + 2
    rain <- read_delim(rain_file, skip = skipnr, col_names = FALSE, delim = " ") %>%
      rename(mins = X1) %>%
      pivot_longer(cols = -mins, names_to = "station", values_to = "p") %>%
      group_by(mins) %>%
      summarise(rain = round(mean(p), digits = 2))
    # select rain source for correct timing of results
    rain_source <- if(str_detect(rain_type, "KNMI")) {
      "P_i_knmi" } else { "rain_tot"}
    
    # collect relevant data from model results
    if (evaluate < 4) {
    file <- paste0(base_dir[j], "/res/hydrographs.csv")
    
    #load model results timeseries
    hy_names <- readLines(file)[2] %>%
      str_split(",", simplify = TRUE) %>%
      str_remove_all(" |#")
    hydrograph <- read_csv(file, skip = 2) %>%
      rename_with(~hy_names) %>%
      mutate(mins = round(Time * 24 * 60, digits = 2)) %>%
      distinct()
    } else {
      hydrograph <- hy_list
    }
    if (evaluate > 1) {
      hydrograph <- hydrograph %>%
        rename_with(~str_replace(., " #", "_"))
    }
    if (evaluate > 2) {
      hydrograph <- hydrograph %>%
        mutate(PCw = PQw / Qall,
               PCs = PQs / Qsall)
    }
    
    # load model totals
    if (evaluate < 4) {
    filetot <- paste0(base_dir[j], "/res/totals.csv")
    
    totals_mod <- read_csv(filetot) %>%
      rename_with( ~ c("var", "val"))
    } else { 
      totals_mod <- totals_list
    }
    event <- Q_dat %>%
      filter(date(timestamp) == event_sel[j])
    
    rainstart <- event %>%
      rename(rain_min = all_of(rain_source)) %>%
      filter(!is.na(rain_min)) %>%
      filter(rain_min > 0) %>%
      summarise(rainstart = min(timestamp))
    
    # load observations for this event
    event <- Q_dat %>%
      filter(date(timestamp) == event_sel[j]) %>%
      mutate(mins = (as.numeric(timestamp) - as.numeric(rainstart$rainstart))/60) %>%
      left_join(hydrograph, by = "mins") %>%
      left_join(rain, by = "mins") %>%
      mutate(Q = Q * 1000)
    
    if (evaluate > 1) {
      event <- event %>%
        left_join(er_dat, by = "timestamp")
    }
    if (evaluate > 2) {
      event <- pest_dat %>%
        filter(date(timestamp) == event_sel[j]) %>%
        mutate(mins = (as.numeric(timestamp) - as.numeric(rainstart$rainstart))/60) %>%
        left_join(hydrograph, by = "mins") %>%
        left_join(rain, by = "mins") %>%
        mutate(Q = Q * 1000,
               PP_load = PP_mg / 60,
               DP_load = DP_mg / 60)
      
    }
    
    ## Evaluation criteria --------------
    
    #mean observed uncertainty
    obs_err <- event %>%
      mutate(Knew = if_else(!is.na(K_pred), K_pred, K),
             rel_err = if_else(Knew < 0.7, 0.03, 0.05),
             obs_un = rel_err * Q) %>%
      summarise(obs_un_wat = mean(obs_un, na.rm = T))
    # totals model & evaluation criteria
    mod_totals[[j]] <- tibble(Qmax_sim = filter(totals_mod, str_detect(var, "Peak discharge"))$val,
                              Qtot_sim = filter(totals_mod, str_detect(var, "Total outflow \\(overland"))$val,
                              Ptot_sim = filter(totals_mod, str_detect(var, "Total Precipitation"))$val,
                              Q_peak_min_sim = filter(totals_mod, str_detect(var, "Peak time discharge"))$val,
                              NSE_wat = NSE(event$Q, event$Qall),
                              KGE_wat = KGE(event$Q, event$Qall),
                              RMSE_wat = rmse(event$Q, event$Qall),
                              obs_un_wat = obs_err$obs_un_wat,
                              date = event_sel[j]) %>%
      mutate(across(1:8, as.numeric)) %>%
      mutate(across(1:8, round, 3)) %>%
      select(NSE_wat, KGE_wat, RMSE_wat, obs_un_wat)
    if (evaluate > 1) {
      # totals model
      mod_totals[[j]] <- mod_totals[[j]] %>%
        mutate(NSE_sed = NSE(event$load, event$Qsall),
               KGE_sed = KGE(event$load, event$Qsall),
               RMSE_sed = rmse(event$load, event$Qsall),
               obs_un_sed = mean(event$load * 0.06, na.rm = T)) %>%
        mutate(across(5:8, as.numeric)) %>%
        mutate(across(5:8, round, 3))
    }
    if (evaluate > 2) {
      # totals model
      mod_totals[[j]] <- mod_totals[[j]] %>%
        mutate(RMSEs_pest = rmse(event$PP_load, event$PQs),
               RMSEw_pest = rmse(event$DP_load, event$PQw),
               NSEs_pest = NSE(event$PP_load, event$PQs),
               NSEw_pest = NSE(event$DP_load, event$PQw),
               KGEs_pest = KGE(event$PP_load, event$PQs),
               KGEw_pest = KGE(event$DP_load, event$PQw),
               obs_un_ps = mean(event$PP_load * 0.18, na.rm = T),
               obs_un_pw = mean(event$DP_load * 0.17, na.rm = T)) %>%
        mutate(across(9:16, as.numeric)) %>%
        mutate(across(9:16, round, 3)) %>%
        select(KGEs_pest, KGEw_pest, NSEs_pest, NSEw_pest, RMSEs_pest, RMSEw_pest, obs_un_ps, obs_un_pw, KGE_wat, KGE_sed)
      
      #mass balance errors
      if (evaluate < 4) {
      file_err_pest <- paste0(base_dir[j], "/res/error_pest.txt")
      err_var <- as.numeric(readLines(file_err_pest)[2])
      err_nms <- readLines(file_err_pest)[3:(err_var+2)]
      error <- read_delim(file_err_pest, skip = err_var+2, col_names = F, delim = " ") %>%
        select(-X1) %>%
        rename_with(~err_nms) %>%
        rename(runstep = 'run step') %>%
        mutate(time = runstep / 60)
      } else {
        error <- error_list
      }
      mb_error[[j]] <- error[nrow(error), ] %>%
        select(-runtime, -time, -runstep, -WMerr, -SMerr) %>%
        mutate(compound = comp)
      mod_totals[[j]] <- bind_cols(mod_totals[[j]], mb_error[[j]])
      if (evaluate == 4 | evaluate == 3) {
        # also add total transport to model
        mod_totals[[j]] <- mod_totals[[j]] %>%
          mutate(DP_tot_sim = as.numeric(filter(totals_mod, str_detect(var, "Total dissolved pesticide transport"))$val),
                 PP_mod_sim = as.numeric(filter(totals_mod, str_detect(var, "Total particulate pesticide transport"))$val))
      }
    }
    
    # Figures -----------------------------------------------------------------
    if (figure > 1) {
    # plot with rainfall data for the event
    P_scale_min <- 30
    P_scale <- if_else(max(event$rain, na.rm = T) < P_scale_min, P_scale_min, 1.1 * max(event$rain, na.rm = T))
    plot_rain[[j]] <- ggplot() +
      geom_bar(data = event, aes(x = mins, y = rain), stat = "identity", width = 4) +
      scale_y_reverse(limits = c(P_scale, 0), expand = margin) +
      theme_classic() +
      scale_x_continuous(position = "top", limits = c(0, xlimits[j])) +
      labs(x = "", y = bquote(P~ mm~ h^-1),
           title = event_sel[j], parse = T) +
      my_theme
    # observed and modeled runoff
    Q_scale_min <- 200
    Q_scale <- if_else(max(event$Q_int, na.rm = T) * 1000 < Q_scale_min, Q_scale_min, 1.2 * max(event$Q, na.rm = T))
    if (evaluate == 1){
      plot_list[[j]] <- ggplot() + 
        geom_point(data = event, aes(mins, Q, color = "Q_obs"), shape = 19) +
        geom_line(data = hydrograph, aes(mins, ls, color = "Q_mod"), linetype = "dashed") +
        theme_classic() +
        scale_color_manual(values = c("Q_obs" = colors1[1], "Q_mod" = colors1[1]),
                           guide = guide_legend(override.aes = list(
                             linetype = c("dashed", "blank"),
                             shape = c(NA, 19)
                           )),
                           name = " ") +
        labs(x = "Time (min)", y = "Q (l sec^-1)", parse = TRUE) + 
        scale_x_continuous(limits = c(0, xlimits[j])) + 
        scale_y_continuous(limits = c(0, Q_scale)) +
        theme(legend.position = c(0.9, 0.9))
      
      tab_eval <- ggplot() +
        theme_bw() +
        annotation_custom(tableGrob(mod_totals[[j]], 
                                    theme = ttheme_minimal(base_size = 10)))
    }
    
    # observed and modeled runoff and erosion
    if (evaluate > 1) {
      
      coeff <- 0.1
      p <- ggplot() +
        geom_point(data = event, aes(mins, Q, color = "Q_obs"), shape = 19) +
        geom_point(data = event, aes(mins, load / coeff, color = "S_obs"), shape = 19) +
        geom_line(data = hydrograph, aes(mins, Qall, color = "Q_sim"), linetype = "dashed") +
        geom_line(data = hydrograph, aes(mins, Qsall / coeff, color = "S_sim"), linetype = "dashed") +
        scale_y_continuous(sec.axis = sec_axis(~.*coeff, name = bquote(Sediment~ load~ (kg~ sec^-1))), limits = c(0, Q_scale),
                           expand = margin) +
        scale_color_manual(values = c("Q_obs" = colors1[2], "S_obs" = colors1[1], "Q_sim" = colors1[2], "S_sim" = colors1[1]),
                           guide = guide_legend(override.aes = list(
                             linetype = c("blank", "dashed", "blank", "dashed"),
                             shape = c(19, NA, 19, NA)
                           )),
                           name = " ") +
        labs(x = "Time (min)", y = bquote(Q~ (l~ sec^-1)), 
             parse = TRUE) +
        theme_classic() +
        xlim(c(0, xlimits[j])) + 
        my_theme
      
      if (figure == 2 & j == 1) {
        plot_list[[j]] <- p +
        theme(legend.position = "none")
      } else { 
        plot_list[[j]] <- p +
        theme(legend.position = c(0.8, 0.9))
      }
      
      tab_eval <- ggplot() +
        theme_bw() +
        annotation_custom(tableGrob(mod_totals[[j]], 
                                    theme = ttheme_minimal(base_size = 10)))
    }
    
    # observed and modeled pesticide transport
    if (evaluate == 3 | evaluate == 5) {
      pest_load_plot <-  ggplot() +
        geom_point(data = event, aes(mins, PP_load, color = "PP_obs"), shape = 19) +
        geom_point(data = event, aes(mins, DP_mg/60, color = "DP_obs"), shape = 19) +
        geom_point(data = event, aes(mins, DP_loq/60, color = "LOQ_DP"), shape = 24) +
        geom_line(data = hydrograph, aes(mins, PQs, color = "PP_sim"), linetype = "dashed") +
        geom_line(data = hydrograph, aes(mins, PQw, color = "DP_sim"), linetype = "dashed") +
        scale_color_manual(values = c("DP_obs" = colors1[2], "PP_obs" = colors1[1], "DP_sim" = colors1[2], "PP_sim" = colors1[1], "LOQ_DP" = colors1[2]),
                           guide = guide_legend(override.aes = list(
                             linetype = c("blank", "dashed", "blank", "blank", "dashed"),
                             shape = c(19, NA, 24, 19, NA)
                           )),
                           name = " ") +
        
        theme_classic() +
        xlim(c(0,xlimits[j])) +
        labs(x = "Time (min)", y = bquote(Load~ (mg~sec^-1))) +
        my_theme
      if (leg == TRUE) {
      pest_load_plot <- pest_load_plot + theme(legend.position = c(0.8,0.6))
      } else {
          pest_load_plot <- pest_load_plot + theme(legend.position = "none")
          }
      
      pest_conc_plot <-  ggplot() +
        geom_point(data = event, aes(mins, conc_S), color = colors1[1]) +
        geom_point(data = event, aes(mins, conc_W), color = colors1[2]) +
        geom_line(data = hydrograph, aes(mins, PCs), color = colors1[1], linetype = "dashed") +
        geom_line(data = hydrograph, aes(mins, PCw), color = colors1[2], linetype = "dashed") +
        theme_classic() +
        xlim(c(0,xlimits[j])) +
        ylim(c(0, 1.5*max(event$conc_S, na.rm = T))) +
        labs(x = "Time (min)", y = bquote(Concentration~ (mg~ kg^-1))) +
        my_theme
      pest_plot <- plot_grid(pest_load_plot, pest_conc_plot, rel_widths = c(1, 1), 
                             ncol = 2, align = "h", axis = "rl")
      
      tab_eval <- ggplot() +
        theme_bw() +
        annotation_custom(tableGrob(mod_totals[[j]], 
                                    theme = ttheme_minimal(base_size = 10)))
      
      tab_mb <- ggplot() +
        theme_bw() +
        annotation_custom(tableGrob(mb_error[[j]],
                                    theme = ttheme_minimal(base_size = 10)))
    }
    
    # 1 rain and water
    # 2 rain and erosion/water
    comb_plot[[j]] <- plot_grid(plot_rain[[j]], plot_list[[j]], 
                                rel_heights = c(1,2) , nrow = 2, 
                                align = "v", axis = "rl")
    
    # 3 all + eval table + err table
    if (evaluate == 3) {
      # figures load / conc 7 x 5 (2x)
      # table error 15 x 4
      title_plot <-ggdraw() +
        draw_label(str_c("Observed and modelled, water & sediment dynamics \n and ", comp, " transport on ", date(event$timestamp[1])),
                   x = 0, hjust = 0)
      
      plot <- plot_grid(title_plot, comb_plot[[j]], tab_eval, pest_plot, tab_mb, 
                        rel_heights = c(0.5, 5, 1, 3, 1) , nrow = 5, ncol = 1, 
                        align = "hv", axis = "rl") +
        theme(plot.background = element_rect(fill = "white"))
      ggsave(plot = plot, filename = paste0("images/lisem_", fig_name, comp, "_", 
                                            str_remove_all(as.character(date(event_sel[j])), "-"), ".tiff"), 
             width = 17, height = 24, units = "cm", dpi = 600, device = "tiff")
      
    }
   
    }
    if (figure ==1) {fig_name <- "initial_"}
    
    
  } # end figure per event loop
  if (figure == 2) {
  # combine comb_plot to 1 figure for paper
  plot <- plot_grid(comb_plot[[1]], comb_plot[[2]], rel_widths = c(1.5, 2), align = "hv", axis = "rl",
                    ncol = 2)#,
                    # labels = c("28-05-2019", "16-08-2020"), label_fontfamily = "Times New Roman",
                    # label_x = 0.6, label_y = 0.8, label_size = 11)
  ggsave(plot = plot, filename = paste0("images/figure3.tif"),
         width = 20, height = 10, units = "cm", dpi = 600, device = "tiff")
  }
  performance <- bind_rows(mod_totals)
  if (evaluate == 5) {
    # make figure for paper
    performance <- pest_plot
  }
  return(performance)
  # write_csv(performance, paste0("lisem_runs/performance_", fig_name, comp, ".csv"))
}

