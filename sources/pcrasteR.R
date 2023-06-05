# functions to integrate PCraster with R

# set_pcraster

# miniconda: set path to your miniconda installation
# env: give the name of the conda environment where pcraster is installed
set_pcraster <- function(miniconda = "~/miniconda3", env = "pcraster") {
  sys_type <- Sys.info()['sysname']
  lib <- ifelse(sys_type == "Windows", "/Library", "")
  dir <- paste0(miniconda, "/envs/", env, lib, "/bin/")
  assign("pcr_dir", dir,envir=parent.env(environment()))
  message <- paste("PCraster directory is set to:", dir)
  return(message)
}

# asc2map

# this function assumes 'set_pcraster()' is already used.
# for explanation of function see pcraster documentation
# https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/app_asc2map.html
#clone: name of clone map
#asc_header: (TRUE) if true the option '-a' will be passed which sets all other options from the ASCII header
#options: a string with all other options separated by " ". With asc_header = TRUE these are neglected.
#map_in: name of input ASCII map
#map_out: name of output map
#sub_dir: if the clone, map_in and map_out are in a subdirectory add here as a string. If only some of these are in a subdirectory,
# add this to the name itself.

asc2map <- function(clone = "mask.map", asc_header = TRUE, options = "", map_in = "in.asc", map_out = "out.map",
                     sub_dir = "") {
  # check if pcr_dir exists
  if (!exists("pcr_dir")) {
    stop("Please set PCraster installation with 'set_pcraster()")
  }
  sys_type <- Sys.info()['sysname']
  exe <- ifelse(sys_type == "Windows", ".exe", "")
  asc_option = ifelse(asc_header == T, "-a ", " ") 
  command <- paste0(pcr_dir, "asc2map", exe, " --clone ", sub_dir, clone, " ", asc_option,
                    options, " ", sub_dir, map_in, " ", sub_dir, map_out)
  system(command)
}

# pcrcalc

# this function assumes 'set_pcraster()' is already used.
# pcrcalc is the general function of PCRaster which allows all kind of 
# map manipulations, for further details see:
# https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/index.html#functional-list-of-applications-and-operations
# the options for a specific pcrcalc command should be given as a character string with 'spaces' in between.
# work_dir: directory where the script should be executed

pcrcalc <- function(options = "", work_dir) {
  # check if pcr_dir exists
  if (!exists("pcr_dir")) {
    stop("Please set PCraster installation with 'set_pcraster()")
  }
  # set working directory to execute command - return to project at end of function
  projwd <- getwd()
  setwd(paste0("./", work_dir))
  
  sys_type <- Sys.info()['sysname']
  exe <- ifelse(sys_type == "Windows", ".exe", "")
  command <- paste0(pcr_dir, "pcrcalc", exe, " ", options)
  system(command)
  setwd(projwd)
}


# execute PCRaster script file

# this function assumes 'set_pcraster()' is already used.
# a script file with pcrcalc commands can be executed from the command line
# for further explanation see
# https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/secdyn.html#secseqbla

#script: name of script file
#script_dir: directory where the script is stored, either relative to project or absolute
#work_dir: directory where the script should be executed
#script_path_rel: (default = TRUE) relative or absolute path for script_dir

pcr_script <- function(script, script_dir, work_dir, script_path_rel = TRUE) {
  # check if pcr_dir exists
  if (!exists("pcr_dir")) {
    stop("Please set PCraster installation with 'set_pcraster()")
  }
  sys_type <- Sys.info()['sysname']
  exe <- ifelse(sys_type == "Windows", ".exe", "")
  # set working directory to execute command - return to project at end of function
  projwd <- getwd()
  setwd(paste0("./", work_dir))
  # make correct path to script
  if (!script_path_rel) {
    script_prefix <- paste0(script_dir, "/")
  } else {
  script_prefix <- paste0(projwd, "/", script_dir, "/")
  }
  command <- paste0(pcr_dir, "pcrcalc", exe, " -f ", script_prefix, script)
  system(command)
  setwd(projwd)
}

# make a new map - newmap
# this function assumes 'set_pcraster()' is already used.
# for explanation of function see pcraster documentation
# https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/app_mapattr.html

#mapname: name of the new map, including .map
#nrows: number of rows of the map
#ncols: number of columns of the map
#cellsize: resolution of the map (meters)
#datatype: datatype that the map should have (B, N, S, O, D or L), check documentation
#   for details. default = scalar (S)
#xulc: x coordinate of upper left corner
#yulc: y coordinate of upper left corner
#dir: directory in which the map will be made
#WARNING: this function overwrites any map with the same name that already exists

newmap <- function(mapname = "mask.map", nrows = NULL, ncols = NULL,
                   datatype = "S", xulc = NULL, yulc = NULL,
                   cellsize = NULL, dir = "") {
  # check if pcr_dir exists
  if (!exists("pcr_dir")) {
    stop("Please set PCraster installation with 'set_pcraster()")
  }
  map <- paste0(dir, mapname)
  if (file.exists(map)) {
    command <- paste0("rm ", map)
    system(command)
  }
  sys_type <- Sys.info()['sysname']
  exe <- ifelse(sys_type == "Windows", ".exe", "")
  command <- paste0(pcr_dir, "mapattr", exe, " -s -R ", nrows, " -C ", ncols, 
                    " -", datatype, " -P yb2t -x ", xulc, " -y ", yulc, " -l ", 
                    cellsize, " ", dir, mapname)
  system(command)
}

# resample - set same location
# this function assumes 'set_pcraster()' is already used.
# for explanation of function see pcraster documentation
# https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/app_resample.html

#clone: name of clone map
#asc_header: (TRUE) if true the option '-a' will be passed which sets all other options from the ASCII header
#options: a string with all other options separated by " ". With asc_header = TRUE these are neglected.
#map_in: name of input ASCII map
#map_out: name of output map
#sub_dir: if the clone, map_in and map_out are in a subdirectory add here as a string. If only some of these are in a subdirectory,
# add this to the name itself.

resample_location <- function(clone = "mask.map", map_in = "in.map", map_out = "out.map",
                    dir = "") {
  # check if pcr_dir exists
  if (!exists("pcr_dir")) {
    stop("Please set PCraster installation with 'set_pcraster()")
  }
  
  sys_type <- Sys.info()['sysname']
  exe <- ifelse(sys_type == "Windows", ".exe", "")
  command <- paste0(pcr_dir, "resample", exe, " --clone ", dir, clone, " ", 
                    dir, map_in, " ", dir, map_out)
  system(command)
}



## test -------------
#set_pcraster(env = "lisem")

#asc2map(map_in = "R1/dem1.asc", map_out = "R1/dem1.map", sub_dir = "RS_data/")

#script <- "RS_data/template_script_cal.txt"
#dir <- "RS_data/shared_maps"
