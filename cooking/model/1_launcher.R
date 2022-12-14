## Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

# set general arguments
user            <- Sys.info()['user']
core_repo            <- file.path('/homes', user, '_code/lbd/lbd_core/')
core_repo            <- file.path('/homes', user, '_code/lbd/hap/')
my_repo            <- file.path('/homes', user, '_code/lbd/hap/')
indicator_group <- 'cooking'
parallel_script <- file.path('model/3_orbit')

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# set covariate arguments
plot_covariates <- TRUE
covariate_plotting_only <- FALSE

# indicate whether to use old run date
use_old_run_date <- T
old_run_date_input <- "2020_10_27_23_38_30"

# set run date
if (use_old_run_date == FALSE) {
  run_date <- make_time_stamp(TRUE)
} else {
  run_date <- old_run_date_input
}

# set config and covariate files
config_par   <- 'hap_sp_fine'
covar_par      <- 'region_specific'

# set whether running for individual countries
individual_countries <- FALSE


# indicate holdout (also need to do so in config)
holdout <- F # only matters if running aggregation, set to TRUE for holdouts

# list all regions or countries
# standard regions
regions <- c('dia_afr_horn', 'dia_cssa', 'dia_wssa', 'dia_name-ESH', 'dia_sssa', 
             'dia_mcaca', 'dia_s_america-GUF', 'dia_central_asia', 'dia_chn_mng', 
             'dia_se_asia', 'dia_malay', 'dia_south_asia', 'dia_mid_east', 'dia_essa')

# custom country-specifics
regions <- c('essa-ERI-DJI-YEM', "ERI+DJI+YEM",
             'sssa-ZAF', 'ZAF',
             'cssa-GNQ',
             'wssa-CPV-NGA', 'NGA',
             'noaf-ESH',
             'caca-CUB',
             'ansa-VEN', 'trsa-GUF',
             'stan-TKM',
             'CHN', 'MNG',
             'ocea-MYS',
             'seas',
             'mide+TKM', 'soas')

regions <- c('CHN', 'noaf-ESH', 'ocea-MYS', 'ansa-VEN')

# list indicators
indics <- 'cooking_fuel_solid'

## Run launch scripts -------------------------------------------------------------------------

for (i in indics) {
   
  # make sure that only selecting a previous run_date intentionally
  if (use_old_run_date == TRUE) {
    prev <- readline('Are you sure you want to use a previous run date? Y or N: ')
    if (prev != 'Y') stop('Set use_old_run_date to FALSE.')
  }
  
  for (reg in regions) {
    
    # set specific arguments
    indicator       <- i
    jname           <- paste('rocket', indicator_group, reg, sep = '_')
    mymem           <- '20G'
    
    # set region specific covariates, if desired
    if (covar_par == 'region_specific') cov_par <- paste0('cooking_', reg)
    else cov_par <- covar_par
    
    # some quick checks for the arguments
    if(use_old_run_date == TRUE & old_run_date_input == '') stop('You indicated using an old run date; please provide an old run date')
    
    # set up qsub
    sys.sub <- paste0('qsub -e /share/temp/sgeoutput/', user,'/errors -o /share/temp/sgeoutput/', user, '/output ', 
                      '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                      '-l fthread=1 -l h_rt=', ifelse(use_geos_nodes, '16:00:00:00', '1:00:00:00'),
                      ' -v sing_image=default -N ', jname, ' -l archive=TRUE ')
    r_shell <- file.path(core_repo, 'mbg_central/share_scripts/shell_sing.sh')
    script <- file.path('/homes', user, '_code/lbd/hap', indicator_group, 'model/2_rocket.R')
    args <- paste(user, core_repo, indicator_group, indicator, config_par, cov_par, reg, parallel_script,
                  plot_covariates, covariate_plotting_only, proj_arg, use_geos_nodes, run_date, my_repo)

    # run launch script
    paste(sys.sub, r_shell, script, args) %>% 
      system
    
  }
  
}