  ##############################################################################
  ##############################################################################

  
  ## Setup -------------------------------------------------------------------------
  
  # clear environment
  rm(list = ls())
  
  # set general arguments
  user            <- Sys.info()['user']
  core_repo            <- file.path('/homes', user, '_code/lbd/hap/')
  my_repo            <- file.path('/homes', user, '_code/lbd/hap/')
  indicator_group <- 'cooking'
  parallel_script <- file.path(indicator_group, 'model/parallel_hap')
  
  # Load MBG packages and functions
  message('Loading in required R packages and MBG functions')
  source(paste0(core_repo, '/mbg_central/setup.R'))
  mbg_setup(package_list = package_list, repos = core_repo)
  
  # set cluster arguments
  use_geos_nodes  <- T
  proj_arg        <- ifelse(use_geos_nodes, 'proj_geo_nodes', 'proj_geospatial')
  proj            <- ifelse(use_geos_nodes, paste0(' -P ', proj_arg, ' -l gn=TRUE '), paste0(' -P ', proj_arg, ' '))
  
  # set script arguments
  skip_entry <- F #use to skip to descent stage
  
  # set covariate arguments
  plot_covariates <- TRUE
  covariate_plotting_only <- FALSE
  
  # indicate whether to use old run date
  use_old_run_date <- T
  old_run_date_input <- '2020_10_27_23_38_30'
  
  # set run date
  if (use_old_run_date == FALSE) {
    run_date <- make_time_stamp(TRUE)
  } else {
    run_date <- old_run_date_input
  }
  
  # list all regions or countries
  # standard regions
  regions <- c('dia_afr_horn', 'dia_cssa', 'dia_wssa', 'dia_name-ESH', 'dia_name', 'dia_sssa', 
               'dia_mcaca', 'dia_s_america-GUF', 'dia_s_america', 'dia_central_asia', 'dia_chn_mng', 
               'dia_se_asia', 'dia_malay', 'dia_south_asia', 'dia_mid_east', 'dia_essa')
  
  # custom region list
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
  
  #regions <- c('ocea-MYS', 'trsa-GUF')
  regions <- c('CHN', 'noaf-ESH', 'ocea-MYS', 'ansa-VEN')
  
  ## Set repo location, indicator group, and some arguments
  user <- 'jfrostad'
  indicator_group <- 'cooking'
  indicator <- 'cooking_fuel_solid'
  config_par   <- 'hap_sp_fine'
  holdout <- 5
  age <- 0
  measure <- 'prev'



if(!skip_entry) {  
  for (reg in regions) {
  
    # set memory based on region
    if (reg %in% c('trsa-GUF')) { mymem <- '750G'
    } else if (reg %in% c('wssa-CPV-NGA', 'CHN', 'soas', 'ansa-VEN', 'ocea-MYS')) { mymem <- '500G'
    } else if (reg %in% c('seas-VNM-THA', 'VNM', 'THA', "ERI+DJI+YEM", 'sssa-ZAF')) { mymem <- '200G'
    } else mymem <- '350G'
    
      for (ho in holdout) {
  
        #name job
        jname           <- paste('eDL', reg, indicator, ho, sep = '_')
  
        #setup covars file
        cov_par <- paste(indicator_group, reg, sep='_')
  
        # set up qsub
        sys.sub <- paste0('qsub -e /share/temp/sgeoutput/', user,'/errors -o /share/temp/sgeoutput/', user, '/output ',
                          '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                          '-l fthread=1 -l h_rt=', ifelse(use_geos_nodes, '16:00:00:00', '3:00:00:00'),
                          ' -v sing_image=default -N ', jname, ' -l archive=TRUE ')
        r_shell <- file.path(core_repo, 'mbg_central/share_scripts/shell_sing.sh')
        script <- file.path(my_repo, indicator_group, 'post/2_entry.R')
        args <- paste(user, core_repo, indicator_group, indicator, config_par, cov_par, reg, run_date, measure, ho, my_repo)
  
        # run launch script
        paste(sys.sub, r_shell, script, args) %>%
          system
      }  
  
  }

} else {
  #use to launch descent instead if necessary
  
  for (reg in regions) {
  
    # set memory based on region
    if (reg %in% c('CHN', 'trsa-GUF', 'soas', 'noaf-ESH', 'wssa-CPV-NGA')) { mymem <- '900G'
    } else if (reg %in% c('ansa-VEN', 'cssa-GNQ', 'ocea-MYS', 'stan-TKM', 'mide+TKM')) { mymem <- '750G'
    } else if (reg %in% c("ERI+DJI+YEM")) { mymem <- '250G'
    } else mymem <- '500G'
  
      #name job
      jname           <- paste('EdL', reg, indicator, sep = '_')
  
      #setup covars file
      cov_par <- paste(indicator_group, reg, sep='_')
      
      # set up qsub
      sys.sub <- paste0('qsub -e /share/temp/sgeoutput/', user,'/errors -o /share/temp/sgeoutput/', user, '/output ',
                        '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                        '-l fthread=1 -l h_rt=', ifelse(use_geos_nodes, '16:00:00:00', '3:00:00:00'),
                        ' -v sing_image=default -N ', jname, ' -l archive=TRUE ')
      r_shell <- file.path(core_repo, 'mbg_central/share_scripts/shell_sing.sh')
      script <- file.path(my_repo, indicator_group, 'post/3_descent.R')
      args <- paste(user, core_repo, indicator_group, indicator, config_par, cov_par, reg, run_date, measure, holdout, my_repo)
  
      # run launch script
      paste(sys.sub, r_shell, script, args) %>%
        system
  
  }
}