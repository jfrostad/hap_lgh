# ----HEADER------------------------------------------------------------------------------------------------------------
#***********************************************************************************************************************

# ----SETUP-------------------------------------------------------------------------------------------------------------
# clear memory
rm(list=ls())
# runtime configuration
if (Sys.info()["sysname"] == "Linux") {
  j_root <- "/home/j/"
  h_root <- file.path("/ihme/homes", Sys.info()["user"])
  
  package_lib    <- file.path(h_root, '_code/_lib/pkg')
  ## Load libraries and  MBG project functions.
  .libPaths(package_lib)
  
  # necessary to set this option in order to read in a non-english character shapefile on a linux system (cluster)
  Sys.setlocale(category = "LC_ALL", locale = "C")
  
} else {
  j_root <- "J:"
  h_root <- "H:"
}

#load external packages
pacman::p_load(assertthat, data.table, dplyr, mgsub, raster, sf, fasterize, fst)

#detect if running interactively
interactive <- F  %>% #manual override
  ifelse(., T, !length(commandArgs())>2) %>%  #check length of arguments being passed in
  ifelse(., T, !(is.na(Sys.getenv("RSTUDIO", unset = NA)))) #check if IDE

if (interactive) {

} else {
  
  ## Set repo location, indicator group, and some arguments
  user            <- commandArgs()[4]
  core_repo       <- commandArgs()[5]
  indicator_group <- commandArgs()[6]
  indicator       <- commandArgs()[7]
  config_par      <- commandArgs()[8]
  cov_par         <- commandArgs()[9]
  reg             <- commandArgs()[10]
  run_date        <- commandArgs()[11]
  measure         <- commandArgs()[12]
  holdout         <- as.numeric(commandArgs()[13])
  my_repo         <- commandArgs()[14]
  age             <- 0

}

#analysis options
new_gbd_estimates <- T #set TRUE if GBD best model has been updated
interval_mo = 12

#***********************************************************************************************************************

# ---FUNCTIONS----------------------------------------------------------------------------------------------------------
#use your own diacritics fx, due to inscrutable error
#note: requires mgsub pkg
fix_diacritics <<- function(x) {
  
  require(mgsub)
  
  #first define replacement patterns as a named list
  defs <-
    list('??'='S', '??'='s', '??'='Z', '??'='z', '??'='A', '??'='A', '??'='A', '??'='A', '??'='A', '??'='A', '??'='A', 
         '??'='C', '??'='E', '??'='E','??'='E', '??'='E', '??'='I', '??'='I', '??'='I', '??'='I', '??'='N', '??'='O', 
         '??'='O', '??'='O', '??'='O', '??'='O', '??'='O', '??'='U','??'='U', '??'='U', '??'='U', '??'='Y', '??'='B', 
         '??'='a', '??'='a', '??'='a', '??'='a', '??'='a', '??'='a', '??'='a', '??'='c','??'='e', '??'='e', '??'='e', 
         '??'='e', '??'='i', '??'='i', '??'='i', '??'='i', '??'='o', '??'='n', '??'='o', '??'='o', '??'='o', '??'='o',
         '??'='o', '??'='o', '??'='u', '??'='u', '??'='u', '??'='y', '??'='y', '??'='b', '??'='y', '??'='Ss')
  
  #then force conversion to UTF-8 and replace with non-diacritic character
  enc2utf8(x) %>% 
    mgsub(., pattern=enc2utf8(names(defs)), replacement = defs) %>% 
    return
  
}

## Load custom post-estimation functions
file.path(my_repo, '_lib', 'post', 'aggregate_inputs.R') %>% source

#***********************************************************************************************************************

# ---PREP CONFIG--------------------------------------------------------------------------------------------------------
## Read config file and save all parameters in memory
config_filepath <- 'cooking/model/configs/'
config <- set_up_config(repo            = my_repo,
                        indicator_group = indicator_group,
                        indicator       = indicator,
                        config_name     = paste0('/model/configs/config_', config_par),
                        covs_name       = paste0('/model/configs/covs_', cov_par),
                        run_tests       = F,
                        post_est_only   = T,
)

# Get the necessary variables out from the config object into global env
rake_countries <- eval(parse(text = config[V1 == 'rake_countries', V2]))
rake_subnational <- eval(parse(text = config[V1 == 'subnational_raking', V2]))
modeling_shapefile_version <- config[V1 == 'modeling_shapefile_version', V2]
raking_shapefile_version <- config[V1 == 'raking_shapefile_version', V2]
countries_not_to_rake <- config[V1 == 'countries_not_to_rake', V2]
countries_not_to_subnat_rake <- config[V1 == 'countries_not_to_subnat_rake', V2]
year_list <- eval(parse(text = config[V1 == 'year_list', V2]))
metric_space <- config[V1 == 'metric_space', V2]
summstats <- eval(parse(text = config[V1 == 'summstats', V2]))

# Print some settings to console
message(indicator)
message(indicator_group)
message(run_date)
message(reg)
message(pop_measure)
message(holdout)
## Define raking parameters ---------------------------------------------------------------------------

## If this is not a solid fuel indicator model run set all countries to not be raked
if (!(indicator %like% 'solid')) {
  
  rake_countries <- F

  # Create function to pull isos for diarrhea custom regs
  get_region_isos <- function(modeling_region) {
    
    # define regions
    region_list <- list(
      'dia_afr_horn' = 'dji+eri+eth+sdn+som+ssd+yem',
      'dia_cssa-cod' = 'ago+caf+cog+gab+gnq+stp',
      'dia_wssa-nga' = 'ben+bfa+civ+cmr+cpv+gha+gin+gmb+gnb+lbr+mli+mrt+ner+sen+sle+tcd+tgo',
      'dia_name' = 'dza+egy+esh+lby+mar+tun',
      'dia_sssa' = 'bwa+nam+zaf',
      'dia_mcaca' = 'blz+cri+cub+dma+dom+grd+gtm+hnd+hti+jam+lca+mex+nic+pan+slv+vct',
      'dia_s_america-bra' = 'bol+col+ecu+guf+guy+per+pry+sur+tto+ven',
      'dia_central_asia' = 'kgz+tjk+tkm+uzb',
      'dia_se_asia' = 'khm+lao+mmr+mys+tha+vnm',
      'dia_malay' = 'idn+phl+png+tls',
      'dia_south_asia-ind-pak' = 'bgd+btn+lka+npl',
      'dia_mid_east' = 'afg+irn+irq+jor+pse+syr',
      'dia_essa-zwe-ken' = 'bdi+com+lso+mdg+moz+mwi+rwa+swz+syc+tza+uga+zmb',
      'PAK' = 'pak',
      'KEN' = 'ken',
      'NGA' = 'nga',
      'COD' = 'cod',
      'IND' = 'ind',
      'ZWE' = 'zwe',
      'MNG' = 'mng',
      'dia_cssa' = 'ago+caf+cod+cog+gab+gnq+stp',
      'dia_wssa' = 'ben+bfa+civ+cmr+cpv+gha+gin+gmb+gnb+lbr+mli+mrt+ner+nga+sen+sle+tcd+tgo',
      'dia_s_america' = 'bol+bra+col+ecu+guf+guy+per+pry+sur+tto+ven',
      'dia_chn_mng' = 'chn+mng',
      'dia_south_asia' = 'bgd+btn+ind+lka+npl+pak',
      'dia_south_asia-ind' = 'bgd+btn+lka+npl+pak',
      'dia_essa' = 'bdi+com+ken+lso+mdg+moz+mwi+rwa+swz+syc+tza+uga+zmb+zwe'
    )
    # return region list
    return(toupper(region_list[[modeling_region]]))
  }
    
  # apply function
  countries_not_to_rake <- get_region_isos(reg)
  countries_not_to_subnat_rake <- get_region_isos(reg)
  rake_subnational <- FALSE
  
} else {
  
  #rake all countries subnationally
  rake_subnational <- TRUE
  countries_not_to_subnat_rake <- NULL
  #countries_not_to_rake <- 'GAB' #problems with GBD estimate due to faulty 2012 DHS
  countries_not_to_rake <- NULL #problems with GBD estimate due to faulty 2012 DHS
  
}

# Determine if a crosswalk is needed
crosswalk <- ifelse(modeling_shapefile_version != raking_shapefile_version, T, F)

# Using logit raking
rake_method <- 'logit'

# Print raking info
message('Metric Space                       : ', metric_space)
message('Subnational raking                 : ', rake_subnational)
message('Countries not to rake at all       : ', countries_not_to_rake)
message('Countries not to rake subnationally: ', countries_not_to_subnat_rake)
#***********************************************************************************************************************

# ---LOAD GBD TARGETS---------------------------------------------------------------------------------------------------
#reload estimates from the central db if they have changed
if(rake_countries) {
  if (new_gbd_estimates) {
    
    source("/ihme/cc_resources/libraries/current/r/get_location_metadata.R")
    source("/ihme/cc_resources/libraries/current/r/get_draws.R")
    locations <- get_location_metadata(location_set_id = 35, gbd_round_id = 6, decomp_step = "step4") %>% as.data.table
    loc_ids <- unique(locations$location_id)
    
    #pull the draws from 
    #got this pull from Sarah Wozniak
    #note that the age/sex IDs are arbitrary since this model doesnt vary over those params
    #note, could also just pull the model from GPR, but the results are the same
    #source("/ihme/code/st_gpr/central/stgpr/r_functions/utilities/utility.r")
    #hap <- model_load(102800,"raked")
    gbd <-
    get_draws("rei_id", 87, source="exposure", status="best", year_id=1990:2019,
              location_id=loc_ids, sex_id=2, age_group_id=11, gbd_round_id=6, decomp_step="step4",
              num_workers=8) %>% 
      .[parameter=='cat1'] #only want exposure, not the inverse

    #summarize
    draw.cols <- paste0('draw_', 0:999) #i hate zero indexing
    gbd[, lower := apply(.SD, 1, quantile, c(.025)), .SDcols=draw.cols]
    gbd[, mean := apply(.SD, 1, mean), .SDcols=draw.cols]
    gbd[, upper := apply(.SD, 1, quantile, c(.975)), .SDcols=draw.cols]
    #gbd[, (draw.cols) := NULL] #no longer need
    
    #format for MBG
    setnames(gbd, c('location_id', 'year_id'), c('name', 'year'))
  
    write.csv(gbd, 
              file=file.path('/share/geospatial/jfrostad', indicator_group, 'data', 
              paste0('gbd_2019_best_', indicator, '_', measure, '.csv')),
              row.names=F)
  
  } else {
  
    #kw_gbd <- as.data.table(read.csv(paste0(rake_to_path, 'gbd_',  measure, '.csv'), stringsAsFactors = FALSE))
    
    gbd <- file.path('/share/geospatial/jfrostad', indicator_group, 'data', 
                              paste0('gbd_2019_best_', indicator, '_', measure, '.csv')) %T>% 
      message('reading GBD best estimates from this path\n', .) %>% 
      fread 
  
  }
}

#***********************************************************************************************************************

# ---PREP INPUTS--------------------------------------------------------------------------------------------------------
# Load cell draws
message("Loading Data...")
load(paste0(outputdir, indicator, "_cell_draws_eb", pathaddin, ".RData"))

# Check if load was successful; stop if not
if (!exists("cell_pred")) {
  message(filename_rds)
  stop("Unable to load cell_pred object! Check to make sure that the relevant object exists.")
}

# Rake estimates
if (rake_countries) {
  if (!exists("gbd")) {
    stop("rake_countries was specified as T in config, gbd raking targets must be provided.")
  }
  
  ## determine if a crosswalk is needed
  crosswalk <- modeling_shapefile_version != raking_shapefile_version
  
  # Assume linear raking unless specified as logit
  if (rake_transform == "logit") rake_method <- "logit" else rake_method <- "linear"
  
  
  ##### Prep input data into raking:
  
  ## Get the simple and new_simple rasters prepped up for us
  print("Getting simple and prepped rasters")
  raster_outputs <- prep_shapes_for_raking(
    reg = reg,
    modeling_shapefile_version = modeling_shapefile_version,
    raking_shapefile_version = raking_shapefile_version,
    field = "loc_id",
    use_sf=T
  )
  
  ## Take out the objects from the list that actually matters to us:
  simple_raster <- raster_outputs[["simple_raster"]]
  new_simple_raster <- raster_outputs[["new_simple_raster"]]
  
  simple_polygon <- raster_outputs[["simple_polygon"]]
  new_simple_polygon <- raster_outputs[["new_simple_polygon"]]
  
  pixel_id <- raster_outputs[["pixel_id"]]

  # Load populations
  message('Loading populations')
  gbd_pops <- prep_gbd_pops_for_fraxrake(pop_measure = pop_measure, reg = reg, year_list = year_list, 
                                         gbd_round_id = 6,
                                         decomp_step='step4')
#***********************************************************************************************************************

# ---RAKE/AGG-----------------------------------------------------------------------------------------------------------
#Using fractional raking
if (metric_space == "rates") {
    print("Using the rates raking and aggregation functions:")

    outputs <- 
      fractionally_rake_rates(
        cell_pred = cell_pred,
        simple_raster = simple_raster,
        simple_polygon = simple_polygon,
        pixel_id = pixel_id,
        shapefile_version = raking_shapefile_version,
        reg = reg,
        pop_measure = pop_measure,
        measure = measure,
        year_list = year_list,
        interval_mo = interval_mo,
        rake_subnational = rake_subnational,
        sharedir = sharedir,
        run_date = run_date,
        indicator = indicator,
        main_dir = outputdir,
        rake_method = 'logit',
        age = age,
        holdout = holdout,
        countries_not_to_rake = countries_not_to_rake,
        countries_not_to_subnat_rake = countries_not_to_subnat_rake,
        return_objects = TRUE
      )
    
    ## Get the necessary outputs and rename the columns
    rf <- data.table(outputs[["rf"]])[, .(loc = location_id, year, start_point = mbg_prev, target = gbd_prev, raking_factor = rf)]
    raked_cell_pred <- outputs[["raked_cell_pred"]]
    
    
    ## Raked simple raster has been made above
    raked_simple_raster <- new_simple_raster
    
  } else if (metric_space == "counts") {
    print("Using the counts raking and aggregation functions:")
    
    ## Rake counts
    outputs <- fractionally_rake_counts(
      count_cell_pred = data.table(cell_pred),
      rake_to = gbd,
      reg = reg,
      year_list = year_list,
      rake_subnational = rake_subnational,
      countries_not_to_subnat_rake = countries_not_to_subnat_rake,
      countries_not_to_rake = countries_not_to_rake,
      simple_raster = simple_raster,
      modeling_shapefile_version = modeling_shapefile_version,
      raking_shapefile_version = raking_shapefile_version
    )
    
    ## Get the necessary outputs
    rf <- outputs$raking_factors
    raked_cell_pred <- outputs$raked_cell_pred
    raked_simple_raster <- new_simple_raster
    
    ## Save out the aggregate files
    raked_frax_counts_save(output_list = outputs, sharedir = sharedir, indicator = indicator, age = age, reg = reg, holdout = holdout)
  }

} else {
  rf <- NULL
  raked_cell_pred <- NULL
  
  # Define simple raster for mask
  simple_polygon_list <- load_simple_polygon(
    gaul_list = get_adm0_codes(reg,
                               shapefile_version = modeling_shapefile_version
    ),
    buffer = 0.4, subset_only = FALSE,
    shapefile_version = modeling_shapefile_version
  )
  subset_shape <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  raster_list <- build_simple_raster_pop(subset_shape)
  simple_raster <- raster_list[["simple_raster"]]
  rm(simple_polygon_list)
}
#***********************************************************************************************************************

# ---SAVE --------------------------------------------------------------------------------------------------------------
message("Saving results...")

## save RF
save_post_est(rf, "csv", paste0(reg, "_rf"))

# make and save summaries

save_cell_pred_summary <- function(summstat, raked, ...) {
  message(paste0("Making summmary raster for: ", summstat, " (", raked, ")"))
  
  if (raked == "unraked") {
    cpred <- "cell_pred"
    mask_raster <- "simple_raster"
  }
  if (raked == "raked") {
    cpred <- "raked_cell_pred"
    mask_raster <- "raked_simple_raster"
  }
  if (raked == "raked_c") {
    cpred <- "raked_cell_pred_c"
    load(paste0(sharedir, "/output/", run_date, "/", indicator, "_raked_c_cell_draws_eb", pathaddin, ".RData" ))
    mask_raster <- "raked_simple_raster"
  }
  ras <- make_cell_pred_summary(
    draw_level_cell_pred = get(cpred),
    mask = get(mask_raster),
    return_as_raster = TRUE,
    summary_stat = summstat,
    ...
  )
  save_post_est(ras,'raster',paste0(reg, ifelse(raked == "raked", "_raked", ifelse(raked == 'raked_c', '_raked_c', '')), '_', summstat, '_raster'))
}

# Do this as lapply to not fill up memory in global env with big obs
if (!exists('gbd')) {
  rake_list <- c("unraked")
} else {
  rake_list <- c("unraked", "raked")
}

summ_list <- expand.grid(summstats[summstats != "p_below"], rake_list)

lapply(1:nrow(summ_list), function(i) {
  summstat <- as.character(summ_list[i, 1])
  raked <- as.character(summ_list[i, 2])
  save_cell_pred_summary(summstat, raked)
})

## Can't pass additional params in the above framework, so will code by hand here
for (r in rake_list) {
  if ("p_below" %in% summstats) {
    save_cell_pred_summary(
      summstat = "p_below",
      raked = r,
      value = 0.8,
      equal_to = F
    )
  }
}

# Write a file to mark done
output_dir <- paste0("/share/geospatial/mbg/", indicator_group, "/", indicator, "/output/", run_date)
pathaddin <- paste0("_bin0_", reg, "_0") # To allow us to use waitformodelstofinish()
write(NULL, file = paste0(output_dir, "/fin_", pathaddin))

# All done
message(paste0("Done with post-estimation and aggregation for ", reg))
#***********************************************************************************************************************

# ---POST_ESTIMATION----------------------------------------------------------------------------------------------------
## Launch post-estimation (TAP calculation)
if(indicator%like%'solid|dirty' & holdout==0) {
  # Define best LRI run_date and cluster specs
  lri_run_date = '2019_10_23_16_13_17'
  proj_arg <- 'proj_geo_nodes'
  use_geos_nodes <- T
  
  # set memory based on region
  if (reg %in% c('CHN', 'trsa-GUF', 'ansa-VEN')) { mymem <- '650G'
  } else if (reg %in% c('wssa-CPV-NGA', 'soas', 'ocea-MYS', 'caca-CUB')) { mymem <- '500G'
  } else mymem <- '350G'
  
  jname           <- paste('EdL', reg, indicator, sep = '_')
  
  # set up qsub
  sys.sub <- paste0('qsub -e ', outputdir, '/errors -o ', outputdir, '/output ', 
                    '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                    '-l fthread=2 -l h_rt=03:00:00:00 -v sing_image=default -N ', jname, ' -l archive=TRUE ')
  r_shell <- file.path(core_repo, 'mbg_central/share_scripts/shell_sing.sh')
  script <- file.path(my_repo, indicator_group, 'post/3_descent.R')
  args <- paste(user, core_repo, indicator_group, indicator, config_par, cov_par, reg, run_date, measure, holdout, my_repo)
  
  
  # submit qsub
  paste(sys.sub, r_shell, script, args) %>% 
    system
  
  ## All done
  message('Post-est submitted, now producing diagnostics for ', reg)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LINEPLOT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Produce diagnostic lineplot for region

## Aggregate data and stackers ------------------------------------------------------
# Aggregate data to admin 0 and 1
dat <- aggregate_input_data(reg,
                            indicator, 
                            indicator_group, 
                            run_date,
                            modeling_shapefile_version,
                            build=T)

# Aggregate stackers to admin 0 and 1
stack <- aggregate_child_stackers(reg,
                                  indicator, 
                                  indicator_group, 
                                  run_date, 
                                  modeling_shapefile_version,
                                  pop_measure=pop_measure,
                                  build=T)
