# ---------------------------------------------------------------------------------------------
# aggregate_input_data()
#
# Written by Ani Desphande
# Modified by Kirsten Wiens
#
# 07/23/2019
# Function that aggregates model input data in a population-weighted manner
#
# Inputs:
# Indicator, indicator group, run date, reg (which can be a vector of regions),
# and shapefile version
#
# Outputs:
# List of admin 0 (ad0) and admin 1 (ad1) aggregated input data
# ---------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# Start function
aggregate_input_data <- function(reg,
                                 indicator,
                                 indicator_group,
                                 run_date,
                                 modeling_shapefile_version, 
                                 build=T, 
                                 debug=F) {
  # -----------------------------------------------------------------------------------
  # Work interactively to build function
  if (debug) browser()
  message(paste0('Aggregating input data for: ', reg))
  # ---------------------------------------------------------------
  # Set-up
  # Set arguments
  output_dir <- file.path('/share/geospatial/mbg', indicator_group, 
                         indicator, 'output', run_date, '/pred_derivatives/admin_summaries/') %T>% 
    dir.create(recursive=T) # create if first run
  output_file <- paste0(output_dir,  '/', reg, "_adm_input_data.fst")
  
  # Save files if not already produced once
  if (build | !(output_file %>% file.exists)) {
    # ------------------------------------------------------------------------------
    # Load necessary objects
    
    # Load input data
    input_data <- paste0('/share/geospatial/mbg/input_data/', indicator, '.csv') %>% fread %>% 
      .[, prop := get(indicator)/N]
    
    # Reading in shapefiles & population raster
    file_dir <- paste0('/home/j/WORK/11_geospatial/admin_shapefiles/', modeling_shapefile_version, '/')
    a0_shp <- readRDS(paste0(file_dir,'lbd_standard_admin_0.rds'))
    a1_shp <- readRDS(paste0(file_dir,'lbd_standard_admin_1.rds'))
    pop_ras <- file.path('/home/j/WORK/11_geospatial/01_covariates/00_MBG_STANDARD/worldpop/',
              'total',
              pop_release,
              '1y',
              'worldpop_total_1y_2010_00_00.tif') %>% raster
    
    # Read in lookup table
    lookup <- fread('/home/j/WORK/11_geospatial/10_mbg/stage_master_list.csv')
    # ------------------------------------------------------------------------------
    
    # ----------------------------------------------------------------------------------------------------------------
    # Aggregate data
    # Format shps as rasters
    a0_shp$ADM0_CODE <- as.numeric(as.character(a0_shp$ADM0_CODE))
    a1_shp$ADM1_CODE <- as.numeric(as.character(a1_shp$ADM1_CODE))
    a0_raster <- fasterize(st_as_sf(a0_shp), pop_ras, field = 'ADM0_CODE')
    a1_raster <- fasterize(st_as_sf(a1_shp), pop_ras, field = 'ADM1_CODE')
    
    # Format input data & extract adminIDs
    gaul_list <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
    input_data <- input_data[country %in% lookup[gadm_geoid %in% gaul_list, iso3]]
    input_data[, ADM1_CODE := raster::extract(a1_raster, input_data[, .(longitude, latitude)])]
    input_data[, ADM0_CODE := raster::extract(a0_raster, input_data[, .(longitude, latitude)])]

    # Summarize input data by admin 0 IDs
    a0_input <- copy(input_data) %>% 
      setkeyv(., c('nid', 'year', 'ADM0_CODE')) %>% 
      .[, input_mean := weighted.mean(prop, sum_of_sample_weights*weight), by = key(.)] %>% 
      .[, input_ss := sum(weight*N), by = key(.)] %>% 
      unique(., by=key(.)) %>% 
      .[, .(ihme_loc_id, nid, year, ADM0_CODE, input_mean, input_ss)] %>% 
      .[, lvl := 'adm0']
    
    # Summarize input data by admin 1 IDs                                  
    a1_input <- copy(input_data) %>% 
      setkeyv(., c('nid', 'year', 'ADM0_CODE', 'ADM1_CODE')) %>% 
      .[, input_mean := weighted.mean(prop, sum_of_sample_weights*weight), by = key(.)] %>% 
      .[, input_ss := sum(weight*N), by = key(.)] %>% 
      unique(., by=key(.)) %>% 
      .[, .(ihme_loc_id, nid, year, ADM0_CODE, ADM1_CODE, input_mean, input_ss)] %>% 
      .[, lvl := 'adm1']
  
  # Save regional data to single DT
    list(a0_input, a1_input) %>% 
      rbindlist(use.names=T, fill=T) %>% 
      #.[, region := reg] %>% 
      write_fst(., path=output_file) %T>%
    return
  # ----------------------------------------------------------------------------------------------------------------
  # otherwise just read them in and return to preserve pipeline
  } else read_fst(output_file, as.data.table=T) %>% return

}
# -----------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# aggregate_child_stackers()
#
# Written by Kirsten Wiens
#
# 07/24/2019
# Function that aggregates child stacker results in a population-weighted manner
#
# Inputs:
# Indicator, indicator group, run date, reg (which can be a vector of regions),
# shapefile version, population measure for population-weighted aggregation, age, and holdout
#
# Outputs:
# List of admin 0 (ad0) and admin 1 (ad1) aggregated child stacker results
#
# TODO:
# - Update so that we don't need to load master shapefile each time
# - Add in option to supply simple polygon or raster separately 
# - Add option to save output (which would be helpful if added to parallel script)
# ---------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# Start function
aggregate_child_stackers <- function(reg,
                                     indicator,
                                     indicator_group,
                                     run_date,
                                     modeling_shapefile_version,
                                     pop_measure,
                                     age = 0,
                                     holdout = 0,
                                     build=T,
                                     debug=F) {
  
  # -----------------------------------------------------------------------------------
  # Work interactively to build function
  if (debug) browser()
  message(paste0('Aggregating stackers for: ', reg))
  
  # Set arguments
  output_dir <- file.path('/share/geospatial/mbg', indicator_group, 
                          indicator, 'output', run_date, '/pred_derivatives/admin_summaries/') %T>% 
    dir.create(recursive=T) # create if first run
  output_file <- paste0(output_dir,  '/', reg, "_adm_input_stackers.fst")
  # --------------------------------------------------------------------------------------
  # Save files if not already produced once
  if (build | !(output_file %>% file.exists)) {
  # ---------------------------------------------------------------------------------------------------------------
    # Load rasters and stackers
    
    # Get the simple and new_simple rasters prepped up for us
    raster_outputs <- prep_shapes_for_raking(
      reg = reg,
      modeling_shapefile_version = modeling_shapefile_version,
      raking_shapefile_version = modeling_shapefile_version,
      field = 'loc_id'
    )
    
    # Take out the objects from the list that actually matters to us:
    simple_raster <- raster_outputs[['simple_raster']]
    simple_polygon <- raster_outputs[['simple_polygon']]
    pixel_id <- raster_outputs[['pixel_id']]
    
    # Add stacking children to this process
    covs = fetch_from_rdata(paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/model_image_history/', 
                                   run_date, '_bin', age,'_', reg, '_', holdout, '.RData'), 'cov_list')
    fes = fetch_from_rdata(paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/model_image_history/', 
                                  run_date, '_bin', age,'_', reg, '_', holdout, '.RData'), 'all_fixed_effects')
    submodels = trimws(strsplit(fes, '+', fixed = T)[[1]])
    covs = covs[submodels]
    
    # Set stacker names
    covnames = names(covs)
    
    # Ensure the dimensions are the same
    for(ccc in covs){
      stopifnot(dim(ccc)[1:2] == dim(simple_raster)[1:2])
    }
    
    # Convert raster to data table and clean
    raster_to_dt = function(ras){
      dt <- crop(ras, extent(simple_raster))
      dt <- data.table(raster::extract(dt, pixel_id))
      dt[, pixel_id := pixel_id]
      dt <- melt(dt, id.vars='pixel_id')
      dt[, c('pixel_id', 'variable') := NULL]
      return(dt)
    }
    
    covdt = lapply(covs, raster_to_dt)
    covdt = do.call(what = cbind, covdt)
    setnames(covdt, names(covs))
    
    # Add pixel_id, but make sure that its recycled explicitly as per data.table 1.12.2 guidelines
    covdt[, pixel_id := rep(pixel_id, times = nrow(covdt) / length(pixel_id))]
    
    # Add year to covdt
    yyy = as.vector(unlist(lapply(min(year_list):max(year_list), function(x) rep.int(x, times = length(pixel_id)))))
    covdt[, year := yyy]
    
    # Free up a bit of space
    rm(covs)
    # ---------------------------------------------------------------------------------------------------------------
    
    
    # ----------------------------------------------------------------------------------------
    # Add in population

    pop <- load_populations_cov(reg, pop_measure=pop_measure, measure = 'count', simple_polygon, 
                                simple_raster, year_list, interval_mo=12, pixel_id = pixel_id)

    # Add to stacker results
    covdt <- merge(covdt, pop, by=c('pixel_id', 'year'))
    # ----------------------------------------------------------------------------------------
    
    
    # -------------------------------------------------------------------------------------------------------
    # Use link table to assign admin codes to stacker results
    
    # Load the cell id to admin units link
    link_table <- get_link_table(simple_raster, shapefile_version = modeling_shapefile_version)
    
    # Prep the cell_pred and link table to be linked by making sure they have the appropriate identifiers
    link <- prep_link_table(
      link_table = link_table,
      simple_raster = simple_raster,
      pixel_id = pixel_id
    )
    
    #added in to be explicit about rep issue:
    covdt[, pixel_id := rep(link_table$pixel_ids, times = nrow(covdt) / length(link_table$pixel_id))]
    
    # Merge identifiers with stacker data frame
    stackers <- merge(covdt, link, by = 'pixel_id', allow.cartesian = T)

    # Free up a bit of space
    rm(covdt, link, link_table, pop)
    # -------------------------------------------------------------------------------------------------------
    
    # -------------------------------------------------------------------
    # Aggregate
    # Get population for calculating weighted average
    stackers[, pop := pop*area_fraction]
    
    # Weighted mean by admin 0
    a0_stackers <- stackers[, lapply(.SD, weighted.mean, w = pop, na.rm = T),
                       by = c('ADM0_CODE', 'year'),
                       .SDcols = covnames] %>% 
      .[, lvl := 'adm0']
    
    # Weighted mean by admin 1
    a1_stackers <- stackers[, lapply(.SD, weighted.mean, w = pop, na.rm = T),
                       by = c('ADM0_CODE', 'ADM1_CODE', 'year'),
                       .SDcols = covnames] %>% 
      .[, lvl := 'adm1']
    # -------------------------------------------------------------------
    
    # Save regional data to single DT
    list(a0_stackers, a1_stackers) %>% 
      rbindlist(use.names=T, fill=T) %>% 
      #.[, region := reg] %>% 
      write_fst(., path=output_file) %T>%
      return
    # ----------------------------------------------------------------------------------------------------------------
    # otherwise just read them in and return to preserve pipeline
  } else read_fst(output_file, as.data.table=T) %>% return
  
}
# -----------------------------------------------------------------------------------