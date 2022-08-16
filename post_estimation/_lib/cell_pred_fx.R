# ---------------------------------------------------------------------------------------------
# Author: JF
# Functions used to create long keyed data tables for doing fast uncertainty interval calculations
# ---------------------------------------------------------------------------------------------

## format_cell_pred ################################################
#TODO write documentation

format_cell_pred <- function(ind_gp,
                             ind,
                             rd,
                             reg,
                             measure,
                             pop_measure,
                             year_start,
                             year_end,
                             var_names = ind, # name using ind by default, but can pass custom name
                             matrix_pred_name = NULL,
                             skip_cols = NULL,
                             rk = T,
                             coastal_fix = T, # if model wasn't run w new coastal rasterization, force old simple raster process 
                             rake_subnational = rk, # TODO is this correct? might need to be defined custom by country
                             shapefile_version = 'current') {
  
  message('loading simple raster & populations')
  
  ## Load simple polygon template to model over
  gaul_list <- get_adm0_codes(reg, shapefile_version = shapefile_version)
  simple_polygon_list <- load_simple_polygon(
    gaul_list = gaul_list, buffer = 1, tolerance = 0.4,
    shapefile_version = shapefile_version
  )
  subset_shape <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  
  ## Load list of raster inputs (pop and simple)
  if (coastal_fix) raster_list <- build_simple_raster_pop(subset_shape, link_table = shapefile_version) #uses new simple_raster 
  else raster_list <- build_simple_raster_pop(subset_shape, link_table = NULL) #uses old rasterize
  
  simple_raster <- raster_list[["simple_raster"]]
  pop_raster <- raster_list[["pop_raster"]]
  pixel_id <- seegSDM:::notMissingIdx(simple_raster)

  ## get number of years
  year_list <- c(year_start:year_end)
  num_yrs <- length(year_list)
  
  message('loading links')
  #####################################################################
  # load the cell id to admin units link
  link_table <- get_link_table(simple_raster, shapefile_version = shapefile_version)
  
  #####################################################################
  # collect and load the population data from the WorldPop rasters
  covdt <- load_populations_cov(reg, pop_measure=pop_measure, measure = measure, simple_polygon, 
                                simple_raster, year_list, interval_mo=12, pixel_id = pixel_id)
  #####################################################################
  # Prepping the cell_pred and link table to be linked by making sure they have the appropriate identifiers.  Also performs a
  # zippering at the region boundary where cells that have a portion of their area outside of the modeling region are reintegrated
  # as a whole cell and considered to be only in one region.  This works becasue a given cell is only modeled in one region.
  link <- prep_link_table(
    link_table = link_table,
    simple_raster = simple_raster,
    pixel_id = pixel_id
  )

  # getting the connector for sub-national or national raking, This connector gets the IHME location_code for our
  # gbd targets and connects that to the ADM0_CODE or ADM1_CODE as nessecary
  connector <- get_gbd_locs(
    rake_subnational, 
    reg = reg,
    shapefile_version = shapefile_version
  )
  
  # merge the connector on to the link table by making sure that each cell fragment gets connected to the appropriate
  # raking geography
  link <- sub_nat_link_merge(
    rake_subnational,
    link,
    connector
  )

  #helper function to load a list of cell preds and then merge them together
  loadCellPreds <- function(i) { 
  
    message('~>loading cell pred for: ', ind[[i]])
    
    # Load the relevant pred object - loads an object named cell_pred
      rdata_file <- paste0('/share/geospatial/mbg/', ind_gp[[i]], '/', ind[[i]], '/output/', rd[[i]], '/',
                           ind[[i]], 
                           ifelse(rk, "_raked_cell_draws_eb_bin0_","_cell_draws_eb_bin0_"), 
                           reg, "_0.RData")
      
    #TODO improve logic
    if (rk) {
      if (file.exists(rdata_file)) {
        load(rdata_file)
        cell_pred <- raked_cell_pred
        rm(raked_cell_pred)
      }
    } else {
      if (file.exists(rdata_file)) {
        load(rdata_file)
      }
    }
    
    # Check to make sure loaded correctly
    if (!exists("cell_pred")) stop("Unable to load raked cell pred object!")
    
    # If extra columns at front of cell_pred, can skip here
    #TODO what is this for
    if(!(is.null(skip_cols))) cell_pred <- as.matric(cell_pred[, (skip_cols+1):ncol(cell_pred)])
    
    # Verify alignment  
    if(nrow(cell_pred)!=nrow(covdt)) stop('Dimensional mismatch between cell pred and simple raster!!')
    
    # set cell pred as a data table, and rename things
    cell_pred <- prep_cell_pred(
      cell_pred = cell_pred,
      cell_ids = link_table[[2]],
      pixel_id = pixel_id,
      covdt = covdt
    )
    
    message('~~>reshape long and formatting')
    #reshape long draws and key on the pixel ID/draw/year
    dt <- melt(cell_pred,
               measure = patterns("V"),
               variable.name = "draw",
               value.name = var_names[[i]]) %>% 
      #convert draw col to int instead of V1-250 as a factor
      .[, draw := substring(as.character(draw), 2) %>% as.integer] %>% #TODO probably can be done in the reshape?
      setkey(., pixel_id, draw, year, cell_pred_id, cell_id, pop) %>% #TODO could we get rid of cell_pred_id earlier?
      return
    
  }
  
  #load/format all the cell preds and then merge them together
  cell_pred <- lapply(1:length(ind), loadCellPreds) %>% 
    Reduce(function(...) merge(..., all = TRUE), .)
  
  #TODO should add a test here to verify that the resulting nrow = original cell_pred/num_yrs/250
  #tested it now and it looks OK but this is important to make more robust

  # merge cell_pred on the link
  # TODO note that pixel_id col in link dt is a duplicate of ID, and causes a merge issue (pixel_id.x % pixel_id.y)
  # eventually should fix this issue upstream but for now removing it pre-merge is sufficient
  cell_pred <- merge(link[, -c('pixel_id')], cell_pred, by.x = "ID", by.y = "cell_id", allow.cartesian = TRUE)
  
  # space
  link <- NULL

  #convert to fractional population 
  #TODO for now i will leave this turned off so user can do this adjustment explicitly
  #cell_pred = cell_pred[,pop := pop * area_fraction]

  #subset to relevant columns and return
  keep_vars <- c('ADM0_CODE', 'ADM1_CODE', 'ADM2_CODE', 
                 'pixel_id', 'year', 'pop', 'area_fraction', 'draw', unlist(var_names))
  cell_pred[, (keep_vars), with=F] %>% 
    return
  
}


#format aggregated results files at admin 0/1/2 lvls
format_admin_results <- function(ind_gp,
                                 ind,
                                 rd,
                                 measure,
                                 suffix,
                                 var_names = ind, # name using ind by default, but can pass custom name
                                 rk) {

  #helper function to load a list of admin results and format them into a single DT
  load_admins <- function(i) { 
    
    message('~>loading results for: ', ind[[i]])

  # def directory, then read in all the admin objects from a single RData file
    combined_file <- file.path('/share/geospatial/mbg', ind_gp[[i]], ind[[i]], 'output', rd[[i]]) %>% 
      paste0(., '/', ind[[i]], 
             ifelse(rk[[i]], '_raked', '_unraked'), 
             measure[[i]], '_admin_draws', suffix[[i]], '.RData')
    
    #helper function to create the combined file if it is not present
    create_combined_results <- function(file, sfx='_0.RData') {

      message('combined file not present...building')

      #find all relevant files
      files <- list.files(gsub(basename(file), '', file), pattern = 'admin_draws', full.names = T) %>% 
        .[grep(sfx, .)]
      
      load_specific_obj <- function(file, obj) {
        
        message('loading ->', obj, ' from: ', file)
        
        #loads an RData file, and returns the requested object
        load(file)
        ls()[ls() == obj] %>% 
          get %>% 
          return
        
      }
      
      #load each of the required objects from all files
      admin_0 <- lapply(files, load_specific_obj, obj='admin_0') %>% rbindlist
      admin_1 <- lapply(files, load_specific_obj, obj='admin_1') %>% rbindlist
      admin_2 <- lapply(files, load_specific_obj, obj='admin_2') %>% rbindlist
      sp_hierarchy_list <- lapply(files, load_specific_obj, obj='sp_hierarchy_list') %>% rbindlist
      
      #return objects in a named list (we can assign this to the parent env using list2env)
      list('admin_0'=admin_0, 'admin_1'=admin_1, 'admin_2'=admin_2, 'sp_hierarchy_list'=sp_hierarchy_list) %>% 
        return()
      
    }

    #load the combined file (create from all files if has not already been created)
    if(combined_file %>% file.exists) load(combined_file, verbose=T)
    else combined_file %>% create_combined_results %>% list2env(., .GlobalEnv) 
    
  # harmonize the hierarchy lists (some are using factor variables which don't merge well)
    factor_vars <- names(sp_hierarchy_list)[vapply(sp_hierarchy_list, is.factor, c(is.factor=FALSE))]
    if (length(factor_vars) > 0) {

      message('formatting spatial hierarchy table to num/chr - input table uses factor variables')
      
      #build helper functions
      format_hierarchy <- function(dt) {
        
      out <- copy(dt)
      
      #helper fx to convert factors without information loss
      #https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information
      facToNum <- function(f) as.numeric(as.character(f))

      #we want to convert the codes to num and the names to chr
      cols_to_num <- factor_vars[factor_vars %like% 'CODE']
      cols_to_chr <- factor_vars[factor_vars %like% 'NAME']
      cols <- list(cols_to_num, cols_to_chr)
      funs <- rep(c(facToNum, as.character), lengths(cols))

      # convert class based on column name
      out <- out[, unlist(cols) := Map(function(f, x) f(x), funs, .SD), .SDcols = unlist(cols)] %>% 
        return
      
      }
      
      #convert vars
      sp_hierarchy_list <- format_hierarchy(sp_hierarchy_list)

      #reassess and test
      #TODO add tests to ensure no information loss?
      if(vapply(sp_hierarchy_list, is.factor, c(is.factor=FALSE)) %>% any) stop('Failed to convert sp hierarchy!')

    }
    
  # format and append
    bind_results <- function(input_dt, info) {

      dt <- copy(input_dt)
      
      #pull out the level of aggregation, rename the variable to harmonize and record the value for later
      lvl_str <- names(dt)[names(dt) %like% 'CODE']
      lvl <- substr(lvl_str, start=1, stop=4) #extract level value
      setnames(dt, lvl_str, 'code') # rename
      dt[, agg_level := lvl] #record the level
      dt[, c('pop', 'region') := NULL] #remove unecessary vars

      #merge on location names while simultaneously formatting them
      out <- merge(dt,
                   info[, .(code=get(lvl_str), 
                            name=paste0(lvl, '_NAME') %>% get)] %>% unique, 
                  by='code',
                  all.x=T)

      #melt and output
      melt(out,
           measure = patterns("V"),
           variable.name = "draw",
           value.name = var_names[[i]]) %>% 
        return
      
    }
    
    dt <- list(admin_0, admin_1, admin_2) %>% 
      lapply(., bind_results, info=sp_hierarchy_list) %>% 
      rbindlist %>% 
      return
    
  }

  #load/format all the admin results and then merge them together
  dt <- lapply(1:length(ind), load_admins) %>% 
    Reduce(function(...) merge(..., all = TRUE), .) %>% 
    return
  
}


# -------------------------------------------------------------------
# Helper function to turn a tif raster file into a dt
raster_to_dt <- function(the_raster,
                         simple_polygon,
                         simple_raster,
                         year_list,
                         interval_mo,
                         outputdir,
                         pixel_id) { #will pull names from raster by default but user can pass in

  
  message(paste0("Prepping the ", names(the_raster), " raster for this region"))
  ## extend and crop pop raster to ensure it matches the simple raster #not convinced this is totally needed
  out  <- extend(the_raster, simple_raster, values = NA)
  out  <- crop(out, extent(simple_raster))
  out  <- setExtent(out, simple_raster)
  out  <- raster::mask(out, simple_raster)
  
  ## check to ensure the pop raster matches the simple raster in extent and resolution
  if (extent(out) != extent(simple_raster)) {
    stop("raster extent does not match simple raster")
  }
  if (any(res(out) != res(simple_raster))) {
    stop("raster resolution does not match simple raster")
  }
  
  #ensure the dimensions are the same
  stopifnot(dim(out)[1:2] == dim(simple_raster)[1:2])
  
  message("converting the raster into a data.table")
  #convert to datables, reshape and stuff
  brick_to_dt = function(bbb, pixel_id = pixel_id){
    dt = setDT(as.data.frame(bbb))
    dt[, pxid := .I] #probably uncessary
    
    #drop rows now in cellIdx
    dt = dt[pixel_id,]
    
    dt = melt(dt, id.vars = 'pxid', variable.factor = F)
    dt = dt[,.(value)]
    return(dt)
  }
  
  dt <- brick_to_dt(bbb = out, pixel_id = pixel_id) %>% 
    setnames(., names(the_raster))
  
  # Add pixel_id, but make sure that its recycled explicitly as per data.table 1.12.2 guidelines
  dt[, pixel_id := rep(pixel_id, times = nrow(dt) / length(pixel_id))]
  
  #add year to covdt
  yyy = as.vector(unlist(lapply(min(year_list):max(year_list), function(x) rep.int(x, times = length(pixel_id)))))
  dt[,year := yyy]
  
  return(dt)
  
}

#TODO write documentation
format_rasters <- function(these_rasters,
                           reg,
                           measure,
                           pop_measure,
                           covs = NULL,
                           cov_measures = NULL,
                           year_start,
                           year_end,
                           var_names = sapply(these_rasters, names), # name using rasters by default, but can pass custom name
                           matrix_pred_name = NULL,
                           skip_cols = NULL,
                           rk = T,
                           coastal_fix = T, # if model wasn't run w new coastal rasterization, force old simple raster process 
                           rake_subnational = rk, # TODO is this correct? might need to be defined custom by country
                           shapefile_version = 'current') {
  
  message('loading simple raster & populations')
  
  ## Load simple polygon template to model over
  gaul_list <- get_adm0_codes(reg, shapefile_version = shapefile_version)
  simple_polygon_list <- load_simple_polygon(
    gaul_list = gaul_list, buffer = 1, tolerance = 0.4,
    shapefile_version = shapefile_version
  )
  subset_shape <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  
  ## Load list of raster inputs (pop and simple)
  if (coastal_fix) raster_list <- build_simple_raster_pop(subset_shape, link_table = shapefile_version) #uses new simple_raster 
  else raster_list <- build_simple_raster_pop(subset_shape, link_table = NULL) #uses old rasterize
  
  simple_raster <- raster_list[["simple_raster"]]
  pop_raster <- raster_list[["pop_raster"]]
  pixel_id <- seegSDM:::notMissingIdx(simple_raster)
  
  ## get number of years
  year_list <- c(year_start:year_end)
  num_yrs <- length(year_list)
  
  message('loading links')
  #####################################################################
  # load the cell id to admin units link
  link_table <- get_link_table(simple_raster, shapefile_version = shapefile_version)

  #####################################################################
  #turn the rasters into a DT and merge them together
  dt <- lapply(these_rasters, raster_to_dt, 
               simple_polygon, simple_raster, year_list, interval_mo=12, pixel_id = pixel_id) %>% 
    Reduce(function(...) merge(..., all = TRUE), .) %>% 
    #force names, auto-extracting from raster can be variable depending on the upload name of the cell pred obj
    setnames(., c('pixel_id', 'year', var_names))
  
  #also merge on population
  dt <- load_populations_cov(reg, pop_measure=pop_measure, measure=measure, simple_polygon, 
                              simple_raster, year_list, interval_mo=12, pixel_id = pixel_id) %>% 
   merge(dt, ., by=c('pixel_id', 'year'), all.x=T) #TODO is it possible to have missing pop values?

  #also load/merge on any user-provided covariates
  if (!is.null(covs)) {
    
    #load the covariates as a raster
    dt <- load_and_crop_covariates_annual(covs = covs,
                                          measures = cov_measures,
                                          simple_polygon = simple_polygon,
                                          start_year  = min(year_list),
                                          end_year    = max(year_list),
                                          interval_mo = 12,
                                          agebin = 1) %>% 
      #convert to DT and combine
      lapply(., raster_to_dt, simple_polygon, simple_raster, year_list, interval_mo=12, pixel_id = pixel_id) %>% 
      Reduce(function(...) merge(..., all = TRUE), .) %>% 
      #force names, auto-extracting from raster can be variable depending on the upload name of the cell pred obj
      setnames(., c('pixel_id', 'year', covs)) %>% 
      #merge to the input rasters DT
      merge(dt, ., by=c('pixel_id', 'year'), all.x=T) #TODO is it possible to have missing covariate values?
  
  }

  #####################################################################
  # Prepping the cell_pred and link table to be linked by making sure they have the appropriate identifiers.  Also performs a
  # zippering at the region boundary where cells that have a portion of their area outside of the modeling region are reintegrated
  # as a whole cell and considered to be only in one region.  This works becasue a given cell is only modeled in one region.
  link <- prep_link_table(
    link_table = link_table,
    simple_raster = simple_raster,
    pixel_id = pixel_id
  )
  
  # getting the connector for sub-national or national raking, This connector gets the IHME location_code for our
  # gbd targets and connects that to the ADM0_CODE or ADM1_CODE as nessecary
  connector <- get_gbd_locs(
    rake_subnational, 
    reg = reg,
    shapefile_version = shapefile_version
  )
  
  # merge the connector on to the link table by making sure that each cell fragment gets connected to the appropriate
  # raking geography
  link <- sub_nat_link_merge(
    rake_subnational,
    link,
    connector
  )
  
  #now create the requisite IDs for merging
  dt[, cell_pred_id := .I] #cell_pred object ID
  dt[,cell_id := rep(link_table[[2]], times = nrow(dt) / length(link_table[[2]]))]  #cell id references the africa map
  dt[,pixel_id := rep(pixel_id, times = nrow(dt) / length(pixel_id))] #pixel id references the regional map  

  #TODO should add a test here to verify that the resulting nrow = original cell_pred/num_yrs/250
  #tested it now and it looks OK but this is important to make more robust
  
  # merge cell_pred on the link
  # TODO note that pixel_id col in link dt is a duplicate of ID, and causes a merge issue (pixel_id.x|pixel_id.y)
  # eventually should fix this issue upstream but for now removing it pre-merge is sufficient
  out <- merge(link[, -c('pixel_id')], dt, by.x = "ID", by.y = "cell_id", allow.cartesian = TRUE)

  #subset to relevant columns and return
  #note that this is why we needed to force the cov names
  keep_vars <- c('ADM0_CODE', 'ADM1_CODE', 'ADM2_CODE', 
                 'pixel_id', 'year', 'pop', 'area_fraction', unlist(var_names), covs)
  out[, (keep_vars), with=F] %>% 
    return
  
}