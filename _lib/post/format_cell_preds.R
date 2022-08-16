# Helper function to turn a tif raster file into a dt
raster_to_dt <- function(the_raster,
                         simple_polygon,
                         simple_raster,
                         year_list,
                         interval_mo,
                         outputdir,
                         pixel_id) { #TODO will pull names from raster by default add user pass in functionality
  
  #some rasters are badly named with multiple slotnames, subset to first name only
  raster_name <- names(the_raster)[1]
  
  message(paste0("Prepping the ", raster_name, " raster for this region"))
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
    setnames(., raster_name)
  
  # Add pixel_id, but make sure that its recycled explicitly as per data.table 1.12.2 guidelines
  dt[, pixel_id := rep(pixel_id, times = nrow(dt) / length(pixel_id))]
  
  #add year to covdt
  yyy = as.vector(unlist(lapply(min(year_list):max(year_list), function(x) rep.int(x, times = length(pixel_id)))))
  dt[,year := yyy]
  
  return(dt)
  
}


link_cell_pred <- function(ind_gp,
                           ind,
                           rd,
                           reg,
                           measure,
                           pop_measure,
                           year_start,
                           year_end,
                           var_names = ind, # name using ind by default, but can pass custom name
                           n_draws = 250, #TODO automatically check the config file for 'samples' value?
                           #matrix_pred_name = NULL,
                           skip_cols = NULL,
                           rk = T,
                           coastal_fix = T, # if model wasn't run w new coastal rasterization, force old simple raster process 
                           rake_subnational = rk, # TODO is this correct? might need to be defined custom by country
                           shapefile_version = 'current',
                           covs=NA, #added fx to bring in covariates and merge
                           debug=F) {
  
  if(debug) browser()
  
  message('Beginning formatting for ', reg)
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
  
  if(pop_measure=='hap') {
    
    #helper function to load multiple worldpop measures
    load_populations_list <- function(subpop) {
      
      message('loading pops for measure=', subpop)
      
      # collect and load the population data from the WorldPop rasters
      out <- load_populations_cov(reg, pop_measure=subpop, measure = 'count', 
                                  simple_polygon=simple_polygon, simple_raster, 
                                  year_list, interval_mo=12, pixel_id = pixel_id) %>% 
        setnames('pop', subpop) %>% 
        return
      
    }
    
    #generate a populations dt with all the relevant hap groupings    
    #TODO test if its worth using a1564t+a65plt * gender_ratio (a1549m:a1549f) instead
    covdt <-
      lapply(c('a1549m', 'a1549f', 'a0514t', 'a0004t'), load_populations_list) %>% 
      Reduce(function(...) merge(..., all = TRUE), .) %>% 
      .[, .(pixel_id, year, #id.vars
            male=a1549m, 
            female=a1549f, 
            child=a0514t,
            under5=a0004t)] %>% 
      melt(id.vars=c('pixel_id', 'year'), value.name='pop', variable.name='grouping')
    
  } else {
    
    # collect and load the population data from the WorldPop rasters
    covdt <- load_populations_cov(reg, pop_measure, measure = 'count', simple_polygon, simple_raster, 
                                  year_list, interval_mo=12, pixel_id = pixel_id)
    
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
  
  #helper function to load a list of cell preds and then merge them together
  loadCellPreds <- function(i, 
                            max_n=n_draws) { 
    
    message('~>loading cell pred for: ', ind[[i]])
    
    # Load the relevant pred object - loads an object named cell_pred
    rdata_file <- paste0('/share/geospatial/mbg/', ind_gp[[i]], '/', ind[[i]], '/output/', rd[[i]], '/',
                         ind[[i]], 
                         ifelse(rk, paste0("_raked_", measure), ""), 
                         "_cell_draws_eb_bin0_", reg, "_0.RData")
    
    #TODO improve logic
    if (rk) {
      if (file.exists(rdata_file)) {
        message('loading raked results')
        load(rdata_file)
        cell_pred <- raked_cell_pred
        rm(raked_cell_pred)
      }
    } else {
      if (file.exists(rdata_file)) {
        message('loading UNraked results')
        load(rdata_file)
      }
    }
    
    # Check to make sure loaded correctly
    if (!exists("cell_pred")) stop("Unable to load raked cell pred object!")
    
    # If extra columns at front of cell_pred, can skip here
    #TODO what is this for
    if(!(is.null(skip_cols))) cell_pred <- as.matric(cell_pred[, (skip_cols+1):ncol(cell_pred)])
    
    # Verify alignment  
    # Noting that if we are working with HAP populations, the covdt will be 4x longer (4 groupings)
    if(nrow(cell_pred)!=nrow(covdt)/ifelse(pop_measure=='hap', 4, 1)) stop('Dims mismatch...cell pred!=simple raster!!')
    
    # Enforce max # of draws
    message('...subsetting to #', max_n, ' draws!')
    cell_pred <- cell_pred[, (1:max_n)]
    
    #set cell pred as a data table, and rename things
    cell_pred <- as.data.table(cell_pred)
    names(cell_pred) <- paste0('V',1:ncol(cell_pred))
    
    cell_pred[, cell_pred_id := .I] #cell_pred object ID
    cell_pred[, cell_id := rep(link_table[[2]], times = nrow(cell_pred) / length(link_table[[2]]))]  #cell id references the africa map
    cell_pred[, pixel_id := rep(pixel_id, times = nrow(cell_pred) / length(pixel_id))] #pixel id references the regional map  
    
    #generate year variable in order to merge on pops
    cell_pred[, year := 1:.N, by='pixel_id']
    cell_pred[, year := (year+min(year_list)-1)] #shift from idx to the actual year using min year in the provided list
    
    #merge population, year and potentially the stackers
    cell_pred <- merge(
      cell_pred, covdt, all.x=T, by=c('pixel_id', 'year'),
      allow.cartesian=pop_measure=='hap' # note that we will expand by the groupings variable if working w hap pops
    )
    
    #make sure it behaved
    stopifnot(any(!(cell_pred[,pixel_id] != rep.int(covdt[,pixel_id], length(year_start:year_end)))))
    
    return(cell_pred)
    
  }
  
  #load/format all the cell preds and then merge them together
  cell_pred <- lapply(1:length(ind), loadCellPreds) %>% 
    Reduce(function(...) merge(..., all = TRUE), .)
  
  #tested it now and it looks OK but this is important to make more robust
  
  # merge cell_pred on the link
  # eventually should fix this issue upstream but for now removing it pre-merge is sufficient
  cell_pred <- merge(link[, -c('pixel_id')], cell_pred, by.x = "ID", by.y = "cell_id", allow.cartesian = TRUE)

  #also load/merge on any user-provided covariates
  if (!is.null(covs)) {
    
    message('adding requested covariates, ', covs)

    message('Grabbing raster covariate layers')
    loader <- MbgStandardCovariateLoader$new(start_year = min(year_list),
                                             end_year = max(year_list),
                                             interval = 12,
                                             covariate_config = covs)

    cell_pred <- loader$get_covariates(simple_polygon) %>% 
      #convert to DT and combine
      lapply(., raster_to_dt, simple_polygon, simple_raster, year_list, interval_mo=12, pixel_id = pixel_id) %>% 
      Reduce(function(...) merge(..., all = TRUE), .) %>% 
      #force names, auto-extracting from raster can be variable depending on the upload name of the cell pred obj
      setnames(., c(covs$covariate %>% unique, 'pixel_id', 'year')) %>% 
      #merge to the input rasters DT
      merge(cell_pred, ., by=c('pixel_id', 'year'), all.x=T) #TODO is it possible to have missing covariate values?

    
  }
  
  #if using any pop measure than total, include total pop as well
  if (pop_measure!='total') {
    
    cell_pred <- load_populations_cov(reg, pop_measure='total', measure = 'count', simple_polygon, 
                                      simple_raster, year_list, interval_mo=12, pixel_id = pixel_id) %>% 
      setnames(., 'pop', 'pop_total') %>% 
      merge(., cell_pred, by=c('pixel_id', 'year'))
    
  }
  
  #drop irrelevant columns and return
  #TODO may be easier to keep columns instead =)
  null_vars <- c('ID', 'ADM_CODE', 'NAME_0', 'NAME_1', 'NAME_2', 'ADM0_NAME', 'ADM1_NAME', 'ADM2_NAME', 
                 'geo_id', 'ad2_id', 'ad0_parent', 'ad1_parent', 'location_id',
                 'start_area', 'end_area', 'total_area',
                 'n', 'reg_id', 'simp_ras_val', 'total_frac', 'old_frac', 'rak_level')
  cell_pred[, -(null_vars), with=F] %>% 
    return
  
}