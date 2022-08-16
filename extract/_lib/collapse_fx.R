# ----HEADER------------------------------------------------------------------------------------------------------------
#***********************************************************************************************************************

# ----Cleaning----------------------------------------------------------------------------------------------------------
#function to do some initial cleaning and prep for the data
initialClean <- function(input.dt, var.fam) {
  
  message("\nBegin Initial Cleaning...")
  message('->Subset to relevant variables')

  if (var.fam == 'cooking') {
    dt <- input.dt[, .(nid, iso3, lat, long, survey_series, hhweight, urban, hh_size, 
                       year_start, int_year, shapefile, location_code,
                       cooking_fuel_mapped, cooking_location_mapped, cooking_type_mapped, cooking_type_chimney_mapped)]
  } else if (var.fam == 'housing') {
    dt <- input.dt[, .(nid, iso3, lat, long, survey_series, hhweight, urban, hh_size, 
                       year_start, int_year, shapefile, location_code,
                       housing_roof, housing_wall, housing_floor,
                       housing_roof_num, housing_wall_num, housing_floor_num)]
  } 

  problem_list <- dt[hh_size <=0]
  ### Standardize iso3s
  setnames(dt, 'iso3', 'ihme_loc_id')
  dt[, iso3 := substr(ihme_loc_id, 1, 3)]
  
  ### Standardize years
  message('-->Create a unique cluster id')
  
  #determine if working with point data
  is.polygon <- all(dt[, lat] %>% unique %>% is.na)

  if (!is.polygon) dt[, cluster_id := .GRP, by=.(iso3, lat, long, nid, year_start)]
  else dt[, cluster_id := .GRP, by=.(iso3, shapefile, location_code, nid, year_start)]
  
  #cluster_id should be unique within NID
  dt[, cluster_id := paste0(nid, '_', cluster_id)]

  dt[, raw_weight := hhweight] #save for later aggregation
  
  ### Standardize point data weighting
  if (!is.polygon) {
    
    message('---->Change weight to 1 for point data')
    message('----->Change shapefile/location_code to NA for point data')
    

    dt[, hhweight := 1]
    dt[, c('shapefile', 'location_code') := NA]
    
  } 
  
  dt[, polygon := is.polygon]
  
  #define a row id for other operations and key on it
  dt[, row_id := .I] 
  setkey(dt, row_id) %>%
    return

}

#***********************************************************************************************************************

# ----Definitions-------------------------------------------------------------------------------------------------------
#function to convert extracted variables based on current definitions
defIndicator <- function(dt, var.fam, definitions, debug=F, clean_up=T) {
  
  message("\nBegin Indicator Definition...")
  
  #allow for interactive debugs
  if (debug) browser()

  #read in the appropriate definitions file
  def.dt <- read_xlsx(path=definitions, sheet=var.fam) %>% as.data.table
  def.dt[, notes := NULL] #cleanup for reshaping
  
  #define a function to merge on definitions for each variable
  mergeDef <- function(data, defs, this.var, impute.vars=NULL) {

    #set the original row count in order to alert if anything goes wrong in remapping merges
    og.count <- nrow(data)
    
    defs <- defs[variable==this.var] #subset defs to working var

    #reshape the definitions and capture the new variable names for output
    new.cols <- names(defs)[!(names(defs) %in% c('variable', 'value'))]
    defs <- dcast(defs, value~variable, value.var=new.cols, sep='_')
    new.cols <- names(defs)[!(names(defs) %in% c('variable', 'value'))] #account for new var names (added stub)

    #remap
    #coerce to character because sometimes NA vars = logical & this cannot merge to the codebook
    data[, (this.var) := get(this.var) %>% as.character] 
    if (this.var %in% impute.vars) data[get(this.var) %>% is.na | get(this.var) == "", c(this.var) := 'unknown' ]
    out <- merge(data, defs, by.x=this.var, by.y='value', all.x=T) #merge onto data using original values
    
    #assert that merge was successful
    if (nrow(out) != og.count) stop(this.var, ': remap is causing data loss, investigate!!')
    else out[, c('row_id', new.cols), with=F] %>% setkey(., row_id) %>% return #keep only the new vars and ID

  }
  
  maps <- lapply(unique(def.dt$variable), mergeDef, data=dt, defs=def.dt, impute.vars=impute.vars)

  #merge the new mapped vars with the original data (note that they all need to be keyed on row ID)
  out <- Reduce(function(...) merge(..., all = T), list(dt, maps))
  
  #generate final indicators based on the intermediate vars
  if (var.fam == 'cooking') {

    #capture unknown/missing values
    unk.vals <- def.dt[is.na(ord), value] %>% 
      c(., '') #add in a missing str value too (different than NA)
        
    #generate categorical cooking fuel
    #pull the fuel types we are interested in from the definitions file
    fuel.types <- def.dt[variable=='cooking_fuel_mapped' , value] %>% unique %>% .[!(. %in% unk.vals)]
    
    #loop over each fuel type and generate the indicator variable
    for (type in fuel.types) {
      
      varname <- paste0('cat_cooking_fuel_', type)
      message('defining -> ', varname)
      #initialize, then fill based on cooking_fuel_mapped
      out[!(is.na(cooking_fuel_mapped) | cooking_fuel_mapped %in% unk.vals), (varname) := 0] 
      out[cooking_fuel_mapped==type, (varname) := 1]
      
    }
    
    #generate binary/ordinal cooking fuel
    message('defining -> cooking_fuel_solid/kerosene/dirty/clean')
    out[!(is.na(cooking_fuel_mapped) | cooking_fuel_mapped %in% unk.vals), 
        `:=` (cooking_fuel_solid=0, 
              cooking_fuel_kerosene=0, 
              cooking_fuel_dirty=0,
              cooking_fuel_clean=0)] #initialize, then fill based on cooking_fuel_mapped
    out[cooking_fuel_mapped %in% c('coal', 'wood', 'crop', 'dung'), cooking_fuel_solid := 1]
    out[cooking_fuel_mapped %in% c('kerosene'), cooking_fuel_kerosene := 1]
    out[cooking_fuel_mapped %in% c('coal', 'wood', 'crop', 'dung', 'kerosene'), cooking_fuel_dirty := 1]
    out[cooking_fuel_mapped %in% c('none', 'electricity', 'gas'), cooking_fuel_clean := 1]
 
    #cleanup intermediate vars
    remove.vars <- names(out)[names(out) %like% "row_id|mapped"]
    if (clean_up==T) out <- out[, (remove.vars) := NULL]
    
  }
  
  #redefine a row id for other operations and key on it
  out[, row_id := .I] 
  setkey(out, row_id) %>%
    return
}

#***********************************************************************************************************************

# ----Missingness-------------------------------------------------------------------------------------------------------
#function to identify missingness in key variables, with option for weighted %
idMissing <- function(input.dt, this.var, criteria=.2, wt.var=NA, check.threshold=F, threshold=NA, debug=F) {
  
  #allow for interactive debugs
  if (debug) browser()
  
  #set as copy so you dont save the new vars
  dt <- copy(input.dt)
  
  #set the original row count in order to print the data loss due to missingness
  og.count <- nrow(dt)
  
  #define the weight variable 
  #make sure that rows where the wtvar is NA are accounted for in this process
  if (wt.var %>% is.na) dt[, wt := 1]
  else dt[, wt := get(wt.var)] %>% .[is.na(wt), wt := 1]

  # Calculate data missingness by cluster for given variable
  # Alternatively check for validity against a minimum threshold (usually 0)
  if (!check.threshold) dt[, miss := lapply(.SD, is.na), .SDcols=this.var]
  else dt[, miss := lapply(.SD, function(x) x <= threshold), .SDcols=this.var]

  #calc pct miss, weight by selected variable (typically hh_size)
  dt[, pct_miss := sum(miss*wt, na.rm=T)/sum(wt, na.rm=T), by=cluster_id] 

  # Return the IDs of any clusters with > specified criteria weighted missingnesss or invalidity
  message("\nidentified #", dt[pct_miss>criteria, cluster_id] %>% uniqueN, " clusters", ' [',
          dt[pct_miss>criteria, nid] %>% uniqueN, ' nids]...or ~', 
          round((nrow(dt[pct_miss>criteria])/og.count)*100), 
          "% of rows \n based on criteria of >", criteria*100, "% ",
          ifelse(!check.threshold, 'missingness', 'invalidity'),
          " of ", this.var)
  
  clusters <- dt[pct_miss>criteria, cluster_id] %>% 
    unique
  
  #save a diagnostic file with the clusters and the type of missingness
  dt[, N := uniqueN(cluster_id), by=nid] %>% # for proportion dropped
    .[, .(nid, ihme_loc_id, int_year, cluster_id, N)] %>%
    setkey(., nid, ihme_loc_id, cluster_id) %>% 
    unique(., by=key(.)) %>% 
    .[, count := as.integer(cluster_id %in% clusters)] %>% 
    .[, count := sum(count), by=nid] %>% 
    .[, prop := count/N, by=nid] %>% 
    .[, var := this.var] %>% 
    .[, type := ifelse(!check.threshold, 'missingness', 'invalidity')] %>% 
    .[, cluster_id := NULL] %>% 
    #save using NID/vartype to uniquely ID the file, we will give a meaningful name to the combined diagnostic later
    unique(., by='nid') %>% 
    write.csv(., file=paste0(doc.dir, '/temp/dropped_clusters_', 
                             dt[, nid] %>% unique %>% sample(2, replace=T) %>% prod, #create a semi-random number
                             '_', this.var, '_', 
                             ifelse(check.threshold, 'invalidity', 'missingness'), '.csv'),
              row.names=F)
  
  if (length(clusters>0)) return(clusters)
  else return(NULL)
  
}

#***********************************************************************************************************************
 
# ----Aggregate---------------------------------------------------------------------------------------------------------
#aggregate the given indicator
aggIndicator <- function(input.dt, var.fam, debug=F) {
  
  #allow for interactive debugs
  if (debug) browser()
  
  #set the variables to work on based on the indicator family
  if (var.fam == 'cooking') {
    these.vars <- names(input.dt)[names(input.dt) %like% 'cooking']
    these.vars <- these.vars[!(these.vars %like% 'risk')] #dont collapse the continuous risk var
  } 
  
  message('collapsing...', paste(these.vars, sep='/'))

  #point data needs to be collapsed using urbanicity
  key.cols <- c('cluster_id', 'nid', 'lat', 'long', 
                'survey_series', 'year_median', 'shapefile', 'location_code', 'urban')

  #set as copy so you dont save the new vars
  dt <- copy(input.dt) %>%
    .[polygon==T, urban := NA] %>%   #for polygon data, urbanicity status is no longer meaningful
    setkeyv(., key.cols)

  #calculations
  dt[, (these.vars) := lapply(.SD, function(x, wt) sum(x*wt, na.rm=T)/sum(wt, na.rm=T), wt=hhweight*hh_size), 
     .SDcols=these.vars, by=key(dt)] #aggregate each variable and modify in place
  
  dt[, N := sum(hhweight*hh_size)^2/sum(hhweight^2*hh_size), by=key(dt)]
  dt[, sum_of_sample_weights := sum(hhweight), by=key(dt)]
  
  dt %>%  #collapse the dt to unique values based on the key and output
    unique(., by=key(dt)) %>% 
    return

}
#***********************************************************************************************************************

# ----Cleanup-----------------------------------------------------------------------------------------------------------
collapseCleanup <- function(var.fam, codebook, test.vars, cleanup=F, debug=F) {
  
  #allow for interactive debugs
  if (debug) browser()
  
  message('creating final diagnostics and cleaning up temp files for\n', var.fam)
  
  #combine the output diagnostics
  dt <- file.path(doc.dir, 'temp') %>% 
    list.files(., pattern='dropped_clusters', full.names = T) %>%
    lapply(., fread) %>% 
    rbindlist
  
  #also read in the final model dataset and collapse it to unique NIDs
  mod <- 
    list.files(model.dir, pattern='.fst', full.names = T) %>% 
    sort(., decreasing=T) %>% 
    .[1] %>% #pull most recent collapsed data 
    read.fst(., as.data.table=T) %>% 
    .[, nid] %>% 
    unique %>% 
    as.data.table %>% 
    setnames(., '.', 'nid') %>% 
    .[, mod_input := 1]
  
  #remove rows that are not codebooked for the test.vars using helper function
  codebook <- codebook[nid %in% unique(dt$nid)] #subset codebook to relevant rows we have collapsed
  checkCodebook <- function(this.var, full.dt, codebook) {
    
    dt <- full.dt[var==this.var] %>% copy #subset to comparison var
    
    #must account for differences in var specification in final model dt
    if (this.var %like% 'cooking_fuel') cb.var <- 'cooking_fuel' 
    else cb.var <- this.var
    
    #extract blank nids (NIDs in which the var has never been codebooked)
    dt[!(nid %in% codebook[is.na(get(cb.var)), unique(nid)])] %>% return
    
  } 
  
  dt <- lapply(test.vars, checkCodebook, full.dt=dt, codebook=codebook) %>% 
    rbindlist %>% 
    list(., dt[!(var %in% test.vars)]) %>% #note that we need to do this to prevent dropping hhweight rows (not in cb)
    rbindlist %>% 
    setkey(., nid, ihme_loc_id, int_year)
  
  #merge model data 
  dt <- merge(dt, mod, by='nid', all.x=T, fill=T) %>% 
    .[is.na(mod_input), mod_input := 0]
  
  #generate total missingness diagnostic and reshape wide
  dt[, newvar := paste0(substr(type, 1, 4), '_', var)] 
  dt <- dcast(dt[, -c('var', 'type'), with=F], ...~newvar, value.var=c('count', 'prop'), fun.aggregate=sum)

  #cleanup
  write.csv(dt, file=paste0(doc.dir, '/', var.fam, '/dropped_clusters.csv'), row.names = F)
  if (cleanup) file.path(doc.dir, 'temp') %>% list.files(., pattern='dropped_clusters', full.names = T) %>% unlink
  
  #output file
  return(dt)
  
}
#***********************************************************************************************************************

