# ----HEADER------------------------------------------------------------------------------------------------------------
#***********************************************************************************************************************

# ----CONFIG------------------------------------------------------------------------------------------------------------
# clear memory
rm(list=ls())

# runtime configuration
if (Sys.info()["sysname"] == "Linux") {
  j_root <- "/home/j/"
  h_root <- "/homes/jfrostad/"
  arg <- commandArgs()[-(1:3)] # First args are for unix use only
  
  if (length(arg)==0) {
    # arg <- c("IND", #current project iteration
    #          "8", #output version
    #          1) #number of cores provided to multicore functions
  }
  
  
  # necessary to set this option in order to read in a non-english character shapefile on a linux system (cluster)
  Sys.setlocale(category = "LC_ALL", locale = "C")
  
} else {
  j_root <- "J:"
  h_root <- "H:"
  # arg <- c("IND", #current project iteration
  #          "4", #output version
  #          1) #number of cores provided to multicore functions
}

#use cairo to render instead of quartz (quartz causes big slowdowns with geom_sf)
if(!identical(getOption("bitmapType"), "cairo") && isTRUE(capabilities()[["cairo"]])){
  options(bitmapType = "cairo")
}

# Use setup.R functions to load common LBD packages and mbg_central "function" scripts
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

#capture date
today <- Sys.Date() %>% gsub("-", "_", .)

#toggle aggregation if data has been updated
agg_results <- T

#options
run_date <- '2020_05_17_11_40_28' #first full run
run_date <- '2020_09_01_11_42_52' #collab submission
run_date <- '2020_10_04_22_20_57' #first submission
#run_date <- '2020_10_05_19_51_35'
lri_run_date <- '2020_06_11_11_19_26'
shapefile <- "2019_09_10"

indicator_group <- 'cooking'
indicator <- 'cooking_fuel_solid'

config_par <- 'hap_sp_fine'
cov_par <- 'cooking_VNM'

type <- 'mean'
raked <- T
start_year <- 2000
end_year <- 2018
cores <- 10

#mapping options
map_ind_gp <- 'hap'

#***********************************************************************************************************************

# ----IN/OUT------------------------------------------------------------------------------------------------------------
#***********************************************************************************************************************

# ---FUNCTIONS----------------------------------------------------------------------------------------------------------
##function lib##
#PE functions#
file.path(my_repo, '_lib', 'viz', 'mapping_fx.R') %>% source
file.path(my_repo, '_lib', 'post', 'landing_gear.R') %>% source

#gbd fx
gbd.shared.function.dir <- '/ihme/cc_resources/libraries/current/r/'
file.path(gbd.shared.function.dir, 'get_location_metadata.R') %>% source
file.path(gbd.shared.function.dir, 'get_age_metadata.R') %>% source
file.path(gbd.shared.function.dir, 'get_draws.R') %>% source
file.path(gbd.shared.function.dir, 'get_population.R') %>% source
#***********************************************************************************************************************

# ---PREP DATA----------------------------------------------------------------------------------------------------------
##read in and prep datasets for analysis##
## Read config file and save all parameters in memory
config <- set_up_config(repo            = my_repo,
                        indicator_group = indicator_group,
                        indicator       = indicator,
                        config_name     = paste0('/model/configs/config_', config_par),
                        covs_name       = paste0('/model/configs/covs_', cov_par),
                        run_tests       = F,
                        post_est_only   = T,
)

#read in link_table
global_link_table <- file.path(global_link_dir, "lbd_full_link.rds") %>% readRDS %>% as.data.table
adm_links <- global_link_table[, .(ADM0_NAME, ADM0_CODE, ADM1_NAME, ADM1_CODE, ADM2_NAME, ADM2_CODE)] %>% unique
admins <- get_sp_hierarchy(shapefile_version = modeling_shapefile_version)

#read in shps
adm2 <- rbind(stage1, stage2)

#read in mbg region info

#read in gbd location info
locs <- get_location_metadata(location_set_id = 35, gbd_round_id = 6) %>% 
  .[, .(iso3=ihme_loc_id, location_name, super_region_id, super_region_name, region_id, region_name, location_id)] #subset to relevant columns

#create file to crosswalk AD0 to iso3
iso3_map <- dplyr::select(adm2, iso3, ADM0_CODE=gadm_geoid) 
iso3_map$geometry <- NULL
iso3_map <- as.data.table(iso3_map) %>% unique
locs <- merge(locs, iso3_map, by='iso3')
adm_links <- merge(adm_links, locs, by=c('ADM0_CODE'), all.x=T)

if(agg_results) {
  
  message('producing aggregated results')

  #combine and save all ad2 level results
  dt <-
    list.files(file.path(data.dir), pattern='ad2_draws.fst', full.names = T) %>% 
    lapply(., read_fst, as.data.table=T) %>% 
    rbindlist(use.names=T, fill=T) 
  
  #generate attributable LRI variables
  #read LRI rates per 1000/counts
  load(file.path(lri_dir, lri_run_date, lri_rate_path), verbose=T)
  lri_dt <- admin_2 %>% 
    melt(.,
         measure = patterns("V"),
         variable.name = "draw",
         value.name = 'rate') %>% 
    .[, `:=` (cause='lri', grouping='under5', pop=NULL, region=NULL)]
  
  load(file.path(lri_dir, lri_run_date, lri_counts_path), verbose=T)
  lri_dt <- admin_2 %>% 
    melt(.,
         measure = patterns("V"),
         variable.name = "draw",
         value.name = 'count') %>% 
    merge(., lri_dt, by=c('year', 'ADM2_CODE', 'draw')) %>% 
    .[, `:=` (cause='lri', grouping='under5', pop=NULL, region=NULL)]
  
  #pull the draws from codcorrect for BRA/CHN
  locations <- get_location_metadata(location_set_id = 35, gbd_round_id = 6, decomp_step = "step4") %>% as.data.table
  loc_ids <- locations[ihme_loc_id %like% 'BRA|CHN' & most_detailed==1, location_id]
  gbd <-
    get_draws("cause_id", 322, 
              source="codcorrect", version=135, year_id=1990:2018,
              location_id=loc_ids, sex_id=3, age_group_id=1, gbd_round_id=6, decomp_step="step5",
              metric_id=c(1,3), measure_id=1,
              num_workers=8) %>% 
    melt(.,
         measure = patterns("draw_"),
         variable.name = "draw",
         value.name = 'state_count') %>% 
    .[, `:=` (cause='lri', grouping='under5', pop=NULL, region=NULL)] %>% 
    .[, draw := gsub('draw_', 'V', draw)]
  
  #calculate rates
  pop_df <- get_population(decomp_step="iterative", gbd_round_id=6,
                           location_id=loc_ids,
                           year_id=1990:2018,
                           age_group_id=1,
                           sex_id=3)
  gbd <- merge(gbd, pop_df[, .(location_id, year_id, population)], by=c('location_id', 'year_id'))
  gbd[, rate := state_count / population]
  
  # add connector object between ADM codes and GBD loc ids
  gbd <- merge(gbd, 
               data.table(get_gbd_locs(rake_subnational = T,
                                       reg = "BRA+CHN",
                                       shapefile_version = raking_shapefile_version)), 
               by='location_id')
  gbd <- merge(gbd, adm_links[ADM0_NAME%in%c('Brazil','China'), .(ADM1_CODE, ADM2_CODE)], by='ADM1_CODE', allow.cartesian=T)
  setnames(gbd, 'year_id', 'year')
  
  #merge LRI results to HAP results
  dt <- merge(dt, lri_dt, by=c('ADM2_CODE', 'year', 'draw', 'grouping', 'cause'), all.x=T) 
  
  #impute BRA/CHN LRI figures from GBD, since LBD has no data and does not produce estimates
  na_dt <- dt[is.na(rate) | is.na(count)] #BRA/CHN
  na_dt[, c('rate', 'count') := NULL]
  na_dt <- merge(na_dt, gbd, by=c('ADM0_CODE', 'ADM2_CODE', 'year', 'draw', 'grouping', 'cause'), all.x=T) 
  na_dt[, pop_share := pop/sum(pop, na.rm=T), by=.(year, ADM1_CODE, draw, type)] #calc dist share of state pop
  na_dt[, count := population * pop_share * rate] #use share of children in ad2 to assign count by GBD u5 pop for ad1
  
  #combine imputations
  dt <- list(
    na_dt[, names(dt), with=F], #imputed BRA/CHN
    dt[!(is.na(rate) | is.na(count))] #all others
  ) %>% rbindlist(use.names=T)
  
  #calculate attrib. using paf and LRI figures
  dt[, atr_rate := rate * paf]
  dt[, atr_count := count * paf]
  
  #merge sr region names/IDs/ADM1_CODES
  dt <- merge(dt, locs, by='ADM0_CODE', all.x=T)

  #save draws for custom draw level calculations
  write.fst(dt[year%in%c(2000, 2018)], file.path(data.dir, 'all_draws.fst'))
  
  #define columns to summarize
  ind_cols <- c('share', 'pm_pc', 'prev', 'paf', 'rate', 'count', 'atr_rate', 'atr_count')
  
  #define list of aggregations
  by_cols <- list('global'=c('year', 'type', 'cause', 'grouping'),
                  'super_region'=c('year', 'type', 'cause', 'grouping', 'super_region_id'),
                  'region'=c('year', 'type', 'cause', 'grouping', 'region_id'),
                  'ad0'=c('year', 'type', 'cause', 'grouping', 'ADM0_CODE'),
                  'ad1'=c('year', 'type', 'cause', 'grouping', 'ADM0_CODE', 'ADM1_CODE'),
                  'ad2'=c('year', 'type', 'cause', 'grouping', 'ADM0_CODE','ADM2_CODE')
                  )
  
  #summarize and produce aggregations
  results <- mclapply(1:length(by_cols), 
                      calcSummaryStats,
                      by_list=by_cols,
                      dt=dt, 
                      ind_cols=ind_cols,
                      year_list=start_year:end_year, 
                      mc.cores=1) %>% 
    rbindlist(use.names=T, fill=T)
  
  #save all aggregations
  write.fst(results, file.path(data.dir, 'all_ad2_summary.fst'))
  
  #helper function to save each aggregation type while cleaning up any irrelevant variables
  saveResults <- function(agg, dt) {
    
    message('saving ', agg, ' results')
    
    dt[dimension==agg] %>% 
      Filter(function(x) !all(is.na(x)), .) %>% 
      write.csv(., paste0(data.dir, '/', agg, '_summary.csv'), row.names = F)
    
    return(NULL)
    
    
  } 
  
  lapply(by_cols %>% names, saveResults, dt=results)
  
} else results <- file.path(data.dir, 'all_summary.fst') %>% read_fst(as.data.table=T)

#subset to ad2 for mapping estimates
dt <- results[dimension=='ad2' & term=='lvl']
dt_d <- results[dimension=='ad2' & term=='change_rate']

#create diagnostic to see if any AD2s are missing
missing_diagnostic <- copy(dt) %>% 
  merge(., adm_links[, .(N_districts=uniqueN(ADM2_CODE)), by=ADM0_CODE], by='ADM0_CODE') %>% 
  .[, .(pct_missing=1-(uniqueN(ADM2_CODE)/N_districts)), by=ADM0_CODE] %>% 
  unique(., by='ADM0_CODE') %>% 
  merge(., iso3_map, by='ADM0_CODE') %>% 
  .[order(pct_missing)] 
#***********************************************************************************************************************

# ---SEV CELL_PREDS-----------------------------------------------------------------------------------------------------
#read in the proper annotations (borders, lakes, mask)
check <- file.exists(annotations_path)
annotations <- ifelse(
  check,
  readRDS(annotations_path),
  load_map_annotations()
)
if(!check) saveRDS(annotations, file=annotations_path)

#append the ADM2 files
sdg_files <-
file.path(data.dir, 'sdg_projections') %>% list.files(pattern='admin_2', full.names = T) %>% 
  lapply(., readRDS)

#extract goal obj to index over
goals <- lapply(1:length(sdg_files), function(x) {
    sdg_files[[x]]$goals
  }) %>% rbindlist %>% unique

#create a dt with all probabilities
probs <- lapply(1:nrow(goals), prepCasts, type='probs', id_dt=goals) %>% 
  rbindlist(fill=T, use.names=T) %>% 
  setnames(c('target_year', 'spatial_idx'), c('year', 'ADM2_CODE')) %>% 
  merge(., adm_links[, .(ADM2_CODE, ADM0_CODE, ADM0_NAME)], by='ADM2_CODE')  

#create a dt with all projections
projs <- 
  lapply(c(2018, seq(2020, 2030, 5)), prepCasts, type='proj', id_var='year') %>% 
  rbindlist(fill=T, use.names=T) %>% 
  setnames('spatial_idx', 'ADM2_CODE') %>% 
  melt(measure = patterns("V"), variable.name = "draw", value.name='sev')

#create a dt with aroc and combine
projs <- prepCasts(2018, type='aroc', id_var='year') %>% 
  melt(measure = patterns("V"), variable.name = "draw", value.name='aroc') %>% 
  merge(., projs, by=c('ADM2_CODE', 'year', 'draw'), all.y=T) %>% 
  merge(., adm_links[, .(ADM2_CODE, ADM0_CODE)], by='ADM2_CODE')

#generate mean/ci
cols <- c('aroc', 'sev')
setkey(projs, year, ADM2_CODE)
projs[, paste0(cols, '_mean') := lapply(.SD, mean, na.rm=T), .SDcols=cols, by=key(projs)]
projs[, paste0(cols, '_lower') := lapply(.SD, quantile, probs=.025, na.rm=T), .SDcols=cols, by=key(projs)]
projs[, paste0(cols, '_upper') := lapply(.SD, quantile, probs=.975, na.rm=T), .SDcols=cols, by=key(projs)]
projs <- unique(projs, by=key(projs)) %>% 
  .[, c(cols, 'draw') := NULL]

#calculate relative uncertainty of SEV
projs[, sev_rel_uncertainty := (sev_upper-sev_lower)/2/sev_mean]
projs[sev_rel_uncertainty>1, sev_rel_uncertainty := 1] #cap at 1

#***********************************************************************************************************************
 
# ---SAVE MAPPING INPUTS------------------------------------------------------------------------------------------------
#output files for Kim to produce key figures
##Figure 1: Selected AD2 Results for 2018##
#A: SFU levels
saveMappingInput(
  dt[type=='HAP'], 
  map_ind='sfu',
  data_ind='prev',
  map_measure='mean',
  data_measure='mean'
)

#B: TAP PC levels
saveMappingInput(
  dt[type=='TAP'], 
  map_ind='tap_pc',
  data_ind='pm_pc',
  map_measure='mean',
  data_measure='mean'
)

#C: HAP PCT levels
saveMappingInput(
  dt[type=='HAP'],  
  map_ind='hap_pct',
  data_ind='share',
  map_measure='mean',
  data_measure='mean'
)

#D: Attributable LRI rates
saveMappingInput(
  dt[type=='TAP'], 
  map_ind='tap_lri',
  data_ind='atr_rate',
  map_measure='mean_rate',
  data_measure='mean'
)

##Figure 2: Selected AD2 Trends 2000-2018##
#A: DFU levels
saveMappingInput(
  dt_d[type=='HAP'&term=='change_rate'], 
  map_ind='sfu',
  data_ind='prev',
  map_measure='change_rate',
  data_measure='mean'
)

#B: Change rate for TAP_PC 2000-2018
saveMappingInput(
  dt_d[type=='TAP'&term=='change_rate'], 
  map_ind='tap_pc',
  data_ind='pm_pc',
  map_measure='change_rate',
  data_measure='mean'
)

#C: Change rate for attributable LRI 2000-2018
saveMappingInput(
  dt_d[type=='TAP'&term=='change_rate'], 
  map_ind='tap_lri',
  data_ind='atr_rate',
  map_measure='change_rate',
  data_measure='mean'
)

#D: SEV in 2018
saveMappingInput(
  projs, 
  map_ind='sev',
  data_ind='sev',
  map_measure='mean',
  data_measure='mean'
)

#E: Relative uncertainty of SEV in 2018
saveMappingInput(
  projs, 
  map_ind='sev',
  data_ind='sev',
  map_measure='relative_uncertainty',
  data_measure='rel_uncertainty'
)

#F: Probability of hitting SDG threshold
threshold <- 0.05
saveMappingInput(
  probs[target==threshold], 
  map_ind='sdg_prob',
  data_ind='absolute_goal',
  map_measure='mean',
  data_measure='prob'
)

#output new files for Kim to produce revised SI uncertainty figures (matt saved the other versions)
#A:M/U/L/SD for sfu
#impute with mean for brazil/china, which came from GBD and don't have SD
tomean <- names(results) %>% .[. %like% 'sd']
dt[, (tomean) := lapply(tomean, function(x) {
  x <- get(x)
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})]


saveMappingInput(
  dt[type=='HAP'], 
  map_ind='sfu',
  data_ind='prev',
  map_measure='sd',
  data_measure='sd'
)


#B: TAP PC levels
saveMappingInput(
  dt[type=='TAP'], 
  map_ind='tap_pc',
  data_ind='pm_pc',
  map_measure='sd',
  data_measure='sd'
)

#C: HAP PCT levels
saveMappingInput(
  dt[type=='HAP'],  
  map_ind='hap_pct',
  data_ind='share',
  map_measure='sd',
  data_measure='sd'
)

#D: Attributable LRI rates
saveMappingInput(
  dt[type=='TAP'], 
  map_ind='tap_lri',
  data_ind='atr_rate',
  map_measure='sd_rate',
  data_measure='sd'
)

#D: Attributable LRI counts
saveMappingInput(
  dt[type=='TAP'], 
  map_ind='tap_lri',
  data_ind='atr_count',
  map_measure='sd_count',
  data_measure='sd'
)
#***********************************************************************************************************************
 
# ---SCRAPS-------------------------------------------------------------------------------------------------------------
# plot(tmp[1], col = "red", main = "Neighbors: st_intersection")
# plot(nc[spare_mtx[[50]], "NAME"], add = TRUE)
# plot(nc[50, "AREA"], add = TRUE)
# 
# FIRSTdegreeNeighbors <- function(x) { 
#   
#   # -- use sf functionality and output to sparse matrix format for minimizing footprint 
#   first.neighbor <- st_touches(x, x, sparse = TRUE) 
#   
#   # -- convert results to data.frame (via sparse matrix)
#   n.ids <- sapply(first.neighbor, length)
#   vals <- unlist(first.neighbor) 
#   out <- sparseMatrix(vals, rep(seq_along(n.ids), n.ids))
# 
#   out.summ <- summary(out)  # -- this is currently only generating row values, need to map to actual obs [next line]
#   data.frame(country = x[out.summ$j,]$NAME, 
#              countryid = x[out.summ$j,]$CNTY_ID, 
#              firstdegreeneighbors = x[out.summ$i,]$NAME, 
#              firstdegreeneighborid = x[out.summ$i,]$CNTY_ID, 
#              stringsAsFactors = FALSE)
# }
# 
# x = st_join(tmp, tmp, join = st_rook)
# 
# sel = unlist(st_relate(tmp, tmp, pattern = "F***0****"))
# sel = st_touches(tmp, tmp) %>% unlist
# plot(st_geometry(tmp))
# plot(st_geometry(tmp[sel,]), add = TRUE, col = 'grey')


# proj_ad0 <- lapply(out, function(x) x[['ad0']]) %>% 
#   .[!sapply(., is.null)] %>% #remove the null tables (missing raster values)
#   rbindlist(use.names=T, fill=T) %>%
#   merge(., stages[, .(ADM0_CODE=gadm_geoid, mbg_reg, spr_reg_id)], by=c('ADM0_CODE')) %>% 
#   merge(., adm_links[, .(ADM0_CODE, ADM0_NAME)] %>% unique, by=c('ADM0_CODE')) %T>% 
#   write.csv(., file.path(data.dir, 'admin_0_sdg_summary.csv'), row.names = F)
# 
# proj_ad2 <- lapply(out, function(x) x[['ad2']]) %>% 
#   .[!sapply(., is.null)] %>% #remove the null tables (missing raster values)
#   rbindlist(use.names=T, fill=T) %>%
#   merge(., stages[, .(ADM0_CODE=gadm_geoid, mbg_reg, spr_reg_id)], by=c('ADM0_CODE')) %>% 
#   merge(., adm_links, by=c('ADM0_CODE', 'ADM2_CODE')) %T>% 
#   write.csv(., file.path(data.dir, 'admin_2_sdg_summary.csv'), row.names = F)

#***********************************************************************************************************************