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

## Set core_repo location and indicator group

#load packages
package_lib    <- sprintf('%s_code/_lib/pkg_lbd',h_root)
## Load libraries and  MBG project functions.
.libPaths(package_lib)
pacman::p_load(data.table, fst, scales, ggplot2, RColorBrewer, viridis, farver, magrittr, sf)
package_list <- c(package_list, 'sf')

# Use setup.R functions to load common LBD packages and mbg_central "function" scripts
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

#capture date
today <- Sys.Date() %>% gsub("-", "_", .)

#options
run_date <- '2020_05_17_11_40_28' #first full run
run_date <- '2020_09_01_11_42_52' #collab submission
run_date <- '2020_10_04_22_20_57' #first submission

lri_run_date <- '2020_06_11_11_19_26'
var_types <- c('lower', 'mean', 'upper')
new_gbd_estimates <- F

indicator_group <- 'cooking'
indicator <- 'cooking_fuel_solid'

config_par <- 'hap_sp_fine'
cov_par <- 'cooking_VNM'
type <- 'mean'
raked <- F
start_year <- 2000
end_year <- 2019
analysis_year <- 2018
cores <- 10
modeling_shapefile_version <- "2019_09_10"
#***********************************************************************************************************************

# ----IN/OUT------------------------------------------------------------------------------------------------------------
###Input###
#***********************************************************************************************************************

# ---FUNCTIONS----------------------------------------------------------------------------------------------------------
##function lib##
#PE functions#
file.path(my_repo, '_lib', 'viz', 'mapping_fx.R') %>% source

#gbd fx
gbd.shared.function.dir <- '/ihme/cc_resources/libraries/current/r/'
file.path(gbd.shared.function.dir, 'get_location_metadata.R') %>% source
file.path(gbd.shared.function.dir, 'get_covariate_estimates.R') %>% source
file.path(gbd.shared.function.dir, 'get_age_metadata.R') %>% source
file.path(gbd.shared.function.dir, 'get_draws.R') %>% source
file.path(gbd.shared.function.dir, 'get_outputs.R') %>% source
file.path(gbd.shared.function.dir, 'get_population.R') %>% source

##custom fx##

#helper fx to pull/prep the appropriate files from our list of SDG projection objects
prepCasts <- function(id, type, list=sdg_files, id_dt=NA, id_var=NA) {
  
  #format ID var if necessary
  if(nchar(id)==4) id <- as.character(id) #if the ID is a year, format as character
  
  #helper function to extract the correct object
  extractObj <- ifelse(type!='aroc',
                       function(x) list[[x]][[type]][[id]] %>% as.data.table,
                       function(x) list[[x]][[type]] %>% as.data.table ) #aroc only has one object

  #do the formatting and extractions
  lapply(1:length(list), extractObj) %>% 
    rbindlist(fill=T, use.names=T) %>% 
    { if(id_var %>% is.na) cbind(., id_dt[id,]) else .[, (id_var) := id] } %>% 
    return
  
}

#spatial functions
#return the max/min districts in a country and the distance between them
exploreRange <- function(explore_country, shp, dt, var, types=var_types, stub=NULL) {

 if (stub %>% is.null) new_vars <- paste(var, types, sep='_')
 else new_vars <- paste(var, types, stub, sep='_')
  
  out <- dt %>% 
    copy %>% 
    setnames(., new_vars, c('var_lower', 'var_mean', 'var_upper')) %>% 
    .[iso3==explore_country & !is.na(var_mean), 
      .(var_lower, var_mean, var_upper, ADM0_NAME, ADM1_NAME,  ADM2_NAME, pop_total)] %>% 
    .[order(var_mean)] %>% 
    .[c(1, nrow(.)), .(var_lower, var_mean, var_upper, ADM0_NAME, ADM1_NAME, ADM2_NAME, pop_total)]

  st_distance(filter(shp, NAME_2==out[1, ADM2_NAME]), 
              filter(shp, NAME_2==out[2, ADM2_NAME])) %>% 
    as.numeric %>% 
    round %>% 
    {if (length(.)>1) message('warning: multipolygon, returning min distance'); min(.)} %>% 
    message('\nDistance is ', ./1e3, ' km')
  
  return(out)
  
}

#function to find polygons that share a border (rook neighbors)
st_rook <- function(a, b = a) st_relate(a, b, pattern = "F***1****") 

tabulateR <- function(results_dt=results, 
                      lvl='ad2', years=analysis_year, types='HAP', terms='lvl',
                      ind='prev', metrics=var_types, stub='',
                      filter=NULL, #should be provided as list(var=,vals=)
                      thresh=NULL, byvar=NULL, popvar='pop_total', inv_threshold=F,
                      cleanup=T,
                      sorted='down') {
  
  #subset DT to requested results
  dt <- results_dt[dimension==lvl & year%in%years & type%in%types & term%in%terms]
  
  #define ind vars in order to remove irrelevant ones
  ind_vars <- names(dt) %>% .[. %like% paste(metrics, collapse='|')]
  irrel_vars <- ind_vars %>% .[!(. %like% ind)]
  rel_vars <- paste(ind, metrics, sep="_")
  mean_var <- rel_vars %>% .[. %like% 'mean']
  name_vars <- names(dt) %>% .[. %like% 'name|NAME']
  
  #remove irrelevant vars and missing rows
  dt[is.na(get(mean_var)), .N] %>% message('missing #', ., ' rows of ', mean_var)
  dt <- dt[!is.na(get(mean_var)), -c(irrel_vars), with=F] 
  setcolorder(dt, neworder=c('year', 'type', rel_vars, name_vars)) #reorder for legibility

  #filter if requested
  if(!is.null(filter)) { 
    message('filtering')
    if(filter$type) dt <- dt[get(filter$var) %in% filter$vals] #inclusive
    else dt <- dt[!(get(filter$var) %in% filter$vals)] #exclusive
  }
  
  #cleanup irrelevant vars if requested
  if(cleanup) dt <- Filter(function(x) !all(is.na(x)), dt)

  #test thresholds, if requested
  if(!is.null(thresh)) {
    
    message('testing threshold of ', ifelse(inv_threshold, 'less', 'greater'), ' than ', thresh)
    
    #helper function to build the thresholds
    testThreshold <- function(x, pop=NULL) {
      if(pop %>% is.null) scalar <- 1 else scalar <- pop
      if(inv_threshold) sum(scalar*(x < thresh))
      else sum(scalar*(x > thresh))
    } 
    
    #get counts
    dt[, paste('count', var_types, sep='_') := lapply(.SD, testThreshold), 
       .SDcols=paste('prev', var_types, sep='_'),
       by=list(byvar %>% get, year)]
    
    #get pop counts
    dt[, paste('pop', var_types, sep='_') := lapply(.SD, testThreshold, pop=get(popvar)), 
       .SDcols=paste('prev', var_types, sep='_'),
       by=list(byvar %>% get, year)]
    
    #get denoms
    dt[, N := .N, by=list(byvar %>% get, year)]
    dt[, pop_N := get(popvar) %>% sum, by=list(byvar %>% get, year)]
    
    #get pcts
    #district pct
    dt[, paste('pct', var_types, sep='_') := lapply(.SD, function(x) x/N), 
       .SDcols=paste('count', var_types, sep='_'),
       by=list(byvar %>% get, year)]
    #pop_pct
    dt[, paste('pop_pct', var_types, sep='_') := lapply(.SD, function(x) x/pop_N), 
       .SDcols=paste('pop', var_types, sep='_'),
       by=list(byvar %>% get, year)]
    
    #take unique
    dt <- unique(dt, by=c(byvar, 'year'))
    
    #sort by the population share
    mean_var <- 'pop_pct_mean'
    
  }
  
  #return table, sorted if requested
  if(sorted=='down') setorderv(dt, cols=c('term', mean_var), order = -1) #descending order
  else if(sorted=='up') setorderv(dt, cols=c('term', mean_var)) #ascending order
  
  return(dt)
  
}

#helper functions
absChange <- function(x) x-data.table::shift(x,n=1)
relChange <- function(x) (x-data.table::shift(x, n=1))/data.table::shift(x,n=1)
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

#read in the proper annotations (borders, lakes, mask)
annotations_path <- file.path(out.dir, 'annotations.RDs')
check <- file.exists(annotations_path)
annotations <- ifelse(
  check,
  readRDS(annotations_path),
  load_map_annotations()
)
if(!check) saveRDS(annotations, file=annotations_path)

#merge sr region names/IDs/SDI quintiles
locs <- get_location_metadata(location_set_id = 35, gbd_round_id = 6) %>% 
  .[, .(iso3=ihme_loc_id, location_id, location_name, super_region_id, super_region_name, region_id, region_name)] #subset to relevant columns

locs_sdi <- get_location_metadata(location_set_id = 40, gbd_round_id = 6) %>% 
  .[, .(location_id, location_name, parent_id, level)] #subset to relevant columns

locs_sdi <- locs_sdi[level!=0] %>% 
  merge(., 
        locs_sdi[level==0, .(location_id, sdi_quintile=location_name)], 
        by.x='parent_id',
        by.y='location_id') %>% 
  .[, `:=` (parent_id=NULL, location_name=NULL, level=NULL)]

#read in link_table
global_link_table <- file.path(global_link_dir, "lbd_full_link.rds") %>% readRDS %>% as.data.table
adm_links <- global_link_table[, .(ADM0_CODE, ADM1_NAME, ADM1_CODE, ADM2_CODE)] %>% unique

#create file to crosswalk AD0 to iso3
iso3_map <- dplyr::select(adm2, iso3, ADM0_CODE=gadm_geoid) 
iso3_map$geometry <- NULL
iso3_map <- as.data.table(iso3_map) %>% unique
locs <- merge(locs, iso3_map, by='iso3')
adm_links <- merge(adm_links, locs, by=c('ADM0_CODE'), all.x=T)

#read in the input data
input_dt <- file.path(share.model.dir, "cooking_fuel_solid.csv") %>% fread
ker_dt <- file.path(share.model.dir, "cooking_fuel_kerosene.csv") %>% fread

#read in results
results <- file.path(data.dir, 'all_summary.fst') %>% read_fst(as.data.table=T)
dt <- results[dimension=='ad2'] %>% Filter(function(x) !all(is.na(x)), .)
dt_d <- results[dimension=='ad2' & term%like%'change'] %>% Filter(function(x) !all(is.na(x)), .)

# #merge sr region names/IDs
dt <- merge(dt, adm_links, by=c('ADM0_CODE', 'ADM2_CODE'), all.x=T)

#also read in the all draws file
all_ad2 <-file.path(data.dir, 'all_draws.fst') %>%  read_fst(as.data.table=T)

#read in input data and prepare it for mapping
data <- load_map_results(indicator, indicator_group, run_date, raked, 
                   year_list=c(2000:2018),
                   custom_path = list('admin2'=dt),
                   geo_levels=c('admin2'),
                   cores=cores)
#define extent of map
zoom.afr <- data.table(x1=-10, x2=50, y1=-20, y2=40)
zoom.global <- data.table(x1=-120, x2=150, y1=-40, y2=55)
#***********************************************************************************************************************
 
# ---GENERAL METRICS----------------------------------------------------------------------------------------------------
#number of LMICs
stages[Stage %in% c('1', '2a', '2b'), uniqueN(iso3)]

#number of countries with data
uniqueN(input_dt$ihme_loc_id)

#number of surveys
uniqueN(input_dt$nid)

#number of people
input_dt[, sum(hh_size, na.rm=T)]
input_dt[source %like% "CENSUS", sum(N, na.rm=T)]
input_dt[!(source %like% "CENSUS"), sum(N, na.rm=T)]

#number of points/polys
table(input_dt$point)

#explore kerosene data
ker_dt <- merge(ker_dt, locs, by.x='ihme_loc_id', by.y='iso3')
ker_dt[, ad0_kerosene := weighted.mean(cooking_fuel_kerosene, w=N), by=.(ihme_loc_id, year)]
ggplot(ker_dt, aes(x=year, y=ad0_kerosene, color=super_region_name)) +
  geom_hline(yintercept=.1, color='grey10', linetype='dashed') +
  geom_point() +
  geom_smooth() +
  scale_color_brewer(type='qual', palette = 'Paired') +
  facet_wrap(~ihme_loc_id) +
  theme_minimal()

#***********************************************************************************************************************

# ---SPATIAL VARIATION--------------------------------------------------------------------------------------------------
##2017 patterns##

#country level patterns
tabulateR(lvl='global', types='HAP', ind='prev', terms=c('lvl'), years=c(2000, 2018))
tabulateR(lvl='super_region', types='HAP', ind='prev', terms=c('lvl'), years=c(2000, 2018))
tabulateR(lvl='ad0', types='HAP', ind='prev', terms=c('lvl'), years=c(2000, 2018))

#district level
#counts based on threshold
threshold <- .95

#what share of population in districts were above .95 in 2018
tabulateR(lvl='ad2', types='HAP', ind='prev', terms=c('lvl'),
          thresh=threshold, byvar='type')
tabulateR(lvl='ad2', types='HAP', ind='prev', terms=c('lvl'),
          filter=list(var='ADM0_CODE', vals=unique(locs[super_region_id==166, ADM0_CODE]), type=T),
          thresh=threshold, byvar='type')

#worst country
tabulateR(lvl='ad2', types='HAP', ind='prev', terms=c('lvl'),
          thresh=threshold, byvar='ADM0_CODE') %>% 
  head(10)

#worst country outside of SSA
tabulateR(lvl='ad2', types='HAP', ind='prev', terms=c('lvl'),
          filter=list(var='ADM0_CODE', vals=unique(locs[super_region_id==166, ADM0_CODE]), type=F),
          thresh=threshold, byvar='ADM0_CODE') %>% 
  head(10)

#worst country in Americas
tabulateR(lvl='ad2', types='HAP', ind='prev', terms=c('lvl'),
          filter=list(var='ADM0_CODE', vals=unique(locs[super_region_id==103, ADM0_CODE]), type=T),
          thresh=threshold, byvar='ADM0_CODE') %>% 
  head(10)

#counts of exposed
dt[term=='lvl' & year==analysis_year & type=='HAP', .(count_lower=sum(prev_lower*pop_total, na.rm=T),
                                                      count=sum(prev_mean*pop_total, na.rm=T),
                                                      count_upper=sum(prev_upper*pop_total, na.rm=T)
), by=.(ADM2_CODE, ADM2_NAME, ADM0_NAME, ADM1_NAME)] %>% 
  .[order(count)]
#***********************************************************************************************************************

# ---INEQUALITY---------------------------------------------------------------------------------------------------------
##produce metrics of inequality for 2017##
#calculate GINI/MAD at country level
dt_ineq <- dt[term=='lvl' & year==analysis_year & type=='HAP',
              .(iso3, year, ADM0_CODE, ADM2_CODE, ADM2_NAME, prev_mean, 
                super_region_id, super_region_name, region_id, region_name, pop_total)]
dt_ineq[, gini := gini(prev_mean, weights = pop_total), by=.(iso3, year)]
dt_ineq[, mean := weighted.mean(prev_mean, weights = pop_total), by=.(iso3, year)]
dt_ineq[, mad := mad(prev_mean, center = mean(prev_mean)), by=.(iso3, year)]
dt_ineq[, mean := mean(prev_mean, na.rm=T), by=.(iso3, year)]
dt_ineq[, max := max(prev_mean, na.rm=T), by=.(iso3, year)]
dt_ineq[, min := min(prev_mean, na.rm=T), by=.(iso3, year)]
dt_ineq[, range := max-min]

#range results
dt_ineq[year==analysis_year, .(range=max-min), by=iso3] %>% 
  unique %>% 
  .[order(range)]

exploreRange('MRT', shp=adm2, dt=dt[term=='lvl' & year==analysis_year & type=='HAP'], var='prev')
exploreRange('GTM', shp=adm2, dt=dt[term=='lvl' & year==analysis_year & type=='HAP'], var='prev')

#AID results
#note that AID = gini * 2 * mean
aid.dt <- dt[cause=='lri' & grouping=='under5' & type=='HAP' & term=='lvl', 
             .(iso3, year, ADM0_CODE, ADM0_NAME, ADM2_CODE, ADM2_NAME, prev_mean, 
               super_region_id, super_region_name, region_id, region_name, pop_total)] %>% 
  na.omit(., cols='prev_mean') %>% 
  .[, mean := weighted.mean(prev_mean, weights = pop_total, na.rm=T), by=.(iso3, year)] %>% 
  .[, .(mean,
        aid=gini(prev_mean, weights=pop_total) * 2 * mean,
        pop=sum(pop_total, na.rm=T)), 
    by=.(super_region_id, region_name, iso3, year)] %>% 
  unique(by=c('iso3', 'year')) %>% 
  .[order(aid),] %>% 
  .[year==min(year), aid_start := aid] %>% 
  .[, aid_start := mean(aid_start, na.rm=T), by=.(iso3)] %>%
  .[, aid_d := (aid-aid_start)] %>% 
  .[, aid_dr := aid_d/aid_start]

#average change in AID
aid.dt[year==max(year)] %>% .[,weighted.mean(aid_dr, weights=pop)]

#regional changes
aid.dt[year==max(year) & region_name=='South Asia'] #south asia
aid.dt[year==max(year) & region_id==167] #central sub-saharan africa
#***********************************************************************************************************************

# ---TEMPORAL VARIATION-------------------------------------------------------------------------------------------------
#global results
#abs/rel change in prevalence
tabulateR(lvl='global', types='HAP', ind='prev', terms=c('change', 'change_rate', 'lvl'), years=c(2000, 2018))

#absolute decrease in exposure
tabulateR(lvl='global', types='HAP', ind='prev', terms=c('lvl'), years=c(2000, 2018)) %>% 
  .[, lapply(.SD, function(x) x * pop_total), .SDcols=paste('prev', var_types, sep='_'), by=year] %T>% 
  print %>% 
  .[, lapply(.SD, absChange), .SDcols=paste('prev', var_types, sep='_')] %>% .[2]

#absolute decrease in proportion at global level
tabulateR(lvl='global', types='HAP', ind='share', terms=c('change', 'change_rate', 'lvl'), years=c(2000, 2018))

#regional/ad0 results
tabulateR(lvl='super_region', types='HAP', ind='prev', terms=c('change_rate'))
tabulateR(lvl='region', types='HAP', ind='prev', terms=c('change', 'change_rate'))
tabulateR(lvl='ad0', types='HAP', ind='prev', terms=c('change', 'change_rate'))

#best in SSA
tabulateR(lvl='ad0', types='HAP', ind='prev', terms=c('change_rate'),
          filter=list(var='ADM0_CODE', vals=unique(locs[super_region_id==166, ADM0_CODE]), type=T),
          sorted='up')

#ad2 results
tabulateR(lvl='ad2', types='HAP', ind='prev', terms=c('change_rate'))

#what percent of people lived in districts where SFU fell by more than 25%
tabulateR(lvl='ad2', types='HAP', ind='prev', terms=c('change'),
          thresh=-.25, inv_threshold = T, byvar='ADM0_CODE') %>% 
  head(30)

#how many people transitioned in china
#absolute decrease in exposure
tabulateR(lvl='ad0', types='HAP', ind='prev', terms=c('lvl'), years=c(2000, 2018),
          filter=list(var='ADM0_CODE', vals=unique(locs[iso3=='CHN', ADM0_CODE]), type=T)) %>% 
  .[, lapply(.SD, function(x) x * pop_total), .SDcols=paste('prev', var_types, sep='_'), by=year] %T>% 
  print %>% 
  .[, lapply(.SD, absChange), .SDcols=paste('prev', var_types, sep='_')] %>% .[2]

#what percent of districts improved in India
tabulateR(lvl='ad2', types='HAP', ind='prev', terms=c('change'),
          thresh=0, inv_threshold = T, byvar='ADM0_CODE') %>% 
  .[ADM0_NAME=='India']

#range results
dt[year==analysis_year & term=='change' & type=='HAP', 
   .(range=max(prev_mean, na.rm=T)-min(prev_mean, na.rm=T)), by=iso3] %>%
  unique %>%
  .[order(range)]

#explore the range to find inequality for the top 3 ranges (all in South Asia)
exploreRange('IND', shp=adm2, 
             dt=dt[term=='change' & year==analysis_year & type=='HAP'],
             var='prev') %T>%
  print %>% 
  .[, lapply(.SD, function(x) x - data.table::shift(x)), .SDcols=paste0('var_', var_types)]

exploreRange('BTN', shp=adm2, 
             dt=dt[term=='change' & year==analysis_year & type=='HAP'],
             var='prev') %T>%
  print %>% 
  .[, lapply(.SD, function(x) x - data.table::shift(x)), .SDcols=paste0('var_', var_types)]

exploreRange('NPL', shp=adm2, 
             dt=dt[term=='change' & year==analysis_year & type=='HAP'],
             var='prev') %T>%
  print %>% 
  .[, lapply(.SD, function(x) x - data.table::shift(x)), .SDcols=paste0('var_', var_types)]

exploreRange('BGD', shp=adm2, 
             dt=dt[term=='change' & year==analysis_year & type=='HAP'],
             var='prev') %T>%
  print %>% 
  .[, lapply(.SD, function(x) x - data.table::shift(x)), .SDcols=paste0('var_', var_types)]
#***********************************************************************************************************************

# ---SDG PROJECTIONS----------------------------------------------------------------------------------------------------
#SDG projection probabilities
#append the ADM0 files
sdg_files_ad0 <-
  file.path(data.dir, 'sdg_projections') %>% list.files(pattern='admin_0', full.names = T) %>% 
  lapply(., readRDS)

#append the ADM2 files
sdg_files <-
  file.path(data.dir, 'sdg_projections') %>% list.files(pattern='admin_2', full.names = T) %>% 
  lapply(., readRDS)

#extract goal obj to index over
goals <- lapply(1:length(sdg_files), function(x) sdg_files[[x]]$goals) %>% rbindlist %>% unique

#ad0 files
#create a dt with all probabilities
probs_ad0 <- lapply(1:nrow(goals), prepCasts, type='probs', id_dt=goals, list=sdg_files_ad0) %>% 
  rbindlist %>% 
  setnames(c('target_year', 'spatial_idx'), c('year', 'ADM0_CODE')) %>% 
  merge(., locs, by='ADM0_CODE') %>% 
  merge(., results[dimension=='ad2' & year==analysis_year & type=='HAP' & term=='lvl', pop_total, by=ADM0_CODE],
        by='ADM0_CODE')

#ad2 files
#create a dt with all probabilities
probs <- lapply(1:nrow(goals), prepCasts, type='probs', id_dt=goals, list=sdg_files) %>% 
  rbindlist %>% 
  setnames(c('target_year', 'spatial_idx'), c('year', 'ADM2_CODE')) %>% 
  merge(., adm_links[, .(ADM2_CODE, ADM0_CODE)], by='ADM2_CODE')  %>% 
  merge(., locs, by='ADM0_CODE') %>% 
  merge(., results[dimension=='ad2' & year==analysis_year & type=='HAP' & term=='lvl', pop_total, by=ADM2_CODE],
        by='ADM2_CODE')

#create a dt with all projections
projs <- 
  lapply(c(2018, seq(2020, 2030, 5)), prepCasts, type='proj', id_var='year') %>% 
  rbindlist %>% 
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
projs[, sev_rel_uncertainty := (sev_upper-sev_lower)/sev_mean]
projs[sev_rel_uncertainty>2, sev_rel_uncertainty := 2] #cap at 2

#district level results
target_threshold <- .05
prob_threshold <- .95

#how many countries will succeed 
ad0_probs <-
  probs[target==target_threshold & year==2030, 
        .(prob=weighted.mean(absolute_goal_prob, w=pop_total), 
          pop=sum(pop_total, na.rm=T)), by=.(super_region_id, region_id, iso3)] %>% 
  .[, global_pop := sum(pop, na.rm=T)] %>% 
  merge(., iso3_map, by='iso3')

ad0_probs[, .(success=sum(prob>=prob_threshold, na.rm=T), fail=sum(prob<=(1-prob_threshold), na.rm=T))]
ad0_probs[prob>=prob_threshold, .(regs=uniqueN(region_id), pop=sum(pop, na.rm=T), global_pop=mean(global_pop))] %>% 
  .[, pop_share := pop/global_pop] %>% print
ad0_probs[prob<=(1-prob_threshold), .(regs=uniqueN(region_id), pop=sum(pop, na.rm=T), global_pop=mean(global_pop))] %>% 
  .[, pop_share := pop/global_pop] %>% print

#how many countries will succeed in every district
ad0_probs <-
probs[target==target_threshold & year==2030] %>% 
  copy %>% 
      .[, .(success=sum(absolute_goal_prob>=prob_threshold, na.rm=T)/.N, 
        fail=sum(absolute_goal_prob<=(1-prob_threshold)/.N, na.rm=T),
        pop=sum(pop_total, na.rm=T)),
      by=.(super_region_id, region_id, iso3)] %>% 
  .[, global_pop := sum(pop, na.rm=T)]


#what is the country with the highest SFU in 2000 that is on track to succeed
tabulateR(lvl='ad0', types='HAP', ind='prev', terms=c('lvl'), years=c(2000)) %>% 
  .[ADM0_CODE %in% ad0_probs[prob>=prob_threshold, unique(ADM0_CODE)]]

#what pct of people live in districts that will meet goal by 2030
probs[target==target_threshold & year==2030 & absolute_goal_prob>=prob_threshold, 
      sum(pop_total, na.rm=T)] / probs[target==target_threshold & year==2030, sum(pop_total, na.rm=T)]

#how many districts have met in 2018
probs[target==target_threshold & year==2018] %>% 
  .[, .(count=sum(absolute_goal_prob>=prob_threshold, na.rm=T), N=.N)] %>% 
  .[, .(count, N, pct=count/N)]
#what share of pop
probs[target==target_threshold & year==2018 & absolute_goal_prob>=prob_threshold, 
      sum(pop_total, na.rm=T)] / probs[target==target_threshold & year==2018, sum(pop_total, na.rm=T)]

#what share of unmet districts will meet between 2018-2030
unmet_districts <- probs[target==target_threshold & year==2018 & absolute_goal_prob<prob_threshold, unique(ADM2_CODE)]
probs[target==target_threshold & year==2030 & ADM2_CODE%in%unmet_districts] %>% 
  .[, .(count=sum(absolute_goal_prob>=prob_threshold, na.rm=T), N=.N)] %>% 
  .[, .(count, N, pct=count/N)]
#what share of pop
probs[target==target_threshold & year==2030 & absolute_goal_prob>=prob_threshold & ADM2_CODE%in%unmet_districts, 
      sum(pop_total, na.rm=T)] / probs[target==target_threshold & year==2018, sum(pop_total, na.rm=T)]

#superregional breakdown
probs[target==target_threshold & year==2030 & ADM2_CODE%in%unmet_districts] %>% 
  .[, .(count=sum(absolute_goal_prob>=prob_threshold, na.rm=T), 
        pop=sum((absolute_goal_prob>=prob_threshold)*pop_total, na.rm=T), 
        pop_N=sum(pop_total, na.rm=T),
        N=.N), by=super_region_name] %>% 
  .[, .(count, N, pct=count/N, pop_share=pop/pop_N), by=super_region_name] %>% 
  .[order(pop_share)]

#regional breakdown
probs[target==target_threshold & year==2030 & ADM2_CODE%in%unmet_districts] %>% 
  .[, .(count=sum(absolute_goal_prob>=prob_threshold, na.rm=T), 
        pop=sum((absolute_goal_prob>=prob_threshold)*pop_total, na.rm=T), 
        pop_N=sum(pop_total, na.rm=T),
        N=.N), by=region_name] %>% 
  .[, .(count, N, pct=count/N, pop_share=pop/pop_N), by=region_name] %>% 
  .[order(pop_share)]

#country breakdown
probs[target==target_threshold & year==2030 & ADM2_CODE%in%unmet_districts] %>% 
  .[, .(count=sum(absolute_goal_prob>=prob_threshold, na.rm=T), 
        pop=sum((absolute_goal_prob>=prob_threshold)*pop_total, na.rm=T), 
        pop_N=sum(pop_total, na.rm=T),
        N=.N), by=iso3] %>% 
  .[, .(count, N, pct=count/N, pop_share=pop/pop_N), by=iso3] %>% 
  .[order(pop_share)]

#country breakdown for SSA
probs[target==target_threshold & year==2030 & ADM2_CODE%in%unmet_districts & region_name %like% 'Latin'] %>% 
  .[, .(count=sum(absolute_goal_prob>=prob_threshold, na.rm=T), 
        pop=sum((absolute_goal_prob>=prob_threshold)*pop_total, na.rm=T), 
        pop_N=sum(pop_total, na.rm=T),
        N=.N), by=iso3] %>% 
  .[, .(count, N, pct=count/N, pop_share=pop/pop_N), by=iso3] %>% 
  .[order(pop_share)]

#country failure breakdown for 2030
probs[target==target_threshold & year==2030] %>% 
  .[, .(count=sum(absolute_goal_prob<=(1-prob_threshold), na.rm=T), 
        pop=sum((absolute_goal_prob<=(1-prob_threshold))*pop_total, na.rm=T), 
        pop_N=sum(pop_total, na.rm=T),
        N=.N), by=iso3] %>% 
  .[, .(count, N, pct=count/N, pop_share=pop/pop_N), by=iso3] %>% 
  .[order(pop_share)]

#country failure breakdown for 2030 (ESSA)
probs[target==target_threshold & year==2030 & region_id %in% c(174)] %>% 
  .[, .(count=sum(absolute_goal_prob<=(1-prob_threshold), na.rm=T), 
        pop=sum((absolute_goal_prob<=(1-prob_threshold))*pop_total, na.rm=T), 
        pop_N=sum(pop_total, na.rm=T),
        N=.N), by=iso3] %>% 
  .[, .(count, N, pct=count/N, pop_share=pop/pop_N), by=iso3] %>% 
  .[order(pop_share)]

#country with the largest divide
probs[target==target_threshold & year==2030] %>% 
  .[, .(count_fail=sum(absolute_goal_prob<.5, na.rm=T), 
        count_success=sum(absolute_goal_prob>.5, na.rm=T),
        pop_fail=sum((absolute_goal_prob<.5)*pop_total, na.rm=T), 
        pop_N=sum(pop_total, na.rm=T),
        N=.N), by=iso3] %>% 
  .[, pop_share := pop_fail/pop_N] %>% 
  .[count_fail!=0&count_success!=0] %>% 
  .[order(pop_share)]

probs[target==target_threshold & year==2030] %>% 
  .[, .(pop_fail=sum((absolute_goal_prob<.5)*pop_total, na.rm=T), 
        pop_N=sum(pop_total, na.rm=T),
        N=.N), by=iso3] %>% 
  .[, pop_share := pop_fail/pop_N] %>% 
  .[count_fail!=0&count_success!=0] %>% 
  .[order(ratio)]
#***********************************************************************************************************************

# ---AIR POLLUTION------------------------------------------------------------------------------------------------------
##analyze relationship to AAP; TAP; HAP_SHARE##

##global levels/change
#current tap_pc/tap_paf global avg
tabulateR(lvl='global', types='TAP', ind='pm_pc', terms=c('lvl', 'change', 'change_rate'))  #highest TAP
tabulateR(lvl='global', types='TAP', ind='paf', terms=c('lvl', 'change', 'change_rate')) #highest TAP PAF
tabulateR(lvl='global', types='HAP', ind='share', terms=c('lvl', 'change', 'change_rate')) #highest HAP share

tabulateR(lvl='super_region', types='TAP', ind='pm_pc') #highest TAP
tabulateR(lvl='super_region', types='TAP', ind='paf', years=c(start_year, end_year)) #highest TAP PAF
tabulateR(lvl='super_region', types='HAP', ind='share', terms=c('lvl', 'change', 'change_rate')) #highest HAP share

tabulateR(lvl='region', types='TAP', ind='pm_pc') #highest TAP
tabulateR(lvl='region', types='TAP', ind='paf', years=c(start_year, end_year)) #highest TAP PAF
tabulateR(lvl='region', types='HAP', ind='share', terms=c('lvl', 'change')) #highest HAP share

#examine south asia
tabulateR(lvl='region', types='HAP', ind='share', terms=c('lvl', 'change'), years=c(start_year, end_year)) %>% 
  .[region_id==159]#change in HAP share in South Asia
tabulateR(lvl='region', types='HAP', ind='prev', terms=c('change_rate')) %>% 
  .[region_id==159]#change in SFU for South Asia
tabulateR(lvl='region', types='AAP', ind='pm_pc', terms=c('lvl', 'change'), years=c(start_year, end_year)) %>% 
  .[region_id==159]#change in AAP dose in South Asia

#deeper dive to countries in south asia
tabulateR(lvl='ad0', types='HAP', ind='prev', terms=c('change_rate'), 
          filter=list(var='ADM0_CODE', vals=unique(locs[super_region_name=='South Asia', ADM0_CODE]), type=T))

##country levels
tabulateR(lvl='ad0', types='TAP', ind='pm_pc') #highest TAP
tabulateR(lvl='ad0', types='TAP', ind='pm_pc', 
          filter=list(var='super_region_id', vals=166, type=F)) #highest TAP excluding SSA

tabulateR(lvl='ad0', types='HAP', ind='share', 
          filter=list(var='super_region_id', vals=166, type=F)) #highest HAP share excluding SSA

tabulateR(lvl='ad0', types='TAP', ind='pm_pc', sorted='up') #lowest TAP
tabulateR(lvl='ad0', types='AAP', ind='share', sorted='up') #lowest AAP share
tabulateR(lvl='ad0', types='AAP', ind='share') #highest AAP share

#country level HAP:AAP ratio, what percent of countries is HAP the main contributor vs AAP
tabulateR(lvl='ad0', types='HAP', ind='share') %>% .[share_mean>.5, .N] #2018
tabulateR(lvl='ad0', types='HAP', ind='share', years=start_year) %>% .[share_mean>.5, .N] #2000

#pct of population that lives below WHO threshold 
threshold <- 5
dt[year==analysis_year & term=='lvl' & type=='TAP', 
   lapply(.SD, function(x) sum((x < threshold)*pop, na.rm=T)/sum(pop, na.rm=T)),
   .SDcols=paste0('pm_pc_', var_types)]

#pct of population that lives below WHO threshold interim-1
threshold <- 35
dt[year==analysis_year & term=='lvl' & type=='TAP', 
   lapply(.SD, function(x) sum((x < threshold)*pop, na.rm=T)/sum(pop, na.rm=T)),
   .SDcols=paste0('pm_pc_', var_types)]
dt[year==min(year) & term=='lvl' & type=='TAP', 
   lapply(.SD, function(x) sum((x < threshold)*pop, na.rm=T)/sum(pop, na.rm=T)),
   .SDcols=paste0('pm_pc_', var_types)]

#***********************************************************************************************************************

# ---ATTRIBUTABLE LRI---------------------------------------------------------------------------------------------------
##analyses of attributable LRI##

#what was the u5 mortality rate of LRI in regions where more than half of LRI was attributed to TAP

#how many children died from LRI attributable to TAP
tabulateR(lvl='global', types='TAP', ind='atr_count', years=c(analysis_year, start_year),
          terms=c('lvl', 'change', 'change_rate')) 
tabulateR(lvl='global', types=c('HAP', 'AAP'), ind='atr_count', years=c(analysis_year, start_year),
          terms=c('lvl', 'change', 'change_rate')) 

tabulateR(lvl='super_region', types='TAP', ind='atr_count') #highest TAP count
tabulateR(lvl='region', types='TAP', ind='atr_count') #highest TAP count
tabulateR(lvl='ad0', types='TAP', ind='atr_count') #highest TAP count

#how did the HAP share change
tabulateR(lvl='global', types='HAP', ind='share', terms=c('lvl', 'change', 'change_rate')) 

#how did the PAFs change
tabulateR(lvl='global', types='TAP', ind='paf', terms=c('lvl', 'change', 'change_rate'), years=c(analysis_year, start_year))

#PAFs by reg/country/district
tabulateR(lvl='super_region', types='TAP', ind='paf') 
tabulateR(lvl='region', types='TAP', ind='paf')
tabulateR(lvl='ad0', types='TAP', ind='paf', sorted='up') 
tabulateR(lvl='ad2', types='TAP', ind='paf')

#How many districts are majority TAP attributable
dt[year==analysis_year & type=='TAP' & term=='lvl', lapply(.SD, function(x) sum(x>.5, na.rm=T)/.N), 
   .SDcols=paste0('paf_', var_types)]

#investigate countries where LRI rate fell but SFU was stable
tabulateR(lvl='ad0', types='TAP', ind='rate', terms='change_rate', years=c(analysis_year, start_year), sorted='up')
tabulateR(lvl='ad0', types='HAP', ind='prev', terms='change_rate', years=c(analysis_year, start_year))

#what is the lowest HAP paf in a country where LRI rate remains above 4/1000
results[year==analysis_year & dimension=='ad0' & type=='HAP' & term=='lvl' & rate_mean>2/1e3] %>% 
  .[order(paf_mean), .(ADM0_NAME, 
                       rate_lower=rate_lower*1e3, rate_mean=rate_mean*1e3, rate_upper=rate_upper*1e3, 
                       prev_lower, prev_mean, prev_upper,
                       paf_lower, paf_mean, paf_upper)]

#proportion of total deaths from HAP
#change in the percent of LRI deaths attributable to HAP
#use the all ad2 results
all_deaths <- all_ad2[, sum(atr_count, na.rm=T), by=.(year, type, draw)] %>% 
  dcast(., year+draw~type, value.var='V1') %>% 
  .[, hap_share := HAP/TAP] %>%
  .[, aap_share := AAP/TAP]

ind_cols <- c('hap_share', 'aap_share')

all_deaths[, paste0(ind_cols, '_lower') := lapply(.SD, quantile, p=.025, na.rm=T), .SDcols=ind_cols, by=year]
all_deaths[, paste0(ind_cols, '_mean') := lapply(.SD, mean, na.rm=T), .SDcols=ind_cols, by=year]
all_deaths[, paste0(ind_cols, '_upper') := lapply(.SD, quantile, p=.975, na.rm=T), .SDcols=ind_cols, by=year]

all_deaths[, c('year',
          names(all_deaths) %>% .[. %like% 'lower|mean|upper'] %>% sort), with=F] %>% 
  unique(by='year')
#***********************************************************************************************************************
