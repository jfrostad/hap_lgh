# ----HEADER------------------------------------------------------------------------------------------------------------
# Author: JF
# Date: 09/11/2018 (never forget)
# Purpose: Exploring data coverage across indicators for hap
# source("/homes/jfrostad/_code/lbd/hap/extract/5_data_coverage.R", echo=T)
#***********************************************************************************************************************

# ----CONFIG------------------------------------------------------------------------------------------------------------
# clear memory
rm(list=ls())

# runtime configuration
if (Sys.info()["sysname"] == "Linux") {
  j_root <- "/home/j/"
  h_root <- file.path("/ihme/homes", Sys.info()["user"])
  arg <- commandArgs()[-(1:3)] # First args are for unix use only
  
  if (length(arg)==0) {
    # arg <- c("IND", #current project iteration
    #          "8", #output version
    #          1) #number of cores provided to multicore functions
  }
  
  package_lib    <- file.path(h_root, '_code/_lib/pkg')
  ## Load libraries and  MBG project functions.
  .libPaths(c( .libPaths(), package_lib))
  .libPaths(package_lib)
  
  # necessary to set this option in order to read in a non-english character shapefile on a linux system (cluster)
  Sys.setlocale(category = "LC_ALL", locale = "C")
  
} else {
  j_root <- "J:"
  h_root <- "H:"
  # arg <- c("IND", #current project iteration
  #          "4", #output version
  #          1) #number of cores provided to multicore functions
}

#load packages
pacman::p_load(data.table, RMySQL, dplyr, feather, ggmosaic, ggplot2, ggrepel, googledrive, gridExtra, maptools, 
               questionr, parallel, raster, RColorBrewer, readxl, rgdal, rgeos, 
               scales, survival, stringr, tm, viridis, ggwordcloud) 

## Set core_repo location and indicator group
user            <- Sys.info()['user']
core_repo       <- file.path(h_root, '_code/lbd/lbd_core')
my_repo         <- file.path(h_root, '_code/lbd/hap')
commondir       <- file.path(my_repo, 'mbg_central/share_scripts/common_inputs')
package_list    <- file.path(commondir, 'package_list.csv') %>% fread %>% t %>% c

#capture date
today <- "2019_03_26" #date of current post-extraction

#options
cores <- 10
new.gbd.results <- F #set T if GBD results have been updated and need to redownload
new.extracts <- F
topic <- "hap"
this.family <- 'cooking'
indicators <- c('cooking_fuel', 'cooking_type', 'cooking_type_chimney', 'cooking_location', 
                'heating_fuel', 'heating_type', 'heating_type_chimney', 'lighting_fuel', 'electricity', 
                'housing_roof', 'housing_wall', 'housing_floor')
#***********************************************************************************************************************

# ----IN/OUT------------------------------------------------------------------------------------------------------------
###Input###
#raw data
doc.dir <- file.path(j_root, 'WORK/11_geospatial/hap/documentation')
def.file <- file.path(doc.dir, 'definitions.xlsx')
raw.dir <- file.path("/ihme/limited_use/LIMITED_USE/LU_GEOSPATIAL/ubCov_extractions/hap/batch")
in.dir <- file.path('/share/limited_use/LIMITED_USE/LU_GEOSPATIAL/collapse/hap/')
temp.dir <- file.path(j_root, 'temp/jfrostad')
share.dir <- file.path('/share/geospatial/jfrostad')
###Output###
out.dir  <- file.path(j_root, 'temp/jfrostad/housing/')
graph.dir <- file.path(j_root, 'WORK/11_geospatial/hap/graphs')
#***********************************************************************************************************************

# ---FUNCTIONS----------------------------------------------------------------------------------------------------------
##function lib##
tabulateIndicators <- function(dt, indicator, byvar) {

  message('tabulating ', indicator)

  out <- dt[, sum(!is.na(get(indicator))), by=byvar]
  setnames(out, 'V1', indicator)
  
  if(byvar=='') out[, merge := 1]
  
  return(out)
  
}

#find values %like% helper
getLikeMe <- function(obj, val, invert=F) {
  
  message('finding values that ', ifelse(invert, 'are NOT', 'are') ,' like: ', val, '\n|', 
          '\no~~> from: ', deparse(match.call()$obj))
  
  #add option to invert
  if (invert) { 
    
    out <- names(obj)[!(names(obj) %like% val)]
    
  } else out <- names(obj)[names(obj) %like% val]
  
  return(out)
  
}

# Load custom functions
source(file.path(my_repo, 'cooking/model/_lib/fx.R'))
hap.function.dir <- file.path(h_root, '_code/lbd/hap/extract/functions')
#this pulls hap collapse helper functions
file.path(hap.function.dir, '/collapse_fx.R') %>% source

#mbg function lib
pkg.list <- c('RMySQL', 'data.table', 'dismo', 'doParallel', 'dplyr', 'foreign', 'gbm', 'ggplot2', 'glmnet', 
              'grid', 'gridExtra', 'gtools', 'magrittr', 'pacman', 'parallel', 'plyr', 'raster', 'rgdal', 'rgeos',
              'seegMBG', 'seegSDM', 'tictoc') #will be loaded by MBG setup
lbd.shared.function.dir <- file.path(my_repo, "mbg_central")
file.path(lbd.shared.function.dir, 'setup.R') %>% source
mbg_setup(repo=lbd.shared.function.dir, package_list=pkg.list) #load mbg functions

##shared functions##
#gbd#
gbd.shared.function.archive <- file.path(j_root,  "temp/central_comp/libraries/2017_archive/r")
gbd.shared.function.dir <- file.path(j_root,  "temp/central_comp/libraries/current/r")
file.path(gbd.shared.function.dir, 'get_covariate_estimates.R') %>% source
file.path(gbd.shared.function.archive, 'get_draws.R') %>% source
file.path(gbd.shared.function.dir, 'get_location_metadata.R') %>% source
file.path(gbd.shared.function.dir, 'get_ids.R') %>% source
file.path(gbd.shared.function.dir, 'get_population.R') %>% source
file.path(gbd.shared.function.dir, 'get_outputs.R') %>% source

stop()
#***********************************************************************************************************************

# ---PREP DATA----------------------------------------------------------------------------------------------------------
codebook <- read_xlsx(file.path(raw.dir, 'hap.xlsx'), sheet='codebook') %>% as.data.table
codebook <- codebook[assigned=='qnguyen1' | assigned == 'jfrostad' | assigned == 'albrja']

##make a table showing the count of each indicator by survey series##
table <- lapply(indicators, tabulateIndicators, dt=codebook, byvar='survey_name') %>%  
  Reduce(function(...) merge(..., all = TRUE, by = "survey_name"), .)

#also add a total row
all <- lapply(indicators, tabulateIndicators, dt=codebook, byvar='') %>%  
  Reduce(function(...) merge(..., all = TRUE, by = "merge"), .)
all[, merge := NULL]
all[, survey_name := '1_Total']

#merge total row with the other survey series rows and return the table
table <- list(table, all) %>% rbindlist(use.names=T)

#also add a total column and sort by
table[, total := rowSums(.SD), .SDcols=indicators]
write.csv(table[order(-total)], file=file.path(doc.dir, 'data_universe.csv'), row.names = F)
#***********************************************************************************************************************

# ---GRAPH CATEGORIES---------------------------------------------------------------------------------------------------
dt <- paste0(in.dir, "/", "data_", this.family, '_', today, ".feather") %>% read_feather %>% as.data.table
dt[, name := tolower(ihme_loc_id) %>% get_adm0_codes, by = ihme_loc_id] #add gaul codes

#add on location hierarchy info
locs <- get_location_hierarchy(41)
dt <- merge(dt, locs[, .(ihme_loc_id, region_name, super_region_name)], by='ihme_loc_id')

#pull hap exposure from gbd2017 - command provided by kate causey
if(new.gbd.results==T) {
  
  hap.exp <- get_draws(gbd_id_type = "rei_id",
                       gbd_id=87,
                       source="exposure",
                       year_id=c(2000:2017),
                       location_id=locs[,location_id],
                       age_group_id=2,
                       sex_id=2,
                       gbd_round_id=5)
  
  hap.exp <- hap.exp[parameter=='cat1'] #cat2 is just the inverse
  
  #produce the mean and CI
  draw.cols <- getLikeMe(hap.exp, 'draw')
  hap.exp <- hap.exp[, gbd_lower := apply(.SD, 1, quantile, probs=.025), .SDcols=draw.cols]
  hap.exp <- hap.exp[, gbd_mean := rowMeans(.SD), .SDcols=draw.cols]
  hap.exp <- hap.exp[, gbd_upper := apply(.SD, 1, quantile, probs=.975), .SDcols=draw.cols]
  hap.exp[, (draw.cols) := NULL] #no longer need
  
  #add location hierarchy to prep for merge to geospatial data
  hap.exp <- merge(hap.exp, 
                   locs[, .(location_id, ihme_loc_id, region_name, region_id, super_region_name)], 
                   by='location_id')
  hap.exp[, reg_iso3 := paste0(region_id, ': ', ihme_loc_id)]

  write.csv(hap.exp, file = file.path(graph.dir, 'gbd_hap_results.csv'))
  
} else hap.exp <- file.path(graph.dir, 'gbd_hap_results.csv') %>% fread

#sort by gbd exposure and then bin years in 5 year groups
plot.dt <- merge(dt, hap.exp, by.x=c('ihme_loc_id', 'int_year'), by.y=c('ihme_loc_id', 'year_id'))
plot.dt <- plot.dt[order(gbd_mean)]
plot.dt[, order := .GRP, by=ihme_loc_id]
plot.dt[, country := paste0(order, "-", ihme_loc_id)]
plot.dt[, year := round_any(int_year, 5)]

#reshape props long
fuel.type.cols <- c('cooking_fuel_')

#categorical: set up order of severity (best to worst)
cat.order <- c('none', 'electricity', 'gas', 'kerosene', 'wood', 'crop', 'coal', 'dung', 'other', 'unknown')
plot.dt[, category := factor(cat_cooking_fuel_mapped, levels=cat.order)]
plot.dt[, N_cat := sum(N), by=.(ihme_loc_id, year, category)]
plot.dt[, prop := N_cat/sum(N), by=.(ihme_loc_id, year)]

#ord good: set up order of severity (best to worst)
# ord.good.order <- c('clean', 'medium', 'dirty')
# plot.dt[, ord_good := factor(ord_good_cooking_fuel_mapped, levels=ord.good.order)]
# plot.dt[, N_ord_good := sum(N), by=.(ihme_loc_id, year, ord_good)]
# plot.dt[, prop_ord_good := N_ord_good/sum(N), by=.(ihme_loc_id, year)]


#set up color scale
colors <- c(plasma(9, direction=-1), "#C0C0C0", "#C0C0C0") #use gray for other and unknown
names(colors) <- cat.order

plot <- ggplot(plot.dt[year>=2000], aes(x=reorder(ihme_loc_id,order,sum))) + 
  geom_bar(aes(fill = category), position = position_stack(reverse = TRUE)) +
  coord_flip() +
  facet_wrap(~year) +
  scale_fill_manual(values = colors) +
  ggtitle("HAP data by number of datapoints", subtitle = "Countries sorted by decreasing use of solid fuels (GBD 2017)") +
  theme(legend.position = "top") +
  theme_minimal()

print(plot)

plot <- ggplot(plot.dt[year>=2000], aes(x=reorder(ihme_loc_id,order,sum))) + 
  geom_bar(aes(fill = category, weight=N), position = position_stack(reverse = TRUE)) +
  coord_flip() +
  facet_wrap(~year) +
  scale_fill_manual(values = colors) +
  ggtitle("HAP data by number of individuals surveyed", subtitle = "Countries sorted by decreasing use of solid fuels (GBD 2017)") +
  theme(legend.position = "top") +
  theme_minimal()

print(plot)

propPlot <- function(region) {
  
  message('plotting: ', region)
  
  plot <- ggplot(plot.dt[year>=2000 & super_region_name==region], aes(x=reorder(ihme_loc_id,order,sum), y=prop, fill=category)) + 
    geom_bar(stat = "identity",position="fill") +
    coord_flip() +
    facet_wrap(~year) +
    scale_fill_manual(values = colors) +
    ggtitle(paste0("HAP data by proportion: ", region), subtitle = "Countries sorted by decreasing use of solid fuels (GBD 2017)") +
    theme(legend.position = "top") +
    theme_minimal()
  
  return(plot)
  
}

pdf(file=file.path(out.dir, 'category_proportions.pdf'), onefile=T, width=11, height=8)
lapply(unique(plot.dt$super_region_name), propPlot) 
dev.off()

#set up color scale
colors <- plasma(3, direction=-1)
names(colors) <- ord.good.order

propPlot <- function(region) {
  
  message('plotting: ', region)
  
  plot <- ggplot(plot.dt[year>=2000 & super_region_name==region], aes(x=reorder(ihme_loc_id,order,sum), y=prop_ord_good, fill=ord_good)) + 
    geom_bar(stat = "identity",position="fill") +
    coord_flip() +
    facet_wrap(~year) +
    scale_fill_manual(values = colors) +
    ggtitle(paste0("HAP data by proportion: ", region), subtitle = "Countries sorted by decreasing use of solid fuels (GBD 2017)") +
    theme(legend.position = "top") +
    theme_minimal()
  
  return(plot)
  
}

pdf(file=file.path(out.dir, 'ord_good_proportions.pdf'), onefile=T, width=11, height=8)
lapply(unique(plot.dt$super_region_name), propPlot) 
dev.off()

#set up color scale
colors <- plasma(3, direction=-1)
names(colors) <- ord.bad.order

propPlot <- function(region) {
  
  message('plotting: ', region)
  
  plot <- ggplot(plot.dt[year>=2000 & super_region_name==region], aes(x=reorder(ihme_loc_id,order,sum), y=prop_ord_bad, fill=ord_bad)) + 
    geom_bar(stat = "identity",position="fill") +
    coord_flip() +
    facet_wrap(~year) +
    scale_fill_manual(values = colors) +
    ggtitle(paste0("HAP data by proportion: ", region), subtitle = "Countries sorted by decreasing use of solid fuels (GBD 2017)") +
    theme(legend.position = "top") +
    theme_minimal()
  
  return(plot)
  
}

pdf(file=file.path(out.dir, 'ord_bad_proportions.pdf'), onefile=T, width=11, height=8)
lapply(unique(plot.dt$super_region_name), propPlot) 
dev.off()
#***********************************************************************************************************************

# ---GRAPH CATS (MOSAIC)------------------------------------------------------------------------------------------------
#compile raw data
compileAndDefine <- function(x, defs) {
  
  dt <- read_feather(x) %>% 
    as.data.table %>% 
    defIndicator(., var.fam='cooking', definitions=defs, debug=F, clean_up=F)
  
}

#not using definitions at this time
cooking <- file.path(share.dir, 'cooking') %>% 
  list.files(full.names = T, pattern='uncollapsed_') %>% 
  mclapply(., function(x) as.data.table(read_feather(x)), mc.cores=15) %>% 
  rbindlist

#cleanup
cooking[cooking_fuel_mapped == 'biomass', cooking_fuel_mapped := 'crop'] #TODO fix biomass in the original coding
cooking[is.na(cooking_fuel_mapped) | cooking_fuel_mapped == '', cooking_fuel_mapped := 'missing']
cooking[is.na(cooking_type_mapped) | cooking_type_mapped == '', cooking_type_mapped := 'missing']
cooking[is.na(cooking_type_chimney_mapped) | cooking_type_chimney_mapped == '', cooking_type_chimney_mapped := 'missing']
cooking[is.na(cooking_location_mapped) | cooking_location_mapped == '', cooking_location_mapped := 'missing']
#remove.vars <- names(cooking)[names(cooking) %like% "row_id|ord_|cat_"]
#cooking <- cooking[, (remove.vars) := NULL]

#merge location info
cooking <- merge(cooking, locs[, .(ihme_loc_id, region_name, super_region_name)], by='ihme_loc_id')

#set up order of fuel types
cat.order <- c('none', 'electricity', 'gas', 'kerosene', 'coal', 'wood', 'crop', 'dung', 'other', 'unknown', 'missing')
cooking[, cooking_fuel := factor(cooking_fuel_mapped, levels=cat.order)]

#set up color scale for fuel types
colors <- c(plasma(8, direction=-1), "#C0C0C0", "#C0C0C0", "#C0C0C0") #use gray for other and unknown
names(colors) <- cat.order

#set up color scale for risk
#TODO investigate assumptions
# cooking[cooking_risk<6, risk := 'clean']
# cooking[cooking_risk==6, risk := 'med']
# cooking[cooking_risk>6, risk := 'dirty']
# risk.order <- c('clean', 'med', 'dirty')
# cooking[, risk := factor(risk, levels=risk.order)]
# risk.colors <- plasma(3, direction=-1)
# names(risk.colors) <- risk.order

#collapse by country to improve speed
setkey(cooking, nid, cluster_id, cooking_fuel, cooking_type_mapped, cooking_type_chimney_mapped, cooking_location_mapped)
cooking[, N := sum(hhweight*hh_size, na.rm=T)^2/sum(hhweight^2*hh_size, na.rm=T), by=key(cooking)]
cooking[, hh_size_sum := sum(hh_size, na.rm=T), by=key(cooking)]
plot.dt <- unique(cooking, by=key(cooking))

#cooking_type
pdf(file=file.path(out.dir, 'global_cooking_fuel_v_cooking_type.pdf'), onefile=T, width=11, height=8)
ggplot(data = plot.dt[cooking_fuel != 'missing']) +
  geom_mosaic(aes(weight = N,
                  x = product(cooking_fuel, cooking_type_mapped), 
                  fill=cooking_fuel), na.rm=T) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 65))
dev.off()

##location
pdf(file=file.path(out.dir, 'global_cooking_fuel_v_cooking_location.pdf'), onefile=T, width=11, height=8)
ggplot(data = plot.dt[!(cooking_fuel %in% c('missing', 'other', 'unknown'))]) +
  geom_mosaic(aes(weight = N,
                  x = product(cooking_fuel, cooking_location_mapped), 
                  fill=cooking_fuel), na.rm=T) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 65))
dev.off()

##chimney
pdf(file=file.path(out.dir, 'global_cooking_fuel_v_cooking_chimney.pdf'), onefile=T, width=11, height=8)
ggplot(data = plot.dt[!(cooking_fuel %in% c('missing', 'other', 'unknown'))]) +
  geom_mosaic(aes(weight = N,
                  x = product(cooking_fuel, cooking_type_chimney_mapped), 
                  fill=cooking_fuel), na.rm=T) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 65))
dev.off()

#build a function to specify the mosaic plots
mosaicPlot <- function(region, var1, var2, dt, facet_var=NA, fill_var, fill_vals, wt, title) {
  
  plot <- ggplot(data = dt[super_region_name == region]) +
    geom_mosaic(aes_string(weight = wt,
                           x = paste0('product(', var1, ',', var2,')'), 
                           fill=fill_var), na.rm=T) +
    scale_fill_manual(values = fill_vals) +
    ggtitle(paste0(title, region)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 65))
  
  #facet if requested
  if (!missing(facet_var)) plot <- plot + facet_grid(.~risk)
  
  return(plot)
  
}

##create plots colored by cooking risk after definition
#super regions
pdf(file=file.path(out.dir, 'regional_cooking_fuel_v_cooking_type.pdf'), onefile=T, width=11, height=8)
mclapply(unique(plot.dt[super_region_name != 'High-income', super_region_name]), mosaicPlot, 
         dt=plot.dt[!(cooking_fuel %in% c('missing', 'other', 'unknown'))],
         var1='cooking_fuel', var2='cooking_type_mapped', fill_var='cooking_fuel', fill_vals=colors, wt='N',
         title="Cooking fuel vs cooking type: ",
         mc.cores=cores) 
dev.off()

pdf(file=file.path(out.dir, 'regional_cooking_fuel_v_cooking_location.pdf'), onefile=T, width=11, height=8)
mclapply(unique(plot.dt[super_region_name != 'High-income', super_region_name]), mosaicPlot, 
         dt=plot.dt[!(cooking_fuel %in% c('missing', 'other', 'unknown'))], 
         var1='cooking_fuel', var2='cooking_location_mapped', fill_var='cooking_fuel', fill_vals=colors, wt='N',
         title="Cooking fuel vs cooking location: ",
         mc.cores=cores) 
dev.off()

pdf(file=file.path(out.dir, 'regional_cooking_fuel_v_cooking_chimney.pdf'), onefile=T, width=11, height=8)
mclapply(unique(plot.dt[super_region_name != 'High-income', super_region_name]), mosaicPlot, 
         dt=plot.dt[!(cooking_fuel %in% c('missing', 'other', 'unknown'))], 
         var1='cooking_fuel', var2='cooking_type_chimney_mapped', fill_var='cooking_fuel', fill_vals=colors, wt='N',
         title="Cooking fuel vs cooking chimney: ",
         mc.cores=cores) 
dev.off()

#subset to MACRO
pdf(file=file.path(out.dir, 'regional_cooking_fuel_v_cooking_type_macro.pdf'), onefile=T, width=11, height=8)
mclapply(unique(plot.dt[super_region_name != 'High-income', super_region_name]), mosaicPlot, 
         dt=plot.dt[!(cooking_fuel %in% c('missing', 'other', 'unknown')) & survey_series %like% 'MACRO'],
         var1='cooking_fuel', var2='cooking_type_mapped', fill_var='cooking_fuel', fill_vals=colors, wt='N',
         title="Cooking fuel vs cooking type: ",
         mc.cores=cores) 
dev.off()

pdf(file=file.path(out.dir, 'regional_cooking_fuel_v_cooking_location_macro.pdf'), onefile=T, width=11, height=8)
mclapply(unique(plot.dt[super_region_name != 'High-income', super_region_name]), mosaicPlot, 
         dt=plot.dt[!(cooking_fuel %in% c('missing', 'other', 'unknown')) & survey_series %like% 'MACRO'],
         var1='cooking_fuel', var2='cooking_location_mapped', fill_var='cooking_fuel', fill_vals=colors, wt='N',
         title="Cooking fuel vs cooking location: ",
         mc.cores=cores) 
dev.off()

pdf(file=file.path(out.dir, 'regional_cooking_fuel_v_cooking_chimney_macro.pdf'), onefile=T, width=11, height=8)
mclapply(unique(plot.dt[super_region_name != 'High-income', super_region_name]), mosaicPlot, 
         dt=plot.dt[!(cooking_fuel %in% c('missing', 'other', 'unknown')) & survey_series %like% 'MACRO'],
         var1='cooking_fuel', var2='cooking_type_chimney_mapped', fill_var='cooking_fuel', fill_vals=colors, wt='N',
         title="Cooking fuel vs cooking chimney: ",
         mc.cores=cores) 
dev.off()

##create plots colored by cooking risk after definition
#super regions
# pdf(file=file.path(out.dir, 'regional_cooking_fuel_v_cooking_type_by_risk.pdf'), onefile=T, width=11, height=8)
# mclapply(unique(plot.dt[super_region_name != 'High-income', super_region_name]), mosaicPlot, 
#          var1='cooking_fuel', var2='risk', fill_var='risk', fill_vals=risk.colors, wt='hh_size',
#          title="Cooking fuel vs cooking type: ",
#          mc.cores=cores) 
# dev.off()
# 
# pdf(file=file.path(out.dir, 'regional_cooking_fuel_v_cooking_location_by_risk.pdf'), onefile=T, width=11, height=8)
# mclapply(unique(plot.dt[super_region_name != 'High-income', super_region_name]), mosaicPlot, 
#          var1='cooking_fuel', var2='risk', fill_var='risk', fill_vals=risk.colors, wt='hh_size',
#          title="Cooking fuel vs cooking location: ",
#          mc.cores=cores) 
# dev.off()
# 
# pdf(file=file.path(out.dir, 'regional_cooking_fuel_v_cooking_chimney_by_risk.pdf'), onefile=T, width=11, height=8)
# mclapply(unique(plot.dt[super_region_name != 'High-income', super_region_name]), mosaicPlot, 
#          var1='cooking_fuel', var2='risk', fill_var='risk', fill_vals=risk.colors, wt='hh_size',
#          title="Cooking fuel vs cooking chimney: ",
#          mc.cores=cores) 
# dev.off()
#***********************************************************************************************************************

# ---WORD CLOUDS FUELTYPE-----------------------------------------------------------------------------------------------
#create function to read in raw data and collapse the number of unique str match combos
readCollapseStrings <- function(file, varlist) {
  
  message(file)
  raw <- fread(file)
  
  #helper fx to loop over varlist
  varLoop <- function(var, input.dt) {
    
    dt <- copy(input.dt)
  
    var.mapped <- paste0(var, '_mapped')
  
    if (any(names(dt)==var.mapped)) { #make sure that var exists in raw data first

      setnames(dt, c(var, var.mapped), c('var_og', 'var_mapped'))
      setkey(dt, var_og, var_mapped)
      dt <- dt[!(is.na(var_og) | var_og == ""),
               .(nid, ihme_loc_id, int_year, survey_name, survey_module, var_og, var_mapped)] 
      dt[, count := .N, by=key(dt)] #generate count of string pairs
      dt[, N := .N] #generate total
      dt[, prop := count/N] #generate proportions
      dt[, var := var]
      unique(dt, by=key(dt)) %>% return
    
    } else return(NULL)
  
  }
  
  lapply(varlist, varLoop, input.dt=raw) %>% 
    rbindlist %>% 
    return
  
}

if (new.extracts==T) {
  
  cooking.vars <- c('cooking_fuel', 'cooking_type', 'cooking_type_chimney', 'cooking_location')
  
  strings <- file.path(raw.dir) %>% 
    list.files(full.names = T, pattern='.csv') %>% 
    mclapply(., readCollapseStrings, varlist=cooking.vars,
             mc.cores=10) %>% 
    rbindlist
  
  write.csv(strings, file=file.path(doc.dir, 'cooking_string_match_tabulations.csv'), row.names = F)

} else strings <- file.path(doc.dir, 'cooking_string_match_tabulations.csv') %>% fread

#cleanup
strings[var=='cooking_fuel' & var_mapped == 'biomass', var_mapped := 'crop'] #TODO fix biomass in the original coding
strings[var=='cooking_fuel' & var_mapped == '', var_mapped := 'missing'] #TODO fix biomass in the original coding

#set up order of fuel types
fuel.order <- c('none', 'electricity', 'gas', 'kerosene', 'coal', 'wood', 'crop', 'dung', 'other', 'unknown', 'missing')
loc.order <- c('outside', 'kitchen', 'inside', 'other')
type.order <- c('open_chimney', 'closed', 'open', 'other')

#set up color scales
fuel.colors <- c(plasma(8, direction=-1), "#C0C0C0", "#C0C0C0", "#C0C0C0") #use gray for other and unknown
loc.colors <- c(plasma(3, direction=-1), "#C0C0C0") #use gray for other 
type.colors <- c(plasma(3, direction=-1), "#C0C0C0") #use gray for other 
names(fuel.colors) <- fuel.order
names(loc.colors) <- loc.order
names(type.colors) <- type.order

#merge locs to use region/super region
strings <- merge(strings, locs[, .(ihme_loc_id, region_name, super_region_name)], by='ihme_loc_id')
strings[, iso3 := substr(ihme_loc_id, start=1, stop=3) %>% as.factor]

#plot word cloud by region
wordCloudRegion <- function(sr, str.dt, this.var, order, colors) {
  
  message(sr)
  
  #subset data
  dt <- str.dt[super_region_name %like% sr & var==this.var]
  dt[, factor_mapped := factor(var_mapped, levels=order)]
  
  #collapse again by region
  setkey(dt, region_name, var_og, var_mapped)
  dt[, count := .N, by=key(dt)] #generate count of string pairs
  dt[, N := .N] #generate total
  dt[, prop := count/N] #generate proportions
  dt <- unique(dt, by=key(dt))
  
  plot <- 
  ggplot(data=dt) +
    aes(x = 1, y = 1, size = log(count), label = var_og,
        color= factor_mapped) +
    geom_text_wordcloud_area()+
    scale_size(range = c(2, 10), guide = FALSE) +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks = NULL) +
    labs(x = '', y = '') +
    facet_grid(~region_name) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    theme(strip.text = element_text(
      color="black", size=16, lineheight=5.0),
      plot.title = element_text(colour = "black",
                                size = 18,
                                hjust = 0.5, vjust = 0.8, angle = 0))
  
  print(plot)
  
  return(NULL)
  
}

pdf(file=file.path(out.dir, 'fueltype_word_cloud.pdf'), onefile=T, width=11, height=8)
lapply(unique(plot.dt[, super_region_name]), wordCloudRegion, str.dt=plot.dt) 
dev.off()
#***********************************************************************************************************************

# ---GRAPH CATS (HOUSING)-----------------------------------------------------------------------------------------------
var.fam <- 'housing'
housing <- file.path(share.dir, var.fam) %>% 
  list.files(full.names = T, pattern='uncollapsed') %>% 
  lapply(., function(x) as.data.table(readRDS(x))) %>% 
  rbindlist

#create row ID for vectorized functions
housing[, index := .I]
setkey(housing, index)

#merge location info
housing <- merge(housing, locs[, .(ihme_loc_id, region_name, super_region_name)], by='ihme_loc_id')

#extra the first digit of floor number as an ordinal indicator
housing[, floor_rank := substr(as.character(housing_floor_num), 0, 1)]

##data prep
#clean the strings a bit
housing[, floor := tolower(housing_floor)]
housing[, floor := str_replace_all(floor, '[.]', '')] #remove common errors
housing[, floor := str_replace_all(floor, '<ff>', '')] #remove common errors
housing[, floor := str_replace_all(floor, '<fb>', '')] #remove common errors
housing[, floor := str_replace_all(floor, '<e0>', '')] #remove common errors
housing[, floor := str_replace_all(floor, '<e1>', '')] #remove common errors
housing[, floor := str_replace_all(floor, '<e3>', '')] #remove common errors
housing[, floor := str_replace_all(floor, '<e9>', '')] #remove common errors
housing[, floor := str_replace_all(floor, '<a4>', '')] #remove common errors
housing[, floor := str_replace_all(floor, '<f1>', '')] #remove common errors
housing[, floor := str_replace_all(floor, ' /', ' ')] #remove common errors
housing[, floor := str_replace_all(floor, '/ ', ' ')] #remove common errors
housing[, floor := str_replace_all(floor, '/', ' ')] #remove common errors
housing[, floor := str_replace_all(floor, ' - ', ' ')] #remove common errors
housing[, floor := str_replace_all(floor, ': ', ' ')] #remove common errors
housing[, floor := str_replace_all(floor, '["]', '')] #remove common errors
housing[, floor := str_replace_all(floor, ',', ' ')] #remove common errors

#tabulate the amount of each floor type using hh_size
setkey(housing, ihme_loc_id, floor)
housing[, floor_N := sum(hh_size, na.rm=T), by=key(housing)]
plot.dt <- unique(housing, by=key(housing))
plot.dt <- plot.dt[floor_N > 350] # keep only the more popular types

#plot word cloud by region
wordCloudRegion <- function(region) {
  
  message(region)
  
ggplot(data=plot.dt[region_name %like% region & floor_rank <4]) +
  aes(x = 1, y = 1, size = floor_N, label = floor,
      color= floor_rank) +
  geom_text_repel(segment.size = 0, force = 100) +
  scale_size(range = c(2, 10), guide = FALSE) +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks = NULL) +
  labs(x = '', y = '') +
  facet_grid(~ihme_loc_id) +
  theme_minimal() +
  theme(strip.text = element_text(
    color="black", size=16, lineheight=5.0),
    plot.title = element_text(colour = "black",
                              size = 18,
                              hjust = 0.5, vjust = 0.8, angle = 0)) +
  ggtitle(paste0("Floor Type,", region))
  
}

pdf(file=file.path(out.dir, 'floor_word_cloud.pdf'), onefile=T, width=11, height=8)
lapply(unique(plot.dt[, region_name]), wordCloudRegion) 
dev.off()

cat.order <- c('0', '1', '2', '3', '9')
plot.dt[, floor_rank := factor(floor_rank, levels=cat.order)]

#set up color scale
colors <- c(plasma(4, direction=1), "#C0C0C0") #use gray for other and unknown
names(colors) <- cat.order

#plot word cloud by iso3
wordCloudRegion <- function(region) {
  
  message(region)
  
  ggplot(data=plot.dt[ihme_loc_id %like% region & floor_rank %in% cat.order]) +
    aes(x = 1, y = 1, size = floor_N, label = floor,
        color= floor_rank) +
    geom_text_repel(segment.size = 0, force = 100) +
    scale_size(range = c(2, 10), guide = FALSE) +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks = NULL) +
    scale_color_manual(values = colors) +
    labs(x = '', y = '') +
    facet_grid(~ihme_loc_id) +
    theme_minimal() +
    theme(strip.text = element_text(
      color="black", size=16, lineheight=5.0),
      plot.title = element_text(colour = "black",
                                size = 18,
                                hjust = 0.5, vjust = 0.8, angle = 0)) +
    ggtitle(paste0("Floor Type,", region))
  
}

pdf(file=file.path(out.dir, 'floor_word_cloud_iso3.pdf'), onefile=T, width=11, height=8)
lapply(unique(plot.dt[!is.na(floor), ihme_loc_id]), wordCloudRegion) 
dev.off()


##SCRAP
plot.dt <- housing[!is.na(hh_size)]

wordcloud(words = plot.dt$housing_floor, freq = plot.dt$hh_size, min.freq = 1, #scale = c(2, 0.2),
          max.words=200, random.order=FALSE, rot.per=0.1, 
          ordered.colors=TRUE,
          colors=brewer.pal(8, "Dark2")[factor(plot.dt$floor_rank)])

ggplot(data = housing[survey_series %like% 'MACRO' & floor_rank <= 3]) +
  geom_mosaic(aes(weight = hh_size,
                  x = product(housing_floor, floor_rank), 
                  fill=floor_rank), na.rm=T)
