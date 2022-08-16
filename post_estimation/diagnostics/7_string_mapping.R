# ----HEADER------------------------------------------------------------------------------------------------------------
# Author: JF
# Date: 09/11/2018 (never forget)
# Purpose: Exploring data coverage across indicators for hap
# source("/homes/jfrostad/_code/lbd/hap/post_estimation/diagnostics/7_string_mapping.R", echo=T)
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
               questionr, parallel, raster, RColorBrewer, readr, readxl, rgdal, rgeos, 
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
new.extracts <- T
topic <- "hap"
this.family <- 'cooking'
indicators <- c('cooking_fuel', 'cooking_type', 'cooking_type_chimney', 'cooking_location', 
                'heating_fuel', 'heating_type', 'heating_type_chimney', 'lighting_fuel', 'electricity', 
                'housing_roof', 'housing_wall', 'housing_floor')
#***********************************************************************************************************************

# ----IN/OUT------------------------------------------------------------------------------------------------------------
###Input###
#raw data
doc.dir <- file.path(j_root, 'WORK/11_geospatial/hap/documentation/str_review')
def.file <- file.path(doc.dir, 'definitions.xlsx')
raw.dir <- file.path("/ihme/limited_use/LIMITED_USE/LU_GEOSPATIAL/ubCov_extractions/hap/")
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
# file.path(gbd.shared.function.dir, 'get_covariate_estimates.R') %>% source
# file.path(gbd.shared.function.archive, 'get_draws.R') %>% source
# file.path(gbd.shared.function.dir, 'get_location_metadata.R') %>% source
# file.path(gbd.shared.function.dir, 'get_ids.R') %>% source
# file.path(gbd.shared.function.dir, 'get_population.R') %>% source
# file.path(gbd.shared.function.dir, 'get_outputs.R') %>% source

#create function to read in raw data and collapse the number of unique str match combos
readCollapseStrings <- function(file, varlist) {
  
  message(file)
  raw <- fread(file, encoding = 'UTF-8')
  
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
#***********************************************************************************************************************

# ---PREP DATA----------------------------------------------------------------------------------------------------------
codebook <- read_xlsx(file.path(raw.dir, 'hap.xlsx'), sheet='codebook') %>% as.data.table
codebook <- codebook[assigned=='qnguyen1' | assigned == 'jfrostad' | assigned == 'albrja']

#add on location hierarchy info
locs <- get_location_hierarchy(41)

#if necessary, compile all the raw data tabulations
if (new.extracts==T) {
  
  cooking.vars <- c('cooking_fuel', 'cooking_type', 'cooking_type_chimney', 'cooking_location')
  
  strings <- file.path(raw.dir) %>% 
    list.files(full.names = T, pattern='.csv') %>% 
    mclapply(., readCollapseStrings, varlist=cooking.vars,
             mc.cores=10) %>% 
    rbindlist

  write_excel_csv(strings, path=file.path(doc.dir, 'cooking_string_match_tabulations.csv'))
  
} else strings <- file.path(doc.dir, 'cooking_string_match_tabulations.csv') %>% fread

#cleanup
#strings[var=='cooking_fuel' & var_mapped == 'biomass', var_mapped := 'crop'] #TODO fix biomass in the original coding
strings[var=='cooking_fuel' & var_mapped == '', var_mapped := 'missing'] #TODO fix biomass in the original coding

#examine these tabulations to see if any strings are unmapped
strings[var=='cooking_fuel' & var_mapped=='missing', table(var_og)]
strings[var=='cooking_fuel' & var_mapped=='missing', ][order(-count)] %>% 
  write_excel_csv(., path=file.path(doc.dir, 'cooking_fuel_missing_strings.csv'))
#***********************************************************************************************************************

# ---WORD CLOUDS FUELTYPE-----------------------------------------------------------------------------------------------
#set up order of fuel types then produce color scales
fuel.order <- c('none', 'electricity', 'gas', 'kerosene', 'coal', 'wood', 'crop', 'dung', 'other', 'unknown', 'missing')
fuel.colors <- c(plasma(8, direction=-1), "#C0C0C0", "#C0C0C0", "#C0C0C0") #use gray for other and unknown
names(fuel.colors) <- fuel.order

#set up order of cooking locations then produce color scales
loc.order <- c('outside', 'kitchen', 'inside', 'other')
loc.colors <- c(plasma(3, direction=-1), "#C0C0C0") #use gray for other
names(loc.colors) <- loc.order

#set up order of cooking types then produce color scales
type.order <- c('open_chimney', 'closed', 'open', 'other')
type.colors <- c(plasma(3, direction=-1), "#C0C0C0") #use gray for other 
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
lapply(unique(plot.dt[, super_region_name]), wordCloudRegion, str.dt=strings,
       this.var='cooking_location', order=fuel.order, colors=fuel.colors) 
dev.off()

pdf(file=file.path(out.dir, 'cooking_location_word_cloud.pdf'), onefile=T, width=11, height=8)
lapply(unique(plot.dt[, super_region_name]), wordCloudRegion, str.dt=strings,
       this.var='cooking_location', order=loc.order, colors=loc.colors) 
dev.off()
#***********************************************************************************************************************