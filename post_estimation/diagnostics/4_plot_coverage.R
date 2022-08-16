########################################################################
###### Short loop to generate data coverage plots for all regions ######
######                Scott Swartz | 8/29/2017                    ######
########################################################################

#source("/homes/jfrostad/_code/lbd/hap/extract/4_plot_coverage.R", echo=T)

# ---CONFIG----------------------------------------------------------------------------------------------------------------------
# clear memory
rm(list=ls())

# disable scientific notation
options(scipen = 999)

library(RMySQL)
library(rgeos)
library(feather)
library(data.table)
library(ggplot2)
library(doParallel)
library(gridExtra)
library(stringr)
library(RColorBrewer)
library(rgdal)
library(raster)
library(magrittr)
library(dplyr)
library(tidyr)

#setup
sing_image = TRUE #whether running on RStudio singularity image
numcores = 3
root <- ifelse(Sys.info()[1]=="Windows", "J:/", "/home/j/")
j_root <- root
repo <- '/share/code/geospatial/kwilson7/lbd_core/'
setwd(repo)
#source('mbg_central/polygon_functions.R')
#source('mbg_central/shapefile_functions.R')

indicator <- 'hap' #water or sani
var <- 'cooking_fuel_solid' #imp, unimp, surface, od, piped

title <- "Cooking Fuel"

date <- "2019_03_26"
data.dir  <- file.path('/share/limited_use/LIMITED_USE/LU_GEOSPATIAL/collapse/hap/')
coverage_data <-  paste0(data.dir, "data_cooking_", date, ".feather") %>% read_feather %>% as.data.table

#import dataset (written during collapse code)
coverage_data <- coverage_data[!(shapefile == 'mombasa' & is.na(lat)),]
coverage_data$point <- ifelse(is.na(coverage_data$lat), 0, 1)

#rename some vars
setnames(coverage_data,
         c('iso3', 'lat', 'long', 'year_start', 'survey_series'),
         c('country', 'latitude', 'longitude', 'year', 'source'))

#create indicator
coverage_data[, cooking_fuel := N * get(var)]


#run loop
source('mbg_central/graph_data_coverage.R')
regions <- c('africa', 'latin_america', 'middle_east','south_asia','se_asia')
#regions <- c('south_asia','se_asia')
#regions <- 'se_asia'

#fix problem with a se asia shapefile (naming issue)
coverage_data <- coverage_data[shapefile == "lg_g2015_2007_1", shapefile := "lf_g2015_2007_1"]
#drop weird shapefiles for now
#TODO investigate these issues
coverage_data <- coverage_data[!(shapefile %like% "2021")] 
#coverage_data <-coverage_data[!(shapefile %like% "gadm_3_4_vnm_adm3")]
coverage_data <-coverage_data[!(shapefile %like% "Chapineiro_COL_CS")]
coverage_data <-coverage_data[!(shapefile %like% "PRY_central_wo_asuncion")]
coverage_data <- coverage_data[!(country %like% "MDV")] 

for (reg in regions[-1]){
  coverage_maps <- graph_data_coverage_values(df = coverage_data,
                                              var = var,
                                              title = title,
                                              legend_title = "Prevalence",
                                              year_min = 1998,
                                              year_max = 2017,
                                              year_var = 'year',
                                              region = reg,

                                              cores = numcores,
                                              indicator = indicator,

                                              extra_file_tag = '',
                                              save_on_share = FALSE,
                                              out_dir = NULL,
                                              core_repo = repo,
                                              log_dir = NULL,

                                              fast_shapefiles = TRUE,
                                              new_data_plots = FALSE,
                                              since_date = NULL,
                                              annual_period_maps = FALSE,
                                              save_period_maps = TRUE,
                                              prep_shiny = FALSE,
                                              return_maps = TRUE,
                                              debug = FALSE,

                                              color_scheme = "classic",
                                              color_scheme_scatter = "brewer",
                                              high_is_bad = TRUE,
                                              cap = 90,
                                              cap_type = "percentile",
                                              legend_min = 0,
                                              legend_max = 1,
                                              endemic_gauls = NULL,
                                              stage_3_gray = TRUE,
                                              simplify_polys = TRUE,
                                              tolerance = 0.03,
                                              base_font_size = 18,
                                              map_point_size = 0.8,
                                              poly_line_width = 0.2,

                                              remove_rank = TRUE
  )
}
