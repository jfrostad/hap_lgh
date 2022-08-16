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
user            <- Sys.info()['user']

#load packages
package_lib    <- sprintf('%s_code/_lib/pkg_slurm',h_root)
## Load libraries and  MBG project functions.
.libPaths(package_lib)
pacman::p_load(data.table, fst, scales, ggplot2, ggpubr, ggridges, ggrepel, gridExtra, isoband, RColorBrewer, 
               sf, viridis, farver, reldist, ggnewscale,
               raster, dismo, gbm, foreign,
               parallel, doParallel, grid, gridExtra, gtools, glmnet, RMySQL, plyr,
               tictoc, dplyr, rgdal) 
package_list <- package_list[!(package_list%like%'seeg|INLA')]


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

indicator_group <- 'cooking'
indicator <- 'hap'
type <- 'mean'
raked <- T
start_year <- 2000
end_year <- 2018
cores <- 10
modeling_shapefile_version <- "2019_09_10"
#***********************************************************************************************************************

# ----IN/OUT------------------------------------------------------------------------------------------------------------
#***********************************************************************************************************************

# ---FUNCTIONS----------------------------------------------------------------------------------------------------------
##function lib##
#PE functions#
file.path(my_repo, '_lib', 'viz', 'mapping_fx.R') %>% source
file.path(my_repo, '_lib', 'viz', 'figure_fx.R') %>% source

#gbd fx
gbd.shared.function.dir <- '/ihme/cc_resources/libraries/current/r/'
file.path(gbd.shared.function.dir, 'get_location_metadata.R') %>% source

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

#scale transformations
sqrt_signed <- scales::trans_new("signed_log",
                                 transform=function(x) sign(x)*sqrt(abs(x)),
                                 inverse=function(x) sign(x)*(abs(x))^2)
#***********************************************************************************************************************

# ---PREP DATA----------------------------------------------------------------------------------------------------------
#read in the proper annotations (borders, lakes, mask)
#read in the proper annotations (borders, lakes, mask)
annotations_path <- file.path(out.dir, 'annotations.RDs')
check <- file.exists(annotations_path)
if(check) {
  annotations <- readRDS(annotations_path)
} else {
  annotations <- load_map_annotations()
  saveRDS(annotations, file=annotations_path)
}

#read in link_table
global_link_table <- file.path(global_link_dir, "lbd_full_link.rds") %>% readRDS %>% as.data.table
adm_links <- global_link_table[, .(ADM0_NAME, ADM0_CODE, ADM1_NAME, ADM1_CODE, ADM2_NAME, ADM2_CODE)] %>% unique

#read in shps
adm2 <- rbind(stage1, stage2)

#read in results
results <- file.path(data.dir, 'all_summary.fst') %>% read_fst(as.data.table=T) 
dt <- results[dimension=='ad2' & term=='lvl'] %>% Filter(function(x) !all(is.na(x)), .)
dt_d <- results[dimension=='ad2' & term%like%'change'] %>% Filter(function(x) !all(is.na(x)), .)

#merge sr region names/IDs
locs <- get_location_metadata(location_set_id = 35, gbd_round_id = 6) %>% 
  .[, .(iso3=ihme_loc_id, location_name, super_region_id, super_region_name, region_id, region_name)] #subset to relevant columns

#create file to crosswalk AD0 to iso3
iso3_map <- dplyr::select(adm2, iso3, ADM0_CODE=gadm_geoid)
iso3_map$geometry <- NULL
iso3_map <- as.data.table(iso3_map) %>% unique
locs <- merge(locs, iso3_map, by='iso3')

# #merge sr region names/IDs
dt <- merge(dt, locs, by='ADM0_CODE', all.x=T)
dt_d <- merge(dt_d, locs, by='ADM0_CODE', all.x=T)

#read in input data and prepare it for mapping
data <- load_map_results(indicator, indicator_group, run_date, raked, 
                         year_list=c(2000:2018),
                         custom_path = list('admin2'=dt),
                         geo_levels=c('admin2'),
                         cores=cores)
data_d <-
  load_map_results(indicator, indicator_group, run_date, raked, 
                   year_list=2018,
                   custom_path = list('admin2'=dt_d),
                   geo_levels=c('admin2'),
                   cores=cores)
#setup the list of top countries
#defined by population in 2018
biggest_countries <- 
  dt[year==max(dt$year) & dimension=='ad2', .(sum=sum(pop_total, na.rm=T)), by=.(iso3)] %>%
  .[order(sum)] %>%
  tail(10) %>%
  .[, unique(iso3)]

#top country per GBD sregion
biggest_countries_sr <- #defined based on pop
  dt[year==min(dt$year), .(sum=sum(pop_total, na.rm=T)), by=.(iso3,region_name)] %>%
  .[, .SD[which.max(sum)], by=region_name] %>% 
  .[order(-sum)]

#defined based on LRI rates/counts
top_countries <- 
  dt[year==min(dt$year), .(mean=weighted.mean(rate_mean, w=pop, na.rm=T)), by=.(iso3)] %>%
  .[order(mean)] %>%
  tail(10) %>%
  .[, unique(iso3)]

top_countries_c <-
  dt[year==min(dt$year), .(sum=sum(count_mean, na.rm=T)), by=.(iso3)] %>%
  .[order(sum)] %>%
  tail(14) %>%
  .[, unique(iso3)]

#top country per GBD region
top_countries_gbdreg <- #defined based on LRI counts
  dt[year==min(dt$year), .(sum=sum(count_mean, na.rm=T)), by=.(iso3,region_name)] %>%
  .[, .SD[which.max(sum)], by=region_name]

#top country per MBG region
regs <- load_adm0_lookup_table() %>% .[,.(mbg_reg_name=reg_name, mbg_reg, iso3=toupper(iso3))] %>% unique
dt <- merge(dt, regs, by='iso3')
top_countries_mbgreg <- #defined based on LRI counts
  dt[year==min(dt$year), .(sum=sum(count_mean, na.rm=T)), by=.(iso3,mbg_reg)] %>%
  .[, .SD[which.max(sum)], by=mbg_reg]

second_countries_reg <- #defined based on LRI counts
  dt[year==min(dt$year), .(sum=sum(count_mean, na.rm=T)), by=.(iso3,mbg_reg_name)] %>% 
  .[!(iso3 %in% unique(top_countries_gbdreg$iso3))] %>% 
  .[, .SD[which.max(sum)], by=mbg_reg_name]

#define extent of map
zoom.afr <- data.table(x1=-10, x2=50, y1=-20, y2=40)
zoom.global <- data.table(x1=-120, x2=150, y1=-40, y2=55)
#***********************************************************************************************************************

# ---SUPERREGION LEGEND-------------------------------------------------------------------------------------------------
#create region colors
reg_colors <- c('4'='#1f78b4',
                '31'='#a65628',
                '103'='#ff7f00',
                '137'='#984ea3',
                '158'='#e31a1c',
                '166'='#4daf4a')

#create data to make regional map color legend
reg_data <- annotations$adm0 %>% 
  merge(., locs[ADM0_CODE %in% unique(results$ADM0_CODE), .(sreg=super_region_id %>% as.factor, ADM0_CODE)], 
        by.y='ADM0_CODE', by.x='geo_id', all.x=T)

annotations_reg <- lapply(annotations, st_crop, 
                          xmin=zoom.global$x1, xmax=zoom.global$x2, ymin=zoom.global$y1, ymax=zoom.global$y2)

#make a map showing all the region colors in place
reg_map <- ggplot() + geom_sf(data = annotations_reg$adm0, lwd=0.1, color = 'black', fill = 'gray60')
reg_map <- reg_map + geom_sf(data = reg_data, aes(fill = sreg), lwd=0) + coord_sf(datum = NA)
reg_map <- reg_map + geom_sf(data = annotations_reg$adm0, lwd=0.1, color = 'black', fill=NA)
reg_map <- reg_map + geom_sf(data = annotations_reg$lakes, lwd=0, color = 'gray60', fill = 'lightblue')
reg_map <- reg_map + geom_sf(data = annotations$stage3, lwd=0, color = 'gray60', fill = 'gray60')
reg_map <- reg_map + scale_fill_manual(values=reg_colors, guide=F,
                                       na.value = "gray60")
reg_map <- reg_map +
  labs(x="", y="", title='') +
  theme_classic() +
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(hjust=0.5), plot.margin=unit(c(0, 0, 0, 0), "in"))
reg_grob <- ggplotGrob(reg_map)

ggsave(filename=file.path(save.dir, 'region_legend.eps'), plot=reg_map, width=12, height=6, units='in', dpi=900)

#***********************************************************************************************************************
 
# ---FIGURE 2-----------------------------------------------------------------------------------------------------------
#SDG projection probabilities
#append the ADM2 files
sdg_files <-
  file.path(data.dir, 'sdg_projections') %>% list.files(pattern='admin_2', full.names = T) %>% 
  lapply(., readRDS)

#extract goal obj to index over
goals <- lapply(1:length(sdg_files), function(x) sdg_files[[x]]$goals) %>% rbindlist %>% unique

#create a dt with all probabilities
probs <- lapply(1:nrow(goals), prepCasts, type='probs', id_dt=goals) %>% 
  rbindlist %>% 
  setnames(c('target_year', 'spatial_idx'), c('year', 'ADM2_CODE')) %>% 
  merge(., adm_links, by='ADM2_CODE')  %>% 
  merge(., locs, by='ADM0_CODE')

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
  merge(., adm_links, by='ADM2_CODE') %>% 
  merge(., locs, by='ADM0_CODE')

#create plot of the divide by country
target_threshold <- .05
prob_threshold <- .95

plot.dt <- 
probs[target==target_threshold & year==2030] %>% 
  merge(., dt[type=='TAP'&year==max(year), .(ADM2_CODE, pop_total)], by='ADM2_CODE') %>% 
  .[absolute_goal_prob<=.5, status := 'may fail'] %>% 
  .[absolute_goal_prob<=(1-prob_threshold), status := 'fail'] %>% 
  .[absolute_goal_prob>.5, status := 'may succeed'] %>% 
  .[absolute_goal_prob>=prob_threshold, status := 'success'] %>% 
  .[, country_pop := sum(pop_total, na.rm=T), by=.(iso3, ADM0_NAME)] %>% 
  .[, pop_fail := sum((absolute_goal_prob<=.5)*pop_total, na.rm=T), by=.(iso3, ADM0_NAME)] %>% 
  .[, pop_success := sum((absolute_goal_prob>.5)*pop_total, na.rm=T), by=.(iso3, ADM0_NAME)] %>% 
  .[, .(pop_share=sum(pop_total, na.rm=T)/country_pop,
        pop_fail,
        pop_success
        ), 
    by=.(iso3, ADM0_NAME, status, super_region_id)] %>% 
  unique(., by=c('ADM0_NAME', 'status', 'super_region_id')) %>% 
  .[, ratio := pop_fail/pop_success] %>% 
  .[status%like%'fail', pop_share := pop_share*-1] %>% #create mirroring effect
  .[order(ratio)] %>% 
  .[, status := factor(status, levels=c('success', 'may succeed', 'fail', 'may fail'))] 

plot <-
plot.dt %>%
  ggplot(aes(x = forcats::fct_reorder(iso3, -ratio), 
             y = pop_share, alpha = status, fill=super_region_id %>% as.factor)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept=0) +
  scale_fill_manual(values=reg_colors, guide=F) +
  #scale_color_manual(values=c('success'='navy', 'may_succeed'='gray', 'may_fail'='gray', 'fail'='red'), guide=F) +
  scale_alpha_manual(values=c('success'=.8, 'may succeed'=.2, 'may fail'=.2, 'fail'=.8), guide=F) +
  scale_x_discrete('') +
  scale_y_continuous('',
                     breaks=c(-1, -.5, 0, .5, 1),
                     labels=c("(100%)", "(50%)", "0%",  "50%", "100%")) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust=1)
  )

#plot <- plot + annotation_custom(grob=reg_grob, ymin=80, ymax=115, xmin=25, xmax=34)

ggsave(filename=file.path(save.dir, 'fig_2_all.png'), plot=plot, 
       width=12, height=8, units='in', dpi=900)

plot <-
  plot.dt[ratio!=0 & !is.infinite(ratio)] %>%
  ggplot(aes(x = forcats::fct_reorder(iso3, -ratio), 
             y = pop_share, alpha = status, fill=super_region_id %>% as.factor)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept=0) +
  scale_fill_manual(values=reg_colors, guide=F) +
  #scale_color_manual(values=c('success'='navy', 'may_succeed'='gray', 'may_fail'='gray', 'fail'='red'), guide=F) +
  scale_alpha_manual(values=c('success'=.8, 'may succeed'=.2, 'may fail'=.2, 'fail'=.8), guide=F) +
  scale_x_discrete('') +
  scale_y_continuous('',
                     breaks=c(-1, -.5, 0, .5, 1),
                     labels=c("(100%)", "(50%)", "0%",  "50%", "100%")) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 55, hjust=1, size=13)
  )

ggsave(filename=file.path(save.dir, 'fig_2.eps'), plot=plot, 
       width=12, height=8, units='in', dpi=900)

pdf(paste0(save.dir, '/fig_3_vector.pdf'),
    height=8, width=12)
plot
dev.off()
#***********************************************************************************************************************

# ---FIGURE 3-----------------------------------------------------------------------------------------------------------
#histogram/density plot of LRI deaths vs tap_pc
#setup plot data
plot.dt <- dt %>% 
  copy %>% 
  .[year %in% c(2000, 2018)] %>% 
  na.omit(., cols=c('prev_mean', 'pop')) %>% 
  .[, .(iso3, ADM0_NAME, ADM2_CODE, year, pop, pop_total, type, prev_mean, share_mean, pm_pc_mean, count_mean, atr_count_mean,
        region_id, super_region_id)] %>% 
  .[, exposed_pop := pop_total*prev_mean]

#generate the plots
pop_plot <-
  makeTapDensityPlot(input_dt=plot.dt, locs=biggest_countries_sr[1:6, iso3], 
                     wt_var='exposed_pop', 
                     tap_cutoff = 750,
                     smoother=.1)
death_plot <-
makeTapDensityPlot(input_dt=plot.dt, locs=biggest_countries_sr[1:6, iso3], 
                   wt_var='atr_count_mean', 
                   tap_cutoff = 750,
                   smoother=.1)

#extract legend and set coordinates
legend <- get_legend(pop_plot)
legend$vp$x <- unit(0, 'npc')
legend$vp$y <- unit(.69, 'npc')
#define layout
lay <- rbind(c(1,1,1,1),
             c(2,2,2,2))
plot <- arrangeGrob(grobs=list(pop_plot + theme(legend.position = 'none'), 
                               death_plot + theme(legend.position = 'none',
                                                  axis.text.x = element_blank())), 
                    layout_matrix=lay) %>% 
  grid.arrange

ggsave(filename=file.path(save.dir, 'fig_3_bare.eps'), plot=plot, width=12, height=8, units='in', dpi=1200)
ggsave(filename=file.path(save.dir, 'fig_3_legend.eps'), plot=legend, width=4, height=4, units='in', dpi=1200)

#draw plot and legend/annotations
png(paste0(out.dir, '/fig_3.png'),
    height=8, width=12, units='in', res=2400)
grid.arrange(plot)
grid.draw(legend)
grid.text('a', x = unit(0.01, "npc"), y = unit(0.95, "npc"), gp = gpar(fontsize = 24, fontface = "bold"))
grid.text('b', x = unit(0.01, "npc"), y = unit(0.05, "npc"), gp = gpar(fontsize = 24, fontface = "bold"))
grid.text('2000', x = unit(.95, "npc"), y = unit(0.7, "npc"), gp = gpar(fontsize = 24, fontface = "bold"))
grid.text('2000', x = unit(.95, "npc"), y = unit(0.2, "npc"), gp = gpar(fontsize = 24, fontface = "bold"))
grid.text('2018', x = unit(.95, "npc"), y = unit(0.76, "npc"), gp = gpar(fontsize = 24, fontface = "bold"))
grid.text('2018', x = unit(.95, "npc"), y = unit(0.29, "npc"), gp = gpar(fontsize = 24, fontface = "bold"))
dev.off()

#draw plot and legend/annotations
pdf(paste0(save.dir, '/fig_4_vector.pdf'),
    height=8, width=12)
grid.arrange(plot)
dev.off()

#also run a version that is not area-scaled to year to provide accurate probabilities for quintiles
pop_plot <-
  makeTapDensityPlot(input_dt=plot.dt, locs=biggest_countries_sr[1:6, iso3], 
                     wt_var='exposed_pop', 
                     tap_cutoff = 750,
                     smoother=.1,
                     normalize_years = F)
death_plot <-
  makeTapDensityPlot(input_dt=plot.dt, locs=biggest_countries_sr[1:6, iso3], 
                     wt_var='atr_count_mean', 
                     tap_cutoff = 750,
                     smoother=.1,
                     normalize_years = F)

#make a version of Figure 3 for each region
#generate the plots
plot.dt <- plot.dt[!(iso3 %in% biggest_countries_sr[1:6, iso3])] #remove the biggest countries
sapply(unique(plot.dt$region_id), makeFig3Loclist, dt=plot.dt) %>%
  flatten2 %>%
    mclapply(., makeRegFigure3s, mc.cores=5)
#***********************************************************************************************************************

# ---FIGURE 4-----------------------------------------------------------------------------------------------------------
#setup plot data using TAP values
plot.dt <- dt %>% 
  copy %>% 
  .[year %in% c(start_year, end_year) & type=='TAP'] %>% 
  na.omit(., cols=c('paf_mean'))

#merge on HAP shares
plot.dt <- dt %>% 
  copy %>% 
  .[year %in% c(start_year, end_year) & type=='HAP', .(year, ADM2_CODE, hap_pct_mean=share_mean)] %>% 
  merge(., plot.dt, by=c('year', 'ADM2_CODE'))

#formatting
plot.dt %>% 
  .[, rate_mean := rate_mean*1e3] %>%  #scale LRI rates to per 1000
  .[iso3=='COD', location_name := 'D.R. Congo'] %>%
  .[, loc_fct := factor(location_name)] %>% 
  .[, loc_fct := forcats::fct_rev(loc_fct)]

#drop this district, outlier without data support
plot.dt <- plot.dt[ADM2_NAME!='Lakshadweep']

#make master plots
master_plot <- makeMasterPlot(custom_countries=c('COD', 'KEN', 'IND', 'MNG', 'THA') ,
                              custom_cols=c('D.R. Congo'='#f781bf', 'Kenya'='#66a61e', 'India'='#a6761d', 
                                            'Mongolia'='#e31a1c', 'Thailand'='#1f78b4'),
                              add_rug=T, rug_limits=c(0,4))
legend <- get_legend(master_plot)

##cloud version with top 14##
#now make a single inset plot for each of the remaining top by region (and add Pakistan to round it out)
inset_countries <- c(top_countries_gbdreg[, iso3]) %>% 
  .[!(. %in% c('COD', 'KEN', 'IND', 'MNG', 'THA'))]

#duplicate one to help guide the reader
#replaced ETH
inset_countries[2] <- 'COD' 

insets <- sort(inset_countries) %>% 
  lapply(1:length(.), makeInset, loclist=., 
         scale_labels=T, type='cloud_contour', contour_bw=10,
         rug_limits=c(0,4))

#arrange into master figure
all_grobs <- copy(insets)
all_grobs[[13]] <- master_plot
lay <- rbind(c(13,13,13,13,1,2,3),
             c(13,13,13,13,4,5,6),
             c(13,13,13,13,7,8,9),
             c(13,13,13,13,10,11,12))
plot <- arrangeGrob(grobs=all_grobs, layout_matrix=lay, 
                    bottom=textGrob("Percent of TAP contributed by ambient sources", 
                                    gp = gpar(fontsize=17))
) %>% 
  grid.arrange

ggsave(plot=plot, filename=file.path(save.dir, 'fig_4.png'),
       width=12, height=8, units='in', dpi=900)
ggsave(plot=legend, filename=file.path(save.dir, 'fig_4_legend.png'),
       width=4, height=8, units='in', dpi=900)

pdf(paste0(save.dir, '/fig_5.pdf'),
    height=8, width=12)
arrangeGrob(grobs=all_grobs, layout_matrix=lay, 
                    bottom=textGrob("Percent of TAP contributed by ambient sources", 
                                    gp = gpar(fontsize=17))
) %>% 
  grid.arrange
dev.off()

#make a scatter version
insets <- sort(inset_countries) %>% 
  lapply(1:length(.), makeInset, loclist=., 
         scale_labels=T, type='scatter', contour_bw=10,
         rug_limits=c(0,4))

#arrange into master figure
all_grobs <- copy(insets)
all_grobs[[13]] <- master_plot
lay <- rbind(c(13,13,13,13,1,2,3),
             c(13,13,13,13,4,5,6),
             c(13,13,13,13,7,8,9),
             c(13,13,13,13,10,11,12))
plot <- arrangeGrob(grobs=all_grobs, layout_matrix=lay, 
                    bottom=textGrob("Percent of TAP contributed by ambient sources", 
                                    gp = gpar(fontsize=17))
) %>% 
  grid.arrange

ggsave(plot=plot, filename=file.path(save.dir, 'fig_4_scatter.png'),
       width=12, height=8, units='in', dpi=900)

#make a version of Figure 4 for each region
#generate the plots
sapply(unique(plot.dt$region_id), makeFig3Loclist, dt=plot.dt) %>% 
  flatten2 %>% 
  mclapply(., makeRegFigure4s, mc.cores=5)
#***********************************************************************************************************************

# ---EXT FIG 1----------------------------------------------------------------------------------------------------------
#set up plot data
dt_ineq <- dt[cause=='lri' & grouping=='under5' & type=='HAP' & year %in% c(start_year, end_year), 
              .(iso3, year, ADM0_CODE, ADM0_NAME, ADM2_CODE, ADM2_NAME, prev_mean, 
                super_region_id, super_region_name, region_id, region_name, pop_total)]

dt_ineq[, mean := weighted.mean(prev_mean, w=pop_total, na.rm=T), by=.(ADM0_CODE, year)]
dt_ineq[, max := max(prev_mean, na.rm=T), by=.(ADM0_CODE, year)]
dt_ineq[, min := min(prev_mean, na.rm=T), by=.(ADM0_CODE, year)]

#reorder regions
dt_ineq[, super_region := factor(super_region_id, levels=c(166, 158, 4, 31, 103, 137) %>% rev)]

#range results
dt_ineq <- unique(dt_ineq, by=c('ADM0_CODE', 'year'))
dt_ineq_end <- dt_ineq[year==max(year)]
dt_ineq_end <- dt_ineq_end[dt_ineq_end[,do.call(order, .SD), .SDcols = c('super_region', 'mean')]]
dt_ineq_end[, country := factor(iso3, levels=unique(iso3))]
#dt_ineq_start[, country := forcats::fct_rev(country)]
dt_ineq_start <- dt_ineq[year==min(year)]
dt_ineq_start[, country := factor(iso3, levels=dt_ineq_end$country)]
dt_ineq_start[, country := forcats::fct_rev(country)]

#combine
dt_ineq <- list(
  dt_ineq_start,
  dt_ineq_end
) %>% 
  rbindlist

#plot absolute inequality
plot <-
  ggplot(dt_ineq_end, aes(x=country, y=mean, ymax=max, ymin=min, fill=super_region)) +
  geom_crossbar(data=dt_ineq_start, color=NA, fill='gray80', width=.4) +
  geom_crossbar(color=NA, position=position_nudge(x=.4), width=.4) +
  geom_point(data=dt_ineq_start, shape=4, color='gray40', size=2.5) +
  geom_point(color='#000000', size=3, position=position_nudge(x=.4)) +
  scale_fill_manual(values=reg_colors, guide=F) +
  scale_y_continuous('Dirty Fuel Use', labels = scales::percent) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        text = element_text(size=16),
        axis.text.x = element_text(size=7.5, angle = 45, hjust = 1))

plot <- plot + annotation_custom(grob=reg_grob, ymin=0, ymax=.25, xmin=0, xmax=20)

ggsave(filename=file.path(save.dir, 'ex_fig_1.png'),
       width=12, height=6, units='in', dpi=900)

#try as a dotplot too
plot <-
  ggplot(dt_ineq, aes(x=country, y=mean, color=super_region, fill=super_region, 
                      shape=year %>% as.character,
                      linetype=year %>% as.character)) +
  geom_errorbar(aes(ymin = min, ymax = max)) +
  geom_point(alpha=1, size=3) +   
  scale_color_manual(values=reg_colors, guide=F) +
  scale_fill_manual(values=reg_colors, guide=F) +
  scale_shape_manual(values=c('2000'=21, '2018'=4), guide=F) +
  scale_linetype_manual(values=c('2000'='dotted', '2018'='solid'), guide=F) +
  scale_y_continuous('Dirty Fuel Use', labels = scales::percent) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        text = element_text(size=16),
        axis.text.x = element_text(size=7, angle = 60, hjust = 1)) 

plot <- plot + annotation_custom(grob=reg_grob, ymin=0, ymax=.25, xmin=0, xmax=20)

ggsave(filename=file.path(save.dir, 'ex_fig_1_lolli_a.png'),
       width=12, height=6, units='in', dpi=900)

#try as a dotplot too with grey
plot <-
  ggplot(dt_ineq, aes(x=country, y=mean, color=super_region, fill=super_region, 
                      shape=year %>% as.character,
                      linetype=year %>% as.character)) +
  geom_errorbar(aes(ymin = min, ymax = max)) +
  geom_point(alpha=1, size=3) +   
  scale_color_manual(values=reg_colors, guide=F) +
  scale_fill_manual(values=reg_colors, guide=F) +
  scale_shape_manual(values=c('2000'=21, '2018'=4), guide=F) +
  scale_linetype_manual(values=c('2000'='dotted', '2018'='solid'), guide=F) +
  scale_y_continuous('Solid Fuel Use', labels = scales::percent) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        text = element_text(size=16),
        axis.text.x = element_text(size=7, angle = 60, hjust = 1)) 

plot <- plot + annotation_custom(grob=reg_grob, ymin=0, ymax=.25, xmin=0, xmax=20)

ggsave(filename=file.path(save.dir, 'figure_2_lolli.eps'),
       width=12, height=6, units='in', dpi=900)

pdf(paste0(save.dir, '/fig_2_lolli.pdf'),
    height=8, width=12)
plot
dev.off()

plot <-
  ggplot(dt_ineq_end, aes(x=forcats::fct_rev(country), y=mean, ymax=max, ymin=min, color=super_region, fill=super_region)) +
  geom_crossbar(width=.4) +
  geom_errorbar(data=dt_ineq_start, color='gray20', linetype='dotted') +
  geom_point(data=dt_ineq_start, shape=15, color='gray30', size=2) +
  geom_point(shape=15, color='black', size=2, position=position_nudge(x=0)) +
  scale_color_manual(values=reg_colors, guide=F) +
  scale_fill_manual(values=reg_colors, guide=F) +
  scale_y_continuous('Dirty Fuel Use', labels = scales::percent) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        text = element_text(size=16),
        axis.text.x = element_text(size=7.5, angle = 60, hjust = 1))

plot <- plot + annotation_custom(grob=reg_grob, ymin=0, ymax=.25, xmin=0, xmax=20)

ggsave(filename=file.path(save.dir, 'ex_fig_1_lolli_c.png'),
       width=12, height=6, units='in', dpi=900)
#***********************************************************************************************************************

# ---EXT FIG 2----------------------------------------------------------------------------------------------------------
#Figure 5: Inequality over time
dt_ineq <- dt[cause=='lri' & grouping=='under5' & type=='HAP', 
              .(iso3, year, ADM0_CODE, ADM0_NAME, ADM2_CODE, ADM2_NAME, prev_mean, 
                super_region_id, super_region_name, region_id, region_name, pop_total)]
#AID results
#note that AID = gini * 2 * mean
aid.dt <-
  na.omit(dt_ineq, cols='prev_mean') %>% 
  .[, mean := weighted.mean(prev_mean, weights = pop_total, na.rm=T), by=.(iso3, year)] %>% 
  .[, .(mean,
        aid=gini(prev_mean, weights=pop_total) * 2 * mean,
        pop=sum(pop_total, na.rm=T)), 
    by=.(super_region_id, region_name, iso3, ADM0_NAME, year)] %>% 
  unique(by=c('iso3', 'year')) %>% 
  .[order(aid),] %>% 
  .[year==min(year), aid_start := aid] %>% 
  .[, aid_start := min(aid_start, na.rm=T), by=.(iso3)] %>%
  .[, aid_d := (aid-aid_start)] %>% 
  .[, aid_dr := aid_d/aid_start] %>% 
  .[year==min(year), sfu_start := mean] %>% 
  .[, sfu_start := min(sfu_start, na.rm=T), by=.(iso3)] %>%
  .[, dfu_d := (mean-sfu_start)] %>% 
  .[, dfu_dr := dfu_d/sfu_start] %>% 
  .[, pred := predict(loess(aid~(1-mean)))] %>% 
  .[, resid := abs(aid-pred)] %>% 
  .[year==max(year)&resid>.075, label := iso3]

plot <-
  ggplot(aid.dt, aes(x=1-mean, y=aid, color=super_region_id %>% as.factor, label=label)) +
  geom_smooth(method='gam', aes(x=1-mean, y=aid, color=super_region_id %>% as.factor), se=T, size=1.5) +
  geom_point(data=aid.dt[year==max(year)], aes(x=1-mean, y=aid, color=super_region_id %>% as.factor)) +
  geom_text_repel(force=3, color='black') +
  scale_color_manual(values=reg_colors, guide=F) +
  scale_y_continuous('Average interpersonal difference') +
  scale_x_continuous('Clean fuel use (%)', limits=c(0, 1)) +
  theme_minimal()
ggsave(filename=file.path(out.dir, 'aid_vs_cfu_sr.png'), plot=plot, 
       width=12, height=8, units='in', dpi=500)

plot <-
  ggplot(aid.dt, aes(x=1-mean, y=aid, label=label)) +
  geom_text_repel(force=2, color='black') +
  geom_smooth(method='gam', aes(x=1-mean, y=aid), color='navy', se=F, size=1, linetype='solid') +
  geom_point(data=aid.dt[year==max(year)], aes(x=1-mean, y=aid, color=super_region_id %>% as.factor), size=2) +
  scale_color_manual(values=reg_colors, guide=F) +
  scale_y_continuous('Average interdistrict difference') +
  scale_x_continuous('Clean fuel use (%)', limits=c(0, 1)) +
  theme_minimal() +
  theme(text = element_text(size=16))

plot <- plot + annotation_custom(grob=reg_grob, ymin=.5, ymax=.7, xmin=.75, xmax=1)

ggsave(filename=file.path(save.dir, 'ex_fig_2.png'), plot=plot, 
       width=12, height=8, units='in', dpi=900)
#***********************************************************************************************************************

# ---LRI CHANGE VS HAP CHANGE-------------------------------------------------------------------------------------------
#scatter the change in LRI rates vs the change in HAP prevalence
plot.dt <- dt_d[term=='change_rate' & type=='HAP']

plot.dt <-  dt %>% 
  copy %>% 
  .[year %in% c(start_year) & type=='HAP', .(ADM2_CODE, rate_start=rate_mean*1e3)] %>% 
  merge(., plot.dt, by=c('ADM2_CODE'))

#custom fx to help
makeChangeScatters <- function(reg, rug_var='rate_start', rug_limits=c(0, 5)) {
  
  
  message('building limits/scale')
  ## Enforce limits & define plot_var for simplicity
  plot.dt$plot_var <- pmax(rug_limits[1], pmin(rug_limits[2], as.data.frame(plot.dt)[, rug_var])) 
  start_range <- range(plot.dt$plot_var, na.rm = T)
  
  legend_colors <- inferno(n=6)
  legend_color_values <- c(rug_limits[1], .25, .5, 1, 3, rug_limits[2])
  
  ## Create breaks
  breaks <- pretty(rug_limits, 5)
  if (rug_limits[1] < 0 & rug_limits[2] > 0) breaks <- sort(unique(c(0, breaks)))
  
  ## Create labels
  labels <- format(breaks, nsmall = 2)
  if (min(rug_limits) >= 0) divider <- "-" else divider <- " to "
  if (start_range[1] < rug_limits[1]) {
    labels[1] <- paste0(format(floor(100*start_range[1])/100, nsmall=2), divider, labels[1])
  }
  if (start_range[2] > rug_limits[2]) {
    labels[length(labels)] <- paste0(labels[length(labels)], divider, format(ceiling(100*start_range[2])/100, nsmall=2))
  }

  plot <-
  ggplot(plot.dt[region_id==reg], aes(x=rate_mean, y=prev_mean, color=plot_var)) +
    geom_point() +
    geom_hline(yintercept=0, linetype='dashed') +
    geom_vline(xintercept=0, linetype='dashed') +
    geom_abline(slope=1) +
    facet_wrap(~ADM0_NAME) +
    scale_x_continuous('Rate of change for LRI Rates', limits=c(-1,1)) +
    scale_y_continuous('Rate of change for SFU%', limits=c(-1,1)) +
    scale_color_gradientn(colors = legend_colors, values = legend_color_values,
                          limits = range(breaks), breaks = breaks, labels = labels, 
                          name = 'LRI deaths / 1,000 children') +
    ggtitle(plot.dt[region_id==reg, unique(region_name)])+
    theme_bw()
  
  print(plot)
  
}

pdf(paste0(out.dir, '/hap_lri_change_scatters.pdf'),
    height=8, width=12)

lapply(unique(dt_d$region_id), makeChangeScatters)

dev.off()
#***********************************************************************************************************************
