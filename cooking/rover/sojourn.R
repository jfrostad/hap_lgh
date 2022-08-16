##############################################################################
## MBG diagnostics functions and plots 
#source("/homes/jfrostad/_code/lbd/hap/cooking/rover/sojourn.R", echo=T)
##############################################################################


## Setup -------------------------------------------------------------------------

## clear environment
rm(list=ls())

# runtime configuration
if (Sys.info()["sysname"] == "Linux") {
  j_root <- "/home/j/"
  h_root <- file.path("/ihme/homes", Sys.info()["user"])
  
  package_lib    <- file.path(h_root, '_code/_lib/pkg')
  ## Load libraries and  MBG project functions.
  .libPaths(package_lib)
  
  # necessary to set this option in order to read in a non-english character shapefile on a linux system (cluster)
  Sys.setlocale(category = "LC_ALL", locale = "C")
  
} else {
  j_root <- "J:"
  h_root <- "H:"
}

#load external packages
pacman::p_load(data.table, fasterize, fst, gargle, googledrive, magrittr, mgsub, sf, farver, mltools)

#detect if running interactively
interactive <- F  %>% #manual override
  ifelse(., T, !length(commandArgs())>2) %>%  #check length of arguments being passed in
  ifelse(., T, !(is.na(Sys.getenv("RSTUDIO", unset = NA)))) #check if IDE

if (interactive) {

  
} else {
  
  ## Set repo location, indicator group, and some arguments
  user            <- commandArgs()[4]
  core_repo       <- commandArgs()[5]
  indicator_group <- commandArgs()[6]
  indicator       <- commandArgs()[7]
  config_par      <- commandArgs()[8]
  cov_par         <- commandArgs()[10]

}

#if new vetting activity has occured, need to refresh the local sheet
new_vetting <- F

#is this a raked model?
raked <- T

message(indicator)

## Load MBG packages
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

## Load custom post-estimation functions
file.path(my_repo, '_lib', 'post', 'aggregate_inputs.R') %>% source
file.path(my_repo, '_lib', 'post', 'plot_pv_results.R') %>% source
file.path(my_repo, '_lib', 'viz', 'stacker_admin_time_series.R') %>% source
file.path(my_repo, '_lib', 'viz', 'plot_hyperparameters.R') %>% source


#use your own diacritics fx, due to inscrutable error
#note: requires mgsub pkg
fix_diacritics <<- function(x) {
  
  #first define replacement patterns as a named list
  defs <-
    list('Š'='S', 'š'='s', 'Ž'='Z', 'ž'='z', 'À'='A', 'Á'='A', 'Â'='A', 'Ã'='A', 'Ä'='A', 'Å'='A', 'Æ'='A', 
         'Ç'='C', 'È'='E', 'É'='E','Ê'='E', 'Ë'='E', 'Ì'='I', 'Í'='I', 'Î'='I', 'Ï'='I', 'Ñ'='N', 'Ò'='O', 
         'Ó'='O', 'Ô'='O', 'Õ'='O', 'Ö'='O', 'Ø'='O', 'Ù'='U','Ú'='U', 'Û'='U', 'Ü'='U', 'Ý'='Y', 'Þ'='B', 
         'à'='a', 'á'='a', 'â'='a', 'ã'='a', 'ä'='a', 'å'='a', 'æ'='a', 'ç'='c','è'='e', 'é'='e', 'ê'='e', 
         'ë'='e', 'ì'='i', 'í'='i', 'î'='i', 'ï'='i', 'ð'='o', 'ñ'='n', 'ò'='o', 'ó'='o', 'ô'='o', 'õ'='o',
         'ö'='o', 'ø'='o', 'ù'='u', 'ú'='u', 'û'='u', 'ý'='y', 'ý'='y', 'þ'='b', 'ÿ'='y', 'ß'='Ss')
  
  #then force conversion to UTF-8 and replace with non-diacritic character
  enc2utf8(x) %>% 
    mgsub(., pattern=enc2utf8(names(defs)), replacement = defs) %>% 
    return
  
}

## Throw a check for things that are going to be needed later
message('Looking for things in the config that will be needed for this script to run properly')

config <- set_up_config(repo            = my_repo,
                        indicator_group = indicator_group,
                        indicator       = indicator,
                        config_name     = paste0('/model/configs/config_', config_par),
                        covs_name       = paste0('/model/configs/covs_', cov_par),
                        run_tests       = F,
                        post_est_only   = T,
)


## Create proper year list object
if (class(year_list) == 'character') year_list <- eval(parse(text=year_list))

## Get regions that have successfully completed through aggregation step
Regions <- list.files(outputdir, pattern = '*_admin_draws')
Regions <- gsub('.*eb_bin0_', '', Regions)
for (r in 1:length(Regions)) Regions[[r]] <- substr(Regions[[r]], start = 1, stop = nchar(Regions)[[r]]-8)
Regions <- unique(Regions)
Regions <- Regions[Regions != '']
message(paste0(Regions, '\n'))

## Combine and summarize aggregated results --------------------------

# combine raked results
message('Combining raked/unraked aggregated results')
if (indicator == 'cooking_fuel_solid') {
  combine_aggregation(rd       = run_date,
                      indic    = indicator,
                      ig       = indicator_group,
                      ages     = 0,
                      regions  = Regions,
                      holdouts = holdouts,
                      raked    = T:F,
                      measure = measure,
                      delete_region_files = F)
  
  # summarize admins
  summarize_admins(ad_levels = c(0,1,2), raked = T:F, measure = measure, metrics = 'rates')
} else {
  
  # combine unraked results
  message('Combining unraked aggregated results')
  combine_aggregation(rd       = run_date,
                      indic    = indicator,
                      ig       = indicator_group,
                      ages     = 0,
                      regions  = Regions,
                      holdouts = holdouts,
                      raked    = F,
                      measure = measure,
                      delete_region_files = F)
  
  # summarize admins
  summarize_admins(ad_levels = c(0,1,2), raked = F, measure = measure, metrics = 'rates')
  
}

if (length(holdouts)==1) {

## Aggregate data and stackers ------------------------------------------------------
# Aggregate data to admin 0 and 1
dat <- mclapply(Regions,function(x) 
  aggregate_input_data(reg=x,
                       indicator, 
                       indicator_group, 
                       run_date,
                       modeling_shapefile_version, 
                       build=F),
                       mc.cores=5
  ) %>% 
  rbindlist

# Aggregate stackers to admin 0 and 1
stack <- mclapply(Regions, function(x) 
  aggregate_child_stackers(reg=x,
                           indicator, 
                           indicator_group, 
                           run_date, 
                           modeling_shapefile_version,
                           pop_measure=pop_measure,
                           build=F),
                           mc.cores=5
  ) %>% 
  rbindlist

## Combine data and stackers with summary results ------------------------------------------

# Load and combine estimates
mbg <- list(
  paste0(outputdir, '/pred_derivatives/admin_summaries/', indicator, '_admin_0_unraked_summary.csv') %>% 
    fread %>%
    .[, lvl := 'adm0'],
  paste0(outputdir, '/pred_derivatives/admin_summaries/', indicator, '_admin_1_unraked_summary.csv') %>% 
    fread %>%
    .[, lvl := 'adm1']
)  %>% 
  rbindlist(use.names=T, fill=T)

# raked results
if (raked) {
  mbg_raked <- list(
    paste0(outputdir, '/pred_derivatives/admin_summaries/', indicator, paste0('_admin_0_raked_', measure, '_summary.csv')) %>% 
      fread %>%
      .[, lvl := 'adm0'],
    paste0(outputdir, '/pred_derivatives/admin_summaries/', indicator, paste0('_admin_1_raked_', measure, '_summary.csv')) %>%
      fread %>%
      .[, lvl := 'adm1']
  )  %>% 
    rbindlist(use.names=T, fill=T)
}

# Combine all
# raked results
if (raked) {
  #modify colnames
  c('mean', 'upper', 'lower', 'cirange') %>% 
    setnames(mbg_raked, ., paste0(., '_raked'))
  
  mbg <- merge(mbg, mbg_raked, 
               by = names(mbg) %>% .[grep('ADM|region|year|pop|lvl', .)],
               all.x = T)
}

# stackers
mbg <- merge(mbg, stack,
             by =  names(mbg) %>% .[grep('CODE|year|lvl', .)],
             all.x = T)

# data
mbg <- merge(mbg, dat,
             by = names(mbg) %>% .[grep('CODE|year|lvl', .)],
             all.x = T)

# save
write.csv(mbg, paste0(outputdir, '/pred_derivatives/admin_summaries/', indicator, '_mbg_data_stackers.csv' ))

##classify datapoints based on HAP vetting
#update vetting sheet if necessary

#read in vetting sheet
vetting <- file.path(doc.dir, 'HAP Tracking Sheet.xlsx') %>% readxl::read_xlsx(sheet='1. Vetting', skip=1) %>% 
  as.data.table %>% 
  .[, .(nid, vetting=`HAP Vetting Status`, svy_iso3=ihme_loc_id)] %>%  #subset to relevant columns
  unique(., by=names(.)) 
#Fix name for bobby =)
vetting[vetting=='Not started', vetting := 'Adequate']

#merge onto data
mbg <- merge(mbg, vetting, by='nid', all.x=T)

#define colorscale
# build color scheme for the vetting sheet values
vetting_colors <- c("Adequate"='grey4', 
                    "Problematic"='darkorange1',
                    "Completed"='forestgreen',
                    "Flagged"='purple1',
                    "Excluded"='gray71',
                    "In progress"='indianred2',
                    "Ready for Review"='indianred2')

## Plot stackers and covariates ------------------------------------------------------
message('Making time series plots for stackers by admin unit')
dir.create(paste0(outputdir, '/diagnostic_plots/'))

if (use_stacking_covs) {

  # plot stackers over time aggregated to admins
  message('Making time series plots for stackers by admin unit')
  mclapply(Regions %>% rev, function(x) 
    stacker_time_series_plots(reg=x,
                              dt=mbg,
                              indicator, 
                              indicator_group, 
                              run_date, 
                              raked=raked,
                              vetting_colorscale=vetting_colors,
                              label='config',
                              debug=F),
                              mc.cores=5
  )
  
}



# Plot priors for spatial hyperparameters --------------------------------------------------

message('Plotting spatial hyperparameters prior and posteriors')

# Plot a la HIV team
plot_hyperparameters(indicator = indicator,
                     indicator_group = indicator_group,
                     run_date = run_date,
                     age = 0,
                     holdout = holdouts,
                     save_file = NULL,
                     regs = Regions,
                     debug=F)

} else {

    # Combine csv files
    csvs <- list.files(outputdir, pattern = 'input_data_(.*).csv', full.names = T)
    csv_master <- rbindlist(lapply(csvs, fread))
    csv_master[, V1 := NULL]
    write.csv(csv_master, file=paste0(outputdir, '/input_data.csv'))

    # Get in and out of sample draws
    run_in_oos <- get_is_oos_draws(ind_gp        = indicator_group,
                                   ind           = indicator,
                                   rd            = run_date,
                                   ind_fm        = 'binomial',
                                   rk            = F,
                                   age           = 0,
                                   yrs           = c(2000, 2005, 2010, 2015, 2018),
                                   get.oos       = as.logical(makeholdouts),
                                   year_col      = 'year',
                                   write.to.file = TRUE,
                                   shapefile_version = modeling_shapefile_version,
                                   reg_list = Regions)
    
    # Set out_dir
    dir.create(paste0(outputdir, '/summary_metrics_unraked/'), recursive = T, showWarnings = F)
    
    # Calculate and save PV summary statistics
    pvtab <- rbindlist(lapply(list(c('oos'), c('oos', 'year'), c('oos', 'region')), function(results_by) {
      pv <- get_pv_table(d               = run_in_oos,
                         yrs             = c(2000, 2005, 2010, 2015, 2018),
                         indicator_group = indicator_group,
                         rd              = run_date,
                         indicator       = indicator,
                         aggregate_on    = c('country', 'ad1', 'ad2'),
                         result_agg_over = results_by,
                         coverage_probs  = c(95),
                         plot_ci         = TRUE,
                         plot_ncol       = 2,
                         draws           = as.numeric(samples),
                         save_csv        = FALSE,
                         out.dir         = paste0(outputdir, '/summary_metrics_unraked/'))
      ldply(pv, .id = 'Level')
    }), fill = T)
    write.csv(pvtab, file = paste0(outputdir, '/summary_metrics_unraked/pv_metrics.csv'), row.names=F)

    #also produce raked metrics
    # Get in and out of sample draws
    run_in_oos_rk <- get_is_oos_draws(ind_gp        = indicator_group,
                                   ind           = indicator,
                                   rd            = run_date,
                                   ind_fm        = 'binomial',
                                   rk            = T,
                                   age           = 0,
                                   yrs           = c(2000, 2005, 2010, 2015, 2018),
                                   get.oos       = as.logical(makeholdouts),
                                   year_col      = 'year',
                                   write.to.file = TRUE,
                                   shapefile_version = modeling_shapefile_version,
                                   reg_list = Regions)
    
    # Set out_dir
    dir.create(paste0(outputdir, '/summary_metrics_raked/'), recursive = T, showWarnings = F)
    
    # Calculate and save PV summary statistics
    #draws.df <- fread(paste0(outputdir, '/output_draws_data.csv'))
    pvtab <- rbindlist(lapply(list(c('oos'), c('oos', 'year'), c('oos', 'region')), function(results_by) {
      pv <- get_pv_table(d               = run_in_oos_rk,
                         indicator_group = indicator_group,
                         rd              = run_date,
                         indicator       = indicator,
                         aggregate_on    = c('country', 'ad1', 'ad2'),
                         result_agg_over = results_by,
                         coverage_probs  = c(95),
                         plot_ci         = TRUE,
                         plot_ncol       = 2,
                         draws           = as.numeric(samples),
                         save_csv        = FALSE,
                         out.dir         = paste0(outputdir, '/summary_metrics_raked/'))
      ldply(pv, .id = 'Level')
    }), fill = T)
    write.csv(pvtab, file = paste0(outputdir, '/summary_metrics_raked/pv_metrics.csv'), row.names=F)
}
  