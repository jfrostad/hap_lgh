# ---------------------------------------------------------------------------------------------
# plot_stackers_by_adm01
#
# Written by Dan Casey
# Modified by Kirsten Wiens
# Modified by JF
#
# 08/28/2018
# Function that pulls in aggregated results, data, and stackers
# and plots them for admin 0 and admin 1 over time
#
# Inputs:
# indicator - modeled indicator
# indicator_group - modeled indicator group
# run_date - run date corresponding to model run
# reg - region modeling (NOTE: the function only takes one at a time)
# measure - the name you'd like for your measure, e.g. 'prevalence' or 'proportion'
# draws - whether or not to use draws
# raked - whether or not to plot raked estimates
# credible_interval - credible interval to plot for the mean unraked estimates
# N_breaks - breaks you'd like to use to plot points in proportion to N
# admin_data - optional
#              if NULL:
#                will aggregate data pulled from /share/geospatial/mbg/input_data/
#              if proivded, must be:
#                list of aggregated admin 0 and admin 1 data with data.table/data.frame 
#                elements named ad0 and ad1, respectively, and the following collumns:
#                ad0: svy_id, ADM0_NAME, ADM0_CODE, year, outcome, N
#                ad1: svy_id, ADM0_NAME, ADM0_CODE, ADM1_CODE, ADM1_CODE, year, outcome, N
#                which can be created using mbg_central function input_aggregate_admin()
#
# Outputs:
# plots of data, stackers, mean unraked estimates with upper and lower credible intervals,
# and mean raked estimates (optionally) over time saved in /share/ mbg output folder
# corresponding to indicator and run date in file /diagnositc_plots/
# ---------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# Start function
stacker_time_series_plots <- function(reg,
                                      dt,
                                      indicator = indicator,
                                      indicator_group = indicator_group,
                                      run_date = run_date,
                                      raked = T,
                                      vetting_colorscale = NULL,
                                      label = NULL,
                                      debug=F) {
  # -----------------------------------------------------------------------------------
  # Work interactively to build function
  if (debug) browser()

  message(paste0('Aggregating input data for: ', reg))
  # Load packages
  pacman::p_load(data.table, ggplot2, ggthemes, ggrepel, magrittr, raster, RColorBrewer, rgeos, rgdal, sp, sf)
  # --------------------------

  # ------------------------------------------------------------------------------------------------------------------------
  # Load and clean data and model information
  
  # set model output directories
  message('Loading and cleaning model information')
  outputdir <- paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date, '/')
  imagedir <- paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/model_image_history/')

  # subset input data
  dt <- dt[region == reg]
  # load submodel names
  config <- fetch_from_rdata(paste0(imagedir, 'pre_run_tempimage_', run_date, '_bin0_', reg, '_0.RData'), 'config')
  
  #load geographic information to help ID out of country points (artifact of frax aggregation)
  lookup <- fread('/home/j/WORK/11_geospatial/10_mbg/stage_master_list.csv')
  
  # check to see if new stacking was used
  contains_object <- function(file, obj){
    load(file)
    return(exists(obj))
  }
  
  if (contains_object(paste0(imagedir, 'pre_run_tempimage_', run_date, '_bin0_', reg, '_0.RData'), 'stacking_models')){
    
    stack <- fetch_from_rdata(paste0(imagedir, 'pre_run_tempimage_', run_date, '_bin0_', reg, '_0.RData'), 'stacking_models')
    submodels <- unlist(lapply(stack, function(x) x$model_name))
    
  } else submodels <- trimws(strsplit(config[V1 == 'stacked_fixed_effects',V2], '+', fixed = T)[[1]])
  
  #check relevant xgboost settings
  if (config[V1 == 'xg_second_best',V2]) submodels <- c(submodels, 'xgboost2')

  # load coefficients
  mod <- lapply(paste0(outputdir, indicator, '_model_eb_bin0_',reg,'_0.RData'), 
               function(x) fetch_from_rdata(x, 'res_fit'))
  if (as.logical(use_inla_country_fes)) {
    # fixed effects
    mod <- mod[[1]]
    coefs <- data.table(child = submodels, coef = invlogit(mod$summary.fixed$mean[2:4]), region = reg)
  } else {
    # random effects
    mod <- lapply(mod, function(x) data.table(child = submodels, coef = x$summary.random$covar$mean))
    coefs <- rbindlist(lapply(seq(reg), function(x) mod[[x]][,region:=reg[x]]))
  }
  
  #linear combination of stacker columns by coefficient
  dt$stack <- 0
  for (stacker in submodels){
    dt$stack<- dt[,..stacker]*coefs[child == stacker,coef] + dt$stack
  }
  submodels <- c(submodels, 'stack')

  # get unraked results and reshape long
  mbg <- melt(dt, variable.factor = F,
              measure.vars = c(submodels, 'mean', 
                               ifelse(raked, 'mean_raked', NA) %>% #account for possibility of raked values
                                 na.omit) 
              )

  #merge on coefficients
  #also rename stacker results, to make the colors match up later (we want mean to be alphabetically first)
  mbg <- merge(mbg, coefs, by.x = c('variable', 'region'), by.y = c('child', 'region'), all.x = T)
  mbg[!is.na(coef), variable := paste('stk:', variable, round(coef, 2))]
  
  #ID out of country points using lookup info
  mbg <- merge(mbg, lookup[, .(iso3, ADM0_CODE=gadm_geoid)], by='ADM0_CODE')
  mbg[, border_data := iso3 != svy_iso3]
  mbg[, iso3_label := '']
  mbg[border_data==T, iso3_label := paste0(svy_iso3, ': ')]
  
  #add NID label to sample sizes
  mbg[!is.na(input_ss), N_label := paste0(iso3_label, nid, ' [#', round(input_ss), ']')]
  # ------------------------------------------------------------------------------------------------------------------------
  
  
  # ------------------------------------------------------------------------------------------------------------------------
  # Make plots
  message('Making plots')
  
  #add label to plot filepath
  if (label %>% is.null) plot_label <- ''
  else if (label=='config') plot_label <- list.files(outputdir, pattern = '.note') %>% gsub('.note|000', '', .)
  else plot_label <- label

  # create pdf to write to
  dir.create(paste0(outputdir, 'diagnostic_plots/'))
  pdf(paste0(outputdir, 'diagnostic_plots/admin_stacker_line_plots_', reg, plot_label, '.pdf'), height = 10, width = 14)

  # test data structure
  if (nrow(dt[lvl=='adm0' & !is.na(input_mean)]) != dt[lvl=='adm0' & !is.na(input_mean), .(nid,ADM0_CODE, year)] %>% uniqueN) {
    warning('The number of years of admin 0 data does not match number of years of admin 1 data.')
  }
  
  # build custom function to create timeseries plot suite
  timePlot <- function(admin_id, dt, this_lvl,
                       year_start=2000, year_end=2017) {

    #subset data
    this_codevar <- this_lvl %>% toupper %>% paste0(., '_CODE')
    this_namevar <- this_lvl %>% toupper %>% paste0(., '_NAME')
    plot_dt <- dt[lvl==this_lvl & get(this_codevar) == admin_id,]
    
    adname <- plot_dt[, get(this_namevar)] %>% unique

    # time series plot
    if (na.omit(plot_dt, cols='value') %>% nrow == 0) message(paste0('Skipping ', adname, ' , no estimates!'))
    else {
      
      message('Plotting ', adname)
      
      # set ylims
      maxval <- (plot_dt[,.(upper, value, input_mean)] %>% max) + .05 #add 5% buffer

      #create plot
      plot <- ggplot(data = plot_dt[variable == 'mean']) + 
        geom_point(aes(x = year, y = input_mean, size = input_ss)) +
        geom_text_repel(aes(x = year, y = input_mean, label = N_label),
                        size=5, segment.colour = 'grey17', segment.size = .2, direction = 'x', nudge_y=.05, force=5) +
        geom_line(data = plot_dt, aes(x = year, y = value, color = variable), size = 1.2) +
        geom_ribbon(data = , aes(x = year, ymin = lower, ymax = upper), na.rm = TRUE, alpha = 0.2, fill = 'indianred3') +
        xlim(year_start, year_end) +
        ylim(0, maxval) +
        scale_color_brewer(palette = 'Set1') +
        scale_size_continuous('N', range=c(2,7)) +
        labs(title=paste("Model Results for", adname, ' (', reg, ')'),
             subtitle = paste0('Indicator=', indicator),
             x='Year',
             y='Proportion') + 
        theme_minimal()
      
      #if raked, add uncertainty ribbon
      if (raked) { 
        
        plot <- plot + 
          geom_ribbon(data = , aes(x = year, ymin = lower_raked, ymax = upper_raked), 
                      na.rm = TRUE, alpha = 0.2, fill = 'dodgerblue1')
      
      }  
      
      #add vetting colors to points if scheme is provided
      if (!is.null(vetting_colorscale)) {
        
        plot <- plot + 
          geom_point(aes(x = year, y = input_mean, size = input_ss, fill=vetting), shape=21, stroke=0) +
          scale_fill_manual('Vetting Status', values = vetting_colorscale)
        
      }

      #print to dev
      print(plot)
      
    }
      

    # make an admin 1 plot by child model
    if (this_lvl=='adm0') {
      
      #subset data
      plot_dt <- mbg[lvl=='adm1' & get(this_codevar) == admin_id,]
      
      #reset ylims
      maxval <- (plot_dt[,.(upper, value, input_mean)] %>% max) + .05 #add 5% buffer

      # build a palette with sufficient colors by sampling from colorbrewer
      admin_N <- uniqueN(plot_dt$ADM1_NAME)
      colors <- brewer.pal.info[brewer.pal.info$category == 'qual' & brewer.pal.info$colorblind==T,]
      colors <- mapply(brewer.pal, colors$maxcolors, rownames(colors)) %>% unlist
      if (admin_N>length(colors)) colors <- sample(colors, admin_N, replace=T)
      else sample(colors, admin_N)

      childplot <- ggplot(plot_dt, aes(x = year, y = value, color = ADM1_NAME)) + 
        geom_point(data = plot_dt, mapping = aes(x = year, y = input_mean, size = input_ss)) +
        #label border datapoints with their svy iso3
        geom_text_repel(data = plot_dt[border_data==T],
                        aes(x = year, y = input_mean, label = svy_iso3),
                        size=5, segment.colour = 'grey17', segment.size = .2, direction = 'x', nudge_y=.05, force=5) +
        geom_line() + 
        facet_wrap(~variable) + 
        xlim(year_start, year_end) +
        ylim(0, maxval) +
        scale_color_manual(values=colors) +
        scale_size_continuous('N', range=c(2,7)) +
        labs(title=paste("Model Results for", adname, 'at admin 1 level'),             
             x='Year',
             y='Proportion') +
        theme_minimal() +
        theme(legend.position="bottom")

        print(childplot)
        
    }
    
    return(NULL) #no need to return anything

  }
  
  # build custom function to loop over each ADM0
  countryLoop <- function(country, dt) {

    #first, make country level dataplots to get a feel for outliers
    plot <- ggplot(dt[variable=='mean' & lvl=='adm0']) + 
      geom_point(mapping=aes(x=year, y=input_mean, color=ADM0_NAME, size=input_ss, shape=border_data)) + 
      geom_line(mapping=aes(x=year, y=value, color=ADM0_NAME)) +
      scale_size_continuous('N', range=c(2,7)) +
      labs(title=paste("Region-level MBG Results/Data"),
           subtitle=reg,
           x='Year',
           y='Proportion') +
      theme_minimal()
    print(plot)
    
    #first, create the country level lineplots
    timePlot(admin_id=country, this_lvl='adm0', dt=dt)
    
    #finally, create the ad1 plots
    # for each ad1 in the country
    dt[ADM0_CODE==country, ADM1_CODE] %>% 
      unique %>% 
      lapply(., timePlot, this_lvl='adm1', dt=dt)
    
    return(NULL) #no need to return anything
    
  }

  # for each country in the region
  mbg[, ADM0_CODE] %>% 
    unique %>% 
    intersect(., get_adm0_codes(reg)) %>%
    lapply(., countryLoop, dt=mbg)

  dev.off()
  message(paste0('Plots saved at ', outputdir, 'diagnostic_plots/admin_stacker_line_plots_', reg, '.pdf'))

  # ------------------------------------------------------------------------------------------------------------------------
  # End function
  return('Stacker line plots complete!')
  
}
# ------------------------------------------------------------------------------------------------------------------------
