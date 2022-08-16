require(INLA)
require(data.table)
require(ggplot2)

plot_hyperparameters <- function(indicator, indicator_group, run_date, age, holdout, save_file = NULL, regs,
                                 debug=F) {
  
  if(debug) browser()
  
  # get regions
  outputdir <- paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date, '/')
  regions <- regs
  
  # extract prior and posterior distributions from INLA model objects
  message('Load models & extract priors and posteriors for hyper-parameters')
  dist <- rbindlist(lapply(regions, function(r) {

    # load model
    message(paste0('...', r))
    result = tryCatch({
      load(paste0(outputdir, indicator, '_model_eb_bin', age, '_', r, '_', holdout, '.RData'))
    }, error = function(e) {
      return()
    })
    
    # extract hyper-priors from INLA (based on plot.inla() code)
    all.hyper <- INLA:::inla.all.hyper.postprocess(res_fit$all.hyper)
    hyper <- res_fit$marginals.hyperpar
    id <- strsplit(sapply(hyper, attr, 'hyperid'), split = '\\|')
    prior <- rbindlist(lapply(names(id), function(x) {
      print(x)
      if (grepl('Theta. for', x)) range <- c(-5, 5)
      if (grepl('GroupRho for', x)) range <- c(-0.999, 0.999)
      if (grepl('Group PACF. for', x)) range <- c(-0.999, 0.999)
      if (grepl('Precision for', x)) range <- c(1, 1000)
      p <- INLA:::inla.get.prior.xy(section = tolower(id[[x]][2]), hyperid = id[[x]][1], all.hyper = all.hyper, range = range, intern = F)
      if (grepl('Precision for', x)) {
        p <- inla.tmarginal(function(x) sqrt(1/x), p, method = 'linear')
        x <- gsub('Precision for', 'SD for', x)
      }
      data.table(region = r, type = 'prior', name = x, x = p$x, y = p$y)
    }))
    
    # extract corresponding posteriors from INLA
    post <- rbindlist(lapply(names(hyper), function(x) {
      p <- hyper[[x]]
      if (grepl('Precision for', x)) {
        p <- inla.tmarginal(function(x) sqrt(1/x), p, method = 'linear')
        x <- gsub('Precision for', 'SD for', x)
      }
      data.table(region = r, type = 'posterior', name = x, x = p[, 'x'], y = p[, 'y'])
    }))
    
    # combine
    all <- rbind(prior, post)
    all[, name := factor(name, unique(name))]
    return(all)
  }))
  
  # make plots
  message('Plotting hyper-parameters')

  #custom function to control the hyper param plots
  distPlot <- function(dist_dt, reg) {
    
    #subset data
    if (reg=='regional') { dt <- dist_dt %>% copy %>% .[nchar(region)>3]; reg_label <- ' (Regional Models)'
    } else if (reg=='country') { dt <- dist_dt %>% copy %>% .[nchar(region)==3]; reg_label <- ' (Ctry Models)'
    } else dt <- dist_dt %>% copy %>% .[region==reg]; reg_label <- paste0(' (', reg, ')')

    gg <- ggplot(dt[y > 1e-8,], aes(x = x, y = y, color = region, linetype = type)) +
      facet_wrap(~ name, scales = 'free') +
      geom_line() +
      scale_color_brewer(palette = 'Set3') +
      labs(x = '', y = '', title = paste0('Hyper-parameter prior and posterior distributions for ', indicator, reg_label)) +
      theme_minimal()
    
    print(gg)
    
    return(NULL)
    
  }

  
  #make series of plots
  
  if (is.null(save_file)) save_file <- paste0(outputdir, '/diagnostic_plots/inla_hyperparameters.pdf')
  pdf(save_file, width = 14, height = 8)
  
  if (length(regions)>1) {
    
    distPlot(dist, reg='regional')
    distPlot(dist, reg='country')
    
  } 
  
  lapply(unique(dist$region), distPlot, dist_dt=dist)

  dev.off()
  
  return(dist)
}