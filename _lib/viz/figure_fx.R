#***********************************************************************************************************************

# ---FIGURE 4-----------------------------------------------------------------------------------------------------------
#helper function to make the master transition plot
makeMasterPlot <- function(dt=plot.dt,
                           custom_countries,
                           custom_cols=NA,
                           add_rug=F, rug_var='rate_mean', rug_limits=c(0, 5),
                           contour_bw=7,
                           debug=F){
  
  if (debug) browser()
  
  #subset to these countries
  dt <- dt[iso3 %in% custom_countries]
  
  #build a color scale using color brewer if not provided
  if(custom_cols %>% is.na) custom_cols <- brewer.pal(uniqueN(plot.dt$loc_fct), "Dark2")
  
  #generate plot
  plot <- 
    ggplot(dt, 
           aes(x=1-hap_pct_mean, y=paf_mean, color=loc_fct, shape=year %>% as.factor, group=ADM2_CODE)) + 
    annotate("rect", xmin = 0, xmax = .5, ymin = .0, ymax = .61, fill='steelblue', alpha = .075) +
    annotate("rect", xmin = .5, xmax = 1, ymin = .0, ymax = .61, fill='tomato', alpha = .075) +
    annotate("text", x = .32, y = .60, label = "majority household sources") +
    annotate("text", x = .66, y = .60, label = "majority ambient sources") +
    annotate(
      geom = "segment",
      x = .81, 
      xend = .85, 
      y = .60,
      yend = .60,
      colour = "black",
      arrow = arrow(length = unit(0.2, "cm"), type = "closed")
    ) +
    annotate(
      geom = "segment",
      x = .15, 
      xend = 0.11,
      y = .60, 
      yend = .60,
      colour = "black",
      arrow = arrow(length = unit(0.2, "cm"), type = "closed")
    ) +
    geom_line(alpha=.1) +
    geom_point(size=1) +
    geom_point(data=dt[year==2018], size=1) + #overdraw the later years
    geom_vline(xintercept=.5, linetype="dashed") +
    #geom_density_2d(aes(linetype=year %>% as.factor, group=paste0(year, loc_fct)), bins=contour_bw, size=.2) +
    scale_color_manual('Country', values=custom_cols, guide = guide_legend(ncol=5, keywidth = .75)) +
    scale_shape_manual('Year', values=c(1, 16), guide=F) +
    scale_linetype_manual('Year', values=c('dotted', 'solid'), guide=F) +
    scale_x_continuous('', limits=c(-.02, 1.02), labels = scales::percent, breaks=c(.25, .5, .75), expand = c(0.01, 0.01)) +
    scale_y_continuous('Percent of LRI attributable to TAP (PAF)', 
                       limits=c(0,.61), labels = scales::percent, breaks=c(0, .2, .4, .6), expand = c(0, 0)) +
    theme_bw(base_size = 16) +
    theme(
      legend.position = c(.03, .2),
      legend.direction="horizontal",
      legend.justification = c("left", "top"),
      legend.box.just = "left",
      legend.text = element_text(size=8),
      legend.title = element_text(size=8),
      legend.margin = margin(6, 6, 6, 6),
      plot.margin = margin(12, 0, 0, 6, "pt"),
      axis.text.y = element_text(angle = 90, hjust=.4),
      axis.title.x = element_blank()
    ) 
  
  #add a rug to show the LRI rates
  if (add_rug) {
    
    message('building limits/scale')
    ## Enforce limits & define plot_var for simplicity
    dt$plot_var <- pmax(rug_limits[1], pmin(rug_limits[2], as.data.frame(dt)[, rug_var])) 

    plot <- plot + 
      new_scale_color() +
      geom_rug(data=dt[year==min(year)], aes(color=plot_var), sides='l', alpha=.5, position='jitter', 
               length=unit(0.035, "npc")) + 
      geom_rug(data=dt[year==max(year)], aes(color=plot_var), sides='r', alpha=.5, position='jitter', 
               length=unit(0.035, "npc")) +
      scale_color_gradientn(colors = plasma(n=8), limits=c(rug_limits[1], rug_limits[2]), 
                           na.value = "#FDE725FF",
                           name = 'LRI mortality per 1,000',
                           guide = guide_colourbar(title.position="top", title.hjust = 0.5,
                                                   keyheight = unit(2, 'npc'))) 
    
  }
  
  return(plot)
  
}


#helper function to create smaller inset versions of the transition plot
makeInset <- function(i, loclist, 
                      colors=c('2000'='#737373', '2018'='#08519c'),
                      type='scatter', scale_labels=F,
                      contour_bw=20,
                      add_rug=T, rug_var='rate_mean', rug_limits=c(0, 4),
                      debug=F) {
  
  if(debug) browser()
  
  loc <- loclist[i]
  n <- length(loclist)
  message('plotting ', loc)
  
  dt <- plot.dt[iso3%in%loc]
  cap <- plot.dt[iso3%in%loc, unique(location_name)]
  if(cap %like% 'Congo') cap <- 'D.R. Congo'
  if(cap %like% 'Tanzania') cap <- 'Tanzania'

  #build either a scatterplot or a contoured cloudplot
  if (type=='scatter') {
    
    plot <- ggplot(dt, aes(x=1-hap_pct_mean, y=paf_mean, group=ADM2_CODE)) + 
      geom_line(alpha=.1) +
      geom_point(data=dt[year==2000], color=colors[[1]], size=.5, alpha=.5) +
      geom_point(data=dt[year==2018], color=colors[[2]], size=.5, alpha=.5)
      
    
  } else if (type %like% 'cloud') {
    
    #Egypt needs to be noised slightly in 2018 because hap is so universally low that it cannot be plotted as a contour
    dt[, noise := runif(.N, min=0, max=0.025)]
    dt[iso3=='EGY' & year==2018, hap_pct_mean := hap_pct_mean + noise]
    
    if (type=='cloud_contour') {
      
      #increase bin width for countries with less variance
      if(sd(dt$hap_pct_mean, na.rm=T)<.05) contour_bw <- contour_bw * 3.5
      
      plot <- ggplot(dt, aes(x=1-hap_pct_mean, y=paf_mean)) + 
        annotate("rect", xmin = 0, xmax = .5, ymin = .05, ymax = .61, fill='steelblue', alpha = .075) +
        annotate("rect", xmin = .5, xmax = 1, ymin = .05, ymax = .61, fill='tomato', alpha = .075) +
        geom_density_2d(data=dt[year==2000], color=colors[[1]], bins=contour_bw, size=.2, linetype='solid') +
        geom_density_2d(data=dt[year==2018], color=colors[[2]], bins=contour_bw, size=.2, linetype='solid')
      
    } else if (type=='cloud_alpha') {
      
      plot <- ggplot(dt, aes(x=1-hap_pct_mean, y=paf_mean)) + 
        stat_density_2d(data=dt[year==2000], aes(alpha = after_stat(level)), geom = "polygon", fill=colors[[1]], contour=TRUE, binwidth=20) +
        stat_density_2d(data=dt[year==2018], aes(alpha = after_stat(level)), geom = "polygon", fill=colors[[2]], contour=TRUE, binwidth=20) +
        scale_alpha(guide=F, range=c(.1, 0.01))
      
    }
    
  }
  
  if (add_rug) {

    message('building limits/scale')
    ## Enforce limits & define plot_var for simplicity
    dt$plot_var <- pmax(rug_limits[1], pmin(rug_limits[2], as.data.frame(dt)[, rug_var])) 
    start_range <- range(dt$plot_var, na.rm = T)

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
    
    plot <- plot + 
      geom_rug(data=dt[year==min(year)], aes(color=plot_var), sides='l', alpha=.75, position='jitter', size=1,
               length = unit(0.075, "npc")) + 
      geom_rug(data=dt[year==max(year)], aes(color=plot_var), sides='r', alpha=.75, position='jitter', size=1,
               length = unit(0.075, "npc")) + 
      scale_color_gradientn(colors = plasma(n=8), limits=c(rug_limits[1], rug_limits[2]), 
                           na.value = "#FDE725FF", 
                           guide=F) 
    
  }
  
  #standard settings
  plot <- plot +
    geom_vline(xintercept=.5, linetype="dashed") +
    scale_x_continuous("", limits=c(-.05, 1.05), labels = scales::percent, breaks=c(.25, .75), expand = c(0.02, 0.02)) +
    scale_y_continuous("", limits=c(.05,.61), labels = scales::percent, breaks=c(.3, .5), expand = c(0, 0),
                       position='right') +
    ggtitle(cap) + #if you want the title on the top
    theme_bw(base_size=16) +
    theme(plot.title=element_text(hjust=0.55, size=12, margin=margin(0,0,0,0)),
          axis.title=element_blank())
  
  #remove the labels for the top left insets or if there are no scale labels requested
  if (!scale_labels | (n==9 & i %in% c(1, 3, 5)) | (n==12 & i %in% c(1:2, 4:5, 7:8))) {
    
    plot <- plot + theme(axis.text=element_blank(), 
                         axis.ticks=element_blank(),
                         plot.margin=margin(2, 4, 12, 8, "pt"))
    
    #remove just the x labels for right insets
  } else if ((n==9 & i %in% c(2,4,6)) | (n==12 & i %in% c(3,6,9))) {
    
    plot <- plot + theme(axis.text.x=element_blank(), 
                         axis.ticks.x=element_blank(),
                         axis.text.y.right=element_text(angle = -90, vjust=0, hjust=.4),
                         plot.margin=margin(2, 0, 12, 0, "pt"))
    
    #remove just th y labels for the bottom left insets  
  } else if ((n-i) %in% (1:2)) {
    
    plot <- plot + theme(axis.text.y=element_blank(), 
                         axis.ticks.y=element_blank(),
                         plot.margin=margin(0, 4, 0, 8, "pt"))
    
    #keep all labels for the bottom right insets  
  } else if (i==n) {
    
    plot <- plot + theme(axis.text.y.right=element_text(angle = -90, vjust=0, hjust=.4),
                         plot.margin=margin(0, 0, 0, 0, "pt"))
    
  }
  
  return(plot)
  
}

makeRegFigure4s <- function(loc_list) {
  
  plot <-
  makeMasterPlot(custom_countries=loc_list,
                 add_rug=T, rug_limits=c(0,4))
  
  ggsave(plot=plot, filename=file.path(save.dir, 'si',  
                                       paste0('fig_4_',
                                              unique(plot.dt[iso3%in%loc_list, region_id]), 
                                              '_', paste0(loc_list, collapse="_"), '.png')
                                       ),
         width=12, height=8, units='in', dpi=600)
  
}

#***********************************************************************************************************************

# ---FIGURE 3-----------------------------------------------------------------------------------------------------------
#helper function to give us the intervals for binwidth after cutting
cut_borders <- function(x){
  pattern <- "(\\(|\\[)(-*[0-9]+\\.*[0-9]*),(-*[0-9]+\\.*[0-9]*)(\\)|\\])"
  
  chr <- as.character(x)
  
  start <- as.numeric(gsub(pattern,"\\2", chr))
  end <- as.numeric(gsub(pattern,"\\3", chr))
  
  list(start, end)
}

#helper function to create a weighted density plot
makeTapDensityPlot <- function(
  input_dt, wt_var='pop_total', #which variables to plot and dataset
  locs, #list of iso3 to plot separately
  start_year=2000, end_year=2018, #time variables
  tap_cutoff=750, #where to cut off the x axis
  log_seq=T, #use a logged sequence to generate bins (else regular bins)
  smoother=.1, #what span to use for the smoother
  nudge=1e-5, #how much to nudge values in order to work in logspace (prevent NaNs for 0)
  normalize_years=T, #normalize the areas by using the max sum to calculate density
  flipped=F, #flip to look like population pyramid
  density_max=NA, #pass in a limit for the density axis
  plot_quartiles=F, #if wanting to display quartiles as dotted lines
  reg_plot=F, #set to subset by super_region isntead of by all LMICs
  debug=F
) {
  
  if(debug) browser()
  
  #subset if regplot
  if (reg_plot) { 
    message('working on locs=', paste0(locs, collapse='+'))
    input_dt <- input_dt[iso3 %in% locs]
  }
  
  #prepare input data
  dt <- input_dt[type!='TAP'] %>% copy #we will recalculate TAP later
  dt[iso3 %in% locs, loc := ADM0_NAME]
  dt[is.na(loc), loc := 'All other LMICs']
  dt[loc=='Democratic Republic of the Congo', loc := 'D.R. Congo']

  #generate logarithmically spaced sequence
  sequence <- ifelse(rep(log_seq, 400),
                     c(exp(seq(log(.01), log(tap_cutoff), length.out = 399)), max(plot.dt$pm_pc_mean, na.rm=T)),
                     c(seq(0, 750, length.out = 399), max(plot.dt$pm_pc_mean, na.rm=T)))
  
  #create data.table from the sequence
  seq.dt <- data.table(tap_bin_id=cut(sequence, sequence, include.lowest = TRUE, right=F, labels=F))
  seq.dt[, c('start', 'end') := cut(sequence, sequence, include.lowest = TRUE, right=F, dig.lab=4) %>% cut_borders]
  seq.dt[1, `:=` (start=1e-5, end=sequence[1], tap_bin_id=0)]
  
  #expand the dt to include all relevant vars
  skeleton <- expand.grid(year=c(start_year, end_year),
                          type=c('AAP', 'HAP'),
                          loc=unique(dt$loc),
                          stringsAsFactors = F)
  seq.dt <- tidyr::crossing(seq.dt, skeleton) %>% as.data.table
  
  #merge TAP values onto the nearest bin
  dt[, tap_pm := sum(pm_pc_mean), by=.(year, ADM2_CODE)] #first, recalculate TAP
  dt[, tap_bin := cut(tap_pm, sequence, include.lowest = TRUE, right=F, dig.lab=4), by=.(loc, year)]
  dt[, tap_bin_id := cut(tap_pm, sequence, include.lowest = TRUE, right=F, labels=F), by=.(loc, year)]
  
  #generate sums
  dt[, sum_var := get(wt_var)]
  dt[, sum := sum(sum_var, na.rm=T), by=.(year)]
  # collapse to bin/types
  dt[, bin_sum := sum(sum_var, na.rm=T), by=.(loc, year, type, tap_bin_id)]
  
  #collapse
  bin.dt <- unique(dt, by=c('loc', 'type', 'year', 'tap_bin')) %>% .[(order(type, year, tap_bin_id))]
  
  #add on the missing bins
  bin.dt <- merge(bin.dt, seq.dt, by=c('loc', 'type', 'year', 'tap_bin_id'), all.y=T, allow.cartesian=T)
  bin.dt[, sum := zoo::na.aggregate(sum), by=.(year)]
  bin.dt[is.na(bin_sum), bin_sum := 0]

  #if we want the areas to be comparable across years (ie not equal), we can use the max sum for both
  if (normalize_years) bin.dt[, sum := max(sum, na.rm=T)]
  
  #calculate density
  bin.dt[, bin_width := end-start]
  bin.dt[, prob := bin_sum/sum]
  bin.dt[, density := prob/bin_width]
  bin.dt[, mid := start+((end-start)/2)]
  
  #make smoothed predictions
  bin.dt[, smooth := loess(log(density+nudge) ~ mid, span=smoother) %>% 
           predict(.SD) %>% 
           exp, 
         by=.(type, loc, year)]

  #calculate cumulative density across sample
  stats <- copy(bin.dt) %>% setkey(., mid, year)
  stats[, sum := sum(prob, na.rm=T), by=key(stats)]
  stats <- stats[, .(sum, mid, year)] %>% unique(., by=key(.))
  stats[, cum := cumsum(sum), by = year]
  
  #define stats to print
  quantiles <- expand.grid(year=c(start_year, end_year),
                            cum=c(.2, 0.25, .4, .5, .6, .75, .8)) %>% as.data.table
  
  quantiles <- stats[quantiles, on=.(year, cum), roll='nearest'] 
  print(quantiles)

  #melt years to make a new variable
  bin.dt <- bin.dt[, .(loc, type, year, mid, smooth)] #drop everything else that varies by year and cleanup dt
  bin.dt[, year := paste0('smooth_', year)]
  bin.dt <- dcast(bin.dt, ...~year, value.var='smooth')

  #make color palette
  bin.dt[, loc := factor(loc)]
  if(uniqueN(bin.dt$loc)>2) loc_colors <- brewer.pal(uniqueN(bin.dt$loc)-ifelse(reg_plot, 0, 1), "Dark2")
  else loc_colors <- c('#a6cee3', '#fdbf6f')
  if (!reg_plot) loc_colors <- c('#bdbdbd', loc_colors) #force other LMICs to grey

  #make plot
  plot <-
    ggplot(data=bin.dt[type!='TAP'], 
           aes(x=mid, y=smooth_2018, fill=loc, color=loc, alpha=type)) +
    #endyear plots
    geom_area() +
    #start year plots
    geom_area(aes(y=smooth_2000*-1)) + #invert for the mirrored comparison
    #thresholds
    geom_vline(xintercept=35, linetype='dashed', color='dodgerblue2') +
    geom_hline(yintercept=0, color='gray10', size=1) +
    #scales
    scale_x_sqrt('', limits=c(5, tap_cutoff),
                 breaks=c(15, 35, 50, 100, 300, 600)) +
    scale_y_continuous(expand = c(0.000075, 0.000075)) +
    scale_color_manual('', values=loc_colors) +
    scale_fill_manual('', values=loc_colors) +
    scale_alpha_manual(values=c('AAP'=.4, 'HAP'=1), guide=F) +
    theme_minimal() +
    #theme
    theme(axis.title=element_blank(),
          text = element_text(size=16),
          legend.position = c(.01, .99),
          legend.justification = c("left", "top"),
          legend.box.just = "left",
          legend.margin = margin(6, 6, 6, 6),
          plot.margin = margin(0, 0, 0, 0, "pt")
          )
  
  if(plot_quartiles) {
    
    plot <- plot +     #quartiles
    geom_segment(aes(y=0, yend=.01, x=quartiles[year==2018 & cum==.25, mid],
                     xend=quartiles[year==2018 & cum==.25, mid]), linetype='dotted', color='gray10') +
    geom_segment(aes(y=0, yend=.01, x=quartiles[year==2018 & cum==.5, mid],
                     xend=quartiles[year==2018 & cum==.5, mid]), linetype='dotted', color='gray10') +
    geom_segment(aes(y=0, yend=.01, x=quartiles[year==2018 & cum==.75, mid],
                     xend=quartiles[year==2018 & cum==.75, mid]), linetype='dotted', color='gray10') +
    geom_segment(aes(y=0, yend=-.01, x=quartiles[year==2000 & cum==.25, mid],
                     xend=quartiles[year==2000 & cum==.25, mid]), linetype='dotted', color='gray10') +
    geom_segment(aes(y=0, yend=-.01, x=quartiles[year==2000 & cum==.5, mid],
                     xend=quartiles[year==2000 & cum==.5, mid]), linetype='dotted', color='gray10') +
    geom_segment(aes(y=0, yend=-.01, x=quartiles[year==2000 & cum==.75, mid],
                     xend=quartiles[year==2000 & cum==.75, mid]), linetype='dotted', color='gray10')
    
  }
  
  if(flipped) plot <- plot + coord_flip() +  theme(axis.text.x=element_blank())
  else plot <- plot <- plot + theme(axis.text.y=element_blank())
  
  if(!is.na(density_max)) plot <- plot + scale_y_continuous(limits=c(-1*density_max, density_max))
  return(plot)
  
}

makeRegFigure3s <- function(loc_list) {

  #generate the plots
  pop_plot <-
    makeTapDensityPlot(input_dt=plot.dt, locs=loc_list, 
                       wt_var='exposed_pop', 
                       tap_cutoff = 750,
                       smoother=.1,
                       reg_plot=T)
  death_plot <-
    makeTapDensityPlot(input_dt=plot.dt, locs=loc_list, 
                       wt_var='atr_count_mean', 
                       tap_cutoff = 750,
                       smoother=.1,
                       reg_plot=T)

  #extract legend and set coordinates
  legend <- get_legend(pop_plot)
  legend$vp$x <- unit(0, 'npc')
  legend$vp$y <- unit(.65, 'npc')
  #define layout
  lay <- rbind(c(1,1,1,1),
               c(2,2,2,2))
  plot <- arrangeGrob(grobs=list(pop_plot + theme(legend.position = 'none'), 
                                 death_plot + theme(legend.position = 'none',
                                                    axis.text.x = element_blank())), 
                      layout_matrix=lay) %>% 
    grid.arrange

  #draw plot and legend/annotations
  png(paste0(save.dir, '/si/', 
             paste0('fig_3_',
                     unique(plot.dt[iso3%in%loc_list, region_id]), '_', paste0(loc_list, collapse="_")),
             '.png'),
      height=8, width=12, units='in', res=2400)
  grid.arrange(plot)
  grid.draw(legend)
  grid.text('a', x = unit(0.01, "npc"), y = unit(0.95, "npc"), gp = gpar(fontsize = 24, fontface = "bold"))
  grid.text('b', x = unit(0.01, "npc"), y = unit(0.05, "npc"), gp = gpar(fontsize = 24, fontface = "bold"))
  grid.text('2000', x = unit(.96, "npc"), y = unit(0.1, "npc"), gp = gpar(fontsize = 24, fontface = "bold"))
  grid.text('2018', x = unit(.96, "npc"), y = unit(0.4, "npc"), gp = gpar(fontsize = 24, fontface = "bold"))
  grid.text('2000', x = unit(.96, "npc"), y = unit(0.6, "npc"), gp = gpar(fontsize = 24, fontface = "bold"))
  grid.text('2018', x = unit(.96, "npc"), y = unit(0.9, "npc"), gp = gpar(fontsize = 24, fontface = "bold"))
  dev.off()
  
}

#make location lists (groups of 6 per super region)
makeFig3Loclist <- function(this.reg, dt=plot.dt) {

  locs <- dt[region_id==this.reg, unique(iso3)]
  if(length(locs)<=6) return(locs)
  else split(locs, cut(seq_along(locs), ceiling(length(locs)/6), labels = FALSE))
  
}

#helper function to flatten the above output
flatten2 <- function(x) {
  len <- sum(rapply(x, function(x) 1L))
  y <- vector('list', len)
  i <- 0L
  rapply(x, function(x) { i <<- i+1L; y[[i]] <<- x })
  y
}
#***********************************************************************************************************************

# ---FIGURE ?-----------------------------------------------------------------------------------------------------------