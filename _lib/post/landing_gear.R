#custom function to save files into the central mapping format
saveMappingInput <- function(dt, 
                             map_ind,
                             data_ind=map_ind, #but it could be different if misnamed
                             map_measure,
                             data_measure=map_measure, #but it could be different if misnamed 
                             raking_label=raked,
                             admin_level=2,
                             parent_dir=map.dir
) {
  
  #define dirs and paths (also create recursively)
  out_dir <- file.path(parent_dir, map_ind, run_date, 'inputs') %T>% dir.create(recursive=T) 
  out_path <- file.path(out_dir, paste0(map_ind, '_', map_measure,
                                        ifelse(raking_label, '_raked_', '_unraked_'), 
                                        'ad', admin_level, '.csv')) %T>%
    message('saving data for ind=', data_ind, ' (', map_measure, ') as\n',
            .)
  
  #setup dt
  out <- copy(dt) %>% 
    .[, c(paste0('ADM', admin_level, '_CODE'), paste(data_ind, data_measure, sep='_'), 'year'), with=F] %>% 
    setnames(paste(data_ind, data_measure, sep='_'), 'value')
  
  #save
  write.csv(out, file=out_path, row.names=F)
  
}

#helper fx to pull/prep the appropriate files from our list of SDG projection objects
prepCasts <- function(id, type, list=sdg_files, id_dt=NA, id_var=NA) {
  
  #format ID var if necessary
  if(nchar(id)==4) id <- as.character(id) #if the ID is a year, format as character
  
  message(id)
  
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


#helper fx to resample  draws for each variable and then do weighted aggregations in order to produce the AD0 stats
resampleDraws <- function(dt, var, by_vars,
                          n_draws=50,
                          wt_var='pop',
                          aggregate=F,
                          make_ci=F,
                          debug=F) {
  
  if (debug) browser()
  
  message('working on ', var)
  
  #define variables
  var_types <- c('lower', 'mean', 'upper')
  id_vars <- names(dt) %>% .[!(. %like% 'pop|lower|mean|upper')]
  ind_vars <- paste(var, var_types, sep="_")
  draw_vars <- paste0('draw_', 1:n_draws)
  
  #define transformation based on var
  logit_space <- !(var %like% 'rate|count|_pc')
  count_space <- var %like% '_count|_c'
  
  #copy dt with relevant vars
  out <- dt[, c(id_vars, ind_vars, wt_var), with=F] %>% 
    copy %>% 
    setnames(., c(ind_vars, wt_var), c(var_types, 'n'))
  
  #resample draws using delta method
  message('resampling')
  out[, row_id := .I]
  out[, se := (upper-lower)/3.92] #median N per cluster
  
  if(logit_space) {
    
    out[, logit_se := sqrt((1/(mean - mean^2))^2 * se^2)]
    out[, logit_mean := logit(mean)]
    out[, (draw_vars) := lapply(1:n_draws, function(x) rnorm(.N, logit_mean, logit_se))]
    out[, (draw_vars) := lapply(.SD, inv.logit), .SDcols=draw_vars]
    
  } else out[, (draw_vars) := lapply(1:n_draws, function(x) msm::rtnorm(.N, mean, se, lower=0))]
  
  if(aggregate) {
    
    message('aggregating')
    if (count_space) out[, (draw_vars) := lapply(.SD, sum, na.rm=T), .SDcols=draw_vars, by=by_vars]
    else out[, (draw_vars) := lapply(.SD,  Hmisc::wtd.mean, weights=n), .SDcols=draw_vars, by=by_vars]
    out <- unique(out, by=by_vars)
    
  }
  
  if(make_ci) {
    
    message('summarizing')
    out[, mean := rowMeans(.SD), .SDcols=draw_vars]
    out[, upper := apply(.SD, 1, quantile, p=.975, na.rm=T), .SDcols=draw_vars]
    out[, lower := apply(.SD, 1, quantile, p=.025, na.rm=T), .SDcols=draw_vars]
    
    out[, (draw_vars) := NULL]
    
    setnames(out, var_types, ind_vars)
    setkeyv(out, id_vars)
    
    out[, c(id_vars, ind_vars), with=F] %>% return
    
  } else {
    
    message('reshaping long')  
    
    out <- out[, c(id_vars, draw_vars), with=F] %>% 
      melt(., id.vars=id_vars, value.name=var) %>% 
      return
    
  }
  
  
}

#calculate mean and CI at different levels of aggregation
calcSummaryStats <- function(i, by_list, dt, labs_dt=adm_links, ind_cols, 
                             year_list=c(start_year:end_year),
                             change_years=list(start=start_year, end=end_year)) {

  by_cols <- by_list[[i]]
  this_dimension <- names(by_list)[i]
  
  message('summarizing ', paste(ind_cols, collapse='|'), ' -- lvl=',  paste(this_dimension, collapse='|'))
  
  agg <- dt[year %in% year_list] %>% copy
  
  #decide if collapse by population is required based on dimensions of bycols
  collapsing <- nrow(unique(agg,  by=c(by_cols, 'draw')))!=nrow(agg)
  
  if(collapsing) {
    
    message('population weighting to aggregate')
    
    setkeyv(agg, c(by_cols, 'draw')) #rekey for speed
    
    #distinguish that count columns should be a weighted sum instead of a weighted mean
    sum_cols <- ind_cols %>% .[. %like% 'count'] 
    mean_cols <- ind_cols %>% .[!(. %in% sum_cols)] 
    
    #check which pop columns were returned (pop_total is produced if using a non-total pop to aggregate)
    #append them to sum_cols, they will be likewise summed in the aggregation step
    pop_cols <- names(agg) %>% .[(. %like% 'pop')]
    sum_cols <- c(sum_cols, pop_cols)

    agg[, paste0(mean_cols) := lapply(.SD, function(x) sum(x*pop, na.rm=T)/sum(pop, na.rm=T)), 
        .SDcols=mean_cols, by=key(agg)] 
    agg[, (sum_cols) := lapply(.SD, sum, na.rm=T), .SDcols=sum_cols, by=key(agg)] 
    agg <- unique(agg, by=key(agg))
    
  }
  
  message('calculating change stats')
  #agg[, term := 'lvl'] #note that vals are currently calculating the lvl for a given year
  
  #add term variable to key
  by_cols <- c(by_cols, 'term')

  dt <- copy(agg)
  agg[, term:='lvl']
  agg[, start_year := start_year]
  
    
    setkeyv(agg, by_cols)  #rekey for speed
    agg[, paste0(ind_cols, '_lower') := lapply(.SD, quantile, p=.025, na.rm=T), .SDcols=ind_cols, by=by_cols]
    agg[, paste0(ind_cols, '_mean') := lapply(.SD, mean, na.rm=T), .SDcols=ind_cols, by=by_cols]
    agg[, paste0(ind_cols, '_upper') := lapply(.SD, quantile, p=.975, na.rm=T), .SDcols=ind_cols, by=by_cols]
    agg[, paste0(ind_cols, '_sd') := lapply(.SD, sd, na.rm=T), .SDcols=ind_cols, by=by_cols]

    
  agg <- unique(agg,  by=by_cols)
  
  #reorder names for legibility
  agg <- agg[, c(by_cols,
                 'start_year',
                 names(agg) %>% .[. %like% 'pop'],
                 names(agg) %>% .[. %like% 'lower|mean|upper|sd'] %>% sort), with=F]
  
  #merge on labels, unless global
  if(this_dimension!='global') { 
    
    message('merging labels')
    loc_cols <- by_cols %>% .[. %like% 'id|CODE']
    lab_cols <- loc_cols %>% 
      str_replace(., 'id', 'name') %>% 
      str_replace(., 'CODE', 'NAME')
    agg <- labs_dt[, c(loc_cols, lab_cols), with=F] %>% 
      unique(by=loc_cols) %>% 
      merge(agg, ., by=loc_cols)
    
  }
  
  agg[, dimension := this_dimension] %>% return
  
}

