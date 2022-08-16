#custom function for calculating the LRI RR based on PM2.5
calcRR <- function(col.n, ratio, exp, exp.cols, curve, curve.n) {
  
  crude <- fobject$eval(exp[, get(exp.cols[col.n])], curve[[curve.n]][col.n, ])
  out <- ratio * crude - ratio + 1
  
  return(out)
  
}

#load the mr brt output
formatRR <- function(i, dir, debug=F) {
  
  if(debug) browser()
  
  message('reading mrBRT results for ', i)
  
  out <- file.path(dir, i, paste0(i, '_y_samples.csv')) %>% 
    fread %>% 
    melt(., id.var='exposure_spline', variable.name = 'draw', value.name='rr') %>% 
    .[, draw := gsub('draw_', 'V', draw)] %>% 
    .[, cause := i]
  
  #do some unique processing for cardiovascular outcomes
  if(i %like% 'cvd') {
    
    out <- out[, age_group_id := str_sub(cause, start=-2, end=-1) %>% as.numeric]
    out <- out[, cause := str_sub(cause, start=1, end=-4)]
    
    return(out)
    
  } else return(out)
  
}

# find and replace values in your data (also gives you the option to change the variable name, or to expand the rows by specifying more than one replacement for a single input)
findAndReplace <- function(table, #data table that you want to replace values in
                           input.vector, # vector of the old values you want to replace
                           output.vector, # vector of the new values you want to replace them with
                           input.variable.name, # the current variable name
                           output.variable.name, # new variable name (same as previous if no change required)
                           expand.option=FALSE) { # set this option TRUE if you are doing a one:many replace (IE expanding rows)
  
  values <- input.vector
  new.values <- output.vector
  
  # Replacement data table
  replacement.table = data.table(variable = values, new.variable = new.values, key = "variable")
  
  # Replace the values (note the data table must be keyed on the value to replace)
  setkeyv(table, input.variable.name) #need to use setkeyv (verbose) in order to pass the varname this way
  table <- setnames(
    replacement.table[table, allow.cartesian=expand.option][is.na(new.variable), new.variable := variable][,variable := NULL],
    'new.variable',
    output.variable.name)
  
  return(table)
  
}

#custom function to expand causes to the appropriate age/sex bins
causeExpansion <- function(input.table, cause.code) {
  
  message('expanding for ', cause.code)
  
  #expansion helper
  expandR <- function(i, var, dt) {
    
    message('->adding ', var, '=', i)
    
    copy(dt) %>% 
      .[, (var) := i] %>% 
      return
    
  }
  
  # Take out this cause
  out <- input.table[cause == cause.code]
  
  if (cause.code %in% c("cvd_ihd", "cvd_stroke")) {
    
    message('->replacing age name with age code')
    
    age.codes <- seq(25, 95, by=5)
    age.ids <- c(10:20, 30:32, 235)
    
    # then pass to your custom function
    out <- findAndReplace(out,
                          age.codes,
                          age.ids,
                          "age_group_id",
                          "age_group_id")
    
  } else {
    
    # Add back in with proper ages
    if (cause.code == "lri") age.ids <- c(2:20, 30:32, 235)
    if (cause.code %in% c("neo_lung", "resp_copd", "t2_dm")) age.ids <- c(10:20, 30:32, 235)
    
    out <- lapply(age.ids, expandR, var='age_group_id', dt=out) %>% rbindlist
    
  }
  
  out <- lapply(1:2, expandR, var='sex_id', dt=out) %>% rbindlist
  
  return(out)
  
}