####################################################################################################
## Description:   Make maps of estimates in Africa for specified years and year pairs.
##
## Inputs:        mean rasters
##                adm0 and adm1 estimates
##                country outlines, lakes, and population mask
##
##                  [run_date]/results_maps/[indicator]_raked_mean.pdf').
####################################################################################################
## load_map_annotations ----------------------------------------------------------------------------

load_map_annotations <- function(use.sf=T, mask_stage3=T) {

  if(use.sf) {
    
    ## Base shapefile (country outlines)
    message('->loading country borders')
    stage1 <- st_read('/home/j/WORK/11_geospatial/09_MBG_maps/misc_files/shps_by_stage/stage1_ad0_gadm.shp')
    stage2 <- st_read('/home/j/WORK/11_geospatial/09_MBG_maps/misc_files/shps_by_stage/stage2_ad0_gadm.shp')
    adm0 <- rbind(stage1, stage2)
    
    ## Stage 3 mask
    message("---->loading stage 3 mask")
    stage3 <- st_read('/home/j/WORK/11_geospatial/09_MBG_maps/misc_files/global_files/stage_3_mask.shp')

    ## Lakes
    message('-->loading lake mask')
    lakes <- raster('/home/j/WORK/11_geospatial/09_MBG_maps/misc_files/global_files/global_lakes.tif') %>% 
      { if (mask_stage3) mask(x=., mask=stage3, inverse=T) %>% crop(., adm0) else . } %>% 
      as(., 'SpatialPolygonsDataFrame') %>% 
      st_as_sf

    
    ## Population mask
    message('--->loading population mask')
    mask <- raster('/home/j/WORK/11_geospatial/09_MBG_maps/misc_files/global_files/global_mask_master.tif') %>% 
      { if (mask_stage3) mask(x=., mask=stage3, inverse=T) %>% crop(., adm0) else . } %>% 
      as(., 'SpatialPolygonsDataFrame') %>% 
      st_as_sf

  } else {
  
    ## Base shapefile (country outlines)
    message('->loading country borders')
    stage1 <- shapefile('/home/j/WORK/11_geospatial/09_MBG_maps/misc_files/shps_by_stage/stage1_ad0_gadm.shp')
    stage2 <- shapefile('/home/j/WORK/11_geospatial/09_MBG_maps/misc_files/shps_by_stage/stage2_ad0_gadm.shp')
    adm0 <- bind(stage1, stage2) %>% fortify
  
    ## Lakes
    message('-->loading lake mask')
    lakes <- raster('/home/j/WORK/11_geospatial/09_MBG_maps/misc_files/global_files/global_lakes.tif') %>% 
      rasterToPoints %>% 
      as.data.table %>% 
      setnames(., c("long", 'lat', 'lakes'))
  
    ## Population mask
    message('--->loading population mask')
    mask <- raster('/home/j/WORK/11_geospatial/09_MBG_maps/misc_files/global_files/global_mask_master.tif') %>% 
    rasterToPoints %>% 
      as.data.table %>% 
      setnames(., c("long", 'lat', 'mask'))
    
    ## Stage 3 mask
    message("---->loading stage 3 mask")
    stage3 <- shapefile('/home/j/WORK/11_geospatial/09_MBG_maps/misc_files/global_files/stage_3_mask.shp') %>% fortify
  
  }

  return(list(adm0 = adm0, lakes = lakes, mask = mask, stage3 = stage3))

}

## load_map_results --------------------------------------------------------------------------------

load_map_results <- function(indicator=NA, indicator_group=NA, run_date=NA, raked=NA, year_list=NA,
                             geo_levels = c("raster", "admin1", "admin2"),
                             custom_path=NULL, subvars=NULL,
                             use.sf=T,
                             cores=1, 
                             debug=F) {
  
  if (debug) browser()

  ## Set the input directory
  maps_path <- paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date)
  
  ## Helper function to prepare the admin level results
  prep_admin_x <- function(x) {
    
    #setup the column name with requested division
    code_var <- paste0('ADM', x, '_CODE')
    
    #read from standard location if not provided with custom path
    if(custom_path %>% is.null) {
      pred <- paste0(maps_path, "/pred_derivatives/admin_summaries/", indicator, "_admin_", x, 
                     ifelse(raked, "_raked", "_unraked"), 
                     "_summary.csv") %>% fread
    } else if (custom_path[[ paste0('admin', x)]] %>% is.data.table) { #or just use the file if in memory
      pred <- custom_path[[ paste0('admin', x)]]
    } else pred <- custom_path[[ paste0('admin', x)]] %>% fread
    
    message('-> admin', x, ' found and fread')
    
    #remove ADM NAMES if already added to pred dt
    col_names <- names(pred) %>% .[!(. %like% 'NAME')]
    
    if (subvars %>% is.null) pred <- pred[year %in% year_list, (col_names), with=F] 
    else pred <- pred[year %in% year_list, c('ADM0_CODE', code_var, 'year', subvars), with = F] #subset vars if req

    if (use.sf) {
      
      shp <- get_admin_shapefile(x) %>% 
        st_read %>% 
        filter(get(code_var) %in% pred[,  get(code_var)]) %>% 
        merge(., pred, by=c('ADM0_CODE', code_var), allow.cartesian=T)
      
      message('--> admin', x, ' results merged to sf')
      
    } 
    return(shp)
    
  }

  ## raster estimates
  if ("raster" %in% geo_levels) {
    message('loading raster data')
    if(missing(custom_path)) {
      raster <- paste0(maps_path, "/", indicator, "_", 
                       ifelse(type == "cirange", "range", type),
                       "_", 
                       ifelse(raked, "raked_", ""), 
                       years, ".tif") %>% brick
    } else raster <- brick(custom_path$raster)
    
    message('-> raster found and bricked')

    raster <- mclapply(year_list, 
                       function(y) {

                         message('--> sending to points (year=', y, ')')
                         df <- rasterToPoints(raster[[y - 1999]]) %>% data.table
                         setnames(df, c("long", 'lat', 'outcome'))
                         df[, year := y]
                         
                         return(df)
                         
                       },
                       mc.cores=cores) %>% 
      rbindlist %>% 
      setkey(., long, lat, year)
      
    message('--> converted to dt and keyed')
    
  }
  
  ## admin1 estimates and shape file
  if ("admin0" %in% geo_levels) admin0 <- prep_admin_x(0)

  ## admin1 estimates and shape file
  if ("admin1" %in% geo_levels) admin1 <- prep_admin_x(1)

  ## admin2 estimates and shape file
  if ("admin2" %in% geo_levels) admin2 <- prep_admin_x(2)
    
  ## combine and return all estimates
  mget(geo_levels) %>% return

}

## calc_diff_map -----------------------------------------------------------------------------------

calc_diff_map <- function(pred, diff_years) {
  diff <- lapply(names(pred), function(g) {
    rbindlist(lapply(diff_years, function(y) {
      temp <- pred[[g]][year %in% y, ]
      temp <- temp[, list(outcome = outcome[year == y[2]] - outcome[year == y[1]]), by=setdiff(names(temp), c("outcome", "year"))]
      temp[, years := paste(y, collapse="-")]
      temp
    }))
  })
  names(diff) <- names(pred)
  return(diff)
}

## plot_map ----------------------------------------------------------------------------------------
plot_map <- function(map_data, annotations, title, limits, this_var='outcome',
                     legend_colors, legend_color_values, legend_breaks, legend_labels, legend_title, 
                     legend_flip=F,
                     pop_mask=T, lake_mask=T, borders=T, stage3_mask=T,
                     subset, zoom,
                     debug=F) {
  
  if (debug) browser()
  
  message('checking arguments')
  #set function arguments based on argument classes
  map_sf <- 'sf' %in% class(map_data)
  annotate_sf <- 'sf' %in% sapply(annotations, class)
  if(map_sf&annotate_sf) message('-> using SF to render plots, vroom!')
  else stop('sorry, geom_polygons have been deprecated, please provide data & annotations as sf objects')
  
  
  ##subset map_data
  if(!missing(subset)) {
    
    if(subset %>% is.list) {
      
      message('subsetting data to ', subset$var, '==', subset$value)
      map_data <- filter(map_data, get(subset$var)==subset$value)
      
    } else stop('subset must be provided as a list with named objects var and value')
    
  }
  
  message('building limits/scale')
  ## Enforce limits & define plot_var for simplicity
  map_data$plot_var <- pmax(limits[1], pmin(limits[2], as.data.frame(map_data)[, this_var])) 
  
  ##Setup legend and scales if not provided by user
  if(missing(legend_colors)) { #Default is a 10 color magma from viridis scale
    
    message('->color scale not provided, using -magma- from viridis scales')
    
    
    legend_colors <- c("#FCFDBFFF", "#FEC98DFF", "#FD9567FF", "#F1605DFF", "#CD4071FF",
                       "#9F2F7FFF", "#721F81FF", "#451077FF", "#180F3EFF", "#000004FF")
    
  }
  
  if(missing(legend_color_values)) { #Default is to use the distribution of data to build a color scale
    
    message('->color values not provided, building from data IQR')
    
    data_min <- map_data$plot_var %>% min(na.rm=T) %>% floor
    data_p25 <- map_data$plot_var %>% quantile(probs=.2, na.rm=T)
    data_p75 <- map_data$plot_var %>% quantile(probs=.8, na.rm=T)
    data_max <- map_data$plot_var %>% max(na.rm=T) %>% ceiling
    
    legend_color_values <- c(seq(data_min, data_p25, length.out = 2), 
                             seq(data_p25, data_p75, length.out = 8), #bring out the variation in the IQR
                             seq(data_p75, data_max, length.out = 2)) %>%
      unique %T>% 
      print %>%
      rescale
      
  }

  if (missing(legend_breaks) & missing(legend_labels)) {
    
    message('->labels not provided, building from color values')

    start_range <- range(map_data$plot_var, na.rm = T)

    ## Create breaks
    breaks <- pretty(limits, 5)
    if (limits[1] < 0 & limits[2] > 0) breaks <- sort(unique(c(0, breaks)))
  
    ## Create labels
    labels <- format(breaks, nsmall = 2)
    if (min(limits) >= 0) divider <- "-" else divider <- " to "
    if (start_range[1] < limits[1]) {
      labels[1] <- paste0(format(floor(100*start_range[1])/100, nsmall=2), divider, labels[1])
    }
    if (start_range[2] > limits[2]) {
      labels[length(labels)] <- paste0(labels[length(labels)], divider, format(ceiling(100*start_range[2])/100, nsmall=2))
    }
    
  } else {breaks <- legend_breaks; labels <- legend_labels}
  
  message('plotting canvas')
  
  ## Crop the annotations for speed
  if (!missing(zoom) & annotate_sf) { 
    message('->cropping canvas to zoom')
    
    if(zoom %>% is.data.table) { 
      
      message('-->using custom zoom')
      annotations <- lapply(annotations, st_crop, xmin=zoom$x1, xmax=zoom$x2, ymin=zoom$y1, ymax=zoom$y2)
    
    } else {
      
      message('-->using zoom from data with 1 degree buffer')
      
      coords <- st_coordinates(map_data) %>% as.data.table
      zoom <- data.table(
        x1=floor(min(coords$X)) - 1,
        x2=ceiling(max(coords$X)) + 1,
        y1=floor(min(coords$Y)) - 1, 
        y2=ceiling(max(coords$Y)) + 1
      )
      
      annotations <- lapply(annotations, st_crop, xmin=zoom$x1, xmax=zoom$x2, ymin=zoom$y1, ymax=zoom$y2)
      
    }
  }  
  
  ## Plot the base map (this is what shows in places with no estimates and no mask)
  canvas <- ggplot() + geom_sf(data = annotations$adm0, lwd=0.1, color = 'black', fill = 'gray90')

  message('plotting outcome')

  ## Plot predictions
  if (map_sf) {
    message('-> outcome is aggregated')
    gg <- canvas + geom_sf(data = map_data, aes(fill = plot_var), lwd=0) + coord_sf(datum = NA)
  } else {
    message('-> outcome is at pixel level')
    gg <- canvas + geom_raster(data = map_data, aes(fill = plot_var, y = lat, x = long)) + 
      coord_equal(ratio = 1)
  }

  message('plotting annotations')
  
  ## Plot mask, lakes, and adm borders using SF
  if (pop_mask) gg <- gg + geom_sf(data = annotations$mask, lwd=0, color = 'gray70', fill = 'gray70')
  if (lake_mask) gg <- gg + geom_sf(data = annotations$lakes, lwd=0, color = 'gray70', fill = 'lightblue')
  if (borders) gg <- gg + geom_sf(data = annotations$adm0, lwd=0.1, color = 'black', fill=NA)
  if (stage3_mask) gg <- gg + geom_sf(data = annotations$stage3, lwd=0, color = 'gray70')

  message('defining aesthetics')
  
  ## Scales
  gg <- gg +
    scale_fill_gradientn(colors = legend_colors, values = legend_color_values,
                         limits = range(breaks), breaks = breaks, labels = labels, name = legend_title)

  ## Labels & aesthetics
  gg <- gg +
    labs(x="", y="", title=title) +
    theme_classic() +
    theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          legend.position = c(ifelse(legend_flip, .8, 0), 0), legend.justification = c(0, 0),
          #legend.text=element_text(size=10),
          plot.title = element_text(hjust=0.5), plot.margin=unit(c(0, 0, 0, 0), "in")) +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 7))

  return(gg)
  
}

#function to make quick map data for ctry diagnostics
quick_data <- function(dt, var, country, year) {

  load_map_results(indicator='x', indicator_group='x', run_date='x', raked='x', year_list=year,
                   custom_path = list('admin2'=dt),
                   geo_levels=c('admin2'), subvars=var,
                   cores=cores) %>%
    .[['admin2']] %>% 
    filter(NAME_0==country) %>% 
    return

} 
