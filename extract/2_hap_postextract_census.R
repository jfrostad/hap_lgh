rm(list=ls())

#Define values
topic <- "hap"
redownload <- F #update the codebook from google drive
cluster <- TRUE #running on cluster true/false
geos <- TRUE #running on geos nodes true/false
cores <- 3
#FOR THE CLUSTER:
#qlogin -now n -pe multi_slot 5 -P proj_geospatial -l geos_node=TRUE

# ####### YOU SHOULDN'T NEED TO CHANGE ANYTHING BELOW THIS LINE. SORRY IF YOU DO ##################################################


#Load packages
package_lib    <- file.path(h, '_code/_lib/pkg')
  .libPaths(package_lib)
pacman::p_load(data.table, dplyr, haven, feather, fst, googledrive, magrittr, parallel, doParallel, readxl, stringr)

#timestamp
today <- Sys.Date() %>% gsub("-", "_", .)

message("Getting common column names")
if (topic == "hap" & geos){
  #get the most recent pt and poly files and parse them for column names

    list.files(pattern="points", full.names=T) %>% 
    grep(pattern=".fst$", value=T) %>% 
    .[length(.)] %>% 
    read_fst(., from = 1, to = 2) %>% 
    names

    list.files(pattern="poly", full.names=T) %>% 
    grep(pattern=".fst$", value=T) %>% 
    .[length(.)] %>% 
    read_fst(., from = 1, to = 2) %>% 
    names
  noms <- c(pt_names, poly_names) %>% unique
} else{
  message("The error you're about to get has to do with the fact that you're not running on geos and/or you're not prepping hap data.")
  stop("I don't know how to parse .Rdata files for column headers in a timely way. Please figure something out and make a pull request.")
}

message("List IPUMS dtas")
extractions <- list.files(folder_in, pattern="IPUMS", full.names=T)

#define function to merge IPUMS files with geographies
ipums_merge <- function(file, geo, folder_out, noms){
  
  dt <- fread(file)

  #get survey info
  nid <- dt$nid[1] %>% as.character
  message(nid)
  iso3 <- dt$ihme_loc_id[1] %>% as.character
  year_start <- dt$year_start[1] %>% as.character
  year_end <- dt$year_end[1] %>% as.character
  survey_module <- dt$survey_module[1] %>% as.character
  
  #use new geocodebook database 
  geo <- get_geocodebooks(nids = nid)
  
  #skip bad data
  has_pweight_not_hhweight <- "pweight" %in% names(dt) & !("hhweight" %in% names(dt))
  missing_all_weights <- !("pweight" %in% names(dt)) & !("hhweight" %in% names(dt))
  missing_gid <- !("geospatial_id" %in% names(dt))
  missing_hap <- !("cooking_fuel" %in% names(dt)) & !(grepl('housing', names(dt)) %>% any)
  
  if (has_pweight_not_hhweight) setnames(dt, "pweight", "hhweight")

  if (missing_all_weights | missing_hap | missing_gid) return(nid)
  
  if (!missing_gid) dt[, geospatial_id := as.character(geospatial_id)] #force geospatial IDs to character to match the sheet
  m <- try(merge(dt, geo, by.x=c("nid", "ihme_loc_id", 'geospatial_id'), by.y=c("nid", 'iso3', 'geospatial_id'), all.x=T))

  if (class(m) == "try-error") {message(paste("Check try error", nid)); return(nid)}
  else{
    outname <- paste("IPUMS_CENSUS", nid, survey_module, iso3, year_start, year_end, sep="_")
    m$survey_series <- m$survey_name
    m$iso3 <- m$ihme_loc_id
    orig_names <- names(m)
    new_names <- noms[!(noms %in% orig_names)]
    for (nam in new_names){
      m[, nam] <- NA
    }
    #m <- m[, c("nid", 'survey_series', "year_start", "year_end", "iso3", "lat", "long", "shapefile", "location_code", "w_source_drink", "t_type", "sewage", "shared_san", "mins_ws", "dist_ws", "hw_station", "hw_water", "hw_soap", "hhweight", "hh_size", 'urban')]
    has_shp <- any(!is.na(m$shapefile)) & any(!is.na(m$location_code))
    has_lat_long <- any(!is.na(m$lat)) & any(!is.na(m$long))
    has_geography <- has_shp | has_lat_long
    if (!has_geography){message(paste("Check geography", nid, year_end)); return(nid)}

    out <- paste0(folder_out, "/", outname, ".fst")
    write.fst(m, out)
    return(NULL)
  }
}

# need_some_love <- "106512"
# files <- grep(need_some_love, files, value=T)

message("starting mclapply")
bad_nids <- mclapply(extractions, ipums_merge, geo=geo, folder_out=folder_out, noms=noms, mc.cores=cores) %>% unlist
write.csv(bad_nids, paste0(folder_out, "/fix_these_nids.csv"), row.names = F, na='')
#write.csv(to_do, to_do_outpath, row.names=F)

#####################################################################
#######Find broken extractions###################
#####################################################################

files <- substr(files, 1, nchar(files)-4)
files <- strsplit(files,'_')
files[6] #TODO jank
for(f in 1:length(files)){
  temp <- unlist(files[f])
  temp <- temp[6]
  files[f] <- temp
}
files <- unlist(files)
#find nids in codebook not in this list ^
codebook.nids <- codebook[!(year_end < 2000 | ihme_loc_id %in% stages[Stage==3, iso3]), nid] %>% unique
bad_ipums <- codebook.nids[!(codebook.nids %in% files)]

#check against bad nids csv
broken_ipums <- bad_nids[!(bad_nids %in% bad_ipums)]

##Ipums failed extractions
write.csv(broken_ipums, paste0(folder_out, '/failed_extractions_ipums.csv'))