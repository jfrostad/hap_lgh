##############################################################################
##############################################################################


## Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

# set general arguments
user            <- Sys.info()['user']
repo            <- file.path('/homes', user, '_code/lbd/hap/')
indicator_group <- 'cooking'

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
mbg_setup(package_list = package_list, repos = repo)

# set cluster arguments
use_geos_nodes  <- T
proj_arg        <- ifelse(use_geos_nodes, 'proj_geo_nodes', 'proj_geospatial_dia')
proj            <- ifelse(use_geos_nodes, paste0(' -P ', proj_arg, ' -l gn=TRUE '), paste0(' -P ', proj_arg, ' '))

# list stages
stages <- c('1', '2a', '2b')

## Run launch scripts -------------------------------------------------------------------------

for (s in stages) {

    
    # set specific arguments
    jname           <- paste(indicator_group, 'post_extract', s, sep = '_')
    mymem           <- '500G'

    # set up qsub
    sys.sub <- paste0('qsub -e /share/temp/sgeoutput/', user,'/errors -o /share/temp/sgeoutput/', user, '/output ', 
                      '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                      '-l fthread=10 -l h_rt=00:16:00:00 -v sing_image=default -N ', jname, ' -l archive=TRUE ')
    r_shell <- file.path(repo, 'mbg_central/share_scripts/shell_sing.sh')
    script <- file.path(repo, 'extract/2b_postextract.R')
    args <- paste(s, repo, indicator_group)

    # run launch script
    paste(sys.sub, r_shell, script, args) %>% 
      system

}