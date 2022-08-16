/***********************************************************************************************************
***********************************************************************************************************/


//////////////////////////////////
// Setup
//////////////////////////////////
clear all
set more off
set obs 1

///////////////////////////////////////////////////////////////////////
/*  Arguments:
		- topic: your research topic
        - array: Ubcov_id. The id of the codebook row
		- outpath_L: output file path for limited drive files
			(GBD: `l'/LIMITED_USE/LU_GBD/ubcov_extractions/{your topic folder})
			(LBD: `l'/LIMITED_USE/LU_GEOSPATIAL/ubCov_extractions/{your topic folder})
		- outpath_J: output file path for general files
    Optional:
        - keep: Keeps both raw data and extracted data, allows manual extraction check before final output.
        - bypass: Skips the extraction check, output the data 
        - run_all: Loops through all ubcov_ids in the codebook.
*/
////////////////////////////////////////////////////////////////////////

local topics hap
local array 1075
local options bypass //leave it blank if you don't want to keep, bypass, or run_all. If you are running on cluster, you have to use bypass as an option otherwise it will error.



///////////////////////////////////////////////////////////////////////
//Extraction
//////////////////////////////////////////////////////////////////////


// Load functions
cd "`central_root'"
do "`central_root'/modules/extract/core/load.do"

// Load the base code for ubCov
local paths  `central_root'/modules/extract/core/ `central_root'/modules/extract/core/addons/
foreach path in `paths' {
    local files : dir "`path'" files "*.do"
    foreach file in `files' {
        if "`file'" != "run.do" do "`path'/`file'"
    }
}

// Make sure you're in central
cd `central_root'

// Initialize the system
ubcov_path
init, topics(`topics')

// Launches extract
foreach number in `array'{
    local uid `number'
    run_extract `uid', `options'
    tostring year_start, gen(year)
    tostring year_end, gen(end_year)
    tostring nid, gen(nid_n)
    local filename = ihme_loc_id + "_" + survey_name + "_" + year + "_" + end_year + "_" + nid_n
    local filename = subinstr("`filename'", "/", "_",.)
    drop year end_year nid_n
	//if (strpos("$file_path", "LIMITED_USE")|strpos("$file_path", "IDENT")){
	//	local outpath = "`outpath_L'"
	//}
	//else{
	//	local outpath = "`outpath_J'"
	//}
	cd  "`outpath'"
    export delimited using "`filename'.csv", replace
}
