## Topeka Shiner CASM Set Up ###################################################
## Author: CR
## Last Edited: 2019-04-03
################################################################################
rm(list = ls())

## filepaths and variable definition ###########################################
# Initial Folder Locations/Names
dir.prj  = 'D:\\Topeka-shiner-model-master' # project directory
dir.run  = 'CASMRun' # scenario name

# CASM Files (note: relative to dir.prj)
dir.exe = 'casm_ibm_tsv42.exe'# CASM exe
dir.cntl = 'casm_IBM_TS_control.dat' # CASM control file
dir.hab = 'casm_TS_bio_parms_RSK_WithTS_08Oct2018.txt' # CASM habitat file
dir.web = 'web_casmTS_HRM_20Jun2018.dat' # CASM web file
dir.env = 'env_casmTS_2010_08May2018.prn' # CASM environmental file
dir.con = 'ATZ_exp_zero_test.dat' # CASM contaminant file
dir.tox = 'TS_toxicity_atrazine_tri_base.dat' # CASM toxicity file

# Model Parameters
casm.yrs = 1 # years to simulate (changes dir.cntl)

## processing ##################################################################
setwd(dir.prj) # change directory
dir.create(dir.run) # create scenario folder

# create meta data file
writeLines(
  as.character(c('Meta Data for TS CASM Model Run:',
    paste('Scenario Name (dir.nam):', dir.run),
    paste('Scenario Location (dir.run):', dir.run),
    paste('Created:', Sys.time()),
    paste('Project Directory (dir.prj):', dir.prj),
    paste('Model Directory (dir.exe):', dir.exe),
    paste('Control File Location (dir.cntl):', dir.cntl),
    paste('Habitat File Location (dir.hab):', dir.hab),
    paste('Web File Location (dir.web):', dir.web),
    paste('Environmental File Location (dir.env):', dir.env),
    paste('Contaminant File Location (dir.con):', dir.con),
    paste('Toxicity File Location (dir.tox):', dir.tox),
    paste('Years Run (casm.yrs):', casm.yrs)
  )), file.path(dir.run, 'meta_scenario.txt'))

## Set Up Simulation ###########################################################
dir.create(file.path(dir.run, 'TS_CASM_Out')) # temp output folder
dir.create(file.path(dir.run, 'TS_CASM_Out', 'plots')) # temp output folder

# copy needed files
for(temp in c(dir.exe, dir.cntl, dir.hab, dir.web, dir.env, dir.con, dir.tox,
  'libgcc_s_seh-1.dll', 'libgfortran-3.dll', 'libquadmath-0.dll')){
  file.copy(temp, file.path(dir.run, basename(temp)))
}
rm(temp)

## change file paths for dir.cntl and other model parameters
temp = readLines(file.path(dir.run, basename(dir.cntl)))
temp[3]  = 0
temp[7]  = paste0('.\\', basename(dir.hab))
temp[9]  = paste0('.\\', basename(dir.web))
temp[11] = paste0('.\\', basename(dir.env))
temp[21] = paste0('.\\', basename(dir.con))
temp[33] = paste0('.\\', basename(dir.tox))
temp[39:44] = paste0('.\\TS_CASM_Out\\', basename(temp[39:44]))
temp[45:72] = paste0('.\\TS_CASM_Out\\plots\\', basename(temp[45:72]))
temp[74] = '.\\IBM_TS_master_day1.out'
temp[75] = '.\\IBM_TS_master_env1.out'
temp[77] = 0 # don't use IBM model
temp[79] = 0 # 0 - ref, 1 - eff
temp[81] = '.\\IBM_TS_transfer.out'
temp[82] = '.\\IBM_TS_transfer_env.out'
temp[84] = '.\\IBM_TS_master_ref.out'
temp[85] = '.\\IBM_TS_master_eff.out'
temp[87] = casm.yrs # run x years
writeLines(temp, file.path(dir.run, basename(dir.cntl)))
rm(temp)

## Clean Up environment ########################################################
rm(dir.prj, dir.run, dir.exe, dir.cntl, dir.hab, dir.web, dir.env, dir.con,
  dir.tox, casm.yrs)
