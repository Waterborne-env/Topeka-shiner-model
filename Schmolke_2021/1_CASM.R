## Topeka Shiner CASM Set Up ###################################################
## Author: CR
## Last Edited: 2020-08-25
################################################################################
rm(list = ls())

## filepaths and variable definition ###########################################
dir.prj = 'D:\\796.173\\draft_github' # project directory
nam.exp = 'Pond19'  # Control, 0VFS, 15VFS, Pond18, Pond19
nam.mul = 10  # multiplication factor (ignored if dir.exp = Control)

## Initial Processing ##########################################################
# Scenario Name and Exposure Profile File
nam = paste0('CASM_', nam.exp)
if(nam.exp == 'Control'){
  dir.exp = 'exp_zero.dat'
} else if (nam.mul == 1) {
  dir.exp = paste0('exp_', nam.exp, '.dat')
} else {
  nam = paste0(nam, 'x', nam.mul)
  dir.exp = paste0('exp_', nam.exp, 'x', nam.mul, '.dat')
}

# Create Needed Folders
dir.create(file.path(dir.prj, nam)) # create scenario folder
setwd(file.path(dir.prj, nam)) # change directory
dir.create('OC_CASM_Out') # temp output folder

# create meta data file
writeLines(
  as.character(c('Meta Data for TS CASM Model Run:',
    paste('Scenario Name (nam):', nam),
    paste('Created:', Sys.time()),
    paste('Project Directory (dir.prj):', dir.prj),
    paste('Exposure Scenario (nam.exp):', nam.exp),
    paste('Exposure Multiplication Factor (nam.mul):', nam.mul)
  )), 'meta_scenario.txt')
rm(nam.exp, nam.mul, nam)

## Set Up Simulation ###########################################################
# copy needed files
ls.casm = c('casm_ibm_ocv12.exe', 'casm_IBM_OC_control.dat',
  'casm_OC_bio_parms_NT_26Jan2020_003_39_365.txt',
  'web_casmOC_12Aug2019.txt', 'casmOC_Env_Data_11Aug2019.prn',
  'OC_toxicity_Solatenol_tri_base.txt', dir.exp)
for(temp in ls.casm){
  file.copy(file.path(dir.prj, 'CASM_in', temp), temp)
}
rm(ls.casm, dir.prj, temp)

## change file paths for control file and other model parameters
temp = readLines('casm_IBM_OC_control.dat')
temp[21] = paste0('.\\', dir.exp)
writeLines(temp, 'casm_IBM_OC_control.dat')
rm(dir.exp, temp)

# Run CASM
system('casm_ibm_ocv12.exe')
