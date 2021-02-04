## Topeka Shiner IBM Only Set Up & Run #########################################
## Author: CR
## Last Edited: 2020-08-25
################################################################################
rm(list = ls())

## filepaths and variable definition ###########################################
dir.prj = 'D:\\796.173\\draft_github' # project directory
nam.exp = 'Pond'  # Control, 0VFS, 15VFS, Pond
nam.mul = 10  # multiplication factor (ignored if dir.exp = Control)
ibm.eff = 'GUTS-SD + Sublethal' # Effects Input Value: 'GUTS-IT', 'GUTS-SD',
  # 'GUTS-IT + Sublethal', 'GUTS-SD + Sublethal', or 'Sublethal'

## needed packages #############################################################
library('RNetLogo')

## Initial Processing ##########################################################
# Scenario Name and Exposure Profile File
ls.nam = c('GUTS-IT' = 'IT', 'GUTS-SD' = 'SD', 'Sublethal' = 'sub',
  'GUTS-IT + Sublethal' = 'ITsub', 'GUTS-SD + Sublethal' = 'SDsub')
casm.nam = paste0('CASM_', nam.exp)
nam = paste0('IBM_', nam.exp)
if(nam.exp == 'Control'){
  dir.exp = 'exp_zero.dat'
} else if (nam.mul == 1) {
  dir.exp = paste0('exp_', nam.exp, '.dat')
  nam = paste0(nam, '_', ls.nam[[ibm.eff]])
} else {
  dir.exp = paste0('exp_', nam.exp, 'x', nam.mul, '.dat')
  casm.nam = paste0(casm.nam, 'x', nam.mul)
  nam = paste0(nam, 'x', nam.mul, '_', ls.nam[[ibm.eff]])
}
rm(ls.nam)
if(nam.exp == 'Pond'){
  dir.exp = c(sub('d', 'd18', dir.exp), sub('d', 'd19', dir.exp))
  casm.nam = c(sub('d', 'd18', casm.nam), sub('d', 'd19', casm.nam))
}

# Create Needed Folders
dir.create(file.path(dir.prj, nam)) # create scenario folder
setwd(file.path(dir.prj, nam)) # change directory

# create meta data file
writeLines(
  as.character(c('Meta Data for TS IBM Model Run:',
    paste('Scenario Name (nam):', nam),
    paste('Created:', Sys.time()),
    paste('Project Directory (dir.prj):', dir.prj),
    paste('Exposure Scenario (nam.exp):', nam.exp),
    paste('Exposure Multiplication Factor (nam.mul):', nam.mul),
    paste('IBM Effects Input (ibm.eff):', ibm.eff)
  )), 'meta_scenario.txt')
rm(nam.mul)

## Set Up Simulation ###########################################################
# copy needed files
ls.ibm = c('TS_IBM_Effects_May2020.nlogo', 'GUTS_sublethal_May2020.nls',
  'InputParameters_casmOC_27Feb2020.txt')
for(temp in ls.ibm){
  file.copy(file.path(dir.prj, 'IBM_in', temp), temp)
}
rm(ls.ibm, temp)

# copy needed files
ls.casm = c('IBM_TS_master_ref.out', 'casmOC_Env_Data_11Aug2019.prn')
if(nam.exp != 'Control'){
  ls.casm = c(ls.casm, dir.exp[1], 'IBM_TS_master_eff.out')
}
file.copy(file.path(dir.prj, 'CASM_Control', 'exp_zero.dat'), 'exp_zero.dat')
for(temp in ls.casm){
  file.copy(file.path(dir.prj, casm.nam[1], temp), temp)
}
if(nam.exp == 'Pond'){
  file.rename('IBM_TS_master_eff.out', 'IBM_TS_master_eff18.out')
  file.copy(file.path(dir.prj, casm.nam[2], 'IBM_TS_master_eff.out'),
    'IBM_TS_master_eff19.out')
  file.copy(file.path(dir.prj, casm.nam[2], dir.exp[2]), dir.exp[2])
}
rm(ls.casm, casm.nam, temp)

## alter CASM outputs for simulation
for(i in list.files(pattern = 'IBM_TS_master_*')){
  temp = readLines(i)
  temp = c('IBM-CASM transfer values:', 
    temp[2], substr(temp[3], 1, 1136), substr(temp[4:368], 1, 1136))
  substr(temp[4], 323, 336) = formatC(0.195, digits = 8, format = 'E')
  writeLines(temp, i)
}
rm(i, temp)

## IBM Yearly files
writeLines(c(12, rep('casmOC_Env_Data_11Aug2019.prn', 12)), 'WCI.txt')
if(nam.exp == 'Control'){
  writeLines(c(12, rep('IBM_TS_master_ref.out', 12)), 'DPB.txt')
  writeLines(c(12, rep('exp_zero.dat', 12)), 'EP.txt')
} else if(nam.exp %in% c('0VFS', '15VFS')){
  writeLines(c(12, rep('IBM_TS_master_ref.out', 5), 'IBM_TS_master_eff.out',
    rep('IBM_TS_master_ref.out', 6)), 'DPB.txt')
  writeLines(c(12, rep('exp_zero.dat', 5), dir.exp, rep('exp_zero.dat', 6)),
    'EP.txt')
} else if(nam.exp == 'Pond'){
  writeLines(c(12, rep('IBM_TS_master_ref.out', 5), 'IBM_TS_master_eff18.out',
    'IBM_TS_master_eff19.out', rep('IBM_TS_master_ref.out', 5)), 'DPB.txt')
  writeLines(c(12, rep('exp_zero.dat', 5), dir.exp, rep('exp_zero.dat', 5)),
    'EP.txt')
}
rm(nam.exp, dir.exp)

## Run TS IBM Stand Alone ######################################################
ibm.rns = c(6,35,10,88,58,87,36,2,79,77,21,26,64,8,31,70,92,39,15,32)
ibm.inp  = list(
  'InputParameters' = 'InputParameters_casmOC_27Feb2020.txt',
  'WaterConditionInput' = 'WCI.txt', 'DailyPreyBiomass' = 'DPB.txt',
  'ExposureProfile' = 'EP.txt', 'DailyIO' = FALSE, 'MultiYearInputs' = TRUE,
  'pondArea' = 100, 'CASM_pool_area' = 1125, 'Sunfish' = TRUE,
  'Predation' = FALSE, 'DD_egg_larva' = TRUE, 'AverageKmin' = 0.68,
  'ScalingSearchArea' = 1.0, 'MaxDetritusDepth' = 1.5,
  'EffectsModule' = TRUE, 'Effects' = ibm.eff,
  'PesticideName' = 'Solatenol', 'ConstantExposure' = FALSE,
  'guts_n' = 1440, 'GUTS-SD-ke' = 1.28, 'GUTS-SD-kk' = 0.42,
  'GUTS-SD-z' = 3.85, 'GUTS-IT-ke' = 0.47, 'GUTS-IT-a' = 3.97,
  'GUTS-IT-b' = 5.54, 'TKTD-sl-ke' = 0.00061, 'TKTD-sl-cT' = 8.8389,
  'TS_max_SL' = 53.2, 'TKTD-NOEC' = 0.95)
rm(ibm.eff)
for(rns in ibm.rns){
  ## Load netlogo model
  setwd('C:/Program Files/NetLogo 6.0.2/app')
  NLStart('C:/Program Files/NetLogo 6.0.2/app',
    nl.jarname = 'netlogo-6.0.2.jar', gui = FALSE)
  NLLoadModel(file.path(dir.prj, nam, 'TS_IBM_Effects_May2020.nlogo'))
  
  ## Netlogo interface inputs
  NLCommand(paste0('set BiomassOutput "output_', rns, '.txt"'))
  NLCommand(paste('set RandomNumberSeed', rns))
  for (input in names(ibm.inp)){
    value = ibm.inp[[input]] # get value
    if(is.character(value)) value = paste0('"', value, '"') # add quotes to char
    if(is.logical(value)) value = tolower(value) # make T/F lowercase
    NLCommand(paste('set', input, value))
  }
  rm(input, value)
  
  ## Run every day
  setwd(file.path(dir.prj, nam)) # change directory
  NLCommand('setup') # set up run
  for(i in 1:(365*12)) NLCommand('go') # run TS
  NLCommand('write_output') # write TS output
  NLQuit() # quit session
}
rm(rns, ibm.inp, ibm.rns, nam, dir.prj)
