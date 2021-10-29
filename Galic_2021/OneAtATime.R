## Topeka Shiner One At A Time Runs ############################################
## Author: CR
## Last Edited: 2021-09-20
################################################################################
rm(list = ls())

## filepaths and variable definition ###########################################
dir.prj = '' # project directory
nam.oat = '' # scenario name
yr.ibm = 50 # years to run simulation
ts.ibm = 0.195 # gC/m2, intial TS biomass for IBM (if NA, uses CASM's value)
rns.ibm = c(6, 35, 10, 88, 58, 87, 36,  2, 79, 77)
var.oat = 'Te_Th' # 'Te_Th', 'N', 'P', 'TIS', 'depth', 'I0'
chg.oat = c(-3.6, -1.8, 1.8, 3.6) # changes

## needed packages #############################################################
library(RNetLogo)

## Initial Processing ##########################################################
# Create Needed Folders
dir.create(file.path(dir.prj, nam.oat)) # create scenario folder
setwd(file.path(dir.prj, nam.oat)) # change directory

# create meta data file
temp = c('Meta Data for TS OAAT Model Run:',
  paste('Scenario Name (nam.oat):', nam.oat),
  paste('Created:', Sys.time()),
  paste('Project Directory (dir.prj):', dir.prj),
  paste('Years to run Scenario (yr.ibm):', yr.ibm),
  paste('Initial TS Biomass Overwrite (ts.ibm):', ts.ibm),
  paste('Random Number Seeds (rns.ibm):', paste(rns.ibm, collapse = ', ')),
  paste('Variable to Change (var.oat):', var.oat),
  paste('Value Changes (chg.oat):', paste(chg.oat, collapse = ', ')))
writeLines(temp, 'meta_scenario.txt')
rm(temp)

## Set Up Simulation ###########################################################
# folders for single run
dir.create('CASM')
setwd('CASM')
dir.create('OC_CASM_Out')

# copy needed files
for(temp in list.files(file.path(dir.prj, 'CASM_in'))){
  file.copy(file.path(dir.prj, 'CASM_in', temp), temp)
}
rm(temp)

system('casm_ibm_ocv11.exe') # Run CASM
setwd('..')

## Prep Raw Files ##############################################################
dir.create('OAAT_inp')
setwd('OAAT_inp')

# copy needed files
for(temp in list.files(file.path(dir.prj, 'CASM_in'))){
  file.copy(file.path(dir.prj, 'CASM_in', temp), temp)
}
for(temp in list.files(file.path(dir.prj, 'IBM_in'))){
  file.copy(file.path(dir.prj, 'IBM_in', temp), temp)
}
file.copy(file.path(dir.prj, nam.oat, 'CASM', 'IBM_TS_master_ref.out'),
  'IBM_TS_master_ref.out')

# alter control file
temp = readLines('casm_IBM_OC_control.dat')
temp[11] = '.\\OAAT_env.prn'
temp[77] = 1 # use IBM model
temp[87] = 1 # run 1 year
writeLines(temp, 'casm_IBM_OC_control.dat')
rm(temp)

## get day 1 files
temp = readLines('IBM_TS_master_ref.out')
writeLines(
  c('IBM-CASM transfer values:', temp[2],
    substr(temp[3], 13, 1218), substr(temp[4], 13, 1218)), 
  'IBM_TS_transfer.out')
writeLines(
  c('IBM-CASM transfer environmental values:', 
    'Input values initial environmental factors (mg/L):',
    paste0(substr(temp[3], 13, 18), substr(temp[3], 1219, nchar(temp[3]))),
    paste0(substr(temp[4], 13, 18), 
      substr(temp[4], 1219, nchar(temp[3])))),
  'IBM_TS_transfer_env.out')
rm(temp)
file.remove('IBM_TS_master_ref.out')
setwd('..')

## Run OAAT Scenarios ##########################################################
for(j in c(NA, chg.oat)){
  ## Scenario set up
  run.fp = file.path(dir.prj, nam.oat, ifelse(is.na(j), 'base', j))
  dir.create(run.fp) # create run folder
  setwd(run.fp) # set working directory
  dir.create('OC_CASM_Out') # extra output folder
  
  ## Copy files over
  for(temp in list.files(file.path(dir.prj, nam.oat, 'OAAT_inp'))){
    file.copy(file.path(dir.prj, nam.oat, 'OAAT_inp', temp), temp)}
  
  ## Alter env inputs
  if(is.na(j)){
    file.rename('casmOC_Env_Data_11Aug2019.prn', 'OAAT_env.prn')
  } else {
    temp = read.table('casmOC_Env_Data_11Aug2019.prn', skip = 7, 
      col.names = c('DOY', 'Te_C', 'Th_C', 'I0_eins.m2', 'N_mg.L', 'P_mg.L',
        'Si_mg.L', 'depth_m', 'velocit_m.s', 'TIS_mg.L', 'POCe_mg.L', 'wind_m.s',
        'salinity_psu'))
    
    if(var.oat == 'Te_Th'){
      temp$Te_C = temp$Te_C + j
      temp$Te_C[temp$Te_C < 0] = 0
      temp$Th_C = temp$Th_C + j
      temp$Th_C[temp$Th_C < 0] = 0
    } else if (var.oat == 'N'){
      temp$N_mg.L = temp$N_mg.L + j
      temp$N_mg.L[temp$N_mg.L < 0] = 0
    } else if (var.oat == 'P'){
      temp$P_mg.L = temp$P_mg.L + j
      temp$P_mg.L[temp$P_mg.L < 0] = 0
    } else if (var.oat == 'TIS'){
      temp$TIS_mg.L = temp$TIS_mg.L + j
      temp$TIS_mg.L[temp$TIS_mg.L < 0] = 0
    } else if (oaat_var == 'depth'){
      temp$depth_m = j
    } else if (var.oat == 'I0'){
      temp$I0_eins.m2 = temp$I0_eins.m2 * j
    }
  
  # format & write inputs
  temp = c('CASM-OC Iowa    oxbow:  2010.000',
    'fktrN  1.0000', 'fktrP  1.0000', 'fktrS  1.0000',
    paste0('                   Te      Th      I0         N       P       Si',
      '      depth   velocit TIS     POCe    wind    salinity'),
    paste0('             DOY   C       C       eins/m2    mg/L    mg/L',
      '    mg/L    m       m/s     mg/L    mg/L    m/s     psu'),
    paste0(rep('-', 116), collapse = ''),
    paste0(
      formatC(temp$DOY, width = 16),
      formatC(sprintf('%.3f', round(temp$Te_C, 3)), width = 8),
      formatC(sprintf('%.3f', round(temp$Th_C, 3)), width = 8),
      formatC(sprintf('%.3f', round(temp$I0_eins.m2, 3)), width = 8),
      formatC(sprintf('%.3f', round(temp$N_mg.L, 3)), width = 11),
      formatC(sprintf('%.3f', round(temp$P_mg.L, 3)), width = 8),
      formatC(sprintf('%.3f', round(temp$Si_mg.L, 3)), width = 8),
      formatC(sprintf('%.3f', round(temp$depth_m, 3)), width = 8),
      formatC(sprintf('%.3f', round(temp$velocit_m.s, 3)), width = 8),
      formatC(sprintf('%.3f', round(temp$TIS_mg.L, 3)), width = 8),
      formatC(sprintf('%.3f', round(temp$POCe_mg.L, 3)), width = 8),
      formatC(sprintf('%.3f', round(temp$wind_m.s, 3)), width = 8),
      formatC(sprintf('%.3f', round(temp$salinity_psu, 3)), width = 8)))
  writeLines(temp, 'OAAT_env.prn')
  file.remove('casmOC_Env_Data_11Aug2019.prn')
  rm(temp)
  }
  
  ## Run CASM
  lns.bio = readLines('IBM_TS_transfer.out')
  lns.bio[3] = paste0(' Yr  AccDay ', substr(lns.bio[3], 1, 1124))
  lns.bio[4] = paste0('   1      1 ', substr(lns.bio[4], 1, 1124))
  substr(lns.bio[4], 323, 336) = 
    formatC(ts.ibm, width = 12, format = 'E', digits = 8)
  
  for(i in 2:365){
    system('casm_ibm_ocv11.exe')
    temp = readLines('IBM_TS_transfer.out')[4]
    temp = paste('   1   ',  formatC(i, width = 3), substr(temp, 1, 1124))
    lns.bio = c(lns.bio, temp)
    rm(temp)
  }
  writeLines(lns.bio, 'IBM_TS_master_ref_LHC.out')
  rm(i, lns.bio)
  
  for(rs in rns.ibm){
    setwd('C:/Program Files/NetLogo 6.0.3/app')
    NLStart('C:/Program Files/NetLogo 6.0.3/app',
      nl.jarname = 'netlogo-6.0.3.jar', gui = FALSE)
    NLLoadModel(file.path(run.fp, 'TS_IBM_Effects_Sep2019_v2.nlogo'))
    
    NLCommand('set InputParameters "InputParameters_casmOC_27Feb2020.txt"')
    NLCommand('set WaterConditionInput "OAAT_env.prn"')
    NLCommand('set DailyPreyBiomass "IBM_TS_master_ref_LHC.out"')
    NLCommand(paste0('set BiomassOutput "output_', rs, '.txt"'))
    NLCommand('set DailyIO false')
    NLCommand(paste('set RandomNumberSeed', rs))
    NLCommand('set pondArea 100')
    NLCommand('set CASM_pool_area 1125')
    NLCommand('set Sunfish true')
    NLCommand('set Predation false')
    NLCommand('set DD_egg_larva true')
    NLCommand('set AverageKmin 0.68')
    NLCommand('set ScalingSearchArea 1.0')
    NLCommand('set MaxDetritusDepth 1.5')
    NLCommand('set YearsToRun 40')
    NLCommand('set EffectsModule false')
    
    ## Run every day
    setwd(run.fp)
    NLCommand('setup')
    for(i in 1:(365*yr.ibm)) NLCommand('go')
    NLCommand('write_output')
    NLQuit()
    rm(i)
  }
  rm(rs, run.fp)
}
rm(j)

## Clean up environment ########################################################
rm(dir.prj, nam.oat, yr.ibm, ts.ibm, var.oat, chg.oat, rns.ibm)
