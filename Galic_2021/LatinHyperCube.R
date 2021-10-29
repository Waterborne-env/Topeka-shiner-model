## Topeka Shiner Latin Hyper Cube Runs #########################################
## Author: CR
## Last Edited: 2021-09-20
################################################################################
rm(list = ls())

## filepaths and variable definition ###########################################
dir.prj = '' # project directory
nam.lhc = '' # scenario name
yr.ibm = 50 # years to run simulation
ts.ibm = 0.195 # gC/m2, intial TS biomass for IBM (if NA, uses CASM's value)
cnt.lhc = 1250 # file path to csv OR number of LHC runs

## needed packages #############################################################
library(RNetLogo)
library(lhs)

## Initial Processing ##########################################################
# Create Needed Folders
dir.create(file.path(dir.prj, nam.lhc)) # create scenario folder
setwd(file.path(dir.prj, nam.lhc)) # change directory

# create meta data file
temp = c('Meta Data for TS LHC Model Run:',
  paste('Scenario Name (nam.lhc):', nam.lhc),
  paste('Created:', Sys.time()),
  paste('Project Directory (dir.prj):', dir.prj),
  paste('Years to run Scenario (yr.ibm):', yr.ibm),
  paste('Initial TS Biomass Overwrite (ts.ibm):', ts.ibm))
if(class(cnt.lhc) == 'numeric'){
  temp = c(temp, paste('LHC Count (cnt.lhc):', cnt.lhc))
} else {
  temp = c(temp, paste('LHC Source (cnt.lhc):', cnt.lhc))
}
writeLines(temp, 'meta_scenario.txt')
rm(temp)

## LHC Data ####################################################################
if(class(cnt.lhc) == 'character'){
  df.lhc = read.csv(cnt.lhc)
} else { # create LHC
  # add baseline
  df.lhc = data.frame(Te_Th = 0, N = 0, P = 0, TIS = 0, depth = 1.2, I0 = 1)
  rownames(df.lhc) = 'base'
  
  while(nrow(df.lhc) < cnt.lhc + 1){ # ensure all unique
    temp = randomLHS(cnt.lhc - nrow(df.lhc) + 1, 6)
    temp = as.data.frame(temp)
    colnames(temp) = c('Te_Th', 'N', 'P', 'TIS', 'depth', 'I0')
    
    # first 5 have 5 options
    temp[1:5][temp[1:5] >= 0.0 & temp[1:5] < 0.2] = 11
    temp[1:5][temp[1:5] >= 0.2 & temp[1:5] < 0.4] = 12
    temp[1:5][temp[1:5] >= 0.4 & temp[1:5] < 0.6] = 13
    temp[1:5][temp[1:5] >= 0.6 & temp[1:5] < 0.8] = 14
    temp[1:5][temp[1:5] >= 0.8 & temp[1:5] < 1.0] = 15
    
    # last 1 has 4 options
    temp[6][temp[6] >= 0.00 & temp[6] < 0.25] = 11
    temp[6][temp[6] >= 0.25 & temp[6] < 0.50] = 12
    temp[6][temp[6] >= 0.50 & temp[6] < 0.75] = 13
    temp[6][temp[6] >= 0.75 & temp[6] < 1.00] = 14
    
    
    ## reassign values
    temp$Te_Th[temp$Te_Th == 11] = -3.6
    temp$Te_Th[temp$Te_Th == 12] = -1.8
    temp$Te_Th[temp$Te_Th == 13] =  0
    temp$Te_Th[temp$Te_Th == 14] =  1.8
    temp$Te_Th[temp$Te_Th == 15] =  3.6
    
    temp$N[temp$N == 11] = -1.2
    temp$N[temp$N == 12] = -0.6
    temp$N[temp$N == 13] =  0
    temp$N[temp$N == 14] =  0.6
    temp$N[temp$N == 15] =  1.2
    
    temp$P[temp$P == 11] = -1
    temp$P[temp$P == 12] = -0.5
    temp$P[temp$P == 13] =  0
    temp$P[temp$P == 14] =  0.5
    temp$P[temp$P == 15] =  1
    
    temp$TIS[temp$TIS == 11] = -7
    temp$TIS[temp$TIS == 12] = -3.5
    temp$TIS[temp$TIS == 13] =  0
    temp$TIS[temp$TIS == 14] =  3.5
    temp$TIS[temp$TIS == 15] =  7
    
    temp$depth[temp$depth == 11] = 0.1
    temp$depth[temp$depth == 12] = 0.4
    temp$depth[temp$depth == 13] = 0.8
    temp$depth[temp$depth == 14] = 1.2
    temp$depth[temp$depth == 15] = 1.6
    
    temp$I0[temp$I0 == 11] = 1
    temp$I0[temp$I0 == 12] = 0.9
    temp$I0[temp$I0 == 13] = 0.75
    temp$I0[temp$I0 == 14] = 0.5
    
    df.lhc = rbind(df.lhc, temp) # add new rows
    df.lhc = unique(df.lhc) # make sure each run unique
    rm(temp)
  }
  
  df.lhc$run = c('base', paste0('r', 1:cnt.lhc)) # add run names
  df.lhc = df.lhc[c('run', 'Te_Th', 'N', 'P', 'TIS', 'depth', 'I0')] # reorder
}

write.csv(df.lhc, 'lhc.csv', row.names = FALSE) # write file

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
dir.create('LHC_inp')
setwd('LHC_inp')

# copy needed files
for(temp in list.files(file.path(dir.prj, 'CASM_in'))){
  file.copy(file.path(dir.prj, 'CASM_in', temp), temp)
}
for(temp in list.files(file.path(dir.prj, 'IBM_in'))){
  file.copy(file.path(dir.prj, 'IBM_in', temp), temp)
}
file.copy(file.path(dir.prj, nam.lhc, 'CASM', 'IBM_TS_master_ref.out'),
  'IBM_TS_master_ref.out')

# alter control file
temp = readLines('casm_IBM_OC_control.dat')
temp[11] = '.\\LHC_env.prn'
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

## Run LHC Scenarios ###########################################################
for(j in 1:nrow(df.lhc)){
  ## Scenario set up
  run.fp = file.path(dir.prj, nam.lhc, df.lhc$run[j])
  dir.create(run.fp) # create run folder
  setwd(run.fp) # set working directory
  dir.create('OC_CASM_Out') # extra output folder
  
  ## Copy files over
  for(temp in list.files(file.path(dir.prj, nam.lhc, 'LHC_inp'))){
    file.copy(file.path(dir.prj, nam.lhc, 'LHC_inp', temp), temp)}
  
  ## Alter env inputs
  temp = read.table('casmOC_Env_Data_11Aug2019.prn', skip = 7, 
    col.names = c('DOY', 'Te_C', 'Th_C', 'I0_eins.m2', 'N_mg.L', 'P_mg.L',
      'Si_mg.L', 'depth_m', 'velocit_m.s', 'TIS_mg.L', 'POCe_mg.L', 'wind_m.s',
      'salinity_psu'))
  
  # alter inputs from LHC data
  temp$Te_C = temp$Te_C + df.lhc$Te_Th[j]
  temp$Te_C[temp$Te_C < 0] = 0
  temp$Th_C = temp$Th_C + df.lhc$Te_Th[j]
  temp$Th_C[temp$Th_C < 0] = 0
  temp$N_mg.L = temp$N_mg.L + df.lhc$N[j]
  temp$N_mg.L[temp$N_mg.L < 0] = 0
  temp$P_mg.L = temp$P_mg.L + df.lhc$P[j]
  temp$P_mg.L[temp$P_mg.L < 0] = 0
  temp$TIS_mg.L = temp$TIS_mg.L + df.lhc$TIS[j]
  temp$TIS_mg.L[temp$TIS_mg.L < 0] = 0
  temp$depth_m = temp$depth[j]
  temp$I0_eins.m2 = temp$I0_eins.m2 * df.lhc$I0[j]
  
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
  writeLines(temp, 'LHC_env.prn')
  file.remove('casmOC_Env_Data_11Aug2019.prn')
  rm(temp)
  
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
  
  ## Run IBM
  setwd('C:/Program Files/NetLogo 6.0.3/app')
  NLStart('C:/Program Files/NetLogo 6.0.3/app',
    nl.jarname = 'netlogo-6.0.3.jar', gui = FALSE)
  NLLoadModel(file.path(run.fp, 'TS_IBM_Effects_Sep2019_v2.nlogo'))
  
  NLCommand('set InputParameters "InputParameters_casmOC_27Feb2020.txt"')
  NLCommand('set WaterConditionInput "LHC_env.prn"')
  NLCommand('set DailyPreyBiomass "IBM_TS_master_ref_LHC.out"')
  NLCommand('set BiomassOutput "TS-IBM_output.txt"')
  NLCommand('set DailyIO false')
  NLCommand('set RandomNumberSeed 6')
  NLCommand('set pondArea 100')
  NLCommand('set CASM_pool_area 1125')
  NLCommand('set Sunfish true')
  NLCommand('set Predation false')
  NLCommand('set DD_egg_larva true')
  NLCommand('set AverageKmin 0.68')
  NLCommand('set ScalingSearchArea 1.0')
  NLCommand('set MaxDetritusDepth 1.5')
  NLCommand('set YearsToRun 50')
  NLCommand('set EffectsModule false')

  setwd(run.fp)
  NLCommand('setup')
  for(i in 1:(365*yr.ibm)) NLCommand('go')
  NLCommand('write_output')
  NLQuit()
  
  rm(i, run.fp)
}
rm(j)

## Clean up environment ########################################################
rm(dir.prj, nam.lhc, cnt.lhc, yr.ibm, ts.ibm, df.lhc)
