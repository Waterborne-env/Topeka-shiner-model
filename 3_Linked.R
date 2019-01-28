## Topeka Shiner Linked Set Up & Run ###########################################
## Author: CR
## Last Edited: 2019-01-22
################################################################################
rm(list = ls())

## filepaths and variable definition ###########################################
# local locations/folder names
dir.prj  = '' # project directory
dir.run  = '' # scenario name

# NetLogo information
dir.nl = 'C:/Program Files/NetLogo 6.0.2/app' # NetLogo
typ.nl = 'netlogo-6.0.2.jar' # NetLogo version

# IBM locations/folder names (all relative to dir.ibm)
dir.ibm = ''
dir.mod = 'TS_Model_NoPesticide_27Nov2018.nlogo' # TS model
dir.inp = 'InputParameters_20Jun2018.txt' # TS input file

# CASM scenario/file names
dir.casm = '' # CASM base scenario
dir.exe = 'casm_ibm_tsv42.exe'# CASM exe
dir.cntl = 'casm_IBM_TS_control.dat' # CASM control file
dir.hab = 'casm_TS_bio_parms_RSK_WithTS_08Oct2018.txt' # CASM habitat file
dir.web = 'web_casmTS_HRM_20Jun2018.dat' # CASM web file
dir.env = 'env_casmTS_2010_08May2018.prn' # CASM environmental file
dir.con = 'ATZ_exp_zero_test.dat' # CASM contaminant file
dir.tox = 'TS_toxicity_atrazine_tri_base.dat' # CASM toxicity file

# Model Parameters
casm.yrs =  # year of day 1 to pull from CASM

# Model Parameters - on interface of IBM
sim.RandomNumberSeed = 6 # random seeds, from sample.int(1000, 20)
sim.pondArea = 100 # m2
sim.CASM_pool_area = 1125 # m2
sim.Sunfish = TRUE
sim.Predation = FALSE
sim.DD_egg_larva = TRUE
sim.AverageKmin = 0.68
sim.ScalingSearchArea = 1.0
sim.MaxDetritusDepth = 1.5 # cm
sim.YearsToRun = 15

## needed packages #############################################################
library('RNetLogo')

## processing ##################################################################
setwd(dir.prj) # change directory

# create overall scenario name and folder
dir.create(dir.run)
setwd(dir.run)

# create meta data file
writeLines(
  as.character(c('Meta Data for TS Linked Model Run:',
    paste('Scenario Name (dir.nam):', dir.run),
    paste('Scenario Location (dir.run):', dir.run),
    paste('Created:', Sys.time()),
    paste('Project Directory (dir.prj):', dir.prj),
    
    paste('NetLogo Directory (dir.nl):', dir.nl),
    paste('NetLogo Version (typ.nl):', typ.nl),
    paste('IBM File Locations (dir.ibm.base):', dir.ibm),
    paste('Model Directory (dir.mod):', dir.mod),
    paste('Input File Location (dir.inp):', dir.inp),
    
    paste('CASM File Locations (dir.casm):', dir.casm),
    paste('CASM Model Directory (dir.exe):', dir.exe),
    paste('CASM Control File Location (dir.cntl):', dir.cntl),
    paste('CASM Habitat File Location (dir.hab):', dir.hab),
    paste('CASM Web File Location (dir.web):', dir.web),
    paste('CASM Environmental File Location (dir.env):', dir.env),
    paste('CASM Contaminant File Location (dir.con):', dir.con),
    paste('CASM Toxicity File Location (dir.tox):', dir.tox),
    paste('CASM Data Year used (casm.yrs):', casm.yrs),
    
    paste('Random Num Seeds (sim.RandomNumberSeed):', sim.RandomNumberSeed),
    paste('Pond Area (sim.pondArea):', sim.pondArea, 'm2'),
    paste('CASM Pond Area (sim.CASM_pool_area):', sim.CASM_pool_area, 'm2'),
    paste('Sunfish Simulation (sim.Sunfish):', sim.Sunfish),
    paste('LMB Predation (sim.Predation):', sim.Predation),
    paste('Eggs/Larvae Dens. Dependence (sim.DD_egg_larva):', sim.DD_egg_larva),
    paste('Average Kmin (sim.AverageKmin):', sim.AverageKmin),
    paste('Search Area (sim.ScalingSearchArea):', sim.ScalingSearchArea),
    paste('Detritus Depth (sim.MaxDetritusDepth):', sim.MaxDetritusDepth, 'cm'),
    paste('Years to Run IBM (sim.YearsToRun):', sim.YearsToRun)
  )), 'meta_scenario.txt')

## Set Up CASM Part of Simulation ##############################################
dir.create('TS_CASM_Out') # temp output folder
dir.create(file.path('TS_CASM_Out', 'plots')) # temp output folder

# copy needed files
for(temp in c(dir.exe, dir.cntl, dir.hab, dir.web, dir.env, dir.con, dir.tox,
  'IBM_TS_master_ref.out')){
  file.copy(file.path(dir.casm, temp), basename(temp))
}
rm(temp)

## change file paths for dir.cntl and other model parameters
temp = readLines(basename(dir.cntl))
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
temp[77] = 1 # don't use IBM model
temp[79] = 0 # 0 - ref, 1 - eff
temp[81] = '.\\IBM_TS_transfer.out'
temp[82] = '.\\IBM_TS_transfer_env.out'
temp[84] = '.\\IBM_TS_master_ref.out'
temp[85] = '.\\IBM_TS_master_eff.out'
temp[87] = 1
writeLines(temp, basename(dir.cntl))
rm(temp)

## Set Up IBM Part of Simulation ###############################################
# copy needed files
for(temp in c(dir.mod, dir.inp)){
  file.copy(file.path(dir.ibm, temp), basename(temp))
}
rm(temp)

# pull CASM day 1 year X for simulation
temp = readLines('IBM_TS_master_ref.out')
temp = c('IBM-CASM transfer values:', 
  temp[2],
  substr(temp[3], 13, 1136), 
  substr(temp[(4 + 365 * (casm.yrs - 1))], 13, 1136))
writeLines(temp, 'IBM_TS_transfer.out') # to keep day 1
rm(temp)

## copy over & format env transfer
temp = readLines('IBM_TS_master_ref.out')
temp = c('IBM-CASM transfer environmental values:', 
  'Input values initial environmental factors (mg/L):',
  paste0(substr(temp[3], 13, 18), substr(temp[3], 1139, nchar(temp[3]))),
  paste0(substr(temp[(4 + 365 * (casm.yrs - 1))], 13, 18), 
    substr(temp[(4 + 365 * (casm.yrs - 1))], 1139, nchar(temp[3]))))
writeLines(temp, 'IBM_TS_transfer_env.out')
rm(temp)

## Run TS Linked ###############################################################
# load netlogo model
setwd(dir.nl)
NLStart(dir.nl, nl.jarname = typ.nl, gui = FALSE)
NLLoadModel(paste(dir.prj, dir.run, dir.mod, sep = '/'))

# netlogo interface inputs
NLCommand(paste0('set InputParameters "', dir.inp, '"'))
NLCommand(paste0('set WaterConditionInput "', dir.env, '"'))
NLCommand('set DailyPreyBiomass "IBM_TS_transfer.out"')
NLCommand('set BiomassOutput "TS-IBM_output.txt"')
NLCommand('set DailyIO true')

NLCommand(paste0('set RandomNumberSeed ', sim.RandomNumberSeed))
NLCommand(paste0('set pondArea ', sim.pondArea))
NLCommand(paste0('set CASM_pool_area ', sim.CASM_pool_area))
NLCommand(paste0('set Sunfish ', tolower(sim.Sunfish)))
NLCommand(paste0('set Predation ', tolower(sim.Predation)))
NLCommand(paste0('set DD_egg_larva ', tolower(sim.DD_egg_larva)))
NLCommand(paste0('set AverageKmin ', sim.AverageKmin))
NLCommand(paste0('set ScalingSearchArea ', sim.ScalingSearchArea))
NLCommand(paste0('set MaxDetritusDepth ', sim.MaxDetritusDepth))
NLCommand(paste0('set YearsToRun ', sim.YearsToRun))

setwd(file.path(dir.prj, dir.run)) # change directory

NLCommand('setup') # set up run

# run each day
for(i in 1:365){
  if(i > 1) system(dir.exe) # run CASM
  # run IBM
  NLCommand('go') # run TS for 1 day
  NLCommand('write_output') # write TS output
  
  ## rewrite transfer file for next day
  tmp.casm = read.table('IBM_TS_transfer.out', header = FALSE, skip = 3)
  tmp.ibm = read.table(paste0('out_rns', sim.rns, '.txt'), 
    header = FALSE, skip = 1)
  if(tmp.casm[1,1] != tmp.ibm[nrow(tmp.ibm),1]){
    print('Days Do Not Match for Transfer!')
    break
  }
  tmp.casm[1,24] = tmp.ibm[nrow(tmp.ibm), 2] # Copepods
  tmp.casm[1,25] = tmp.ibm[nrow(tmp.ibm), 3] # Cladocerans
  tmp.casm[1,26] = tmp.ibm[nrow(tmp.ibm), 4] # Rotifers
  tmp.casm[1,27] = tmp.ibm[nrow(tmp.ibm), 5] # MicroZplankton
  tmp.casm[1,36] = tmp.ibm[nrow(tmp.ibm), 6] # Ephemeroptera
  tmp.casm[1,37] = tmp.ibm[nrow(tmp.ibm), 7] # Trichoptera
  tmp.casm[1,38] = tmp.ibm[nrow(tmp.ibm), 8] # Oligochaetes
  tmp.casm[1,39] = tmp.ibm[nrow(tmp.ibm), 9] # Chironomids
  tmp.casm[1,49] = tmp.ibm[nrow(tmp.ibm),10] # SedimentPOC
  tmp.casm[1,51] = tmp.ibm[nrow(tmp.ibm),11] # TotalPeriphyt
  tmp.casm[1,31] = tmp.ibm[nrow(tmp.ibm),12] # TopShiner(juv)
  tmp.casm[1,32] = tmp.ibm[nrow(tmp.ibm),13] # TopShiner(adt)
  tmp.casm[1,62] = tmp.ibm[nrow(tmp.ibm),15] # numTopShiner(juv)
  tmp.casm[1,63] = tmp.ibm[nrow(tmp.ibm),16] # numTopShiner(adt)
  tmp.casm[1,70] = tmp.ibm[nrow(tmp.ibm),18] # TSnum-eaten
  tmp.casm[1,71] = tmp.ibm[nrow(tmp.ibm),19] # TSgww-eaten
  rm(tmp.ibm)
  tmp.casm.2 = paste0(formatC(tmp.casm[1,1], width = 4),
    formatC(tmp.casm[1,2], format = 'E', width = 12))
  for(j in 3:ncol(tmp.casm)){
    tmp.casm.2 = paste0(tmp.casm.2,
      formatC(tmp.casm[1,j], format = 'E', width = 16))
  }
  rm(tmp.casm, j)
  tmp.casm.2 = c(readLines('IBM_TS_transfer.out', n = 3), tmp.casm.2)
  writeLines(tmp.casm.2, 'IBM_TS_transfer.out')
  rm(tmp.casm.2)
  
  file.remove(paste0('out_rns', sim.rns, '.txt'))
  
}
rm(i)

NLCommand('write_output') # write TS output

NLQuit() # quit session

## Clean Up environment ########################################################
rm(dir.prj, dir.run, dir.nl, typ.nl, dir.ibm, dir.mod, dir.inp, dir.casm,
  dir.env, casm.yrs, sim.RandomNumberSeed, sim.pondArea, sim.CASM_pool_area,
  sim.Sunfish, sim.Predation, sim.DD_egg_larva, sim.AverageKmin,
  sim.ScalingSearchArea, sim.MaxDetritusDepth, sim.YearsToRun)
