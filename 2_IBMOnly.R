## Topeka Shiner IBM Only Set Up & Run #########################################
## Author: CR
## Last Edited: 2019-01-22
################################################################################
rm(list = ls())

## filepaths and variable definition ###########################################
# local locations/folder names
dir.prj = '' # project directory
dir.run = '' # scenario name

# NetLogo information
dir.nl = 'C:/Program Files/NetLogo 6.0.2/app' # NetLogo
typ.nl = 'netlogo-6.0.2.jar' # NetLogo version

# IBM locations/folder names (all relative to dir.ibm)
dir.ibm = ''
dir.mod = 'TS_Model_NoPesticide_27Nov2018.nlogo' # TS model
dir.inp = 'InputParameters_20Jun2018.txt' # TS input file

# CASM scenario/file names (all relative to dir.casm)
dir.casm = '' # CASM base scenario
dir.env = 'env_casmTS_2010_08May2018.prn' # CASM environmental file
casm.yrs = 1 # year of data to pull from CASM

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
  as.character(c('Meta Data for TS IBM Only Model Run:',
    paste('Scenario Name (dir.run):', dir.run),
    paste('Created:', Sys.time()),
    paste('Project Directory (dir.prj):', dir.prj),
    
    paste('NetLogo Directory (dir.nl):', dir.nl),
    paste('NetLogo Version (typ.nl):', typ.nl),
    paste('IBM File Locations (dir.ibm.base):', dir.ibm),
    paste('Model Directory (dir.mod):', dir.mod),
    paste('Input File Location (dir.inp):', dir.inp),
    
    paste('CASM Scenario Folder (dir.casm):', dir.casm),
    paste('CASM Environmental File Location (dir.env):', dir.env),
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

## Set Up Simulation ###########################################################
# copy needed files
for(temp in c(file.path(dir.ibm, c(dir.mod, dir.inp)), 
  file.path(dir.prj, dir.casm, c(dir.env, 'IBM_TS_master_ref.out')))){
  file.copy(temp, basename(temp))
}
rm(temp)

# pull CASM year for simulation
temp = readLines('IBM_TS_master_ref.out')
temp = c('IBM-CASM transfer values:', 
  temp[2],
  substr(temp[3], 1, 1136), 
  substr(temp[(4 + 365 * (casm.yrs - 1)):(3 + 365 * (casm.yrs))], 1, 1136))
writeLines(temp, 'IBM_TS_master_ref.out')
rm(temp)

## Run TS IBM Stand Alone ######################################################
# load netlogo model
setwd(dir.nl)
NLStart(dir.nl, nl.jarname = typ.nl, gui = FALSE)
NLLoadModel(paste(dir.prj, dir.run, dir.mod, sep = '/'))

# netlogo interface inputs
NLCommand(paste0('set InputParameters "', dir.inp, '"'))
NLCommand(paste0('set WaterConditionInput "', dir.env, '"'))
NLCommand('set DailyPreyBiomass "IBM_TS_master_ref.out"')
NLCommand('set BiomassOutput "TS-IBM_output.txt"')
NLCommand('set DailyIO false')

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
for(i in 1:(365*sim.YearsToRun)) NLCommand('go') # run TS
NLCommand('write_output') # write TS output

NLQuit() # quit session

## Clean Up environment ########################################################
rm(dir.prj, dir.run, dir.nl, typ.nl, dir.ibm, dir.mod, dir.inp, dir.casm,
  dir.env, casm.yrs, sim.RandomNumberSeed, sim.pondArea, sim.CASM_pool_area,
  sim.Sunfish, sim.Predation, sim.DD_egg_larva, sim.AverageKmin,
  sim.ScalingSearchArea, sim.MaxDetritusDepth, sim.YearsToRun)
  