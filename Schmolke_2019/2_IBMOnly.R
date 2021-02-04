## Topeka Shiner IBM Only Set Up & Run #########################################
## Author: CR
## Last Edited: 2019-01-22
################################################################################
rm(list = ls())

## filepaths and variable definition ###########################################
# Initial Folder Locations/Names
dir.prj = 'D:\\Topeka-shiner-model-master' # project directory
dir.run = 'IBMOnlyRun' # scenario name

# NetLogo Location/Files
dir.nl = 'C:/Program Files/NetLogo 6.0.2/app' # NetLogo
typ.nl = 'netlogo-6.0.2.jar' # NetLogo version

# IBM Files (relative to dir.prj)
dir.mod = 'TS_IBM_March2019.nlogo' # TS model
dir.inp = 'InputParameters_20Jun2018.txt' # TS input file

# CASM Files (relative to dir.prj)
dir.casm = 'CASMRun' # CASM scenario name from 1_CASM.R
dir.env = 'env_casmTS_2010_08May2018.prn' # CASM environmental file

# CASM Model Parameters
casm.yr = 1 # year of data to pull from CASM

# IBM Model Parameters
sim.RandomNumberSeed = 6 # random seed
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
dir.create(dir.run) # create scenario folder

# create meta data file
writeLines(
  as.character(c('Meta Data for TS IBM Only Model Run:',
    paste('Scenario Name (dir.run):', dir.run),
    paste('Created:', Sys.time()),
    paste('Project Directory (dir.prj):', dir.prj),
    
    paste('NetLogo Directory (dir.nl):', dir.nl),
    paste('NetLogo Version (typ.nl):', typ.nl),
    
    paste('Model Directory (dir.mod):', dir.mod),
    paste('Input File Location (dir.inp):', dir.inp),
    
    paste('CASM Scenario Folder (dir.casm):', dir.casm),
    paste('CASM Environmental File Location (dir.env):', dir.env),
    paste('CASM Data Year used (casm.yr):', casm.yr),
    
    paste('Random Num Seed (sim.RandomNumberSeed):', sim.RandomNumberSeed),
    paste('Pond Area (sim.pondArea):', sim.pondArea, 'm2'),
    paste('CASM Pond Area (sim.CASM_pool_area):', sim.CASM_pool_area, 'm2'),
    paste('Sunfish Simulation (sim.Sunfish):', sim.Sunfish),
    paste('LMB Predation (sim.Predation):', sim.Predation),
    paste('Eggs/Larvae Dens. Dependence (sim.DD_egg_larva):', sim.DD_egg_larva),
    paste('Average Kmin (sim.AverageKmin):', sim.AverageKmin),
    paste('Search Area (sim.ScalingSearchArea):', sim.ScalingSearchArea),
    paste('Detritus Depth (sim.MaxDetritusDepth):', sim.MaxDetritusDepth, 'cm'),
    paste('Years to Run IBM (sim.YearsToRun):', sim.YearsToRun)
  )), file.path(dir.run, 'meta_scenario.txt'))

## Set Up Simulation ###########################################################
# copy needed files
for(temp in c(dir.mod, dir.inp, 
  file.path(dir.casm, c(dir.env, 'IBM_TS_master_ref.out')))){
  file.copy(temp, file.path(dir.run, basename(temp)))
}
rm(temp)

# pull CASM year for simulation
temp = readLines(file.path(dir.run, 'IBM_TS_master_ref.out'))
temp = c('IBM-CASM transfer values:', 
  temp[2],
  substr(temp[3], 1, 1136), 
  substr(temp[(4 + 365 * (casm.yr - 1)):(3 + 365 * (casm.yr))], 1, 1136))
writeLines(temp, file.path(dir.run, 'IBM_TS_master_ref.out'))
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
rm(i)

NLCommand('write_output') # write TS output

NLQuit() # quit session

## Clean Up environment ########################################################
rm(dir.prj, dir.run, dir.nl, typ.nl, dir.mod, dir.inp, dir.casm, dir.env,
  casm.yr, sim.RandomNumberSeed, sim.pondArea, sim.CASM_pool_area,
  sim.Sunfish, sim.Predation, sim.DD_egg_larva, sim.AverageKmin,
  sim.ScalingSearchArea, sim.MaxDetritusDepth, sim.YearsToRun)
  
