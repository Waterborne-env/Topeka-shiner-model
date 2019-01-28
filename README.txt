## Topeka Shiner Hybrid Model: How To ##########################################
## Author: Colleen Roy
## Last Edited: 2019-01-25
################################################################################

## Disclaimer ##################################################################
The software and associated files uploaded in this repository were used to 
generate the results to be published in: 
Schmolke A, Bartell SM, Roy C, Green N, Galic N, Brain R. Species-specific 
population dynamics and their link to an aquatic food web: a hybrid modeling 
approach.
The manuscript has been submitted for review to a scientific journal. 

This software and associated files are provided "as is" with the sole purpose 
to allow the reproduction of the published results without any warranties of 
performance or fitness for any other purpose. 

## Overview ####################################################################
This how to documents how to run CASM, the TS-IBM, and the linked CASM + TS-IBM
in R.

To run either the TS-IBM OR the linked CASM + TS-IBM you must run CASM first

## Run CASM ####################################################################
File: 1_CASM.R

Change the following lines:
	Line 9: Directory of project (must already exist)
	Line 10: Name of scenario
		Note: this will become a folder in the project directory
		
	Line 13: Folder where all CASM files are located
	Line 14: CASM executable file (ex. casm_ibm_tsv42.exe)
	Line 15: CASM control file (ex. casm_IBM_TS_control.dat)
	Line 16: CASM habitat parameters 
		(ex. casm_TS_bio_parms_RSK_WithTS_08Oct2018.txt)
	Line 17: CASM food web parameters (ex. web_casmTS_HRM_20Jun2018.dat)
	Line 18: CASM environmental time series (ex. env_casmTS_2010_08May2018.prn)
	Line 19: CASM contaminant time series (ex. ATZ_exp_zero_test.dat)
	Line 20: CASM toxicity parameters (ex. TS_toxicity_atrazine_tri_base.dat)
		Note: toxicity scenarios are not considered in this study BUT these
			files need to be present for CASM to run correctly)
	Line 23: number of years to simulate

Run the R script.

Go to the scenario folder created by the R script and run the CASM executable

## Run TS-IBM ##################################################################
Before starting: 
	NetLogo 6 must be installed
		Note: for runs in this study, NetLogo 6.0.2 was used.
	R package RNetLogo must be installed
		
File: 2_IBMOnly.R

Change the following lines:
	Line 9: Directory of project (must already exist)
	Line 10: Name of scenario
		Note: this will become a folder in the project directory
		
	Line 13: NetLogo app folder 
		(ex. C:/Program Files/NetLogo 6.0.2/app)
	Line 14: NetLogo jar file (in NetLogo app folder) (ex. netlogo-6.0.2.jar)
	
	Line 17: Folder where all TS-IBM files are located
	Line 18: TS-IBM NetLogo Model (ex. TS_Model_NoPesticide_27Nov2018.nlogo)
	Line 19: TS-IBM Input Parameter Files (ex. InputParameters_20Jun2018.txt)
	
	Line 22: Folder with CASM scenario (from 1_CASM.R)
	Line 23: CASM environmental time series (ex. env_casmTS_2010_08May2018.prn)
	Line 24: CASM simulation year to use in TS-IBM

	Line 27: TS-IBM RandomNumberSeed input, numeric
	Line 28: TS-IBM pondArea input, numeric
	Line 29: TS-IBM CASM_pool_area input, numeric
	Line 30: TS-IBM Sunfish input, logical
	Line 31: TS-IBM Predation input, logical
	Line 32: TS-IBM DD_egg_larva input, logical
	Line 33: TS-IBM AverageKmin input, numeric
	Line 34: TS-IBM ScalingSearchArea input, numeric
	Line 35: TS-IBM MaxDetritusDepth input, numeric
	Line 36: TS-IBM YearsToRun input, numeric
	
Run the R script.

## Run Linked CASM TS-IBM ######################################################
Before starting: 
	NetLogo 6 must be installed
		Note: for runs in this study, NetLogo 6.0.2 was used.
	R package RNetLogo must be installed
		
File: 3_Linked.R

Change the following lines:
	Line 9: Directory of project (must already exist)
	Line 10: Name of scenario
		Note: this will become a folder in the project directory
		
	Line 13: NetLogo app folder 
		(ex. C:/Program Files/NetLogo 6.0.2/app)
	Line 14: NetLogo jar file (in NetLogo app folder) (ex. netlogo-6.0.2.jar)
	
	Line 17: Folder where all TS-IBM files are located
	Line 18: TS-IBM NetLogo Model (ex. TS_Model_NoPesticide_17Jul2018.nlogo)
	Line 19: TS-IBM Input Parameter Files (ex. InputParameters_20Jun2018.txt)
	
	Line 22: Folder with CASM scenario (from 1_CASM.R)
	Line 23: CASM executable file (ex. casm_IBM_TSv32.exe)
	Line 24: CASM control file (ex. casm_IBM_TS_control_v32.dat)
	Line 25: CASM habitat parameters 
		(ex. casm_TS_bio_parms_HRM_WithTS_29Aug2018.dat)
	Line 26: CASM food web parameters (ex. web_casmTS_HRM_20Jun2018.dat)
	Line 27: CASM environmental time series (ex. env_casmTS_2010_08May2018.prn)
	Line 28: CASM contaminant time series (ex. ATZ_exp_zero_test.dat)
	Line 29: CASM toxicity parameters (ex. TS_toxicity_atrazine_tri_base.dat)
		Note: toxicity scenarios are not considered in this study BUT these
			files need to be present for CASM to run correctly)
	Line 30: Year where Day 1 is pulled from CASM scenario
	
	Line 33: TS-IBM RandomNumberSeed input, numeric
	Line 34: TS-IBM pondArea input, numeric
	Line 35: TS-IBM CASM_pool_area input, numeric
	Line 36: TS-IBM Sunfish input, logical
	Line 37: TS-IBM Predation input, logical
	Line 38: TS-IBM DD_egg_larva input, logical
	Line 39: TS-IBM AverageKmin input, numeric
	Line 40: TS-IBM ScalingSearchArea input, numeric
	Line 41: TS-IBM MaxDetritusDepth input, numeric
	Line 42: TS-IBM YearsToRun input, numeric
	
Run the R script.
