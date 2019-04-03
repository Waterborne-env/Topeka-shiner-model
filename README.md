# Topeka Shiner Hybrid Model: How To
Author: Colleen Roy

Last Edited: 2019-04-03


## Disclaimer
The software and associated files uploaded in this repository were used to generate the results to be published in:

Schmolke A, Bartell SM, Roy C, Green N, Galic N, Brain R. Species-specific population dynamics and their link to an aquatic food web: a hybrid modeling approach.

The manuscript has been submitted for review to a scientific journal.

*This software and associated files are provided "as is" with the sole purpose to allow the reproduction of the published results without any warranties of performance or fitness for any other purpose.*

## Overview
This how to documents how to run CASM, the TS-IBM, and the linked CASM + TS-IBM in R.

*In order to run either the TS-IBM or the Linked CASM + TS-IBM you must:*
1. Have NetLogo 6 installed (for runs in this study, NetLogo 6.0.2 was used)
2. Have the R package "RNetLogo" installed
3. Run CASM first (see below for instructions)

## CASM
**File:** 1_CASM.R

### Change the following lines:
**_Initial Folder Locations/Names:_**

	Line 9: Directory of project (must already exist)
	Line 10: Name of scenario (note: this will become a folder in the project directory)
	
**_CASM Files:_**
	
	Line 13: CASM executable file (ex. casm_ibm_tsv42.exe)
	Line 14: CASM control file (ex. casm_IBM_TS_control.dat)
	Line 15: CASM habitat parameters (ex. casm_TS_bio_parms_RSK_WithTS_08Oct2018.txt)
	Line 16: CASM food web parameters (ex. web_casmTS_HRM_20Jun2018.dat)
	Line 17: CASM environmental time series (ex. env_casmTS_2010_08May2018.prn)
	Line 18: CASM contaminant time series (ex. ATZ_exp_zero_test.dat)
	Line 19: CASM toxicity parameters (ex. TS_toxicity_atrazine_tri_base.dat)
		(note: toxicity scenarios are not considered in this study BUT these files 
		       need to be present for CASM to run correctly)
Note: Lines 13-19 refer to file paths relative to the project directory (defined in line 9). If these files are not located in the project directory, full file paths must be provided instead of just the file names

**_Model Parameters:_**

	Line 22: Years to Run CASM (numeric)

### Run the Model

Run the R script.

Run CASM: Run the CASM executable located in the scenario folder (named in line 10 and created by the R script)

The output file containing daily biomasses is named `IBM_TS_master_ref.out` and is located in the scenario folder

## TS-IBM
**File:** 2_IBMOnly.R

### Change the following lines:
**_Initial Folder Locations/Names:_**

	Line 9: Directory of project (must already exist)
	Line 10: Name of scenario (note: this will become a folder in the project directory)

**_NetLogo Location/Files:_**

	Line 13: NetLogo app folder (ex. C:/Program Files/NetLogo 6.0.2/app)
	Line 14: NetLogo jar file (in NetLogo app folder) (ex. netlogo-6.0.2.jar)

**_IBM Files:_**

	Line 17: TS-IBM NetLogo Model (ex. TS_IBM_March2019.nlogo)
	Line 18: TS-IBM Input Parameter Files (ex. InputParameters_20Jun2018.txt)
Note: Lines 17-18 refer to file paths relative to the project directory (defined in line 9). If these files are not located in the project directory, full file paths must be provided instead of just the file names

**_CASM Files:_**
	
	Line 21: Name of CASM scenario created/ran previously (same value as line 9 from 1_CASM.R)
	Line 22: CASM environmental time series used in CASM scenario (ex. env_casmTS_2010_08May2018.prn)
Note: Line 22 refers to a file within the CASM scenario folder named in line 21

**_CASM Model Parameters:_**

	Line 25: CASM simulation year to use in TS-IBM (numeric)

**_IBM Model Parameters:_**

	Line 28: TS-IBM RandomNumberSeed input, numeric
	Line 29: TS-IBM pondArea input, numeric
	Line 30: TS-IBM CASM_pool_area input, numeric
	Line 31: TS-IBM Sunfish input, logical
	Line 32: TS-IBM Predation input, logical
	Line 33: TS-IBM DD_egg_larva input, logical
	Line 34: TS-IBM AverageKmin input, numeric
	Line 35: TS-IBM ScalingSearchArea input, numeric
	Line 36: TS-IBM MaxDetritusDepth input, numeric
	Line 37: TS-IBM YearsToRun input, numeric

### Run the Model

Run the R script.

The output file containing daily biomasses is named `TS-IBM_output.txt` and is located in the scenario folder

## Linked CASM + TS-IBM
**File:** 3_Linked.R

### Change the following lines:
**_Initial Folder Locations/Names:_**

	Line 9: Directory of project (must already exist)
	Line 10: Name of scenario (note: this will become a folder in the project directory)

**_NetLogo Location/Files:_**

	Line 13: NetLogo app folder (ex. C:/Program Files/NetLogo 6.0.2/app)
	Line 14: NetLogo jar file (in NetLogo app folder) (ex. netlogo-6.0.2.jar)

**_IBM Files:_**

	Line 17: TS-IBM NetLogo Model (ex. TS_IBM_March2019.nlogo)
	Line 18: TS-IBM Input Parameter Files (ex. InputParameters_20Jun2018.txt)
Note: Lines 17-18 refer to file paths relative to the project directory (defined in line 9). If these files are not located in the project directory, full file paths must be provided instead of just the file names

**_CASM Files:_**
	
	Line 21: Name of CASM scenario created/ran previously (same value as line 9 from 1_CASM.R)
	Line 22: CASM executable file (ex. casm_ibm_tsv42.exe)
	Line 23: CASM control file (ex. casm_IBM_TS_control.dat)
	Line 24: CASM habitat parameters (ex. casm_TS_bio_parms_RSK_WithTS_08Oct2018.txt)
	Line 25: CASM food web parameters (ex. web_casmTS_HRM_20Jun2018.dat)
	Line 26: CASM environmental time series (ex. env_casmTS_2010_08May2018.prn)
	Line 27: CASM contaminant time series (ex. ATZ_exp_zero_test.dat)
	Line 28: CASM toxicity parameters (ex. TS_toxicity_atrazine_tri_base.dat)
		(note: toxicity scenarios are not considered in this study BUT these files 
		       need to be present for CASM to run correctly)
Note: Line 22-28 refers to files within the CASM scenario folder named in line 21

**_CASM Model Parameters:_**

	Line 31: CASM simulation year to use as day 1 numbers in CASM + TS-IBM (numeric)

**_IBM Model Parameters:_**

	Line 34: TS-IBM RandomNumberSeed input, numeric
	Line 35: TS-IBM pondArea input, numeric, units: m2
	Line 36: TS-IBM CASM_pool_area input, numeric, units: m2
	Line 37: TS-IBM Sunfish input, logical
	Line 38: TS-IBM Predation input, logical
	Line 39: TS-IBM DD_egg_larva input, logical
	Line 40: TS-IBM AverageKmin input, numeric
	Line 41: TS-IBM ScalingSearchArea input, numeric
	Line 42: TS-IBM MaxDetritusDepth input, numeric, units: cm

### Run the Model

Run the R script. *for the linked scenario, R automatically runs CASM, no other steps are needed*

The output file containing daily biomasses is named `TS-IBM_output.txt` and is located in the scenario folder

