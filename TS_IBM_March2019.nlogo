globals [
; lists of waterbody conditions read from input file
; waterbody input files are expected to have 365 entries
  TemperatureList ; list of daily tempartures of waterbody (pond) in Celsius: from input file
  WaterDepthList  ; list of daily water depths in m: from input file
  ; variables for daily values
  temperature ; current tempartures of waterbody (pond) in Celsius: from input file
  pondDepth   ; depth of waterbody in m; from input file
  pondVolume  ; water volume; calculated from pondArea and pondDepth

  limit_detritus ; maximum detritus biomass (in gC/m2) in the pond that is available for consumption by TS
  ;; 27 Nov 2018: flag included to make sure that detritus biomass value returned to CASM / written to output is not truncated by limit_detritus
  flag_limit_detritus ; TRUE if detritus in the pond exceeds limit_detritus
; TS biomass from CASM input/transfer (only used during setup)
  ini_w_TS_juv
  ini_w_TS_adt
; TS eaten by predator
  list_TS_eaten
  w_TS_eaten          ; wet weight [g] of TS eaten by predator (largemouth bass); read from CASM biomass file
  TS_eaten_cumulative ; wet weight [g] of TS eaten; actual in model (sum of weights of eaten individuals)
  TS_eaten_num        ; number of TS eaten
; prey parameters
  ; lists with prey biomasses read from CASM output file in gC (only used for yearly input runs)
  ; NOTE: order of prey biomasses is essential when inputs are read in and outputs are produced; the order in which they are listed here reflects the order in the CASM file
  list_w_copepods     ; column name in CASM files: Copepods
  list_w_cladocerans  ; column name in CASM files: Cladocerans
  list_w_rotifers     ; column name in CASM files: Rotifers
  list_w_microzoo     ; column name in CASM files: MicroZplankton
  list_w_mayflies     ; column name in CASM files: Ephemeroptera
  list_w_caddisflies  ; column name in CASM files: Trichoptera
  list_w_oligochaetes ; column name in CASM files: Oligochaetes
  list_w_chironomids  ; column name in CASM files: Chironomids
  list_w_detritus     ; column name in CASM files: SedimentPOC
  list_w_periphyton   ; column name in CASM files: TotalPeriphyt
  list_w_TS_juv       ; column name in CASM files: TopShiner(juv) ; only used for setup
  list_w_TS_adt       ; column name in CASM files: TopShiner(adt) ; only used for setup
  w_copepods     ; wet weight [g] of prey guild copepods
  w_cladocerans  ; wet weight [g] of prey guild cladocerans
  w_rotifers     ; wet weight [g] of prey guild rotifers
  w_microzoo     ; wet weight [g] of prey guild microzooplankton
  w_mayflies     ; wet weight [g] of prey guild mayflies (Ephemeroptera)
  w_caddisflies  ; wet weight [g] of prey guild caddisflies (Trichoptera)
  w_oligochaetes ; wet weight [g] of prey guild oligochaetes (Oligochaetes)
  w_chironomids  ; wet weight [g] of prey guild chironomids
  w_detritus     ; wet weight [g] of sediment POC
  w_periphyton   ; wet weight [g] of algae guilds summarized as total periphyton
  preyatthestartoftheday ; biomass in g C per m2 of prey guild at the start of the day (before consumption by TS)
  preyattheendoftheday   ; biomass in g C per m2 of prey guild at the end of the day (after consumption by TS)
; functional response parameters according to CASM (Bartell et al. 2013)
; juvenile TS
  juv_pref_copepods    ; preference of juvenile TS for copepods
  juv_pref_cladocerans ; preference of juvenile TS for cladocerans
  juv_pref_rotifers    ; preference of juvenile TS for rotifers
  juv_pref_microzoo    ; preference of juvenile TS for microzooplankton
  juv_pref_detritus    ; preference of juvenile TS for detritus (SedimentPOC)
  juv_ass_copepods    ; assimiltation efficiency of juvenile TS for copepods
  juv_ass_cladocerans ; assimiltation efficiency of juvenile TS for cladocerans
  juv_ass_rotifers    ; assimiltation efficiency of juvenile TS for rotifers
  juv_ass_microzoo    ; assimiltation efficiency of juvenile TS for microzooplankton
  juv_ass_detritus    ; assimiltation efficiency of juvenile TS for detritus (SedimentPOC)
  juv_h_copepods    ; handling efficiency of juvenile TS for copepods
  juv_h_cladocerans ; handling efficiency of juvenile TS for cladocerans
  juv_h_rotifers    ; handling efficiency of juvenile TS for rotifers
  juv_h_microzoo    ; handling efficiency of juvenile TS for microzooplankton
  juv_h_detritus    ; handling efficiency of juvenile TS for detritus (SedimentPOC)
; adult TS
  adt_pref_copepods     ; preference of adult TS for copepods
  adt_pref_cladocerans  ; preference of adult TS for cladocerans
  adt_pref_mayflies     ; preference of adult TS for mayflies (Ephemeroptera)
  adt_pref_caddisflies  ; preference of adult TS for caddisflies (Trichoptera)
  adt_pref_oligochaetes ; preference of adult TS for oligochaetes
  adt_pref_chironomids  ; preference of adult TS for chironomids
  adt_pref_detritus     ; preference of adult TS for detritus (SedimentPOC)
  adt_pref_periphyton   ; preference of adult TS for periphyton
  adt_ass_copepods     ; assimiltation efficiency of adult TS for copepods
  adt_ass_cladocerans  ; assimiltation efficiency of adult TS for cladocerans
  adt_ass_mayflies     ; assimiltation efficiency of adult TS for mayflies (Ephemeroptera)
  adt_ass_caddisflies  ; assimiltation efficiency of adult TS for caddisflies (Trichoptera)
  adt_ass_oligochaetes ; assimiltation efficiency of adult TS for oligochaetes
  adt_ass_chironomids  ; assimiltation efficiency of adult TS for chironomids
  adt_ass_detritus     ; assimiltation efficiency of adult TS for detritus (SedimentPOC)
  adt_ass_periphyton   ; assimiltation efficiency of adult TS for periphyton
  adt_h_copepods     ; handling efficiency of adult TS for copepods
  adt_h_cladocerans  ; handling efficiency of adult TS for cladocerans
  adt_h_mayflies     ; handling efficiency of adult TS for mayflies (Ephemeroptera)
  adt_h_caddisflies  ; handling efficiency of adult TS for caddisflies (Trichoptera)
  adt_h_oligochaetes ; handling efficiency of adult TS for oligochaetes
  adt_h_chironomids  ; handling efficiency of adult TS for chironomids
  adt_h_detritus     ; handling efficiency of adult TS for detritus (SedimentPOC)
  adt_h_periphyton   ; handling efficiency of adult TS for periphyton
; conversion factor between carbon equivalent (gC) and wet weight (in g)
  gC_ww_copepods
  gC_ww_cladocerans
  gC_ww_rotifers
  gC_ww_microzoo
  gC_ww_mayflies ;  (Ephemeroptera)
  gC_ww_caddisflies ; (Trichoptera)
  gC_ww_oligochaetes
  gC_ww_chironomids
  gC_ww_detritus ; (SedimentPOC)
  gC_ww_periphyton
  gC_ww_fish ; applied to TS
; energy density of prey guilds and TS; in [J / g wet weight]
  ED_copepods
  ED_cladocerans
  ED_rotifers
  ED_microzoo
  ED_mayflies ; (Ephemeroptera)
  ED_caddisflies ; (Trichoptera)
  ED_oligochaetes
  ED_chironomids
  ED_detritus ; (SedimentPOC)
  ED_periphyton
  ED_TS_juv
  ED_TS_adt
; output
  Prey_W_output_list ; output of daily consumption by prey group ; list needs to be compiled in correct order!
  TS_W_output_list   ; output of TS weights (juveniles, subadults, adults, all fish)
; Bioenergetics parameters (read in from parameter input file)
  a_cmax ; parameter for relationship between Cmax and fish weight
  b_cmax ; parameter for relationship between Cmax and fish weight
  te1 ; Thornton & Lessem (1978) temperature dependence of consumption [degree C]
  te2 ; Thornton & Lessem (1978) temperature dependence of consumption [degree C]
  te3 ; Thornton & Lessem (1978) temperature dependence of consumption [degree C]
  te4 ; Thornton & Lessem (1978) temperature dependence of consumption [degree C]
  xk1 ; Thornton & Lessem (1978) temperature dependence of consumption
  xk2 ; Thornton & Lessem (1978) temperature dependence of consumption
  xk3 ; Thornton & Lessem (1978) temperature dependence of consumption
  xk4 ; Thornton & Lessem (1978) temperature dependence of consumption
  SDA ; specific dynamic action (metabolism)
  RA  ; respiration parameter
  RB  ; respiration parameter
  RTO ; Kitchell et al. (1977) temperature dependence of respiration [degree C]
  RTM ; Kitchell et al. (1977) temperature dependence of respiration [degree C]
  RQ  ; Kitchell et al. (1977) temperature dependence of respiration [degree C^(-1)]
  FA  ; Egestion loss constant
  UA  ; Excretion loss constant
; Bioenergetics parameters calculated from parameters and temperature
  C_temp_function_list ; list of temperature function for consumption calculated for each day using parameters from input file during setup
  fC_T ; daily temperature factor for consumption from list
  R_temp_function_list ; list of temperature function for respiration calculated for each day using parameters from input file during setup
  fR_T ; daily temperature factor for resipration from list
  OxycalCoeff ; oxycalorific coefficient [J / g O2]
; Parameters for calculation of egg numbers per female (read in from parameter input file)
  beta_0
  beta_1
; Parameters for calculation of length / weight conversion of TS (read in from paramter input file)
  l2w_a
  l2w_b
; Presence or absence of sunfish: presence of sunfish increases egg and larva survival rate
;  Sunfish ; on/off ; default: on
; Fish life history parameters (read in from parameter input file)
  DailySurvival_Egg_Larva    ; daily survival rate applied to TS eggs and larvae
  DailySurvival_Egg_Larva_0  ; maximum daily survival rate applied to TS eggs and larvae for density dependent rate (alternative to fixed DailySurvival_Egg_Larva)
  beta_dd_egg_survival       ; beta in Ricker function for density dependence of daily egg and larva survival rate (alternative to fixed DailySurvival_Egg_Larva)
  SunfishFactor            ; defines difference of daily TS egg and larva mortality rate due to absence of sunfish
  daily_surv_rate_juv      ; daily survival rate of juvenile TS
  daily_surv_rate_subadult ; daily survival rate of subadult TS
  daily_surv_rate_yr2      ; daily survival rate of adult TS, 12-23 months old
  daily_surv_rate_yr3      ; daily survival rate of adult TS, >24 months old
  StartSpawning_early      ; earliest start date (day of year) of spawning
  StartSpawning_late       ; latest start date (day of year) of spawning
  PeakSpawning_early       ; earliest peak date (day of year) of spawning
  PeakSpawning_late        ; latest peak date (day of year) of spawning
  StopSpawning_early       ; earliest end date (day of year) of spawning
  StopSpawning_late        ; latest end date (day of year) of spawning
  EggWeight                ; weight of egg (in g)
  EggDevTime               ; time between spawning and hatching (in days)
  YolkLasting              ; time between hatching and onset of feeding (days fish remain in larva stage)
  MinSizeSubadult          ; size threshold for subadult stage (in mm SL)
  MinSizeMat               ; size threshold for adult stage (in mm SL)
; Yearly spawning period defined stochastically for every simulated year:
  Yr_startSpawning ; start date (day of year) of spawning in given year
  Yr_peakSpawning  ; peak date (day of year) of spawning in given year
  Yr_stopSpawning  ; end date (day of year) of spawning in given year
]

turtles-own [
  age             ; fish age in days (age 0: day of spawning)
  stage           ; fish life stage: egg, larva, juvenile, subadult, adult
  sex             ; fish sex: female, male
  W               ; fish wet weight in g
  SL              ; fish standard length in mm
  K               ; fish condition (actual weight / weight expected for SL)
  Kmin            ; fish-specific minimum condition for survival (if K < Kmin, fish dies)
  day_consumption ; total daily consumption in J per fish
  R_T             ; daily losses from respiration in J / g
  S               ; specific dynamic action losses in J / g
  E_C             ; daily waste losses in J / g
  d_W             ; daily weight change in g
  spawn_date      ; spawning date (doy) for the given spawning season (set to -1 when fish is not spawning)
]

patches-own [ ]

to setup
  clear-all
  file-close
  reset-ticks

  random-seed RandomNumberSeed

  set Prey_W_output_list []; output of daily consumption by prey group
  set TS_W_output_list []; output of TS weights (juveniles, subadults, adults, all fish)

  readInputParameters ; reading parameter input file containing parameters of bionenegetics functions etc.
  readWaterConditions ; reading daily temperature and water depth

  set pondDepth item 0 WaterDepthList
  set pondVolume pondArea * pondDepth ; m^3; pond depth comes from input file, pond area is set on interface

  set limit_detritus MaxDetritusDepth * 10000 * 1.25 / gC_ww_detritus ; in gC/m2 ; assumption that the density of the organic matter in the detritus (Sediment POC) is 1.25 g/ml
  ;; 27 Nov 2018: flag included to make sure that detritus biomass value returned to CASM / written to output is not truncated by limit_detritus
  set flag_limit_detritus FALSE

  ; calculation of temperature-dependence of consumption rates, stored as look-up list for simulation run:
  set C_temp_function_list []
  let cnt 0
  while [ cnt < length TemperatureList ] [
   set C_temp_function_list lput (consumption_temperature_function (item cnt TemperatureList)) C_temp_function_list
   set cnt cnt + 1
  ]

  ; calculation of temperature-dependence of respiration rates, stored as look-up list for simulation run:
  set R_temp_function_list []
  set cnt 0
  while [ cnt < length TemperatureList ] [
   set R_temp_function_list lput (respiration_temperature_function (item cnt TemperatureList)) R_temp_function_list
   set cnt cnt + 1
  ]

  ; if DailyIO is on, IBM reads in prey biomass input at the start of each time step
  ; otherwise, IBM reads in prey biomass for a whole year (365 days)
  ifelse (DailyIO) [
    readDailyPreyBiomass_day
  ] [
    readDailyPreyBiomass_year
    set ini_w_TS_juv item 0 list_w_TS_juv
    set ini_w_TS_adt item 0 list_w_TS_adt
  ]

  ; initial weight of TS is pooled for juveniles and adults (initial weights in CASM):
  ; assumption that the initial date is 1 January when no juveniles are present in IBM
  ; biomass is split up amongst TS stages/ages according to reported stage distribution (Kerns and Bonneau 2002)
  ; average weights of age classes are applied as initial conditions of fish

  ; create subadults
  let num_subadults round (((ini_w_TS_juv + ini_w_TS_adt) * 0.508) / (calc_weight 25.6))
  create-turtles num_subadults [
    setxy random-xcor random-ycor
    set shape "fish"
    set stage "subadult"
    set age 169 + (random 44)
    set SL 25.6
    set W (calc_weight SL )
    set K 1
    set Kmin random-normal AverageKmin (AverageKmin / 10)
    set day_consumption 0
    set spawn_date -1
    set size 3
  ]
  let num_adt_2ndyr round (((ini_w_TS_juv + ini_w_TS_adt) * 0.39) / (calc_weight 40.1))
  create-turtles num_adt_2ndyr [
    setxy random-xcor random-ycor
    set shape "fish"
    set stage "adult"
    set age 534  + (random 44)
    set SL 40.1
    set W (calc_weight SL)
    set K 1
    set Kmin random-normal AverageKmin (AverageKmin / 10)
    set day_consumption 0
    set spawn_date -1
    set size 4
  ]
  let num_adt_3rdyr round (((ini_w_TS_juv + ini_w_TS_adt) * 0.102) / (calc_weight 42.5))
  create-turtles num_adt_3rdyr [
    setxy random-xcor random-ycor
    set shape "fish"
    set stage "adult"
    set age 899 + (random 44)
    set SL 42.5
    set W (calc_weight SL)
    set K 1
    set Kmin random-normal AverageKmin (AverageKmin / 10)
    set day_consumption 0
    set spawn_date -1
    set size 4
  ]

  ; sex ratio of 1:1 assumed; sex assigned randomly
  ask turtles
  [
    let rnd_num random 2
    ifelse rnd_num = 0 [
      set sex "female"
    ] [
      set sex "male"
    ]
  ]

  ; only used for testing: output of individual fish weight etc. for calibration of fish growth; written out for single fish
  if WriteGrowthTestOutput [
    if file-exists? FishGrowthTest [
      if user-yes-or-no? "FishGrowthTest file already exists. Delete?" [
        file-delete FishGrowthTest
      ]
    ]
    file-open FishGrowthTest
    file-print "Version: TS_Model_21May2018.nlogo"
    file-write "Setup fish file:"
    file-write FishStatusIn
    file-print ""
    file-write "day"
    file-write "age"
    file-write "stage"
    file-write "sex"
    file-write "W"
    file-write "SL"
    file-write "day_consumption"
    file-write "respiration"
    file-write "spec dyn action"
    file-write "waste losses"
    file-write "weight change"
    file-write "w_copepods"
    file-write "w_cladocerans"
    file-write "w_rotifers"
    file-write "w_microzoo"
    file-write "w_mayflies" ;  (Ephemeroptera)
    file-write "w_caddisflies" ; (Trichoptera)
    file-write "w_oligochaetes"
    file-write "w_chironomids"
    file-write "w_detritus"  ; (SedimentPOC)
    file-write "w_periphyton"
    file-write "temperature"
    file-print ""
    file-close
  ]

  reset-ticks
end

to SetupDay
  ; alternative setup for growth calibration:
  ; single fish is initialized from input file on a chosen day of year
  clear-all
  file-close
  reset-ticks

  random-seed RandomNumberSeed

  set Prey_W_output_list [] ; output of daily consumption by prey group in gC/m2
  set TS_W_output_list []   ; output of TS weights in gC/m2 (juveniles, subadults, adults, all fish, actual TS eaten)

  readInputParameters ; parameters of bionenegetics functions etc.
  readWaterConditions ; reads daily temperature and water depth

  set limit_detritus MaxDetritusDepth * 10000 * 1.25 / gC_ww_detritus ; in gC/m2 ; assumption that the density of the organic matter in the detritus (Sediment POC) is 1.25 g/ml
  ;; 27 Nov 2018: flag included to make sure that detritus biomass value returned to CASM / written to output is not truncated by limit_detritus
  set flag_limit_detritus FALSE

  ; calculation of temperature-dependence of consumption rates, stored as look-up list for simulation run:
  set C_temp_function_list []
  let cnt 0
  while [ cnt < length TemperatureList ] [
   set C_temp_function_list lput (consumption_temperature_function (item cnt TemperatureList)) C_temp_function_list
   set cnt cnt + 1
  ]

  ; calculation of temperature-dependence of respiration rates, stored as look-up list for simulation run:
  set R_temp_function_list []
  set cnt 0
  while [ cnt < length TemperatureList ] [
   set R_temp_function_list lput (respiration_temperature_function (item cnt TemperatureList)) R_temp_function_list
   set cnt cnt + 1
  ]

  ; if DailyIO is on, IBM reads in prey biomass input at the start of each time step
  ; otherwise, IBM reads in prey biomass for a whole year (365 days)
  if (not DailyIO) [ readDailyPreyBiomass_year ]

  ; reading parameters of single fish from file; this alternative setup does not work without this input file
  ifelse file-exists? FishStatusIn [
    file-open FishStatusIn
    let header file-read-line
;    print header ; for testing
    while [not file-at-end?] [
      create-turtles 1 [
        setxy random-xcor random-ycor
        set shape "fish"
        set age file-read
;        type "read-in age: " print age ; for testing
        set stage file-read
        set sex file-read
        set W file-read
        set SL file-read
        set K 1
        set Kmin file-read
        set day_consumption 0
        if stage = "egg" [
          set shape "dot"
          set size 1
        ]
        if stage = "juvenile" [ set size 2 ]
        if stage = "subadult" [ set size 3 ]
        if stage = "adult" [ set size 4 ]
      ]
    ]
    file-close
  ] [ user-message "FishStatusIn file does not exist in current directory!" ]

  ; daily fish size values are written to output file; if the file with the name already exists in the folder, it is overwritten
  if WriteGrowthTestOutput [
    if file-exists? FishGrowthTest [
      if user-yes-or-no? "FishGrowthTest file already exists. Delete?" [
        file-delete FishGrowthTest
      ]
    ]
    file-open FishGrowthTest
    file-print "Version: TS_Model_21May2018.nlogo"
    file-write "Setup fish file:"
    file-write FishStatusIn
    file-print ""
    file-write "day"
    file-write "age"
    file-write "stage"
    file-write "sex"
    file-write "W"
    file-write "SL"
    file-write "K"
    file-write "day_consumption"
    file-write "respiration"
    file-write "spec dyn action"
    file-write "waste losses"
    file-write "weight change"
    file-write "w_copepods"
    file-write "w_cladocerans"
    file-write "w_rotifers"
    file-write "w_microzoo"
    file-write "w_mayflies" ;  (Ephemeroptera)
    file-write "w_caddisflies" ; (Trichoptera)
    file-write "w_oligochaetes"
    file-write "w_chironomids"
    file-write "w_detritus"  ; (SedimentPOC)
    file-write "w_periphyton"
    file-write "temperature"
    file-print ""
    file-close
  ]

  tick-advance StartDay ; setting simulation run to chosen start day of year

end

to go
  ; main control procedure: called for every time step (day)
    if not any? turtles [ stop ] ; if no TS are left, simulation run is stopped

  ; for each new simulated year, new starting, peak and end dates (day of year) for spawning are chosen from the ranges defined in the parameter input file
  if ticks mod 365 = 0 [
    if ticks / 365 = YearsToRun [ stop ] ; maximum run time in years can be set on the interface
    set Yr_startSpawning StartSpawning_early + random (StartSpawning_late - StartSpawning_early)
    set Yr_peakSpawning PeakSpawning_early + random (PeakSpawning_late - PeakSpawning_early)
    set Yr_stopSpawning StopSpawning_early + random (StopSpawning_late - StopSpawning_early)
  ]

  ; after spawning season of the year has ended:
  if (ticks mod 365) > Yr_stopSpawning [
    ask turtles with [ age >= (3 * 365) ] [ die ] ; fish die after spawning season if they are 3 years old
    ask turtles [ set spawn_date -1 ]
  ]

  ; waterbody conditions for the day of the year are set
  set temperature item (ticks mod 365) TemperatureList
  set pondDepth item (ticks mod 365) WaterDepthList
  set pondVolume pondArea * pondDepth ; m^3; pond depth needs to come from input file

  ; temperature-dependence of consumption and respiration for the day of year
  set fC_T item (ticks mod 365) C_temp_function_list
  set fR_T item (ticks mod 365) R_temp_function_list

  ; bioenergetics variables are reset at the start of each day (new calculation conducted each day)
  ask turtles [
    set d_W 0
    set day_consumption 0
    set R_T 0
    set E_C 0
  ]

  ; extract daily wet weights [g] of prey guilds from data read from biomass input file (from CASM); weight is calculated from gC to g wet weight
  ifelse DailyIO
  [ readDailyPreyBiomass_day ]
  [ setDailyPreyBiomass ]

  ; variables keeping track of food consumption across all TS and predation on TS reset:
  set TS_eaten_cumulative 0
  set TS_eaten_num 0

  ; predation on TS
  if predation = true [ ; mortality due to predation by largemouth bass (assumption that only subadults and adults are affected)
                        ; adjust according to input received from CASM: assumed here that we have the biomass of TS consumed by LB per day
;    type "weight of TS eaten today: " print w_TS_eaten ; for testing
    let flag_TS_pred TRUE
    while [ flag_TS_pred and count turtles with [stage = "subadult" or stage = "adult"] > 0] [
      let focalTS one-of turtles with [stage = "subadult" or stage = "adult"] ; random predation assumed
;      type "fish " print focalTS ; for testing
      ; If biomass of TS eaten is less than the currently chosen fish, the fish has a chance of dying
      ifelse (w_TS_eaten - TS_eaten_cumulative) < 0 [
        set flag_TS_pred FALSE
      ] [
        ifelse ((w_TS_eaten - TS_eaten_cumulative) >= ([W] of focalTS)) [
          set TS_eaten_cumulative TS_eaten_cumulative + ([W] of focalTS)
          ask focalTS [
;            type "fish " type focalTS type " eaten; weight: " print [W] of focalTS ; for testing
            die ]
          set TS_eaten_num TS_eaten_num + 1
        ] [
          let prob random-float 1.0
          let tmp_ratio (w_TS_eaten - TS_eaten_cumulative) / ([W] of focalTS) ; predation always occurs if the currently called fish is smaller than the weight eaten
          ifelse prob < tmp_ratio [
            set TS_eaten_cumulative TS_eaten_cumulative + ([W] of focalTS)
            ask focalTS [
 ;             type "fish " type focalTS type " eaten; weight: " print [W] of focalTS ; for testing
              die ]
            set TS_eaten_num TS_eaten_num + 1
          ] [
            set flag_TS_pred FALSE
          ]
        ]
      ]
    ]
;    type "Fish eaten by largemouth bass, weight: " type TS_eaten_cumulative type ", number: " print TS_eaten_num ; for testing
  ]

  ; background mortality applied: derived from demographic rates
  ask turtles [
    let prob random-float 1.0
    if (K < Kmin) [ die ] ; if fish condition is < Kmin, fish dies; Kmin determination independent of fish stage/age
    if stage = "juvenile" [
      if (prob > daily_surv_rate_juv) [ die ]
    ]
    if stage = "subadult" [
      if (prob > daily_surv_rate_subadult) [ die ]
    ]
    if (age > 365 and age <= 730) [
      if (prob > daily_surv_rate_yr2) [ die ]
    ]
    if (age > 730) [
      if (prob > daily_surv_rate_yr3) [ die ]
    ]
  ]

  ; TS growth resulting from daily consumption, respiration and waste losses
  ; for testing: consumption set to Cmax (independent of prey availability)
  ifelse CmaxConsumption [
    ask turtles with [stage = "juvenile" or stage = "subadult" or stage = "adult"] [
      set day_consumption ConstCmax
      grow
    ]
  ] [
  ; consumption according to functional response from CASM
    ask turtles with [stage = "juvenile" or stage = "subadult" or stage = "adult"] [
      functional_response_CASM
      grow
    ]
  ]

  ; spawning if current day is within spawning period
  ; each adult female is assigned with actual spawn date(s) at onset of spawning period
  if (ticks mod 365) = Yr_startSpawning [
    ask turtles with [ stage = "adult" and sex = "female" ] [
        set spawn_date triang_dist Yr_startSpawning Yr_peakSpawning Yr_stopSpawning
    ]
  ]

  ; fish prompted to spawn on chosen spawn date
  if (ticks mod 365) >= Yr_startSpawning and (ticks mod 365) < Yr_stopSpawning [
    ask turtles with [ stage = "adult" and sex = "female" ] [
      if ((ticks mod 365) = spawn_date) [
          spawn
      ]
    ]
  ]

  ; prey and TS biomasses converted back to gC/m2 for output
  ; note that biomasses are given as biomass per area (not water volume) for inputs and outputs
  ; IMPORTANT: order of biomasses in output list must be identical with input/order in CASM transfer files!!!
  set preyattheendoftheday []
  set preyattheendoftheday lput (w_copepods / (pondArea * gC_ww_copepods)) preyattheendoftheday
  set preyattheendoftheday lput (w_cladocerans / (pondArea * gC_ww_cladocerans)) preyattheendoftheday
  set preyattheendoftheday lput (w_rotifers / (pondArea * gC_ww_rotifers)) preyattheendoftheday
  set preyattheendoftheday lput (w_microzoo / (pondArea * gC_ww_microzoo)) preyattheendoftheday
  set preyattheendoftheday lput (w_mayflies / (pondArea * gC_ww_mayflies)) preyattheendoftheday ;  (Ephemeroptera)
  set preyattheendoftheday lput (w_caddisflies / (pondArea * gC_ww_caddisflies)) preyattheendoftheday ; (Trichoptera)
  set preyattheendoftheday lput (w_oligochaetes / (pondArea * gC_ww_oligochaetes)) preyattheendoftheday
  set preyattheendoftheday lput (w_chironomids / (pondArea * gC_ww_chironomids)) preyattheendoftheday
;; 27 Nov 2018: accounting for detritus limit
  ifelse flag_limit_detritus  [
    let detritus_tmp w_detritus / (pondArea * gC_ww_detritus)
;    type "detritus after consumption in gC/m2: " print detritus_tmp
    let detritus_cons limit_detritus - detritus_tmp ; detritus consumed during the day in gC/m2
;    type "detritus consumption in gC/m2: " print detritus_cons
    if detritus_cons < 0 [ user-message "detritus consumption negative! (from 'go' procedure)" ]
;    type "detritus at the start of the day in gC/m2: " print item 8 preyatthestartoftheday
    let detritus_act (item 8 preyatthestartoftheday) - detritus_cons
;    type "detritus at the end of the day in gC/m2: " print detritus_act
    set preyattheendoftheday lput detritus_act preyattheendoftheday
  ] ; (SedimentPOC)
  [ set preyattheendoftheday lput (w_detritus / (pondArea * gC_ww_detritus)) preyattheendoftheday ] ; (SedimentPOC)
  set preyattheendoftheday lput (w_periphyton / (pondArea * gC_ww_periphyton)) preyattheendoftheday
  set Prey_W_output_list lput preyattheendoftheday Prey_W_output_list
  ; values for TS biomass written out to the output file after the prey biomasses
  let TS_W []
  ; juvenile and subadult + adult TS biomasses reported as gC/m2
  set TS_W lput ((sum [W] of turtles with [stage = "juvenile"]) / (pondArea * gC_ww_fish)) TS_W
  set TS_W lput ((sum [W] of turtles with [stage = "subadult" or stage = "adult"]) / (pondArea * gC_ww_fish)) TS_W
  set TS_W lput ((sum [W] of turtles) / (pondArea * gC_ww_fish)) TS_W ; note that the sum includes eggs and larvae in addition to the three stages above
  ; TS numbers (all TS; stage-specific numbers per m2)
  set TS_W lput ((count turtles with [stage = "juvenile"]) / pondArea) TS_W
  set TS_W lput ((count turtles with [stage = "subadult" or stage = "adult"]) / pondArea) TS_W
  set TS_W lput ((count turtles) / pondArea) TS_W ; note that the sum includes eggs and larvae in addition to the three stages above
  ; TS weight and numbers eaten (by largemouth bass) in IBM (eaten biomass diverges from value read in from file because it is translated into specific individual fish eaten)
  set TS_W lput (TS_eaten_num / pondArea) TS_W
  set TS_W lput (TS_eaten_cumulative / pondArea) TS_W
  ;; 3 December 2018: lines below added to output table
  set TS_W lput ((sum [W] of turtles with [stage = "egg" and sex = "female"]) / (pondArea * gC_ww_fish)) TS_W ; note that using fish conversion factor might not be accurate to apply to egg weight
  set TS_W lput ((sum [W] of turtles with [stage = "egg" and sex = "male"]) / (pondArea * gC_ww_fish)) TS_W ; note that using fish conversion factor might not be accurate to apply to egg weight
  with-local-randomness [ ;; local randomness added here because larvae output was added after initial runs: calling sum of turtles uses random number generator
  set TS_W lput ((sum [W] of turtles with [stage = "larva" and sex = "female"]) / (pondArea * gC_ww_fish)) TS_W ; note that using fish conversion factor might not be accurate to apply to egg weight
  set TS_W lput ((sum [W] of turtles with [stage = "larva" and sex = "male"]) / (pondArea * gC_ww_fish)) TS_W ; note that using fish conversion factor might not be accurate to apply to egg weight
  ]
  set TS_W lput ((sum [W] of turtles with [stage = "juvenile" and sex = "female"]) / (pondArea * gC_ww_fish)) TS_W ; note that using fish conversion factor might not be accurate to apply to egg weight
  set TS_W lput ((sum [W] of turtles with [stage = "juvenile" and sex = "male"]) / (pondArea * gC_ww_fish)) TS_W ; note that using fish conversion factor might not be accurate to apply to egg weight
  set TS_W lput ((sum [W] of turtles with [stage = "subadult" and sex = "female"]) / (pondArea * gC_ww_fish)) TS_W
  set TS_W lput ((sum [W] of turtles with [stage = "subadult" and sex = "male"]) / (pondArea * gC_ww_fish)) TS_W
  set TS_W lput ((sum [W] of turtles with [stage = "adult" and age < 730 and sex = "female"]) / (pondArea * gC_ww_fish)) TS_W ; adults 1 yr old
  set TS_W lput ((sum [W] of turtles with [stage = "adult" and age < 730 and sex = "male"]) / (pondArea * gC_ww_fish)) TS_W ; adults 1 yr old
  set TS_W lput ((sum [W] of turtles with [stage = "adult" and age >= 730 and age < 1095 and sex = "female"]) / (pondArea * gC_ww_fish)) TS_W ; adults 2 yrs old
  set TS_W lput ((sum [W] of turtles with [stage = "adult" and age >= 730 and age < 1095 and sex = "male"]) / (pondArea * gC_ww_fish)) TS_W ; adults 2 yrs old
  set TS_W lput ((sum [W] of turtles with [stage = "adult" and age >= 1095 and sex = "female"]) / (pondArea * gC_ww_fish)) TS_W ; adults 3 yrs old
  set TS_W lput ((sum [W] of turtles with [stage = "adult" and age >= 1095 and sex = "male"]) / (pondArea * gC_ww_fish)) TS_W ; adults 3 yrs old
  set TS_W lput ((count turtles with [stage = "egg" and sex = "female"]) / pondArea) TS_W
  set TS_W lput ((count turtles with [stage = "egg" and sex = "male"]) / pondArea) TS_W
  set TS_W lput ((count turtles with [stage = "larva" and sex = "female"]) / pondArea) TS_W
  set TS_W lput ((count turtles with [stage = "larva" and sex = "male"]) / pondArea) TS_W
  set TS_W lput ((count turtles with [stage = "juvenile" and sex = "female"]) / pondArea) TS_W
  set TS_W lput ((count turtles with [stage = "juvenile" and sex = "male"]) / pondArea) TS_W
  set TS_W lput ((count turtles with [stage = "subadult" and sex = "female"]) / pondArea) TS_W
  set TS_W lput ((count turtles with [stage = "subadult" and sex = "male"]) / pondArea) TS_W
  set TS_W lput ((count turtles with [stage = "adult" and age < 730 and sex = "female"]) / pondArea) TS_W ; adults 1 yr old
  set TS_W lput ((count turtles with [stage = "adult" and age < 730 and sex = "male"]) / pondArea) TS_W ; adults 1 yr old
  set TS_W lput ((count turtles with [stage = "adult" and age >= 730 and age < 1095 and sex = "female"]) / pondArea) TS_W ; adults 2 yrs old
  set TS_W lput ((count turtles with [stage = "adult" and age >= 730 and age < 1095 and sex = "male"]) / pondArea) TS_W ; adults 2 yrs old
  set TS_W lput ((count turtles with [stage = "adult" and age >= 1095 and sex = "female"]) / pondArea) TS_W ; adults 3 yrs old
  set TS_W lput ((count turtles with [stage = "adult" and age >= 1095 and sex = "male"]) / pondArea) TS_W ; adults 3 yrs old
  set TS_W_output_list lput TS_W TS_W_output_list

  ; development of eggs and larvae (neither are consuming prey)
  ask turtles with [stage = "larva"] [ larvaDevelopment ]
  ask turtles with [stage = "egg"] [ eggDevelopment ]

  ; growth output only produced for model testing/calibration; output for single fish
  if WriteGrowthTestOutput [ writeFishGrowthTest ]

  ;; 27 Nov 2018: resetting detritus flag after each simulated day
  set flag_limit_detritus FALSE

  tick
end

to eggDevelopment
  ; egg weight assumed to stay constant during egg development period
  ; applying daily egg background mortality
  ; if egg/larva background mortality is density dependent, Ricker function is applied
  let bkgr_survival 0
  ifelse DD_egg_larva = TRUE ; egg survival dependent on density (biomass) of subadult and adult TS present
  [
    set bkgr_survival DailySurvival_Egg_Larva_0 * exp ((-1) * beta_dd_egg_survival * (sum [W] of turtles with [stage = "subadult" or stage = "adult"]) / pondArea)
;    type "dd survival: " type bkgr_survival type ", TS biomass: " print (sum [W] of turtles with [stage = "subadult" or stage = "adult"]) ; for testing
  ] [
    set bkgr_survival DailySurvival_Egg_Larva
  ]

  ; if sunfish are absent, eggs experience decreased survival due to lack of guarding
  ; decreased egg survival in the absence of sunfish is assumed to occur in addition to density dependence
  ifelse Sunfish = TRUE
  [ if ((random-float 1.0) > bkgr_survival) [ die ] ]
  [ if ((random-float 1.0) > (bkgr_survival * SunfishFactor)) [ die ] ] ; background survival rate decreased in the absence of sunfish hosts

  ; eggs 'hatch' to become larvae after defined time of egg development time
  set age age + 1
  if age > EggDevTime [
    set stage "larva"
    set shape "fish"
    set size 1
  ]
end

to larvaDevelopment
  ; larva do not eat, but are assumed to use up their yolk: their biomass declines due to respiration
  ; daily background survival rate of larvae is identical to eggs

  ; if egg/larva background mortality is density dependent, Ricker function is applied
  let bkgr_survival 0
  ifelse DD_egg_larva = TRUE ; larva survival dependent on density (biomass) of subadult and adult TS present
  [
    set bkgr_survival DailySurvival_Egg_Larva_0 * exp ((-1) * beta_dd_egg_survival * (sum [W] of turtles with [stage = "subadult" or stage = "adult"]) / pondArea)
  ] [
    set bkgr_survival DailySurvival_Egg_Larva
  ]
  if ((random-float 1.0) > bkgr_survival) [ die ]

  set day_consumption 0
  grow ; larvae experience negative growth due to lack of consumption

  ; new 8 June ; I was here
;  let target_SL (calc_length preStarvW)
;  type "preStarvW: " type preStarvW type ", actual W: " print W
;  type "target SL: " type target_SL type ", actual SL: " print SL
;  if d_W < 0 [ set starvation TRUE ]
;  ifelse starvation = TRUE [ ; fish is assumed to be starving until target weight for SL is reached (fish do not shrink in length)
;    type "preStarvW: " print preStarvW
;    if W < 0.68 * preStarvW [ die ] ; apply threshold: if weight drops below 68% to pre-starvation weight, the fish dies
;    if SL >= target_SL [
;      print "reset starvation"
;      set starvation FALSE
;      set preStarvW W
;      set SL (calc_length W)
;    ]
;  ] [
;    set SL (calc_length W)
;    set preStarvW W ; weight prior to starvation
;  ]

  ; YolkLasting defines the age when larvae transition to juvenile stage (i.e. used up their yolk and start feeding)
  if W <= 0 [ die ]
  if age > (EggDevTime + YolkLasting)  [
    set stage "juvenile"
    set SL (calc_length W)
    set K 1
    set size 2
  ]
end

to functional_response_CASM
  ; consumption is based on wet weight (of prey)
  ; Subadults and adults are assumed to have identical diet preferences; juveniles are assumed to eat a different set of food web guilds

  set day_consumption 0 ; in J
  let Cmax_t 0 ; in J / g wet weight of TS
  let ww_Cmax_juv 0 ; in g wet weight of prey / g wet weight of TS
  let ww_Cmax_adt 0 ; in g wet weight of prey / g wet weight of TS

  ; scaling factors of search area/volume dependent of fish SL in functional response
  ; factors are divided by pondArea/Volume because prey biomasses in functional response implementation below are stated as biomasses in the whole waterbody
  let cfa (ScalingSearchArea * (SL / 1000) ^ (3 / 4)) / pondArea ; daily area of TS foraging (dependent on fish length): area used for ground-dwelling prey and detritus
  let cfv (ScalingSearchArea * (SL / 1000) ^ (3 / 4)) / pondVolume ; daily volume of TS foraging (dependent on fish length): volume used for prey found in water column

  ; temperature dependence of consumption
  set Cmax_t (a_cmax * (W ^ b_cmax) * fC_T) ; Cmax per g fish weight
;    type "Cmax_t: " type Cmax_t type ", fC_T: " print fC_T; for testing

;  type "fish stage: " type stage type ", Cmax: " print Cmax_t ; for testing
  ifelse stage = "juvenile" [ ; consumption (in J) of juvenile TS
    set ww_Cmax_juv Cmax_t

    ; Consumption per prey type dependent on prey biomass
    let sum_rel_C juv_pref_copepods * juv_ass_copepods * juv_h_copepods * w_copepods * cfv
                + juv_pref_cladocerans * juv_ass_cladocerans * juv_h_cladocerans * w_cladocerans * cfv
                + juv_pref_rotifers * juv_ass_rotifers * juv_h_rotifers * w_rotifers * cfv
                + juv_pref_microzoo * juv_ass_microzoo * juv_h_microzoo * w_microzoo * cfv
                + juv_pref_detritus * juv_ass_detritus * juv_h_detritus * w_detritus * cfa; (SedimentPOC) ; assumption that fish can access less area of detritus for daily consumption

   ; consumption rates in g prey guild per fish and day
   let C_copepods W * ww_Cmax_juv * (juv_pref_copepods * juv_ass_copepods * juv_h_copepods * w_copepods * cfv) / (1 + sum_rel_C)
   if C_copepods > W * ww_Cmax_juv [
     user-message "CASM functional response: copepod consumption rate too high!"
     stop
   ]
   let C_cladocerans W * ww_Cmax_juv * (juv_pref_cladocerans * juv_ass_cladocerans * juv_h_cladocerans * w_cladocerans * cfv) / (1 + sum_rel_C)
   if C_cladocerans > W * ww_Cmax_juv [
     user-message "CASM functional response: cladoceran consumption rate too high!"
     stop
   ]
   let C_rotifers W * ww_Cmax_juv * (juv_pref_rotifers * juv_ass_rotifers * juv_h_rotifers * w_rotifers * cfv) / (1 + sum_rel_C)
   if C_rotifers > W * ww_Cmax_juv [
     user-message "CASM functional response: rotifers consumption rate too high!"
     stop
   ]
   let C_microzoo W * ww_Cmax_juv * (juv_pref_microzoo * juv_ass_microzoo * juv_h_microzoo * w_microzoo * cfv) / (1 + sum_rel_C)
   if C_microzoo > W * ww_Cmax_juv [
     user-message "CASM functional response: microzooplankton consumption rate too high!"
     stop
   ]
   let C_detritus W * ww_Cmax_juv * (juv_pref_detritus * juv_ass_detritus * juv_h_detritus * w_detritus * cfa) / (1 + sum_rel_C) ; (SedimentPOC)
   if C_detritus > W * ww_Cmax_juv [
     user-message "CASM functional response: detritus consumption rate too high!"
     stop
   ]

   if (C_cladocerans < 0 or C_copepods < 0 or C_rotifers < 0 or C_microzoo < 0 or C_detritus < 0) [ ; (SedimentPOC)
     user-message "WARNING: negative consumption by juvenile"
   ]
   ; calculate actual consumption (by a fish), and reduce prey guild biomasses accordingly
   ifelse C_copepods > w_copepods [
     set C_copepods w_copepods
     set w_copepods 0
   ] [
     set w_copepods w_copepods - C_copepods
   ]
   ifelse C_cladocerans > w_cladocerans [
     set C_cladocerans w_cladocerans
     set w_cladocerans 0
   ] [
     set w_cladocerans w_cladocerans - C_cladocerans
   ]
   ifelse C_rotifers > w_rotifers [
     set C_rotifers w_rotifers
     set w_rotifers 0
   ] [
     set w_rotifers w_rotifers - C_rotifers
   ]
   ifelse C_microzoo > w_microzoo [
     set C_microzoo w_microzoo
     set w_microzoo 0
   ] [
     set w_microzoo w_microzoo - C_microzoo
   ]
   ifelse C_detritus > w_detritus [ ; (SedimentPOC)
     set C_detritus w_detritus
     set w_detritus 0
   ] [
     set w_detritus w_detritus - C_detritus
   ]

  ; daily consumption of the juvenile (in J):
  set day_consumption (C_copepods * ED_copepods) + (C_cladocerans * ED_cladocerans) + (C_rotifers * ED_rotifers)
                    + (C_microzoo * ED_microzoo) + (C_detritus * ED_detritus)

  ] [ ; consumption of adult TS
    set ww_Cmax_adt Cmax_t

  ; Consumption per prey type dependent on prey biomass
    let sum_rel_C adt_pref_copepods * adt_ass_copepods * adt_h_copepods * w_copepods * cfv
                + adt_pref_cladocerans * adt_ass_cladocerans * adt_h_cladocerans * w_cladocerans * cfv
                + adt_pref_mayflies * adt_ass_mayflies * adt_h_mayflies * w_mayflies * cfa ;  (Ephemeroptera)
                + adt_pref_caddisflies * adt_ass_caddisflies * adt_h_caddisflies * w_caddisflies * cfa ; (Trichoptera)
                + adt_pref_oligochaetes * adt_ass_oligochaetes * adt_h_oligochaetes * w_oligochaetes * cfa
                + adt_pref_chironomids * adt_ass_chironomids * adt_h_chironomids * w_chironomids * cfa
                + adt_pref_detritus * adt_ass_detritus * adt_h_detritus * w_detritus * cfa ; (SedimentPOC) ; assumption that fish can access less area of detritus for daily consumption
                + adt_pref_periphyton * adt_ass_periphyton * adt_h_periphyton * w_periphyton * cfv

  ; consumption rates in g prey/ g fish and day
    let C_copepods W * ww_Cmax_adt * (adt_pref_copepods * adt_ass_copepods * adt_h_copepods * w_copepods * cfv) / (1 + sum_rel_C)
    if C_copepods > W * ww_Cmax_adt [
      user-message "CASM functional response: copepod consumption rate too high!"
      stop
    ]
    let C_cladocerans W * ww_Cmax_adt * (adt_pref_cladocerans * adt_ass_cladocerans * adt_h_cladocerans * w_cladocerans * cfv) / (1 + sum_rel_C)
    if C_cladocerans > W * ww_Cmax_adt [
      user-message "CASM functional response: cladoceran consumption rate too high!"
      stop
    ]
    let C_mayflies W * ww_Cmax_adt * (adt_pref_mayflies * adt_ass_mayflies * adt_h_mayflies * w_mayflies * cfa) / (1 + sum_rel_C) ;  (Ephemeroptera)
    if C_mayflies > W * ww_Cmax_adt [
      user-message "CASM functional response: mayfly consumption rate too high!"
      stop
    ]
    let C_caddisflies W * ww_Cmax_adt * (adt_pref_caddisflies * adt_ass_caddisflies * adt_h_caddisflies * w_caddisflies * cfa) / (1 + sum_rel_C) ; (Trichoptera)
    if C_caddisflies > W * ww_Cmax_adt [
    user-message "CASM functional response: caddisfly consumption rate too high!"
      stop
    ]
    let C_oligochaetes W * ww_Cmax_adt * (adt_pref_oligochaetes * adt_ass_oligochaetes * adt_h_oligochaetes * w_oligochaetes * cfa) / (1 + sum_rel_C)
    if C_oligochaetes > W * ww_Cmax_adt [
      user-message "CASM functional response: oligochaete consumption rate too high!"
      stop
    ]
    let C_chironomids W * ww_Cmax_adt * (adt_pref_chironomids * adt_ass_chironomids * adt_h_chironomids * w_chironomids * cfa) / (1 + sum_rel_C)
    if C_chironomids > W * ww_Cmax_adt [
      user-message "CASM functional response: chironomid consumption rate too high!"
    stop
    ]
    let C_detritus W * ww_Cmax_adt * (adt_pref_detritus * adt_ass_detritus * adt_h_detritus * w_detritus * cfa) / (1 + sum_rel_C) ; (SedimentPOC)
    if C_detritus > W * ww_Cmax_adt [
      user-message "CASM functional response: detritus consumption rate too high!"
      stop
    ]
    let C_periphyton W * ww_Cmax_adt * (adt_pref_periphyton * adt_ass_periphyton * adt_h_periphyton * w_periphyton * cfv) / (1 + sum_rel_C)
    if C_periphyton > W * ww_Cmax_adt [
      user-message "CASM functional response: periphyton consumption rate too high!"
      stop
    ]

  if (C_copepods < 0 or C_cladocerans < 0 or C_mayflies < 0 or C_caddisflies < 0 or C_oligochaetes < 0 or C_chironomids < 0 or C_detritus < 0 or C_periphyton < 0) [
    user-message "WARNING: negative consumption by juveniles"
  ]
  ; calculate actual consumption (by a fish), and reduce prey guild biomasses accordingly
    ifelse C_copepods > w_copepods [
      set C_copepods w_copepods
      set w_copepods 0
    ] [
      set w_copepods w_copepods - C_copepods
    ]
    ifelse C_cladocerans > w_cladocerans [
      set C_cladocerans w_cladocerans
      set w_cladocerans 0
    ] [
      set w_cladocerans w_cladocerans - C_cladocerans
    ]
    ifelse C_mayflies > w_mayflies [ ;  (Ephemeroptera)
      set C_mayflies w_mayflies
      set w_mayflies 0
    ] [
      set w_mayflies w_mayflies - C_mayflies
    ]
    ifelse C_caddisflies > w_caddisflies [ ; (Trichoptera)
      set C_caddisflies w_caddisflies
      set w_caddisflies 0
    ] [
      set w_caddisflies w_caddisflies - C_caddisflies
    ]
    ifelse C_oligochaetes > w_oligochaetes [
      set C_oligochaetes w_oligochaetes
      set w_oligochaetes 0
    ] [
      set w_oligochaetes w_oligochaetes - C_oligochaetes
    ]
    ifelse C_chironomids > w_chironomids [
      set C_chironomids w_chironomids
      set w_chironomids 0
    ] [
      set w_chironomids w_chironomids - C_chironomids
    ]
    ifelse C_detritus > w_detritus [ ; (SedimentPOC)
      set C_detritus w_detritus
      set w_detritus 0
    ] [
      set w_detritus w_detritus - C_detritus
    ]
    ifelse C_periphyton > w_periphyton [
      set C_periphyton w_periphyton
      set w_periphyton 0
    ] [
      set w_periphyton w_periphyton - C_periphyton
    ]

    ; daily consumption by adult TS (in J)
    set day_consumption (C_copepods * ED_copepods) + (C_cladocerans * ED_cladocerans) + (C_mayflies * ED_mayflies) + (C_caddisflies * ED_caddisflies)
                      + (C_oligochaetes * ED_oligochaetes) + (C_chironomids * ED_chironomids) + (C_detritus * ED_detritus) + (C_periphyton * ED_periphyton)
  ]

end

to grow
  ; procedure implements growth according to fish bioenergetics
  ; consumption has been determined in functional response procedure

  ; respiration in J / g wet weight TS (applying oxycalorific coefficient)
  set R_T (RA * (W ^ RB) * fR_T) * OxycalCoeff
  ; specific dynamic action in J / g wet weight TS (note: day_consumption is in J per TS)
  set S SDA * (day_consumption / W) * (1 - FA)
  ; waste losses in J / g wet weight TS
  set E_C (FA * (day_consumption / W)) + (UA * ((day_consumption / W) - (FA * (day_consumption / W))))

  ; Resulting growth
  ifelse stage = "larva" or stage = "juvenile" [
    set d_W W * (((day_consumption / W) - (R_T + S) - E_C) / ED_TS_juv)
  ] [
    set d_W W * (((day_consumption / W) - (R_T + S) - E_C) / ED_TS_adt)
  ]

  set W W + d_W ; apply weight change

  if stage != "larva" [
  let tmp_SL (calc_length W)
  if tmp_SL > SL [ set SL tmp_SL ] ; fish can only increase in length, not decrease
  set K W / (calc_weight SL)       ; calculate fish condition
  ]
  ; stage transition according to fish size
  ask turtles with [ stage = "juvenile" ] [
    if SL >= MinSizeSubadult [ ; fish attain subadult pigmentation at SL ~ 15 mm and ~ 2 months old
      set stage "subadult"
      set size 3
    ]
  ]
  ask turtles with [ stage = "subadult" ] [
    if age > 335 and SL >= MinSizeMat [ ; fish reach maturity at 11 month or later as soon as they reach size of maturity
      set stage "adult"
      set size 4
    ]
  ]
  ask turtles with [ stage = "adult" ] [
   if age / 365 > 2 [ set size 4 ]
  ]

  set age age + 1 ; fish age in days

end

to spawn
  ; calculate clutch size dependent on female weight
  let n_eggs (exp (beta_0 + beta_1 * W))

  hatch n_eggs [ ;; new TS individuals are created
    set age 0
    set stage "egg"
    ; assuming sex ratio of 1:1, and randomly assign sex:
    let rnd_num random 2
    ifelse rnd_num = 0
      [ set sex "female" ]
      [ set sex "male" ]
    set W EggWeight
    set SL 0
;    set starvation FALSE
;    set preStarvW W
    set K 1
    set Kmin random-normal AverageKmin (AverageKmin / 10)
    set day_consumption 0
    set R_T 0
    set E_C 0
    set d_W 0
    set spawn_date -1
    setxy random-xcor random-ycor
    set shape "dot"
    set size 1
  ]

  ; reduce female's weight according to clutch produced:
;  set starvation TRUE ; assuming that egg laying results in a weight loss that has to be made up for
;  set preStarvW W
  set W W - (n_eggs * EggWeight)
  set K W / (calc_weight SL)
end

to-report consumption_temperature_function [ tempC ]
  ; temperature dependence of consumption
  ; Thornton and Lessem 1978; Fish Bioenegetics 3.0 (Hanson et al. 1997), consumption equation 3
  let tt5 1 / (te2 - te1)
  let tt7 1 / (te4 - te3)
  let t5 tt5 * ln (xk2 * (1 - xk1) / (0.02 * xk1))
  let t4 exp (t5 * (tempC - te1))
  let t7 tt7 * ln (xk3 * (1 - xk4) / (0.02 * xk4))
  let t6 exp (t7 * (te4 - tempC))
  let gcta (xk1 * t4) / (1 + xk1 * (t4 - 1))
  let gctb (xk4 * t6) / (1 + xk4 * (t6 - 1))
  report gcta * gctb
end

to-report respiration_temperature_function [ tempC ]
  ; temperature dependence of respiration (from fish bioenergetics)
  let V (RTM - tempC) / (RTM - RTO)
  let Z (ln RQ) * (RTM - RTO)
  let Y (ln RQ) * (RTM - RTO + 2)
  let X ((Z ^ 2) * ((1 + sqrt (1 + (40 / Y))) ^ 2)) / 400
  report (V ^ X) * (exp (X * (1 - V)))
end

to-report calc_weight [ fish_SL ]
  ; calculate TS weight from standard length (SL) (from Kerns and Bonneau 2002)
  let log_weight l2w_a * (log fish_SL 10) + l2w_b
  report 10 ^ log_weight
end

to-report calc_length [ fish_weight ]
  ; calculate TS standard length (SL) from weight (from Kerns and Bonneau 2002)
  let log_length ((log fish_weight 10) - l2w_b) / l2w_a
  report 10 ^ log_length
end

to-report triang_dist [a b c]
  ; triangular distribution applied to spawning dates (day of year):
  ; a = minimum (start spawning)
  ; b = maximum (stop spawning)
  ; c = mode (most common value; peak spawning)
	let x random-float 1
	let fc (c - a)/(b - a)
  ifelse x < fc [
    set x (a + sqrt(x * (b - a) * (c - a)))
  ] [
    set	x (b - sqrt((1 - x) * (b - a) * (b - c)))
  ]
  report (round x)
end

to setDailyPreyBiomass
  ; guild biomasses given per waterbody area in CASM in gC/m2
  ; keeping track of the initial biomasses at the start of the day
  set preyatthestartoftheday []
  set preyatthestartoftheday lput (item (ticks mod 365) list_w_copepods) preyatthestartoftheday
  set preyatthestartoftheday lput (item (ticks mod 365) list_w_cladocerans) preyatthestartoftheday
  set preyatthestartoftheday lput (item (ticks mod 365) list_w_rotifers) preyatthestartoftheday
  set preyatthestartoftheday lput (item (ticks mod 365) list_w_microzoo) preyatthestartoftheday
  set preyatthestartoftheday lput (item (ticks mod 365) list_w_mayflies) preyatthestartoftheday ;  (Ephemeroptera)
  set preyatthestartoftheday lput (item (ticks mod 365) list_w_caddisflies) preyatthestartoftheday ; (Trichoptera)
  set preyatthestartoftheday lput (item (ticks mod 365) list_w_oligochaetes) preyatthestartoftheday
  set preyatthestartoftheday lput (item (ticks mod 365) list_w_chironomids) preyatthestartoftheday
  set preyatthestartoftheday lput (item (ticks mod 365) list_w_detritus) preyatthestartoftheday ; (SedimentPOC)
  set preyatthestartoftheday lput (item (ticks mod 365) list_w_periphyton) preyatthestartoftheday

  ; setting prey biomasses (in g) to initial biomasses of the day (these values will be altered by TS consumption)
  set w_copepods ((item (ticks mod 365) list_w_copepods) * gC_ww_copepods * pondArea)
  set w_cladocerans ((item (ticks mod 365) list_w_cladocerans) * gC_ww_cladocerans * pondArea)
  set w_rotifers ((item (ticks mod 365) list_w_rotifers) * gC_ww_rotifers * pondArea)
  set w_microzoo ((item (ticks mod 365) list_w_microzoo) * gC_ww_microzoo * pondArea)
  set w_mayflies ((item (ticks mod 365) list_w_mayflies) * gC_ww_mayflies * pondArea) ;  (Ephemeroptera)
  set w_caddisflies ((item (ticks mod 365) list_w_caddisflies) * gC_ww_caddisflies * pondArea) ; (Trichoptera)
  set w_oligochaetes ((item (ticks mod 365) list_w_oligochaetes) * gC_ww_oligochaetes * pondArea)
  set w_chironomids ((item (ticks mod 365) list_w_chironomids) * gC_ww_chironomids * pondArea)
  ;; 27 Nov 2018: detritus limit applied here rather than upon read-in of file
;  set w_detritus ((item (ticks mod 365) list_w_detritus) * gC_ww_detritus * pondArea) ; (SedimentPOC)
  let w_detritus_tmp (item (ticks mod 365) list_w_detritus) ; (SedimentPOC)
  if w_detritus_tmp > limit_detritus [
        set w_detritus_tmp limit_detritus
       ;; 27 Nov 2018: flag included to make sure that detritus biomass value returned to CASM / written to output is not truncated by limit_detritus
        set flag_limit_detritus TRUE
      ]
  set w_detritus w_detritus_tmp * gC_ww_detritus * pondArea
  set w_periphyton ((item (ticks mod 365) list_w_periphyton) * gC_ww_periphyton * pondArea)
  ; TS eaten is given as g (wet weight)/m3 in CASM ouput file:
  set w_TS_eaten (((item (ticks mod 365) list_TS_eaten) / CASM_pool_area ) * pondVolume ); biomass of TS eaten from the whole population (in the pond); note that TS eaten in listed in g/m3 in CASM input file, not gC/m3 like the other biomasses

end

to  readInputParameters
  ; read constant model parameters from input file
  ; NOTE: parameters have to be listed in the given order in the input file !!!
  ifelse file-exists? InputParameters [
    file-open InputParameters ; file with input parameters
    let temp file-read-line ; first line is a header
    set temp file-read
    set juv_pref_copepods 	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read juv_pref_copepods" ]
    set temp file-read
    set juv_pref_cladocerans	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read juv_pref_cladocerans" ]
    set temp file-read
    set juv_pref_rotifers	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read juv_pref_rotifers" ]
    set temp file-read
    set juv_pref_microzoo	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read juv_pref_microzoo" ]
    set temp file-read
    set juv_pref_detritus	file-read ; (SedimentPOC)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read juv_pref_detritus" ]
    set temp file-read
    set juv_ass_copepods	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read juv_ass_copepods" ]
    set temp file-read
    set juv_ass_cladocerans	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read juv_ass_cladocerans" ]
    set temp file-read
    set juv_ass_rotifers	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read juv_ass_rotifers" ]
    set temp file-read
    set juv_ass_microzoo	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read juv_ass_microzoo" ]
    set temp file-read
    set juv_ass_detritus	file-read ; (SedimentPOC)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read juv_ass_detritus" ]
    set temp file-read
    set juv_h_copepods	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read juv_h_copepods" ]
    set temp file-read
    set juv_h_cladocerans	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read juv_h_cladocerans" ]
    set temp file-read
    set juv_h_rotifers	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read juv_h_rotifers" ]
    set temp file-read
    set juv_h_microzoo	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read juv_h_microzoo" ]
    set temp file-read
    set juv_h_detritus	file-read ; (SedimentPOC)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read juv_h_detritus" ]
    set temp file-read
    set adt_pref_copepods 	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_pref_copepods" ]
    set temp file-read
    set adt_pref_cladocerans	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_pref_cladocerans	" ]
    set temp file-read
    set adt_pref_mayflies	file-read  ;  (Ephemeroptera)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_pref_mayflies" ]
    set temp file-read
    set adt_pref_caddisflies	file-read ; (Trichoptera)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_pref_caddisflies" ]
    set temp file-read
    set adt_pref_oligochaetes	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_pref_oligochaetes" ]
    set temp file-read
    set adt_pref_chironomids	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_pref_chironomids" ]
    set temp file-read
    set adt_pref_detritus	file-read ; (SedimentPOC)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_pref_detritus" ]
    set temp file-read
    set adt_pref_periphyton	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_pref_periphyton" ]
    set temp file-read
    set adt_ass_copepods	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_ass_copepods" ]
    set temp file-read
    set adt_ass_cladocerans	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_ass_cladocerans" ]
    set temp file-read
    set adt_ass_mayflies	file-read  ;  (Ephemeroptera)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_ass_mayflies" ]
    set temp file-read
    set adt_ass_caddisflies file-read ; (Trichoptera)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_ass_caddisflies" ]
    set temp file-read
    set adt_ass_oligochaetes	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_ass_oligochaetes" ]
    set temp file-read
    set adt_ass_chironomids	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_ass_chironomids" ]
    set temp file-read
    set adt_ass_detritus	file-read ; (SedimentPOC)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_ass_detritus" ]
    set temp file-read
    set adt_ass_periphyton	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_ass_periphyton	" ]
    set temp file-read
    set adt_h_copepods	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_h_copepods" ]
    set temp file-read
    set adt_h_cladocerans	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_h_cladocerans" ]
    set temp file-read
    set adt_h_mayflies	file-read ;  (Ephemeroptera)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_h_mayflies" ]
    set temp file-read
    set adt_h_caddisflies	file-read ; (Trichoptera)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_h_caddisflies" ]
    set temp file-read
    set adt_h_oligochaetes	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_h_oligochaetes" ]
    set temp file-read
    set adt_h_chironomids	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_h_chironomids" ]
    set temp file-read
    set adt_h_detritus	file-read ; (SedimentPOC)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_h_detritus" ]
    set temp file-read
    set adt_h_periphyton	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read adt_h_periphyton" ]
    set temp file-read
    set a_cmax	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read a_cmax" ]
    set temp file-read
    set b_cmax	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read b_cmax" ]
    set temp file-read
    set te1	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read te1" ]
    set temp file-read
    set te2	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read te2" ]
    set temp file-read
    set te3	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read te3" ]
    set temp file-read
    set te4	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read te4" ]
    set temp file-read
    set xk1	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read xk1" ]
    set temp file-read
    set xk2	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read xk2" ]
    set temp file-read
    set xk3	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read xk3" ]
    set temp file-read
    set xk4	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read xk4" ]
    set temp file-read
    set SDA	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read SDA" ]
    set temp file-read
    set RA	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read RA" ]
    set temp file-read
    set RB	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read RB" ]
    set temp file-read
    set RTO	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read RTO" ]
    set temp file-read
    set RTM	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read RTM" ]
    set temp file-read
    set RQ	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read RQ" ]
    set temp file-read
    set FA	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read FA" ]
    set temp file-read
    set UA	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read UA" ]
    set temp file-read
    set beta_0 file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read X1" ]
    set temp file-read
    set beta_1 file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read X2" ]
    set temp file-read
    set l2w_a	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read l2w_a" ]
    set temp file-read
    set l2w_b	file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read l2w_b" ]
    set temp file-read
    set gC_ww_copepods file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read gC_ww_copepods" ]
    set temp file-read
    set gC_ww_cladocerans file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read gC_ww_cladocerans" ]
    set temp file-read
    set gC_ww_rotifers file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read gC_ww_rotifers" ]
    set temp file-read
    set gC_ww_microzoo file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read gC_ww_microzoo" ]
    set temp file-read
    set gC_ww_mayflies file-read ;  (Ephemeroptera)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read gC_ww_mayflies" ]
    set temp file-read
    set gC_ww_caddisflies file-read ; (Trichoptera)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read gC_ww_caddisflies" ]
    set temp file-read
    set gC_ww_oligochaetes file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read gC_ww_oligochaetes" ]
    set temp file-read
    set gC_ww_chironomids file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read gC_ww_chironomids" ]
    set temp file-read
    set gC_ww_detritus file-read ; (SedimentPOC)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read gC_ww_detritus" ]
    set temp file-read
    set gC_ww_periphyton file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read gC_ww_periphyton" ]
    set temp file-read
    set gC_ww_fish file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read gC_ww_fish" ]
    set temp file-read
    set ED_copepods file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read ED_copepods" ]
    set temp file-read
    set ED_cladocerans file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read ED_cladocerans" ]
    set temp file-read
    set ED_rotifers file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read ED_rotifers" ]
    set temp file-read
    set ED_microzoo file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read ED_microzoo" ]
    set temp file-read
    set ED_mayflies file-read ;  (Ephemeroptera)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read ED_mayflies" ]
    set temp file-read
    set ED_caddisflies file-read ; (Trichoptera)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read ED_caddisflies" ]
    set temp file-read
    set ED_oligochaetes file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read ED_oligochaetes" ]
    set temp file-read
    set ED_chironomids file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read ED_chironomids" ]
    set temp file-read
    set ED_detritus file-read ; (SedimentPOC)
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read ED_detritus" ]
    set temp file-read
    set ED_periphyton file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read ED_periphyton" ]
    set temp file-read
    set ED_TS_juv file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read ED_TS_juv" ]
    set temp file-read
    set ED_TS_adt file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read ED_TS_adt" ]
    set temp file-read
    set OxycalCoeff file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read OxycalCoeff" ]
    set temp file-read
    set DailySurvival_Egg_Larva file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read DailySurvival_Egg_Larva" ]
    set temp file-read
    set DailySurvival_Egg_Larva_0 file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read DailySurvival_Egg_Larva_0" ]
    set temp file-read
    set beta_dd_egg_survival file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read beta_dd_egg_survival" ]
    set temp file-read
    set SunfishFactor file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read SunfishFactor" ]
    set temp file-read
    set daily_surv_rate_juv file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read daily_surv_rate_juv" ]
    set temp file-read
    set daily_surv_rate_subadult file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read daily_surv_rate_subadult" ]
    set temp file-read
    set daily_surv_rate_yr2 file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read daily_surv_rate_yr2" ]
    set temp file-read
    set daily_surv_rate_yr3 file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read daily_surv_rate_yr3" ]
    set temp file-read
    set StartSpawning_early file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read StartSpawning_early" ]
    set temp file-read
    set StartSpawning_late file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read StartSpawning_late" ]
    set temp file-read
    set PeakSpawning_early file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read PeakSpawning_early" ]
    set temp file-read
    set PeakSpawning_late file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read PeakSpawning_late" ]
    set temp file-read
    set StopSpawning_early file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read StopSpawning_early" ]
    set temp file-read
    set StopSpawning_late file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read StopSpawning_late" ]
    set temp file-read
    set EggWeight file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read EggWeight" ]
    set temp file-read
    set EggDevTime file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read EggDevTime" ]
    set temp file-read
    set YolkLasting file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read YolkLasting" ]
    set temp file-read
    set MinSizeSubadult file-read
    if file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read MinSizeSubadult" ]
    set temp file-read
    set MinSizeMat file-read
    if not file-at-end? [ user-message "WARNING: parameter input file does not have the right length! Read MinSizeMat" ]
    file-close
  ] [ user-message "Parameter input file does not exist in current directory!" ]
end

to readWaterConditions
  ; reading waterbody conditions (temperature and water depth) from yearly input file (input file is expected to have 365 data rows)
  set TemperatureList []
  set WaterDepthList []

  ifelse file-exists? WaterConditionInput [
    file-open WaterConditionInput
    let entry ""
    let cnt 0
    while [ cnt < 7 ] [ ; top 7 lines of file are headers
      set entry file-read-line
;      print entry ; for testing
      set cnt cnt + 1
    ]
    set cnt 0
    let year 365
    while [ cnt < 365 ] [ ; note that the file ends in two lines reporting averages, not daily data
      set entry file-read ; temperature is listed in the second column in the file
      set TemperatureList lput file-read TemperatureList
      let item_cnt 0
      while [ item_cnt < 5 ] [ ; water depth is listed in the 8th column (5 columns skipped inbetween)
        set entry file-read
        set item_cnt item_cnt + 1
      ]
      set WaterDepthList lput file-read WaterDepthList
      set entry file-read-line
      set cnt cnt + 1
    ]
    file-close

  ] [ user-message "Waterbody condition input file does not exist in current directory!" ]
  if (length TemperatureList) > 365 [
    user-message "Temperature input larger than 365 entries: only first 365 entries used."
    set TemperatureList (sublist TemperatureList 0 365)
  ]
  if length TemperatureList < 365 [
    user-message "Temperature input less than 365 entries: input file cannot be used."
    stop
  ]
  if (length WaterDepthList) > 365 [
    user-message "Water depth input larger than 365 entries: only first 365 entries used."
    set WaterDepthList (sublist WaterDepthList 0 365)
  ]
  if length WaterDepthList < 365 [
    user-message "Water depth input less than 365 entries: input file cannot be used."
    stop
  ]
end

to readDailyPreyBiomass_year
  ; read in yearly prey biomasses (input file is expected to have 365 data rows)

  set list_w_copepods []
  set list_w_cladocerans []
  set list_w_rotifers []
  set list_w_microzoo []
  set list_w_mayflies []  ;  (Ephemeroptera)
  set list_w_caddisflies [] ; (Trichoptera)
  set list_w_oligochaetes []
  set list_w_chironomids []
  set list_w_detritus [] ; (SedimentPOC)
  set list_w_periphyton []
  set list_w_TS_juv []
  set list_w_TS_adt []
  set list_TS_eaten []

  ifelse file-exists? DailyPreyBiomass [
    file-open DailyPreyBiomass
    let entry ""
    let cnt 0
    while [ cnt < 2 ] [ ; top 2 lines of file are file info
      set entry file-read-line
;      print entry ; for testing
      set cnt cnt + 1
    ]
    ;; the following section checks whether the needed columns are present in the file
    ;; NOTE: columns are subsequently read in by column number, not name in the header (see below)!!!
    set entry file-read-line
    let flag_prey false
    let flag_TS_eaten false
;    type "TSgww-eaten in file header: " print member? "TSgww-eaten" entry ; for testing
    if (not member? "Copepods" entry) [ set flag_prey true ]
;    type "copepods - flag_prey: " print flag_prey ; for testing
    if (not member? "Cladocerans" entry) [ set flag_prey true ] ; it has this typo in the CASM out file
;    type "Cladocerans - flag_prey: " print flag_prey ; for testing
    if (not member? "Rotifers" entry) [ set flag_prey true ]
;    type "Rotifers - flag_prey: " print flag_prey ; for testing
    if (not member? "MicroZplankton" entry) [ set flag_prey true ]
;    type "MicroZplankton - flag_prey: " print flag_prey ; for testing
    if (not member? "Ephemeroptera" entry) [ set flag_prey true ]
;    type "Ephemeroptera - flag_prey: " print flag_prey ; for testing
    if (not member? "Trichoptera" entry) [ set flag_prey true ]
;    type "Trichoptera - flag_prey: " print flag_prey ; for testing
    if (not member? "Oligochaetes" entry) [ set flag_prey true ]
 ;   type "Oligochaetes - flag_prey: " print flag_prey ; for testing
    if (not member? "Chironomids" entry) [ set flag_prey true ]
;    type "Chironomids - flag_prey: " print flag_prey ; for testing
    if (not member? "SedimentPOC" entry) [ set flag_prey true ]
 ;   type "SedimentPOC - flag_prey: " print flag_prey ; for testing
    if (not member? "TotalPeriphyt" entry) [ set flag_prey true ]
;    type "TotalPeriphyt - flag_prey: " print flag_prey ; for testing
    if (not member? "TSgww-eaten" entry) [
      set flag_TS_eaten true
      if predation [
        ifelse (user-yes-or-no? "WARNING: no TS - eaten (by predator) listed in daily prey biomass file: continue?")
        [ user-message "Continuing without predation." ]
        [ file-close
          user-message "Run 'setup' again."
          stop
        ]
      ]
    ]
    if (flag_prey) [ user-message "WARNING: prey biomass input from file incomplete!" ]
;    type "flag_TS_eaten: " print flag_TS_eaten ; for testing
    set cnt 0
    while [not file-at-end?] [
      let item_cnt 0
      while [ item_cnt < 25 ] [
        set entry file-read
;        type "entry read from file: " print entry ; for testing
        set item_cnt item_cnt + 1
      ]
      set list_w_copepods lput (file-read) list_w_copepods ; column 26 in input file
      set item_cnt item_cnt + 1
      set list_w_cladocerans lput (file-read) list_w_cladocerans ; column 27 in input file
      set item_cnt item_cnt + 1
      set list_w_rotifers lput (file-read) list_w_rotifers ; column 28 in input file
      set item_cnt item_cnt + 1
      set list_w_microzoo lput (file-read) list_w_microzoo ; column 29 in input file
      set item_cnt item_cnt + 1
      while [ item_cnt < 32 ] [
        set entry file-read
        set item_cnt item_cnt + 1
      ]
      set list_w_TS_juv lput (file-read * gC_ww_fish * pondArea) list_w_TS_juv ; column 33 in input file
      set item_cnt item_cnt + 1
      set list_w_TS_adt lput (file-read * gC_ww_fish * pondArea) list_w_TS_adt ; column 34 in input file
      set item_cnt item_cnt + 1
      while [ item_cnt < 37 ] [
        set entry file-read
        set item_cnt item_cnt + 1
      ]
      set list_w_mayflies lput (file-read) list_w_mayflies ; column 38 in input file ;  (Ephemeroptera)
      set item_cnt item_cnt + 1
      set list_w_caddisflies lput (file-read) list_w_caddisflies ; column 39 in input file ; (Trichoptera)
      set item_cnt item_cnt + 1
      set list_w_oligochaetes lput (file-read) list_w_oligochaetes ; column 40 in input file
      set item_cnt item_cnt + 1
      set list_w_chironomids lput (file-read) list_w_chironomids ; column 41 in input file
      set item_cnt item_cnt + 1
      while [ item_cnt < 50 ] [
        set entry file-read
        set item_cnt item_cnt + 1
      ]
      let detritus_tmp file-read
;      if detritus_tmp > limit_detritus [
;        set detritus_tmp limit_detritus
       ;; 27 Nov 2018: flag included to make sure that detritus biomass value returned to CASM / written to output is not truncated by limit_detritus
;        set flag_limit_detritus TRUE
;      ]
      set list_w_detritus lput detritus_tmp list_w_detritus ; column 51 in input file ("Sediment POC")
      set item_cnt item_cnt + 1
      set entry file-read
      set item_cnt item_cnt + 1
      set list_w_periphyton lput (file-read) list_w_periphyton ; column 53 in input file ("Total periphyton")
      ifelse flag_TS_eaten [
        set list_TS_eaten lput 0 list_TS_eaten ; if TS eaten is not listed in the prey biomass input file, it is assumed that no predation occurs
        set entry file-read-line
      ] [
        while [ item_cnt < 71 ] [ ; it is unclear to me what's happening here, but I get the right column with this (even though TSgww-eaten is in copumn 73)
          set entry file-read
          set item_cnt item_cnt + 1
        ]
        set list_TS_eaten lput (file-read) list_TS_eaten ; column 73 in input file (TSgww - eaten)
      ]
      set cnt cnt + 1
    ]
;    print cnt ; for testing
  ] [ user-message "Daily prey biomass input file does not exist in current directory!" ]
  file-close
  if ( (length list_w_copepods) != 365 or (length list_w_cladocerans) != 365
       or (length list_w_rotifers) != 365 or (length list_w_microzoo) != 365
       or (length list_w_mayflies) != 365 or (length list_w_caddisflies) != 365
       or (length list_w_oligochaetes) != 365 or (length list_w_chironomids) != 365
       or (length list_w_detritus) != 365 or (length list_w_periphyton) != 365) [
    user-message "WARNING: daily prey biomass input does not have 365 entries!"
    stop
  ]
end

to readDailyPreyBiomass_day
  ; keeping track of the initial biomasses at the start of the day
  set preyatthestartoftheday []

  ifelse file-exists? DailyPreyBiomass [
    file-open DailyPreyBiomass
    let entry ""
    let cnt 0
    while [ cnt < 2 ] [ ; top 2 lines of file are file info
      set entry file-read-line
;      print entry
      set cnt cnt + 1
    ]
    ;; the following section checks whether the needed columns are present in the file
    ;; NOTE: columns are subsequently read in by column number, not name in the header (see below)!!!
    set entry file-read-line
    let flag_prey false
    let flag_TS_eaten false
    if (not member? "Copepods" entry) [ set flag_prey true ]
    if (not member? "Cladocerans" entry) [ set flag_prey true ] ; it has this typo in the CASM out file
    if (not member? "Rotifers" entry) [ set flag_prey true ]
    if (not member? "MicroZplankton" entry) [ set flag_prey true ]
    if (not member? "Ephemeroptera" entry) [ set flag_prey true ]
    if (not member? "Trichoptera" entry) [ set flag_prey true ]
    if (not member? "Oligochaetes" entry) [ set flag_prey true ]
    if (not member? "Chironomids" entry) [ set flag_prey true ]
    if (not member? "SedimentPOC" entry) [ set flag_prey true ]
    if (not member? "TotalPeriphyt" entry) [ set flag_prey true ]
    if (not member? "TSgww-eaten" entry) [
      set flag_TS_eaten true
      if predation [
        ifelse (user-yes-or-no? "WARNING: no TS - eaten (by predator) listed in daily prey biomass file: continue?")
        [ user-message "Continuing without predation." ]
        [ file-close
          user-message "Run 'setup' again."
          stop
        ]
      ]
    ]
    if (flag_prey) [ user-message "WARNING: prey biomass input from file incomplete!" ]
;    type "flag_TS_eaten: " print flag_TS_eaten ; for testing

    let item_cnt 0
    let currentDay file-read
    if currentDay != ticks + 1 [
      print "Day listed in daily prey biomass input file does not match with current tick: tick is reset to input."
      reset-ticks
      tick-advance currentDay - 1
    ]
    set item_cnt item_cnt + 1
    while [ item_cnt < 23 ] [
      set entry file-read
      set item_cnt item_cnt + 1
    ]
    set entry file-read
    set preyatthestartoftheday lput entry preyatthestartoftheday
    set w_copepods (entry * gC_ww_copepods * pondArea) ; column 24 in input file
    set item_cnt item_cnt + 1
    set entry file-read
    set preyatthestartoftheday lput entry preyatthestartoftheday
    set w_cladocerans (entry * gC_ww_cladocerans * pondArea) ; column 25 in input file
    set item_cnt item_cnt + 1
    set entry file-read
    set preyatthestartoftheday lput entry preyatthestartoftheday
    set w_rotifers (entry * gC_ww_rotifers * pondArea) ; column 26 in input file
    set item_cnt item_cnt + 1
    set entry file-read
    set preyatthestartoftheday lput entry preyatthestartoftheday
    set w_microzoo (entry * gC_ww_microzoo * pondArea) ; column 27 in input file
    set item_cnt item_cnt + 1
    while [ item_cnt < 30 ] [
      set entry file-read
      set item_cnt item_cnt + 1
    ]
    set ini_w_TS_juv (file-read * gC_ww_fish * pondArea) ; column 31 in input file
    set item_cnt item_cnt + 1
    set ini_w_TS_adt (file-read * gC_ww_fish * pondArea) ; column 32 in input file
    set item_cnt item_cnt + 1
    while [ item_cnt < 35 ] [
      set entry file-read
      set item_cnt item_cnt + 1
    ]
    set entry file-read
    set preyatthestartoftheday lput entry preyatthestartoftheday
    set w_mayflies (entry * gC_ww_mayflies * pondArea) ; column 36 in input file ;  (Ephemeroptera)
    set item_cnt item_cnt + 1
    set entry file-read
    set preyatthestartoftheday lput entry preyatthestartoftheday
    set w_caddisflies (entry * gC_ww_caddisflies * pondArea) ; column 37 in input file ; (Trichoptera)
    set item_cnt item_cnt + 1
    set entry file-read
    set preyatthestartoftheday lput entry preyatthestartoftheday
    set w_oligochaetes (entry * gC_ww_oligochaetes * pondArea) ; column 38 in input file
    set item_cnt item_cnt + 1
    set entry file-read
    set preyatthestartoftheday lput entry preyatthestartoftheday
    set w_chironomids (entry * gC_ww_chironomids * pondArea); column 39 in input file
    set item_cnt item_cnt + 1
    while [ item_cnt < 48 ] [
      set entry file-read
      set item_cnt item_cnt + 1
    ]
    set entry file-read
;    type "Sediment POC entry: " type entry type ", limit: " print limit_detritus
;    if entry > limit_detritus [ set entry limit_detritus ]
;    type "Sediment POC entry: " print entry
    set preyatthestartoftheday lput entry preyatthestartoftheday
;    type "Sediment POC entry: " print entry
               ;; 27 Nov 2018: flag included to make sure that detritus biomass value returned to CASM / written to output is not truncated by limit_detritus
    if entry > limit_detritus [
      set entry limit_detritus
      set flag_limit_detritus TRUE
    ]
    set w_detritus (entry * gC_ww_detritus * pondArea) ; column 49 in input file ("Sediment POC")
;    type "w_detritus from list, wet weight: " type w_detritus type ", limit_detritus: " print limit_detritus ; for testing
;; 27 Nov 2018: this line was included wrongly: the limit is already applied correctly above (comparison to gC rather than wet weight
    ;    if w_detritus > limit_detritus [ set w_detritus limit_detritus ] ; new June 5; detritus is only accessible for consumption by TS up to a limit
;    type "w_detritus after threshold application: " print w_detritus ; for testing
    set item_cnt item_cnt + 1
    set entry file-read
    set item_cnt item_cnt + 1
    set entry file-read
    set preyatthestartoftheday lput entry preyatthestartoftheday
    set w_periphyton (entry * gC_ww_periphyton * pondArea) ; column 51 in input file ("Total periphyton")
    ifelse flag_TS_eaten [
      set w_TS_eaten 0 ; if TS eaten is not listed in the prey biomass input file, it is assumed that no predation occurs
      set entry file-read-line
    ] [
      while [ item_cnt < 69 ] [ ; it is unclear to me what's happening here, but I get the right column with this (even though TSgww-eaten is in copumn 73)
        set entry file-read
        set item_cnt item_cnt + 1
      ]
      set w_TS_eaten (file-read * pondVolume) ; column 71 in input file (TSgww - eaten); note: biomass in g wet weight/m3 (not gC/m3)
    ]
  ] [ user-message "Daily prey biomass input file does not exist in current directory!" ]
  file-close
end

to write_output
  ; daily biomass output
  ; note that output file needs to be in the same order as input file and the list "preyattheendoftheday" (see 'go' procedure)

  if file-exists? BiomassOutput [
   if user-yes-or-no? "BiomassOutput file already exists. Delete?" [
      file-delete BiomassOutput
    ]
  ]
  file-open BiomassOutput
  file-write "Day"
  file-write "Copepods"
  file-write "Cladocerans"
  file-write "Rotifers"
  file-write "MicroZplankton"
  file-write "Ephemeroptera"
  file-write "Trichoptera"
  file-write "Oligochaetes"
  file-write "Chironomids"
  file-write "SedimentPOC"
  file-write "TotalPeriphyt"
  file-write "TopShiner(juv)"
  file-write "TopShiner(adt)"
  file-write "TopShiner(all)"
  file-write "numTopShiner(juv)"
  file-write "numTopShiner(adt)"
  file-write "numTopShiner(all)"
  file-write "TSnum-eaten"
  file-write "TSgww-eaten"
  ;; 3 December 2018: culomns below added to output file
  file-write "TopShiner(eggs_f)"
  file-write "TopShiner(eggs_m)"
  file-write "TopShiner(larvae_f)"
  file-write "TopShiner(larvae_m)"
  file-write "TopShiner(juv_f)"
  file-write "TopShiner(juv_m)"
  file-write "TopShiner(subadult_f)"
  file-write "TopShiner(subadult_m)"
  file-write "TopShiner(adt_1yr_f)"
  file-write "TopShiner(adt_1yr_m)"
  file-write "TopShiner(adt_2yr_f)"
  file-write "TopShiner(adt_2yr_m)"
  file-write "TopShiner(adt_3yr_f)"
  file-write "TopShiner(adt_3yr_m)"
  file-write "numTopShiner(eggs_f)"
  file-write "numTopShiner(eggs_m)"
  file-write "numTopShiner(larvae_f)"
  file-write "numTopShiner(larvae_m)"
  file-write "numTopShiner(juv_f)"
  file-write "numTopShiner(juv_m)"
  file-write "numTopShiner(subadult_f)"
  file-write "numTopShiner(subadult_m)"
  file-write "numTopShiner(adt_1yr_f)"
  file-write "numTopShiner(adt_1yr_m)"
  file-write "numTopShiner(adt_2yr_f)"
  file-write "numTopShiner(adt_2yr_m)"
  file-write "numTopShiner(adt_3yr_f)"
  file-write "numTopShiner(adt_3yr_m)"
  file-print " "

  let day 0
  while [ day < length Prey_W_output_list ] [
    file-write (day + 1)
    let i 0
    while [ i < length (item day Prey_W_output_list) ] [
      file-write item i (item day Prey_W_output_list)
      set i i + 1
    ]
    set i 0
    while [ i < length (item day TS_W_output_list) ] [
      file-write item i (item day TS_W_output_list)
      set i i + 1
    ]
    file-print ""
    set day day + 1
  ]
  file-close
end

to writeFishStatus
  ; procedure for testing/calibration only
  ; output of status of single fish
  if file-exists? FishStatusOut [
   if user-yes-or-no? "FishStatusOut file already exists. Delete?" [
      file-delete FishStatusOut
    ]
  ]
  file-open FishStatusOut
  file-write "age" file-write "stage" file-write "sex" file-write "W" file-write "SL" file-print ""
  ask turtles [
    file-write age
    file-write stage ; fish life stage: egg, larva, juvenile, subadult, adult
    file-write sex ; fish sex: female, male
    file-write W ; fish weight in g
    file-write SL ; fish standard length in mm
    file-write K ; fish condition
    file-print ""
  ]
  file-close
end

to writeFishGrowthTest
  ; procedure for testing/calibration only
  ; output of status of single fish
  if (turtle 0) = NOBODY [
    user-message "Focal fish died: run ended."
    stop
  ]
  file-open FishGrowthTest
  ask turtle 0 [
    file-write ticks
    file-write age
    file-write stage ; fish life stage: egg, larva, juvenile, subadult, adult
    file-write sex ; fish sex: female, male
    file-write W ; fish weight in g
    file-write SL ; fish standard length in mm
    file-write K ; fish condition
    file-write day_consumption
    file-write R_T
    file-write S
    file-write E_C
    file-write d_W
    file-write w_copepods ; weight of prey guild copepods
    file-write w_cladocerans ; weight of prey guild cladocerans
    file-write w_rotifers ; weight of prey guild rotifers
    file-write w_microzoo ; weight of prey guild microzooplankton
    file-write w_mayflies ; weight of prey guild ephemeroptera
    file-write w_caddisflies ; weight of prey guild trichoptera
    file-write w_oligochaetes ; weight of prey guild oligochaetes
    file-write w_chironomids ; weight of prey guild chironomids
    file-write w_detritus ; weight of sediment POC
    file-write w_periphyton ; weight of prey guilds summarized as total periphyton
    file-write temperature
    file-print ""
  ]
  file-close
end

to-report DateREP ; converts day ('tick') to date (day-month-year format) ; copied from BEEHAVE (Becher et al. 2014)
  let month-names (list "January" "February" "March" "April" "May" "June" "July" "August" "September" "October" "November" "December")
  let days-in-months (list 31 28 31 30 31 30 31 31 30 31 30 31)


  let year floor (ticks / 365.01) + 1
  let month 0
  let dayOfYear remainder ticks 365
  if dayOfYear = 0 [ set dayOfYear 365 ]
  let dayOfMonth 0
  let sumDaysInMonths 0
  while [ sumDaysInMonths < dayOfYear ]
  [
    set month month + 1
    set sumDaysInMonths sumDaysInMonths + item (month - 1) days-in-months
    set dayOfMonth dayOfYear - sumDaysInMonths + item (month - 1) days-in-months
  ]

  report (word dayOfMonth "  " (item (month - 1) month-names) " " year )

end
@#$#@#$#@
GRAPHICS-WINDOW
350
10
817
478
-1
-1
9.0
1
14
1
1
1
0
1
1
1
-25
25
-25
25
1
1
1
ticks
30.0

BUTTON
11
13
78
46
setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
11
51
78
84
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
11
212
333
245
pondArea
pondArea
1
2000
100.0
1
1
m2
HORIZONTAL

INPUTBOX
11
417
223
492
InputParameters
InputParameters_20Jun2018.txt
1
0
String

INPUTBOX
11
559
223
632
DailyPreyBiomass
IBM_TS_master_ref.out
1
0
String

BUTTON
84
51
147
84
Day
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

BUTTON
11
92
113
125
Write results
write_output
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

INPUTBOX
11
125
281
185
BiomassOutput
biomass_output_yearly_default_13Mar2019_2.txt
1
0
String

INPUTBOX
11
491
223
560
WaterConditionInput
env_casmTS_2010_08May2018.prn
1
0
String

BUTTON
24
779
201
812
NIL
readWaterConditions
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
26
675
176
693
Block for testing only:
11
0.0
1

BUTTON
24
812
201
845
NIL
readDailyPreyBiomass_year
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
712
478
817
523
Date
DateRep
17
1
11

INPUTBOX
156
696
236
756
StartDay
152.0
1
0
Number

BUTTON
24
696
158
729
Setup on StartDay
SetupDay
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
284
729
476
789
FishStatusOut
fish_status.txt
1
0
String

BUTTON
24
738
147
771
Write fish status
writeFishStatus
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
283
789
476
849
FishStatusIn
fish_egg_male.txt
1
0
String

INPUTBOX
282
888
496
957
FishGrowthTest
test_fish_growth_allCASM_6Aug2018_pref0_male_1.txt
1
0
String

SWITCH
283
856
497
889
WriteGrowthTestOutput
WriteGrowthTestOutput
1
1
-1000

PLOT
24
964
498
1183
Bioenergetics - fish 0
Day
Bioenergetics [J/d]
0.0
10.0
0.0
1.0E-4
true
true
"" ""
PENS
"consumption" 1.0 0 -10899396 true "" "plot [day_consumption / W] of turtle 0"
"respiration" 1.0 0 -2674135 true "" "plot (-1) * ([R_T + S] of turtle 0) "
"waste loss" 1.0 0 -955883 true "" "plot (-1) * ([E_C] of turtle 0) "
"0" 1.0 0 -7500403 false "" "plot 0"

PLOT
825
10
1338
232
Bioenergetics - all fish
Day
Energy [J]
0.0
10.0
0.0
1.0E-4
true
true
"" ""
PENS
"consumption" 1.0 0 -10899396 true "" "plot sum [day_consumption] of turtles"
"respiration" 1.0 0 -2674135 true "" "plot (-1) * sum [(R_T + S) * W] of turtles"
"waste loss" 1.0 0 -955883 true "" "plot (-1) * sum [E_C * W] of turtles"
"0" 1.0 0 -7500403 false "" "plot 0"

PLOT
825
231
1340
449
Fish numbers
Day
Fish number
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Larvae/10" 1.0 0 -10899396 true "" "plot 0.1 * (count turtles with [stage = \"larva\"])"
"Juveniles" 1.0 0 -11221820 true "" "plot count turtles with [stage = \"juvenile\"]"
"Subadults" 1.0 0 -13345367 true "" "plot count turtles with [stage = \"subadult\"]"
"Adults" 1.0 0 -5825686 true "" "plot count turtles with [stage = \"adult\"]"

SWITCH
11
279
134
312
Sunfish
Sunfish
0
1
-1000

SWITCH
11
312
134
345
Predation
Predation
1
1
-1000

SWITCH
24
889
182
922
CmaxConsumption
CmaxConsumption
1
1
-1000

SLIDER
24
921
182
954
ConstCmax
ConstCmax
0
100
11.0
1
1
NIL
HORIZONTAL

BUTTON
24
845
201
878
NIL
readDailyPreyBiomass_day
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
84
13
187
46
DailyIO
DailyIO
1
1
-1000

TEXTBOX
192
11
311
50
Switch DailyIO on for new prey and TS biomass inputs in every time step
10
0.0
1

TEXTBOX
13
395
163
413
Input files:
11
0.0
1

PLOT
825
701
1340
954
Prey biomass
Day
Prey guild biomass in gC/m2
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"Copepods" 1.0 0 -955883 true "" "plot item 0 preyatthestartoftheday"
"Cladocerans" 1.0 0 -5298144 true "" "plot item 1 preyatthestartoftheday"
"Rotifers" 1.0 0 -3508570 true "" "plot item 2 preyatthestartoftheday"
"Microzooplankton" 1.0 0 -8630108 true "" "plot item 3 preyatthestartoftheday"
"Mayflies" 1.0 0 -14454117 true "" "plot item 4 preyatthestartoftheday"
"Caddisflies" 1.0 0 -8990512 true "" "plot item 5 preyatthestartoftheday"
"Oligochaetes" 1.0 0 -7500403 true "" "plot item 6 preyatthestartoftheday"
"Chironomids" 1.0 0 -14070903 true "" "plot item 7 preyatthestartoftheday"
"Detritus" 1.0 0 -6459832 true "" "plot item 8 preyatthestartoftheday"
"Periphyton" 1.0 0 -13840069 true "" "plot item 9 preyatthestartoftheday"
"Chiro - after TS consumpt" 1.0 2 -14070903 true "" "plot item 7 preyattheendoftheday"

SWITCH
11
344
134
377
DD_egg_larva
DD_egg_larva
0
1
-1000

INPUTBOX
157
55
281
115
RandomNumberSeed
6.0
1
0
Number

TEXTBOX
13
654
517
682
___________________________________________________________________________________
11
0.0
1

TEXTBOX
508
663
523
1204
|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|
11
0.0
1

TEXTBOX
12
664
27
1196
|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|
11
0.0
1

TEXTBOX
13
1183
513
1207
___________________________________________________________________________________
11
0.0
1

PLOT
825
448
1340
701
Consumed biomass of prey
Day
Biomass consumed by TS in gC/m2
0.0
10.0
0.0
0.01
true
true
"" ""
PENS
"Copepods" 1.0 0 -817084 true "" "plot ((item 0 preyatthestartoftheday) - (item 0 preyattheendoftheday))"
"Cladocerans" 1.0 0 -5298144 true "" "plot (item 1 preyatthestartoftheday) - (item 1 preyattheendoftheday)"
"Rotifers" 1.0 0 -3508570 true "" "plot (item 2 preyatthestartoftheday) - (item 2 preyattheendoftheday)"
"Microzooplankton" 1.0 0 -8630108 true "" "plot (item 3 preyatthestartoftheday) - (item 3 preyattheendoftheday)"
"Mayflies" 1.0 0 -13791810 true "" "plot (item 4 preyatthestartoftheday) - (item 4 preyattheendoftheday)"
"Caddisflies" 1.0 0 -8990512 true "" "plot (item 5 preyatthestartoftheday) - (item 5 preyattheendoftheday)"
"Oligochaetes" 1.0 0 -7500403 true "" "plot (item 6 preyatthestartoftheday) - (item 6 preyattheendoftheday)"
"Chironomids" 1.0 0 -14070903 true "" "plot (item 7 preyatthestartoftheday) - (item 7 preyattheendoftheday)"
"Detritus" 1.0 0 -6459832 true "" "plot (item 8 preyatthestartoftheday) - (item 8 preyattheendoftheday)"
"Periphyton" 1.0 0 -13840069 true "" "plot (item 9 preyatthestartoftheday) - (item 9 preyattheendoftheday)"

SLIDER
185
311
333
344
ScalingSearchArea
ScalingSearchArea
0
10
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
185
344
333
377
MaxDetritusDepth
MaxDetritusDepth
0.1
1.5
1.5
0.1
1
cm
HORIZONTAL

SLIDER
185
278
333
311
AverageKmin
AverageKmin
0
1
0.68
0.01
1
NIL
HORIZONTAL

INPUTBOX
350
477
464
537
YearsToRun
1.0
1
0
Number

TEXTBOX
471
481
621
537
Stops the simulation run after the entered number of years.\nIf set to '-1', simulation goes on forever.
11
0.0
1

SLIDER
11
245
333
278
CASM_pool_area
CASM_pool_area
1
2000
1125.0
1
1
m2
HORIZONTAL

@#$#@#$#@
## Individual-based population model for the Topeka shiner (TS-IBM)

The model is part of the hybrid model approach for the Topeka shiner. 
The model description is available in Appenix B of the publication:

Schmolke A, Bartell SM, Roy C, Green N, Galic N, Brain R. 2019. Species-specific population dynamics and their link to an aquatic food web: a hybrid modeling approach. Ecological Modelling, in press

The hybrid model (including the current NetLogo model, CASM-TS executable, input files and scripts) is available from:
https://github.com/Waterborne-env/Topeka-shiner-model
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.1
@#$#@#$#@
set grass? true
setup
repeat 75 [ go ]
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
