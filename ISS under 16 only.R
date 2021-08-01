#install.packages("devtools")

library(devtools)
library(MASS)
library(parallel)
library(doParallel)

numCores <- (detectCores() -1) #Leave 1 Core for OS

#Global variables
PSA_switch <- 1                 #1=run PSA, 0=deterministic
PSA_numb <- 2000                 #number of PSA runs
pat_numb <- 25000               #number of patients
days_to_discharge <- 30         #number of days to discharge from hospital
days_in_year <- 365.25          #number of days in a year
time_horizon <- 100             #time horizon, years
discount_rate_QALYs <- 0.035    #discount rate, QALYs
discount_rate_costs <- 0.035    #discount rate, costs
Param_export <- 1               #1=save a copy of PSA parameters

#The proportion of clinical benefit that patients with an ISS of over 8 and under 16 receive
#compared to people with an ISS of 16 or more (0 to 1). Default is no benefit (0)
Proportion_RR_MTC_ISS_o8_u16_hosp <- 0
Proportion_RR_MTC_ISS_o8_u16_1yr <- 0

#The proportion of clinical benefit that patients who are initially sent to a non-MTC receive
#compared to people sent straight to an MTC (0 to 1). Default is full benefit (1)
Proportion_RR_MTC_transfer_hosp <- 1
Proportion_RR_MTC_ISS_transfer_1yr <- 1


TARN_mort_eq <- "Old"           #options are new or old. Default is old
MTCs_in_mort_risk <- "No"       #options are Yes or no. Relates to whether the mort eq is a composite risk score for a 

#population who has / has not been to an MTC or a population who hasn't gone to an MTC. Default is no, as the
#default for the  mortality equation is the Old TARN equation.
percent_TARN_cases_reported_ISS_o16 <- 1
percent_TARN_cases_reported_ISS_o9_u16 <- 1

population_source <- "Dutch" #Source of simulated population.Options are UK, Dutch, Dutch_simp. Dutch is the default

population_ISS_over16_only <- "No"  #Option for resampling to produce a population with ISS >= 16. Options are "Yes" or "No". Default is no. 
population_ISS_under16_only <- "No" #Option for resampling to produce a population with ISS < 16. Can be "Yes" or "No". Default is no.

efficent_life_expectancy <- "Yes"   #Options are Yes or No. Default is yes

test_pat_chars <- "No"              #Change this to Yes if you only want to run the base case analysis with patient level results

PSA_rand_no <-  1477                 #random number to determine PSA parameters either -99 (to not reset the seed) or any positive number

date <- "_1_75_ben_ISS_u9"          #name to append to saved files 

#read in files / save files from the X drive (note not on Git due to confidentiality reasons)
file_location <- "\\\\uosfstore.shefuniad.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\Model\\"

param_data <- read.csv("parameters.csv", row.names=1)
life_tabs <- read.csv("ONSlifetables.csv")
future_costs <- read.csv("lifetime-healthcare-costs.csv")

if(population_source=="UK"){
  means <- as.matrix(read.csv(paste(file_location,"means.csv", sep=""),row.names=1))
  covariance <- as.matrix(read.csv(paste(file_location,"covariance.csv", sep=""), row.names=1))
  age_tab <- read.csv(paste(file_location,"age_tab.csv", sep=""),row.names=1)
  gen_tab <- read.csv(paste(file_location,"gen_tab.csv", sep=""),row.names=1)
  ISS_tab <- read.csv(paste(file_location,"ISS_tab.csv", sep=""),row.names=1)
  GCS_tab <- read.csv(paste(file_location,"GCS_tab.csv", sep=""),row.names=1)
}else if (population_source== "Dutch_simp"){
  means <- as.matrix(read.csv("means_dutch_v2.csv",row.names=1))
  covariance <- as.matrix(read.csv("covariance_dutch_v2.csv", row.names=1))
  age_tab <- read.csv("age_tab_dutch_v2.csv",row.names=1)
  gen_tab <- read.csv("male_tab_dutch_v2.csv",row.names=1)
  ISS_tab <- read.csv("ISS_tab_dutch_v2.csv",row.names=1)
  GCS_tab <- read.csv("GCS_tab_dutch_v2.csv",row.names=1)
  blunt_tab <- read.csv("blunt_tab_dutch_v2.csv",row.names=1)
}else{
  means <- as.matrix(read.csv(paste(file_location,"means_dutch_v2.csv", sep=""), row.names=1))
  covariance <- as.matrix(read.csv(paste(file_location,"covariance_dutch_v2.csv", sep=""), row.names=1))
  age_tab <- read.csv(paste(file_location,"age_tab_dutch_v2.csv", sep=""),row.names=1)
  gen_tab <- read.csv(paste(file_location,"male_tab_dutch_v2.csv", sep=""),row.names=1)
  ISS_tab <- read.csv(paste(file_location,"ISS_tab_dutch_v2.csv", sep=""),row.names=1)
  GCS_tab <- read.csv(paste(file_location,"GCS_tab_dutch_v2.csv", sep=""),row.names=1)
  blunt_tab <- read.csv(paste(file_location,"blunt_tab_dutch_v2.csv", sep=""),row.names=1)
}

#Call in all functions
source("Functions.R")

#Analysis###################
param_data_bc <- param_data

##########################################################

#### add in analysis run here
##example sens 99.8%, spec 2.5%, 1000 PSA runs
All <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 1, 0,1)
write.csv(All, "\\\\uosfstore.shefuniad.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\July 21\\All_ISSU16.csv")

WMAS_Step4 <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 1, 0.91,0.25)
write.csv(WMAS_Step4, "\\\\uosfstore.shefuniad.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\July 21\\WMAS_Step4_ISSU16.csv")

MATTS_sens <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 1, 0.70,0.74)
write.csv(MATTS_sens, "\\\\uosfstore.shefuniad.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\July 21\\MATTS_sens_ISSU16.csv")

MATTS_bal <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 1, 0.57, 0.85)
write.csv(MATTS_bal, "\\\\uosfstore.shefuniad.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\July 21\\MATTS_bal_ISSU16.csv")

field_triage <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 1, 0.48, 0.86)
write.csv(field_triage, "\\\\uosfstore.shefuniad.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\July 21\\field_triage_ISSU16.csv")

oregon <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 1, 0.44, 0.90)
write.csv(oregon, "\\\\uosfstore.shefuniad.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\July 21\\oregon_ISSU16.csv")

MATTS_spec <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 1, 0.38, 0.94)
write.csv(MATTS_spec, "\\\\uosfstore.shefuniad.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\July 21\\MATTS_spec_ISSU16.csv")

SWAST <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 1, 0.34, 0.94)
write.csv(SWAST, "\\\\uosfstore.shefuniad.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\July 21\\SWAST_ISSU16.csv")

TTR <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 1, 0.27, 0.95)
write.csv(TTR, "\\\\uosfstore.shefuniad.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\July 21\\TTR_ISSU16.csv")

pre_hosp_index <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 1, 0.20, 0.96)
write.csv(pre_hosp_index, "\\\\uosfstore.shefuniad.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\July 21\\pre_hosp_index_ISSU16.csv")

None <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 1, 0, 1)
write.csv(None, "\\\\uosfstore.shefuniad.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\July 21\\None_ISSU16.csv")

