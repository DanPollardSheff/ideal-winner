#install.packages("devtools")

library(devtools)
library(MASS)

#Global variables
PSA_switch <- 1                 #1=run PSA, 0=deterministic
PSA_numb <- 750                 #number of PSA runs
pat_numb <- 25000               #number of patients
days_to_discharge <- 30         #number of days to discharge from hospital
days_in_year <- 365.25          #number of days in a year
time_horizon <- 100             #time horizon, years
discount_rate_QALYs <- 0.035    #discount rate, QALYs
discount_rate_costs <- 0.035    #discount rate, costs
Param_export <- 1               #1=save a copy of PSA parameters

#The proportion of clinical benefit that patients with an ISS of over 8 and under 16 receive
#compared to people with an ISS of 16 or more. Default is no benefit
Proportion_RR_MTC_ISS_o8_u16_hosp <- 0
Proportion_RR_MTC_ISS_o8_u16_1yr <- 0

#The proportion of clinical benefit that patients who are initially sent to a non-MTC receive
#compared to people sent straight to an MTC. Default is full benefit
Proportion_RR_MTC_transfer_hosp <- 0
Proportion_RR_MTC_ISS_transfer_1yr <- 0


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

PSA_rand_no <-  1346                 #random number to determine PSA parameters either -99 (to not reset the seed) or any positive number

date <- "_2_75_ben_trans"          #name to append to saved files 

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
sens_75_spec_66_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.746, 0.657,1)
write.csv(sens_75_spec_66_PSA, paste(file_location,"PSA results\\sens_75_spec_66",date,".csv", sep=""))
  