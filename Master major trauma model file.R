#install.packages("devtools")

library(devtools)
library(MASS)

#Global variables
PSA_switch <- 1
PSA_numb <- 500
pat_numb <- 25000
days_to_discharge <- 30
days_in_year <- 365.25
time_horizon <- 100
discount_rate_QALYs <- 0.035
discount_rate_costs <- 0.035
Param_export <- 1

#The proportion of clinical benefit that patients with an ISS of over 8 and under 16 receive
#compared to people with an ISS of 16 or more
#Default is no benefit
Proportion_RR_MTC_ISS_o8_u16_hosp <- 0
Proportion_RR_MTC_ISS_o8_u16_1yr <- 0

#The proportion of clinical benefit that patients who are initially sent to a non-MTC receive
#compared to people sent staright to an MTC
#Default is full benefit
Proportion_RR_MTC_transfer_hosp <- 1
Proportion_RR_MTC_ISS_transfer_1yr <- 1


TARN_mort_eq <- "Old" # options are new or old. Default is old
MTCs_in_mort_risk <- "No" #options are Yes or no. Relates to whether the mort eq is a composite risk score for a 
#population who has / has not been to an MTC or a population who hasn't gone to an MTC. Default is no, as the
#default for the  mortality equation is the Old TARN equation.
percent_TARN_cases_reported_ISS_o16 <- 1
percent_TARN_cases_reported_ISS_o9_u16 <- 1
population_source <- "Dutch_simp" # Options are UK, Dutch, Dutch_simp. Dutch is the default
population_ISS_over16_only <- "No" # Options are yes or no. Default is no. 
efficent_life_expectancy <- "Yes" #Options are Yes or No. Default is yes

test_pat_chars <- "No" #Change this to Yes if you only want to run the base case analysis with patient level results

PSA_strat <- "S28_S1" #Option to make sure that each instance only runs one set of PSAs, as it is computationally intensive
#Options are: S100, S95, S90, S88, S75, S70, S64, S57, S28, MTC, nMTC, S100_S1, S95_S1, S90_S1, S88_S1, S75_S1, S70_S1, S64_S1, S57_S1, S28_S1

PSA_rand_no <-  -99 #random number to determine PSA parameters. #if -99 this will not change the seed after randomly determining the number of patients to run through the model. 
#settings for MATTS phase 1 where first 500 runs 26090100 (after generating pat chars), next 1000 runs (ten diagnostic strategies only) 1346

date <- "_1_75_ben_ISS_u9" #name to append to saved files 

#read in files from the X drive (note not on Git due to confidentiality reasons)
file_location <- "\\\\uosfstore.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\Model\\"

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
sens_100_spec_3_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.998, 0.025,1)
  