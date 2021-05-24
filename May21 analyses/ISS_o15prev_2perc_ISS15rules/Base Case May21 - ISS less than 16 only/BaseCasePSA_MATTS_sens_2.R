#install.packages("devtools")

library(devtools)
library(MASS)

#Global variables

PSA_switch <- 1
PSA_numb <- 750

pat_numb <- 25000
days_to_discharge <- 30
days_in_year <- 365.25
time_horizon <- 100
discount_rate_QALYs <- 0.035
discount_rate_costs <- 0.035
Param_export <- 1

#benefit of MTC care for patients with an ISS between over 8 and under 16
Proportion_RR_MTC_ISS_o8_u16_hosp <- 0
Proportion_RR_MTC_ISS_o8_u16_1yr <- 0

#The proportion of clinical benefit that patients who are initially sent to a non-MTC receive
#compared to people sent straight to an MTC
#Default is full benefit
Proportion_RR_MTC_transfer_hosp <- 1
Proportion_RR_MTC_ISS_transfer_1yr <- 1

TARN_mort_eq <- "Old" # options are new or old. Default is old
MTCs_in_mort_risk <- "No" #options are Yes or no. Relates to whether the mort eq is a composite risk score for a 
#population who has / has not been to an MTC or a population who hasn't gone to an MTC. Default is no, as the
#default for the  mortality equation is the Old TARN equation.
percent_TARN_cases_reported_ISS_o16 <- 1
percent_TARN_cases_reported_ISS_o9_u16 <- 1
population_source <- "Dutch" # Options are UK and Dutch. Dutch is the default
population_ISS_over16_only <- "Yes" # Options are yes or no. Default is no. 
population_ISS_under16_only <- "No" # Options are yes or no. Default is no. only this or the above can be "Yes". Function will not work otherwise
efficent_life_expectancy <- "Yes" #Options are Yes or No. Default is yes

test_pat_chars <- "No" #Change this to Yes if you only want to run the base case analysis with patient level results

PSA_rand_no <-  1346 #random number to determine PSA parameters. #if -99 this will not change the seed after randomly determining the number of patients to run through the model. 
#settings for MATTS phase 1 where first 500 runs 26090100 (after generating pat chars), next 1000 runs (ten diagnostic strategies only) 1346

date <- "_2" #name to append to saved files 

#read in files from the X drive (note not on Git due to confidentiality reasons)
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
}else{
  means <- as.matrix(read.csv(paste(file_location,"means_dutch_v2.csv", sep=""),row.names=1))
  covariance <- as.matrix(read.csv(paste(file_location,"covariance_dutch_v2.csv", sep=""), row.names=1))
  age_tab <- read.csv(paste(file_location,"age_tab_dutch_v2.csv", sep=""),row.names=1)
  gen_tab <- read.csv(paste(file_location,"male_tab_dutch_v2.csv", sep=""),row.names=1)
  ISS_tab <- read.csv(paste(file_location,"ISS_tab_dutch_v2.csv", sep=""),row.names=1)
  GCS_tab <- read.csv(paste(file_location,"GCS_tab_dutch_v2.csv", sep=""),row.names=1)
  blunt_tab <- read.csv(paste(file_location,"blunt_tab_dutch_v2.csv", sep=""),row.names=1)
}




#Call in all functions
source("Functions.R")

#Do you want to the use pre-simluated population and PSA?
predefined_pop_PSA <- "No"  # Option to use the pre-simulated population and PSA parameters
#Set to "Yes" if using the publicly shared version of the model
# In the predefined population we have merged some ISS and age categories for potential
#identifiability reasons



#Analysis###################
param_data_bc <- param_data

##########################################################

MATTS_sens_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual",0.702, 0.738,1)
write.csv(MATTS_sens_PSA, paste(file_location,"PSA results\\May21\\MATTSsens",date,".csv", sep=""))

