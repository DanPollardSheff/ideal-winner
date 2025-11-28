  
#install.packages("MASS")
#install.packages("doParallel")

library(MASS)
library(parallel)
library(doParallel)

numCores <- (detectCores() -1)  #Number of cores available minus 1, to
                                #leave 1 Core for OS

#read in the r script that sets the global variables
source("Set global variables.R")


#read in files / save files from the X drive (note not on Git due to confidentiality reasons)
file_location <- "X:/ScHARR/PR_MATTS/General/Health Economics/Phase 1  & 2/Model/"

#read in global data files
param_data <- read.csv("Parameters/parameters.csv", row.names=1)
triage_rules_params <- read.csv("Parameters/MATTSPhase3Rules.csv")
tarn_22_means <- read.csv("Parameters/New TARN Means.csv", row.names = 1)
tarn_22_vcov <- read.csv("Parameters/New TARN vcov matrix.csv", row.names=1)
life_tables <- read.csv("Parameters/ONSlifetables.csv")
future_costs <- read.csv("Parameters/lifetime-healthcare-costs.csv")

if(population_source=="UK"){
  means <- as.matrix(read.csv(paste(file_location,"means.csv", sep=""),row.names=1))
  covariance <- as.matrix(read.csv(paste(file_location,"covariance.csv", sep=""), row.names=1))
  age_tab <- read.csv(paste(file_location,"age_tab.csv", sep=""),row.names=1)
  gen_tab <- read.csv(paste(file_location,"gen_tab.csv", sep=""),row.names=1)
  ISS_tab <- read.csv(paste(file_location,"ISS_tab.csv", sep=""),row.names=1)
  GCS_tab <- read.csv(paste(file_location,"GCS_tab.csv", sep=""),row.names=1)
}else if (population_source== "Dutch_simp"){
  means <- as.matrix(read.csv("Population/means_dutch_v2.csv",row.names=1))
  covariance <- as.matrix(read.csv("Population/covariance_dutch_v2.csv", row.names=1))
  age_tab <- read.csv("Population/age_tab_dutch_v2.csv",row.names=1)
  gen_tab <- read.csv("Population/male_tab_dutch_v2.csv",row.names=1)
  ISS_tab <- read.csv("Population/ISS_tab_dutch_v2.csv",row.names=1)
  GCS_tab <- read.csv("Population/GCS_tab_dutch_v2.csv",row.names=1)
  blunt_tab <- read.csv("Population/blunt_tab_dutch_v2.csv",row.names=1)
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

####Generate patient characteristics
#set the random number seed
set.seed(26090100)
pat_chars <- gen_pat_chars(pat_numb, means, covariance, age_tab, gen_tab, ISS_tab, GCS_tab)

###Generate parameters
if(PSA_rand_no != -99){
  set.seed(PSA_rand_no)
}
#generate the parameters
parameters <- gen_parameters(PSA_switch,PSA_numb, param_data_bc)

#export the parameters, if required (bug checking / SAVI)
if(Param_export==1){
  write.csv(parameters, file = "parameter_outputs.csv")
}

#Generate a random number table for life table analysis
random_numbs_LE <- array(data = runif(pat_numb*101*2), dim = c(pat_numb,101,2))
#Reset the first columns to be patient ID's
random_numbs_LE[,1,1] <- pat_chars[,"ID"]
random_numbs_LE[,1,2] <- pat_chars[,"ID"]

#Proportion Elderly
ISS_o15 <- sum(ifelse(pat_chars[,"ISS"]>15,1,0))
ISS_u15 <- length(pat_chars[,"ISS"])-ISS_o15
Elderly_ISSo15 <- sum(ifelse(pat_chars[,"Age"]>=65&pat_chars[,"ISS"]>15,1,0))
Elderly_ISSu15 <- sum(ifelse(pat_chars[,"Age"]>=65&pat_chars[,"ISS"]<15,1,0))
Elderly <- sum(ifelse(pat_chars[,"Age"]>=65,1,0))

Elderly_ISSo15/ISS_o15
Elderly_ISSu15/ISS_u15
##########################################################

#### add in analysis run here
##example sens 99.8%, spec 2.5%, 1000 PSA runs
start_time <- Sys.time()
P2_WMAS <- run_simulation(pat_chars, parameters, PSA_numb, "Phase2_WMAS", NA, NA,1,random_numbs_LE)
end_time <- Sys.time()
end_time - start_time
model_runtime <- end_time - start_time
write.csv(P2_WMAS, "Results/Phase2_WMAS_PSA.csv")
P3_WMAS <- run_simulation(pat_chars, parameters, PSA_numb, "Phase3_WMAS", NA, NA,1,random_numbs_LE)
write.csv(P3_WMAS, "Results/Phase3_WMAS_PSA.csv")

SWAST <- run_simulation(pat_chars, parameters, PSA_numb, "Phase2_SWAST", NA, NA,1,random_numbs_LE)
write.csv(SWAST, "Results/SWAST_PSA.csv")

LAS <- run_simulation(pat_chars, parameters, PSA_numb, "Phase2_LAS", NA, NA,1,random_numbs_LE)
write.csv(LAS, "Results/LAS_PSA.csv")

P2_YAS <- run_simulation(pat_chars, parameters, PSA_numb, "Phase2_YAS", NA, NA,1,random_numbs_LE)
write.csv(P2_YAS, "Results/Phase2_YAS_PSA.csv")

P3_YAS <- run_simulation(pat_chars, parameters, PSA_numb, "Phase3_YAS", NA, NA,1,random_numbs_LE)
write.csv(P3_YAS, "Results/Phase3_YAS_PSA.csv")
