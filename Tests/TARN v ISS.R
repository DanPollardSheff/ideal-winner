
#install.packages("MASS")
#install.packages("doParallel")

library(MASS)
library(parallel)
library(doParallel)

numCores <- (detectCores() -1)  #Number of cores available minus 1, to
#leave 1 Core for OS

#read in the r script that sets the global variables
source("Set global variables.R")

#Make sure that this is deterministic
PSA_switch <- 0

#read in files / save files from the X drive (note not on Git due to confidentiality reasons)
file_location <- "\\\\uosfstore.shefuniad.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\Phase 1  & 2\\Model\\"

#read in global data files
param_data <- read.csv("Parameters/parameters.csv", row.names=1)
triage_rules_params <- read.csv("Parameters/MATTSPhase3Rules.csv")
tarn_22_means <- read.csv("Parameters/New TARN Means.csv", row.names = 1)
tarn_22_vcov <- read.csv("Parameters/New TARN vcov matrix.csv", row.names=1)
life_tabs <- read.csv("Parameters/ONSlifetables.csv")
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

##########################################################

#### add in analysis run here
#oldTARN
oldTARN <- pat_chars
oldTARN[,"p_death_hosp"] <- TARN_old_mort_pred(oldTARN, parameters,1)
plot(oldTARN[,"ISS"], oldTARN[,"p_death_hosp"])
Mean_ISS_oldTARN <- aggregate(p_death_hosp~ISS, mean, data=oldTARN)
plot(Mean_ISS_oldTARN[,"ISS"], Mean_ISS_oldTARN[,"p_death_hosp"])

#New TARN
NewTARN <- pat_chars
NewTARN[,"p_death_hosp"] <- TARN_mort_pred(NewTARN, parameters,1)
plot(NewTARN[,"ISS"], NewTARN[,"p_death_hosp"])
Mean_ISS_NewTARN <- aggregate(p_death_hosp~ISS, mean, data=NewTARN)
plot(Mean_ISS_NewTARN[,"ISS"], Mean_ISS_NewTARN[,"p_death_hosp"])

#TARN 2015
#Set the functions to use old TARN and regenereate the parameters
TARN_22_params <- F
parameters2 <- gen_parameters(PSA_switch,PSA_numb, param_data_bc)

TARN2015 <- pat_chars
TARN2015[,"p_death_hosp"] <- TARN_mort_pred(TARN2015, parameters,1)
plot(TARN2015[,"ISS"], TARN2015[,"p_death_hosp"])
Mean_ISS_TARN2015 <- aggregate(p_death_hosp~ISS, mean, data=TARN2015)
plot(Mean_ISS_TARN2015[,"ISS"], Mean_ISS_TARN2015[,"p_death_hosp"])
