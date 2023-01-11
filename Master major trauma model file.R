
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
file_location <- "\\\\uosfstore.shefuniad.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\Phase 1  & 2\\Model\\"

#read in global data files
param_data <- read.csv("Parameters/parameters.csv", row.names=1)
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
##example sens 99.8%, spec 2.5%, 1000 PSA runs
All5 <- run_simulation(pat_chars, parameters, PSA_numb, "manual", 0.8, 0.8,1)
col
Util_source             <- "Kruithof"
All <- run_simulation(pat_chars, parameters, PSA_numb, "manual", 0.8, 0.8,1)