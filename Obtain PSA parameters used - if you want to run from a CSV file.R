#Obtain PSA parameters used


library(devtools)
library(MASS)

#call in the functions
source("Functions.R")

#Call in the data
file_location <- "\\\\uosfstore.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\Model\\"


means <- as.matrix(read.csv(paste(file_location,"means_dutch_v2.csv", sep=""),row.names=1))
covariance <- as.matrix(read.csv(paste(file_location,"covariance_dutch_v2.csv", sep=""), row.names=1))
age_tab <- read.csv(paste(file_location,"age_tab_dutch_v2.csv", sep=""),row.names=1)
gen_tab <- read.csv(paste(file_location,"male_tab_dutch_v2.csv", sep=""),row.names=1)
ISS_tab <- read.csv(paste(file_location,"ISS_tab_dutch_v2.csv", sep=""),row.names=1)
GCS_tab <- read.csv(paste(file_location,"GCS_tab_dutch_v2.csv", sep=""),row.names=1)
blunt_tab <- read.csv(paste(file_location,"blunt_tab_dutch_v2.csv", sep=""),row.names=1)
param_inputs <- read.csv(paste(file_location,"parameters.csv", sep=""), row.names=1)

#record the number of patients to generate
pat_numb <- 25000
population_ISS_over16_only <- 0 
population_source <- "Dutch"
#First set of parameters
PSA_switch <- 1
PSA_numb <- 500
PSA_rand_no <-  -99

#Produce the first set of PSA runs, note I need to gen pat chars first as for the 
#first set of runs the RN seed was not reset after generating pat chars

#set the random number seed
set.seed(26090100)
#Generate pat chars to be 
pat_chars <- gen_pat_chars(pat_numb, means, covariance, age_tab, gen_tab, ISS_tab, GCS_tab)

#Generate PSA values
PSA_1 <- gen_parameters(PSA_switch,PSA_numb, param_inputs)

#Second set of PSA runs (note I do not need to generate pat chars as for the 2nd and 3rd set of runs
#an independent RN seed was set to generate the PSA runs)
PSA_numb <- 750
PSA_rand_no <-  1346 
set.seed(PSA_rand_no)
PSA_2 <- gen_parameters(PSA_switch,PSA_numb, param_inputs)

#Second set of PSA runs (note I do not need to generate pat chars as for the 2nd and 3rd set of runs
#an independent RN seed was set to generate the PSA runs)
PSA_rand_no <-  330413
set.seed(PSA_rand_no)
PSA_3 <- gen_parameters(PSA_switch,PSA_numb, param_inputs)

PSA <- rbind(PSA_1,PSA_2,PSA_3)
write.csv(PSA, "PSA_params.csv")