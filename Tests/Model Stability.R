
#install.packages("MASS")
#install.packages("doParallel")

library(MASS)
library(parallel)
library(doParallel)

numCores <- (detectCores() -1)  #Number of cores available minus 1, to
                                #leave 1 Core for OS

#read in the r script that sets the global variables
source("Set global variables.R")
#set number of patients to a large number
pat_numb <- 250000
#set to deterministic
PSA_switch <- 0
PSA_numb <- 1

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
##example sens 99.8%, spec 2.5%, 1000 PSA runs
start_time <- Sys.time()
MATTSP3 <- run_simulation(pat_chars, parameters, PSA_numb, "MATTSP3", NA, NA,0)
end_time <- Sys.time()
end_time - start_time
model_runtime <- end_time - start_time
LAS <- run_simulation(pat_chars, parameters, PSA_numb, "LAS", NA, NA,0)
SWAS <- run_simulation(pat_chars, parameters, PSA_numb, "SWAS", NA, NA,0)
WMAS <- run_simulation(pat_chars, parameters, PSA_numb, "WMAS", NA, NA,0)
YAS <- run_simulation(pat_chars, parameters, PSA_numb, "YAS", NA, NA,0)

#####Store results in a matrix
stability_res <- matrix(data=NA, nrow = pat_numb, ncol = 11)
colnames(stability_res) <- c("ID", "MATTSP3Cost", "LASCost", "SWASCost", "WMASCost", "YASCost",
                             "MATTSP3QALY", "LASQALY", "SWASQALY", "WMASQALY", "YASQALY")
#Record ID
stability_res[,"ID"] <- 1:pat_numb
#Get cumulative costs for each strategy
stability_res[,"MATTSP3Cost"] <- ave(MATTSP3[,"DCosts"],FUN=cumsum)
stability_res[,"MATTSP3Cost"] <- stability_res[,"MATTSP3Cost"]/stability_res[,"ID"]

stability_res[,"LASCost"] <- ave(LAS[,"DCosts"],FUN=cumsum)
stability_res[,"LASCost"] <- stability_res[,"LASCost"]/stability_res[,"ID"]

stability_res[,"SWASCost"] <- ave(SWAS[,"DCosts"],FUN=cumsum)
stability_res[,"SWASCost"] <- stability_res[,"SWASCost"]/stability_res[,"ID"]

stability_res[,"WMASCost"] <- ave(WMAS[,"DCosts"],FUN=cumsum)
stability_res[,"WMASCost"] <- stability_res[,"WMASCost"]/stability_res[,"ID"]

stability_res[,"YASCost"] <- ave(YAS[,"DCosts"],FUN=cumsum)
stability_res[,"YASCost"] <- stability_res[,"YASCost"]/stability_res[,"ID"]

#Get cumulative QALYs for each strategy
stability_res[,"MATTSP3QALY"] <- ave(MATTSP3[,"dQALYS"],FUN=cumsum)
stability_res[,"MATTSP3QALY"] <- stability_res[,"MATTSP3QALY"]/stability_res[,"ID"]

stability_res[,"LASQALY"] <- ave(LAS[,"dQALYS"],FUN=cumsum)
stability_res[,"LASQALY"] <- stability_res[,"LASQALY"]/stability_res[,"ID"]

stability_res[,"SWASQALY"] <- ave(SWAS[,"dQALYS"],FUN=cumsum)
stability_res[,"SWASQALY"] <- stability_res[,"SWASQALY"]/stability_res[,"ID"]

stability_res[,"WMASQALY"] <- ave(WMAS[,"dQALYS"],FUN=cumsum)
stability_res[,"WMASQALY"] <- stability_res[,"WMASQALY"]/stability_res[,"ID"]

stability_res[,"YASQALY"] <- ave(YAS[,"dQALYS"],FUN=cumsum)
stability_res[,"YASQALY"] <- stability_res[,"YASQALY"]/stability_res[,"ID"]

#trun stability res into a dataframe for ggplot 2
stability_res <- as.data.frame(stability_res)

#plots
install.packages("ggplot2")
library(ggplot2)

CostGraph <- ggplot(stability_res, aes(x=ID))+
  geom_line(aes(y = MATTSP3Cost), color = "red")+
  geom_line(aes(y = LASCost), color = "blue", linetype = 2)+
  geom_line(aes(y = SWASCost), color = "yellow", linetype = 3)+
  geom_line(aes(y = WMASCost), color = "purple", linetype = 4)+
  geom_line(aes(y = YASCost), color = "orange", linetype = 5)+
  ylim(30000,35000)

CostGraph
ggsave("Results/StabilityCostGraph.png", plot = CostGraph)

QALYGraph <- ggplot(stability_res, aes(x=ID))+
  geom_line(aes(y = MATTSP3QALY), color = "red")+
  geom_line(aes(y = LASQALY), color = "blue", linetype = 2)+
  geom_line(aes(y = SWASQALY), color = "yellow", linetype = 3)+
  geom_line(aes(y = WMASQALY), color = "purple", linetype = 4)+
  geom_line(aes(y = YASQALY), color = "orange", linetype = 5)+
  ylim(11,14)

QALYGraph
ggsave("Results/StabilityQALYGraph.png", plot = QALYGraph)
