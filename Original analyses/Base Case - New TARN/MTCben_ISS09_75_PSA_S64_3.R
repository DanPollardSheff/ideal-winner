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
Proportion_RR_MTC_ISS_o8_u16_hosp <- 0
Proportion_RR_MTC_ISS_o8_u16_1yr <- 0
TARN_mort_eq <- "Old" # options are new or old. Default is old
MTCs_in_mort_risk <- "No" #options are Yes or no. Relates to whether the mort eq is a composite risk score for a 
#population who has / has not been to an MTC or a population who hasn't gone to an MTC. Default is no, as the
#default for the  mortality equation is the Old TARN equation.
percent_TARN_cases_reported_ISS_o16 <- 1
percent_TARN_cases_reported_ISS_o9_u16 <- 1
population_source <- "Dutch" # Options are UK and Dutch. Dutch is the default
population_ISS_over16_only <- "No" # Options are yes or no. Default is no. 
efficent_life_expectancy <- "Yes" #Options are Yes or No. Default is yes

test_pat_chars <- "No" #Change this to Yes if you only want to run the base case analysis with patient level results

PSA_strat <- "S64_S1" #Option to make sure that each instance only runs one set of PSAs, as it is computationally intensive
#Options are: S100, S95, S90, S88, S75, S70, S64, S57, S28, MTC, nMTC, S100_S1, S95_S1, S90_S1, S88_S1, S75_S1, S70_S1, S64_S1, S57_S1, S28_S1

PSA_rand_no <-  330413 #random number to determine PSA parameters. #if -99 this will not change the seed after randomly determining the number of patients to run through the model. 
#settings for MATTS phase 1 where first 500 runs 26090100 (after generating pat chars), next 1000 runs (ten diagnostic strategies only) 1346

date <- "_3_75_ben_ISS_u9" #name to append to saved files 

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
  
  
#with 20,000 patients the results are stable in the base case
if(PSA_switch==0){
sens_100_spec_3 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.998, 0.025,1)
sens_95_spec_19 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.948, 0.187,1)
sens_90_spec_58 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.904, 0.584,1)
sens_88_spec_63 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.875, 0.628,1)
sens_75_spec_66 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.746, 0.657,1)
sens_70_spec_70 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.698, 0.701,1)
sens_64_spec_76 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.642, 0.761,1)
sens_57_spec_80 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.57, 0.8,1)
sens_28_spec_89 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.284, 0.886,1)

#create a matrix to store all runs
det_analyses <- matrix (nrow = 9, ncol =12)
#name the columns to make analysis easier
colnames(det_analyses) <- c("Sens_DR","Spec_DR", "Number_recieving_MTC_care","proportion_died_before_discharge","proportion_died_between_discharge_and_1_year", "Years_lived",
                             "undiscounted_QALYs", "discounted_QALYs", "undiscounted_Costs", "discounted_Costs", "proportion_ISS_over_16", "proportion_ISS_over_8_under_16")
#name the rows with the appropiate strategy
rownames(det_analyses) <- c("sens_100_spec_3", "sens_95_spec_19", "sens_90_spec_58", "sens_88_spec_63", "sens_75_spec_66",
"sens_70_spec_70", "sens_64_spec_76", "sens_57_spec_80", "sens_28_spec_89")

det_analyses["sens_100_spec_3", ]<- sens_100_spec_3
det_analyses["sens_95_spec_19", ]<- sens_95_spec_19
det_analyses["sens_90_spec_58", ]<- sens_90_spec_58
det_analyses["sens_88_spec_63", ]<- sens_88_spec_63
det_analyses["sens_75_spec_66", ]<- sens_75_spec_66
det_analyses["sens_70_spec_70", ]<- sens_70_spec_70
det_analyses["sens_64_spec_76", ]<- sens_64_spec_76
det_analyses["sens_57_spec_80", ]<- sens_57_spec_80
det_analyses["sens_28_spec_89", ]<- sens_28_spec_89

write.csv(det_analyses,"base case.csv")
}
if(PSA_switch==1){
  if(PSA_strat == "S100"){
  sens_100_spec_3_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.998, 0.025,1)
  write.csv(sens_100_spec_3_PSA, paste(file_location,"PSA results\\sens_100_spec_3_PSA",date,".csv", sep=""))
  use_params_sens_100_spec_3_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_100_spec_3_PSA, "PSA results\\sens_100_spec_3_PSA_params.csv")
  }
  if(PSA_strat == "S95"){
  sens_95_spec_19_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.948, 0.187,1)
  write.csv(sens_95_spec_19_PSA, paste(file_location,"PSA results\\sens_95_spec_19_PSA",date,".csv", sep=""))
  use_params_sens_95_spec_19_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_95_spec_19_PSA, "PSA results\\sens_95_spec_19_PSA_params.csv")
  }
  if(PSA_strat == "S90"){
  sens_90_spec_58_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.904, 0.584,1)
  write.csv(sens_90_spec_58_PSA, paste(file_location,"PSA results\\sens_90_spec_58_PSA",date,".csv", sep=""))
  use_params_sens_90_spec_58_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_90_spec_58_PSA, "PSA results\\sens_90_spec_58_PSA_params.csv")
  }
  if(PSA_strat == "S88"){
  sens_88_spec_63_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.875, 0.628,1)
  write.csv(sens_88_spec_63_PSA, paste(file_location,"PSA results\\sens_88_spec_63_PSA",date,".csv", sep=""))
  use_params_sens_88_spec_63_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_88_spec_63_PSA, "PSA results\\sens_88_spec_63_PSA_params.csv")
  }
  if(PSA_strat == "S75"){
  sens_75_spec_66_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.746, 0.657,1)
  write.csv(sens_75_spec_66_PSA, paste(file_location,"PSA results\\sens_75_spec_66",date,".csv", sep=""))
  use_params_sens_75_spec_66_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_75_spec_66_PSA, "sens_75_spec_66_PSA_params.csv")
  }
  if(PSA_strat == "S70"){
  sens_70_spec_70_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.698, 0.701,1)
  write.csv(sens_70_spec_70_PSA, paste(file_location,"PSA results\\sens_70_spec_70",date,".csv", sep=""))
  use_params_sens_70_spec_70_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_70_spec_70_PSA, "sens_70_spec_70_PSA_params.csv")
  }
  if(PSA_strat == "S64"){
  sens_64_spec_76_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.642, 0.761,1)
  write.csv(sens_64_spec_76_PSA, paste(file_location,"PSA results\\sens_64_spec_76",date,".csv", sep=""))
  use_params_sens_64_spec_76 <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_64_spec_76, "PSA results\\sens_90_spec_58_PSA_params.csv")
  }
  if(PSA_strat == "S57"){
  sens_57_spec_80_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.57, 0.8,1)
  write.csv(sens_57_spec_80_PSA, paste(file_location,"PSA results\\sens_57_spec_80",date,".csv", sep=""))
  use_params_sens_57_spec_80_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_57_spec_80_PSA, "PSA results\\sens_90_spec_58_PSA_params.csv")
  }
  if(PSA_strat == "S28"){
  sens_28_spec_89_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.284, 0.886,1)
  write.csv(sens_28_spec_89_PSA, paste(file_location,"PSA results\\sens_28_spec_89",date,".csv", sep=""))
  use_params_sens_28_spec_89_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_28_spec_89_PSA, "PSA results\\sens_28_spec_89_PSA_params.csv")
  }
#Use the newer TARN mortality equation
TARN_mort_eq <- "New" 
MTCs_in_mort_risk <- "Yes"
if(PSA_switch==1){
  if(PSA_strat == "S100_S1"){
    sens_100_spec_3_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.998, 0.025,1)
    write.csv(sens_100_spec_3_PSA, paste(file_location,"PSA results\\sens_100_spec_3_PSA",date,".csv", sep=""))
    use_params_sens_100_spec_3_PSA <- read.csv("parameter_outputs.csv")
    write.csv(use_params_sens_100_spec_3_PSA, "PSA results\\sens_100_spec_3_PSA_params.csv")
  }
  if(PSA_strat == "S95_S1"){
    sens_95_spec_19_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.948, 0.187,1)
    write.csv(sens_95_spec_19_PSA, paste(file_location,"PSA results\\sens_95_spec_19_PSA",date,".csv", sep=""))
    use_params_sens_95_spec_19_PSA <- read.csv("parameter_outputs.csv")
    write.csv(use_params_sens_95_spec_19_PSA, "PSA results\\sens_95_spec_19_PSA_params.csv")
  }
  if(PSA_strat == "S90_S1"){
    sens_90_spec_58_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.904, 0.584,1)
    write.csv(sens_90_spec_58_PSA, paste(file_location,"PSA results\\sens_90_spec_58_PSA",date,".csv", sep=""))
    use_params_sens_90_spec_58_PSA <- read.csv("parameter_outputs.csv")
    write.csv(use_params_sens_90_spec_58_PSA, "PSA results\\sens_90_spec_58_PSA_params.csv")
  }
  if(PSA_strat == "S88_S1"){
    sens_88_spec_63_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.875, 0.628,1)
    write.csv(sens_88_spec_63_PSA, paste(file_location,"PSA results\\sens_88_spec_63_PSA",date,".csv", sep=""))
    use_params_sens_88_spec_63_PSA <- read.csv("parameter_outputs.csv")
    write.csv(use_params_sens_88_spec_63_PSA, "PSA results\\sens_88_spec_63_PSA_params.csv")
  }
  if(PSA_strat == "S75_S1"){
    sens_75_spec_66_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.746, 0.657,1)
    write.csv(sens_75_spec_66_PSA, paste(file_location,"PSA results\\sens_75_spec_66",date,".csv", sep=""))
    use_params_sens_75_spec_66_PSA <- read.csv("parameter_outputs.csv")
    write.csv(use_params_sens_75_spec_66_PSA, "sens_75_spec_66_PSA_params.csv")
  }
  if(PSA_strat == "S70_S1"){
    sens_70_spec_70_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.698, 0.701,1)
    write.csv(sens_70_spec_70_PSA, paste(file_location,"PSA results\\sens_70_spec_70",date,".csv", sep=""))
    use_params_sens_70_spec_70_PSA <- read.csv("parameter_outputs.csv")
    write.csv(use_params_sens_70_spec_70_PSA, "sens_70_spec_70_PSA_params.csv")
  }
  if(PSA_strat == "S64_S1"){
    sens_64_spec_76_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.642, 0.761,1)
    write.csv(sens_64_spec_76_PSA, paste(file_location,"PSA results\\sens_64_spec_76",date,".csv", sep=""))
    use_params_sens_64_spec_76 <- read.csv("parameter_outputs.csv")
    write.csv(use_params_sens_64_spec_76, "PSA results\\sens_90_spec_58_PSA_params.csv")
  }
  if(PSA_strat == "S57_S1"){
    sens_57_spec_80_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.57, 0.8,1)
    write.csv(sens_57_spec_80_PSA, paste(file_location,"PSA results\\sens_57_spec_80",date,".csv", sep=""))
    use_params_sens_57_spec_80_PSA <- read.csv("parameter_outputs.csv")
    write.csv(use_params_sens_57_spec_80_PSA, "PSA results\\sens_90_spec_58_PSA_params.csv")
  }
  if(PSA_strat == "S28_S1"){
    sens_28_spec_89_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.284, 0.886,1)
    write.csv(sens_28_spec_89_PSA, paste(file_location,"PSA results\\sens_28_spec_89",date,".csv", sep=""))
    use_params_sens_28_spec_89_PSA <- read.csv("parameter_outputs.csv")
    write.csv(use_params_sens_28_spec_89_PSA, "PSA results\\sens_28_spec_89_PSA_params.csv")
  }
}


param_data_MTCs <- param_data
#Change the variables so everyone with a positive rule goes to the MTC
#everyone with a negative rule goes to an nMTC
#set the costs of MTCs to 0
param_data_MTCs["P_MTC_Tri_pos_ISS_o15",1] <- 1
param_data_MTCs["P_MTC_Tri_pos_ISS_o15",3] <- "Fixed"

param_data_MTCs["P_MTC_Tri_neg_ISS_o15",1] <- 0
param_data_MTCs["P_MTC_Tri_neg_ISS_o15",3] <- "Fixed"

param_data_MTCs["Transfer_nMTC_to_MTC_ISSo15_TN",1] <- 0
param_data_MTCs["Transfer_nMTC_to_MTC_ISSo15_TN",3] <- "Fixed"

param_data_MTCs["C_MTC_ISS_o15",1] <- 0
param_data_MTCs["C_MTC_ISS_o15",3] <- "Fixed"

#Change the population matrix to only include people with an ISS of 16 or more
#reset other options to their defaults
TARN_mort_eq <- "Old" 
MTCs_in_mort_risk <- "No"
population_ISS_over16_only <- "Yes"
if(PSA_switch ==0) {

sens_100_spec_10 <- run_simulation(param_data_MTCs, 0, 1, pat_numb, "manual", 1, 0.1,1)
sens_0_spec_90 <- run_simulation(param_data_MTCs, 0, 1, pat_numb, "manual", 0, 0.9,1)

#create a matrix to store all runs
det_analyses <- matrix (nrow = 2, ncol =12)
#name the columns to make analysis easier
colnames(det_analyses) <- c("Sens_DR","Spec_DR", "Number_recieving_MTC_care","proportion_died_before_discharge","proportion_died_between_discharge_and_1_year", "Years_lived",
                            "undiscounted_QALYs", "discounted_QALYs", "undiscounted_Costs", "discounted_Costs", "proportion_ISS_over_16", "proportion_ISS_over_8_under_16")
#name the rows with the appropiate strategy
rownames(det_analyses) <- c("All_MTC", "No_MTC")

det_analyses["All_MTC", ]<- sens_100_spec_10
det_analyses["No_MTC", ]<- sens_0_spec_90

write.csv(det_analyses, "MTC v no MTC.csv")
}

if(PSA_switch==1){
  if(PSA_strat == "MTC"){
  sens_100_spec_10_PSA <- run_simulation(param_data_MTCs, 1, PSA_numb, pat_numb, "manual", 1, 0.1,1)
  write.csv(sens_100_spec_10_PSA, paste(file_location,"PSA results\\sens_100_spec_10_PSA",date,".csv", sep=""))
  use_params_sens_100_spec_10_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_100_spec_10_PSA, "PSA results\\sens_100_spec_10_PSA_params.csv")
  }
  if(PSA_strat == "nMTC"){
  sens_0_spec_90_PSA <- run_simulation(param_data_MTCs, 1, PSA_numb, pat_numb, "manual", 0, 0.9,1)
  write.csv(sens_0_spec_90_PSA, paste(file_location,"PSA results\\sens_0_spec_90_PSA_",date,".csv",  sep=""))
  use_params_sens_0_spec_90_PSA<- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_0_spec_90_PSA, "PSA results\\sens_0_spec_90_PSA_params.csv")
}
}
}