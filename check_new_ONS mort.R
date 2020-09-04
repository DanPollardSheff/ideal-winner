#install.packages("devtools")

library(devtools)
library(MASS)

#Global variables
PSA_switch <- 0
PSA_numb <- 500
pat_numb <- 1000000
days_to_discharge <- 30
days_in_year <- 365.25
time_horizon <- 100
discount_rate_QALYs <- 0.035
discount_rate_costs <- 0.035
Param_export <- 1
Proportion_RR_MTC_ISS_o8_u16 <- 0
TARN_mort_eq <- "Old" # options are new or old. Default is old
MTCs_in_mort_risk <- "No" #options are Yes or no. Relates to whether the mort eq is a composite risk score for a 
#population who has / has not been to an MTC or a population who hasn't gone to an MTC. Default is no, as the
#default for the  mortality equation is the Old TARN equation.
percent_TARN_cases_reported_ISS_o16 <- 1
percent_TARN_cases_reported_ISS_o9_u16 <- 1
population_source <- "Dutch" # Options are UK and Dutch. Dutch is the default
population_ISS_over16_only <- "No" # Options are yes or no. Default is no. 
efficent_life_expectancy <- "No" #Options are Yes or No. Default is yes

test_pat_chars <- "No" #Change this to Yes if you only want to run the base case analysis with patient level results

PSA_strat <- "S100_S1" #Option to make sure that each instance only runs one set of PSAs, as it is computationally intensive
#Options are: S100, S95, S90, S88, S75, S70, S64, S57, S28, MTC, nMTC, S100_S1, S95_S1, S90_S1, S88_S1, S75_S1, S70_S1, S64_S1, S57_S1, S28_S1

PSA_rand_no <-  -99 #random number to determine PSA parameters. #if -99 this will not change the seed after randomly determining the number of patients to run through the model. 
#settings for MATTS phase 1 where first 500 runs 26090100 (after generating pat chars), next 750 runs (ten diagnostic strategies only) 1346, next 750 runs (ten diagnostic strategies only) 330413 

date <- "20200616" #date to append to saved files 

#read in files from the X drive (note not on Git due to confidentiality reasons)
file_location <- "\\\\uosfstore.shef.ac.uk\\shared\\ScHARR\\PR_MATTS\\General\\Health Economics\\Model\\"

param_data <- read.csv(paste(file_location,"parameters.csv", sep=""), row.names=1)
life_tabs <- read.csv(paste(file_location,"ONSlifetables.csv", sep=""))
future_costs <- read.csv(paste(file_location,"lifetime-healthcare-costs.csv", sep=""))

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




#Analysis###################


#do the base casse analysis
#No benefit of MTCs for patients with an ISS between 9 and 15 (inclusive), Old TARN mort equation, 
#no complex adjustment, Dutch population, MTC bext practice tarriff costs
#100% of potentially eligible cases get best practice tarriffs.
param_data_bc <- param_data
Proportion_RR_MTC_ISS_o8_u16 <- 0
TARN_mort_eq <- "Old" 
MTCs_in_mort_risk <- "No" 
population_source <- "Dutch" 
percent_TARN_cases_reported_ISS_o16 <- 1
percent_TARN_cases_reported_ISS_o9_u16 <- 1

#Set a seed, so we always get the same population and parameters
set.seed(1)

pat_chars <- gen_pat_chars(pat_numb, means, covariance, age_tab, gen_tab, ISS_tab, GCS_tab)
parameters <- gen_parameters(PSA_switch, PSA_numb, param_data)

##Aply a random strat for the next step, and change data to match variable names
strat_name <- "manual"
sensitivity <- 0.5
specificity <- 0.5
SOUR <- 1
life_tables <- life_tabs 

###Run through the outcome function code until, I use the ONS functions
#determine triage rule status
#code to be added, this will be specific to our decision rules
pat_chars <- triage_strategies(pat_chars,strat_name, sensitivity, specificity)
#determine initial transport site

#create a vector of positive and negative triage rules
triage_pos <- pat_chars[,"Triage_rule"]==1
triage_neg <- pat_chars[,"Triage_rule"]==0
#create a vector of MTC and non-MTC injuries
ISS_o15 <- pat_chars[,"ISS"] > 15
ISS_u16 <-  pat_chars[,"ISS"] < 16
ISS_u16_o8 <- pat_chars[,"ISS"] < 16 & pat_chars[,"ISS"] > 8 
ISS_u9 <- pat_chars[,"ISS"] < 9

#code here has been commented out, as we are applying sens and spec values

#generate a vector of probabilities of being sent to a major trauma centre (initally)
#Dummy values used, this is adherence to triage rule results and we currently dont know what we will do with this
#prob_MTC <- ifelse (triage_pos==TRUE & ISS_o15 ==TRUE, parameters[SOUR,"P_MTC_Tri_pos_ISS_o15"] ,ifelse (triage_pos==FALSE & ISS_o15 ==TRUE, parameters[SOUR,"P_MTC_Tri_neg_ISS_o15"], ifelse (triage_pos==TRUE & ISS_o15 ==FALSE, parameters[SOUR,"P_MTC_Tri_pos_ISS_u16"], parameters[SOUR,"P_MTC_Tri_neg_ISS_u16"])))
#Determine if there is an event
#create a vector of random numbers equal in length to the prob_MTC vector
#rand_vect <- runif(length(prob_MTC))
#pat_chars[,"MTC"] <- ifelse (rand_vect < prob_MTC,1,0)
#do people receive a transfer to an MTC?
MTC <- pat_chars[,"MTC"]==1
#create a vector of people who went to an nMTC
nMTC <- pat_chars[,"MTC"]==0
#create a vector of probabilities of transfer, irrespective of the site they have been sent to
prob_transfer <- ifelse (triage_pos==TRUE & ISS_o15 ==TRUE, parameters[SOUR,"Transfer_nMTC_to_MTC_ISSo15_TP"] ,ifelse (triage_pos==FALSE & ISS_o15 ==TRUE, parameters[SOUR,"Transfer_nMTC_to_MTC_ISSo15_TN"], ifelse (triage_pos==TRUE & ISS_o15 ==FALSE, parameters[SOUR,"Transfer_nMTC_to_MTC_ISSu16_TP"], parameters[SOUR,"Transfer_nMTC_to_MTC_ISSu16_TN"])))
#create a vector of random numbers equal in length to the prob_transfer vector
rand_vect <- runif(length(prob_transfer))
#Determine whether they would have an event, regardless of MTC status
trans_MTC <- ifelse (rand_vect < prob_transfer,1,0)
#update the patient characteristics, for those patients who where sent to an nMTC
pat_chars[,"MTC_transfer"] <- trans_MTC*nMTC
#update MTC characteristic for those patients who where transfered
pat_chars[,"MTC"] <- ifelse(pat_chars[,"MTC"]==0&pat_chars[,"MTC_transfer"]==1,1,pat_chars[,"MTC"])
#Recalculate whether the patient went the MTC or nMTC
MTC <- pat_chars[,"MTC"]==1
nMTC <- pat_chars[,"MTC"]==0

#in hospital deaths
#Calculate the people who are seriously/moderately injured, have gone to an nMTC, and have not been transfered
ISSo15_nMTC_ntrans <- ISS_o15 & nMTC 
ISSo8_u16_nMTC_ntrans <- ISS_u16_o8 & nMTC 

#record the probability of in hospital death, for severly injured and treatmeted at MTC, severely injured and not treated at MTC,
#step 1: estimate the probability of death from the TARN risk equation
p_death_TARN <- if(TARN_mort_eq == "Old"){
  TARN_old_mort_pred(pat_chars, parameters, SOUR)
}else{TARN_mort_pred(pat_chars, parameters, SOUR)}
#Make sure these adjustments only happen if a composite risk score is used
if(MTCs_in_mort_risk == "Yes"){
  #For patients with an ISS 16 or over
  p_death_hosp_ISSo15_MTC <- p_death_TARN/(parameters[SOUR,"p_MTC_ISS_o15_UK"]+(1-parameters[SOUR,"p_MTC_ISS_o15_UK"])*parameters[SOUR,"RR_p_death_hosp_ISSo15_nMTC"])
  p_death_hosp_ISSo15_nMTC <- p_death_hosp_ISSo15_MTC * as.numeric(parameters[SOUR,"RR_p_death_hosp_ISSo15_nMTC"])
  
  #For patients with an ISS between 9 and 15 inclusive
  #Step 1: Calculate modfied RR (this will be 1 in the base case)
  mod_RR_MTC_ISS_o8_u16 <- as.numeric(1 + (Proportion_RR_MTC_ISS_o8_u16*(parameters[SOUR,"RR_p_death_hosp_ISSo15_nMTC"]-1)))
  
  p_death_hosp_ISSo8_u16_MTC <- p_death_TARN/(parameters[SOUR,"p_MTC_ISS_o15_UK"]+(1-parameters[SOUR,"p_MTC_ISS_o15_UK"])*mod_RR_MTC_ISS_o8_u16)
  p_death_hosp_ISSo8_u16_nMTC <- p_death_hosp_ISSo8_u16_MTC * mod_RR_MTC_ISS_o8_u16
}else{
  #For patients with an ISS 16 or over
  #RR nMTC v MTV = 1/RR is for MTC v nMTC, 
  p_death_hosp_ISSo15_MTC <- p_death_TARN*(1/parameters[SOUR,"RR_p_death_hosp_ISSo15_nMTC"])
  p_death_hosp_ISSo15_nMTC <- p_death_TARN
  
  #For patients with an ISS between 9 and 15 inclusive
  #Step 1: Calculate modfied RR (this will be 1 in the base case)
  mod_RR_MTC_ISS_o8_u16 <- as.numeric(1 + (Proportion_RR_MTC_ISS_o8_u16*(parameters[SOUR,"RR_p_death_hosp_ISSo15_nMTC"]-1)))
  #RR nMTC v MTV = 1/RR is for MTC v nMTC, 
  p_death_hosp_ISSo8_u16_MTC <- p_death_TARN*(1/mod_RR_MTC_ISS_o8_u16)
  p_death_hosp_ISSo8_u16_nMTC <- p_death_TARN
  
}



p_death_hosp_ISSu9 <- p_death_TARN

#create the vector of the probability of death within hospital
p_death_hosp <- ifelse ((ISS_o15==TRUE & MTC == TRUE), p_death_hosp_ISSo15_MTC, ifelse ((ISS_o15==TRUE & MTC == FALSE), p_death_hosp_ISSo15_nMTC,ifelse ((ISS_u16_o8==TRUE & MTC == TRUE), p_death_hosp_ISSo8_u16_MTC, ifelse ((ISS_u16_o8==TRUE & MTC == FALSE), p_death_hosp_ISSo8_u16_nMTC,p_death_hosp_ISSu9))))

#record the p_death in hospital
pat_chars[,"p_death_hosp"] <- p_death_hosp

#create a vector of random numbers, equal in length to the probability of dying in hospital
rand_vect <- runif(length(p_death_hosp))
#determine if the patient has died in hospital
death_hosp <- ifelse (rand_vect < p_death_hosp, 1,0)
#record whether or not the patient has died in hospital 
pat_chars[,"D_bl_disch"] <- death_hosp


#create a vector of people who survived their hospitilisation
alive_disch <- pat_chars[,"D_bl_disch"] == 0



#record the probability of death 
#step 1: estimate the probability of death between hospital discharge and one year post-hospitilisation using US data
p_death_disch_1yr_ISSo15_MTC <- parameters[SOUR,"p_death_y1_ISSo15_MTC"]
p_death_disch_1yr_ISSo15_nMTC <- parameters[SOUR,"p_death_y1_ISSo15_MTC"] * parameters[SOUR,"RR_p_death_y1_nMTC"]

#alter the relative risk by global proportion in the model
mod_RR_MTC_ISS_o8_u16 <- 1 + (Proportion_RR_MTC_ISS_o8_u16*(parameters[SOUR,"RR_p_death_y1_nMTC"]-1))

p_death_disch_1yr_ISSo8_u16_MTC <- parameters[SOUR,"p_death_y1_ISSu16"]
p_death_disch_1yr_ISSo8_u16_nMTC <- parameters[SOUR,"p_death_y1_ISSu16"] * mod_RR_MTC_ISS_o8_u16

p_death_disch_1yr_ISSu9 <- parameters[SOUR,"p_death_y1_ISSu16"]

#create a vector of the values 
p_death_disch_1yr <- alive_disch*ifelse ((ISS_o15==TRUE & MTC == TRUE), p_death_disch_1yr_ISSo15_MTC, ifelse ((ISS_o15==TRUE & MTC == FALSE), p_death_disch_1yr_ISSo15_nMTC,ifelse ((ISS_u16_o8==TRUE & MTC == TRUE), p_death_disch_1yr_ISSo8_u16_MTC, ifelse ((ISS_u16_o8==TRUE & MTC == FALSE), p_death_disch_1yr_ISSo8_u16_nMTC,p_death_disch_1yr_ISSu9))))

#record these probabilities
pat_chars[,"p_death_disch_1yr"] <- p_death_disch_1yr

#create a random vector, which is the length of the probability of dying between discharge and year 1
rand_vect <- runif(length(p_death_disch_1yr))
#Compare the random numbers to the probability of dying between discharge and death
death_disch_1yr <- ifelse(rand_vect < p_death_disch_1yr ,1,0)
#Store these results in the patient characteristics matrix
pat_chars[,"D_disch_1yr"] <- death_disch_1yr

#calculate the lifetables for people who've experienced major trauma. Increase the general population lifetables by hazard ratios

#Step 1 : produce a table of instantaneous rates of general population mortality
#note for all transformations in this section of code, the unit of time is one year
rates_m <- -log(1-life_tables[,2])/1
rates_f <- -log(1-life_tables[,3])/1

#Step 2: apply a hazard ratio to calculate the instanteneous rates, for the population with an ISS > 15
rates_m_ISS_o_15 <- rates_m *as.numeric(parameters[SOUR,"HR_p_death_lm_ISSo15"])
rates_f_ISS_o_15 <- rates_f *as.numeric(parameters[SOUR,"HR_p_death_lm_ISSo15"])

#Step3: create a new life table, based on these instantaneous rates
Life_table_ISS_o_15 <- life_tables
Life_table_ISS_o_15[,2] <- 1 - exp(-rates_m_ISS_o_15*1)
Life_table_ISS_o_15[,3] <- 1 - exp(-rates_f_ISS_o_15*1)

#Step 4: apply a hazard ratio to calculate the instanteneous rates, for the population with an ISS < 15
rates_m_ISS_u_16 <- rates_m *as.numeric(parameters[SOUR,"HR_p_death_lm_ISSu15"])
rates_f_ISS_u_16 <- rates_f *as.numeric(parameters[SOUR,"HR_p_death_lm_ISSu15"])

#Step5: create a new life table, based on these instantaneous rates
Life_table_ISS_u_16 <- life_tables
Life_table_ISS_u_16[,2] <- 1 - exp(-rates_m_ISS_u_16*1)
Life_table_ISS_u_16[,3] <- 1 - exp(-rates_f_ISS_u_16*1)

#Finsihed producing the life tables

#Step 1: produce a logical vector for patients who are alive after one year
alive_1yr <- pat_chars[,"D_bl_disch"]==0 & pat_chars[,"D_disch_1yr"]==0

#produce the names of pat_chars matrix to remind me what to do next
colnames(pat_chars[,])

set.seed(1)
test1 <- life_expectancy_ONS((pat_chars[,"Age"]+1), pat_chars[,"Gender"],life_tabs)
set.seed(1)
test2 <- life_expectancy_ONS2(pat_chars,life_tabs)

summary(test1)
summary(test2)

length(pat_chars[,"Age"])
