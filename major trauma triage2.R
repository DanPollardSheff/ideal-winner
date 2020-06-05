#install.packages("devtools")

library(devtools)
library(MASS)

#Global variables
PSA_switch <- 0
PSA_numb <- 500
pat_numb <- 25000
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

test_pat_chars <- "No" #Change this to Yes if you only want to run the base case analysis with patient level results

PSA_strat <- "S100" #Option to make sure that each instance only runs one set of PSAs, as it is computationally intensive
#Options are: S100, S95, S90, S88, S75, S70, S64, S57, S28, MTC, nMTC

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




#function to sample from the correct distribution, given the selected distribution for PSA parameters

###Define all functions used in the model##########################################################
#Function to sample from PSA distributions for parmeters or if a PSA is not being conducted return the deterministic value. 
#Uses three arguments: parameter 1 (mean/alpha), parameter 2 (SE / beta), parameter 3 (distribution type)
#distributions which are currently implemented are: Beta, log normal, normal, gamma and multivariate normal.

value_selector <- function(a,b,d, PSA_switch, PSA_numb,f,g){
  #set temp to NA, so if there is a bug, it should be obvious
  temp <- NA
  #Sample/calculate the mean of a beta distributed parameter
  if (d == "Beta"){
    if(PSA_switch == 1){
      #If the PSA switch is one produce all samples of this parameter for use in the PSA
      temp <- qbeta(runif(PSA_numb),a,b)
    }else{
      #Otherwise, return the mean value for a deterministic analysis
     temp <- a/(a+b)
    }
  }
  #Sample/calculate the mean of a log normal distributed parameter
  else if (d == "Normal_ln"){
    if(PSA_switch == 1){
    temp <- exp(qnorm(runif(PSA_numb),a,b))
    }else{
    temp <-  exp(a)
    }
 
  }
  #Sample/calculate the mean of a normal distributed parameter  
  else if (d == "Normal"){
    if(PSA_switch == 1){
      temp <- qnorm(runif(PSA_numb),a,b)
    }else{
      temp <-  a
    }
  }
  #Sample/calculate the mean of a normal distributed parameter 
  else if (d == "Gamma"){
    if(PSA_switch == 1){
      temp <- qgamma(runif(PSA_numb),a,rate = 1/b)
    }else{
      temp <-  a*b
    }
  }
  #Sample/calculate the mean of a multi variate normal distributed parameter  
  else if (d == "mult_norm"){
    if(PSA_switch == 1){
      temp <- mvrnorm(runif(PSA_numb),a,b, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    }else{
      temp <-  a
    }
  }
  #else if (d == "Beta_diff"){
    #if(PSA_switch == 1){
      #sample independently from beta distributions
      #temp <-  qbeta(runif(PSA_numb),a,b)
      #temp2 <- qbeta(runif(PSA_numb),f,g)
      #logit transform the samples
      #temp <- log(temp/(1-temp))
      #temp2 <- log(temp2/(1-temp2))
      #calculate means and variances of the logit transformed values
      #m_v1 <- mean(temp)
      #v_v1 <- var(temp)
      #m_v2 <- mean(temp2)
      #v_v2 <- var(temp2)
      #Calculate the mean, variance, and parameters of a gamma distribution
      #m_d <- m_v2 - m_v1
      #v_d <- abs(v_v1 -v_v2)
      #a_d <- m_d^2/v_d
      #b_d <- v_d / m_d
      #sample the differences
      #temp3 <- qgamma(runif(PSA_numb),a_d,1/b_d)
      #if(v_v2 > v_v1){
        #temp4 <- temp
        #temp5 <- temp4 + temp3
      #}else{
        #temp5 <- temp2
        #temp4 <- temp5 - temp3
      #}
      #undo the logit transformation
      #temp <- exp(temp4)/(1+exp(temp4))
      #temp2 <- exp(temp5)/(1+exp(temp5))
      
      #produce a 2xPSA number matrix
      #res <- matrix(nrow = PSA_numb, ncol = 2)
      #res[,1] <- temp
      #res[,2] <- temp2
      #replace temp with the results matrix
      #temp <- res
      
      #else{
        #res <- matrix(nrow = 1, ncol = 2)
        #res[,1] <- a/(a+b)
        #res[,2] <- f/(f+g)
        #replace temp with the results matrix
        #temp <- res
        
      #}
    #}
  else if (d == "Fixed"){
    if(PSA_switch == 1){
      temp <- rep(c(a), times = PSA_numb)
    }else{
      temp <-  a
    }
  }
  #If none of these distributions have been selected, produce an obvious error message
  else{
    temp <- "Error, distribution has not been selected"
  }
  return(temp)
}

##

###Function to calculate life expectancy, based on the probability of dying between age x and x + 1############
#the first argument is a matrix and patient ages
#the second argument is the patient genders
# the third argument is the following three columns of the ONS lifetables: age; qx for men; and, qx for women (see lifetables for defs)

life_expectancy_ONS <- function(age, gender, life_tables){
  #Check that the length of the two vectors (age and gender) is the same
  
  #Generate a temporary vector of the same length to store the results, with all values defaulting to -99.
  #a value of -99 indicates that the age of death has not been determined in the loop
  all_cause_death <- rep(-99, length.out = length(age))
  
  #Loop so that the process is repeated for each patient in the model on the basis of their age
  for (i in 1:length(age)){
    #If the patient is 101 when this function is run, assume they have the remaining life expectancy of a 100 year old for their gender
    #Life expectancies are 103.02 for men and 103.32 for women
    if(age[i]==101){
            #using solver to calibrate the %die before 101*100.5 + %that don't die before 101*unknown life expectancy 
      # is equal to the life expectancy at age 100, in the UK ONS lifetables
      
      all_cause_death[i] = age[i] + ifelse(gender[i]==1,2.02,2.32)
      
    }else{
    
    
    #Loop down so all patients have an all cause death age assigned 
    
    #apply the loop from the baseline age + 1 (the short term model determines outcomes in the first year)
    temp <- as.numeric(age[i])
   
    for(y in temp:100){
      #store a random number to compbl_ageare against the prob of death in each year
      rand <- runif(1)
      if(rand < life_tables[y+1, 3 - gender[i]*1]){
        #if the random number is less than the probability of death from the life, then set their age of death to the 
        #current age plus a sample from a uniform distribution
        #assumes that if they die between age x and x+1, the deaths will happen uniformly accross the year
        all_cause_death[i] <- y+runif(1)
        #stop the inner loop, as the individual has been simulated as dying
        break
      }
    }
  #If the all cause death date is still NA, set the age of death to 100
  if(all_cause_death[i] == -99){
    #If they don't die between age 100 and 101, assume they die at age 103.02 for men and 103.32 for women
    #using solver to calibrate the %die before 101*100.5 + %that don't die before 101*unknown life expectancy 
    # is equal to the life expectancy at age 100
    all_cause_death[i] = 101 + ifelse(gender[i]==1,2.02,2.32)
  }
  }
  }
    return(all_cause_death)
  }

###
  
############Continuous discounting function##################################;
#1st argument is the yearly value, you want to discount (e.g. yearly costs or utilities)
#2nd argument is the time you want the function to start (in years)
#3rd argument is the time you want the function to end (in years)
#4th argument is the periodic discount rate (e.g. 0.035 in the UK for costs and QALYs)

cont_disc <- function(yearly_val, time_start, time_end, discount_rate){
  #from the discount rate, calulcate the instantaenous discount rate
  inst_dr <- log(1+discount_rate)
  
  #calculate the AUC, given the discount rate
  temp<- ((1/-inst_dr)*exp(-inst_dr*time_end))-((1/-inst_dr)*exp(-inst_dr*time_start))
  
  temp2 <- yearly_val*temp
  
  return(temp2)
  
}
  
##

## Function to set all baseline patient characteristics
gen_pat_chars <- function(pat_numb, means, covariance,age_tab, gen_tab, ISS_tab, GCS_tab){
  
  #Create a matrix for all patient characteristics
  test2 = matrix(nrow = pat_numb, ncol = 21)  
  colnames(test2) <- c("ID", "ISS", "Gender", "Age", "Triage_rule", "GCS", "CCI", "Blunt_trauma", "MTC", "MTC_transfer", "D_bl_disch", "D_disch_1yr", "D_1yr_plus", "Age_death", "Life_years", "QALYS", "dQALYS", "Costs", "DCosts", "p_death_hosp", "p_death_disch_1yr")  
  
  #test sampling
  test <- mvrnorm (n = as.numeric(pat_numb), means, covariance)
  
  #if we are checking whether each patient has an ISS of 16 or above, resample any patients with an ISS value that 
  #would be equal to 16
  #change the lookup value dynamically based on whether it is UK or Dutch
  if(population_ISS_over16_only=="Yes"){
  
   if(population_source=="UK"){
    lookupvalue <- ISS_tab[13,4]
  } else{
    lookupvalue <- ISS_tab[14,4]
  }
  
  for (i in 1:length(test[,1])){
    while (test[i,"ISS"] < lookupvalue){
      test[i,] <- mvrnorm(1,means, covariance)
    }
  }
  }
  #Clean the dataset
  #Age
  #for loop to clean to age data
  #do this one as a loop as it will otherwise be horrendous
  #Outerloop does this across each patient
  for (i in 1:length(test[,1])){
    #Inner loop goes between years 
    for (y in 1:99){
      #If the sampled value is smaller than the cutoff value, and we are only assessing the first
      #cutoff, then change the recorded age and break the inner loop
      #Otherwise the sampled value needs to be smaller than the current value been assessed, by larger
      # than the precedding value
      if(y==1){
        if(test[i,1]<=age_tab[y,4]){
          test[i,1] = age_tab[y,1]
          break
        }
      }else{
        if(test[i,1]<=age_tab[y,4]&test[i,1]>age_tab[y-1,4]){
          test[i,1] = age_tab[y,1]
          break
        }
      }
      
    }
    
  }

  
  #Gender
  test[,2] <- ifelse(test[,2]>= gen_tab[1,4],1,0)
  
  #ISS
  #ISS
  #stash the vector somwhere safe
  #ISS
  #stash the vector somwhere safe
if(population_source=="UK"){    
    ISS <- test[,3]
  
  test[,3] <- ifelse(test[,3]<= ISS_tab[1,4], 1, ifelse(test[,3]<= ISS_tab[2,4], 2, ifelse(test[,3]<= ISS_tab[3,4], 4, ifelse(test[,3]<= ISS_tab[4,4],5, ifelse (test[,3]<= ISS_tab[5,4], 6,
                                                                                                                                                                 ifelse(test[,3]<= ISS_tab[6,4],8, ifelse (test[,3]<= ISS_tab[7,4],9, ifelse(test[,3]<= ISS_tab[8,4],10, ifelse(test[,3]<= ISS_tab[9,4], 11, ifelse(test[,3]<= ISS_tab[10,4], 12,
                                                                                                                                                                                                                                                                                                                    ifelse(test[,3]<= ISS_tab[11,4],13, ifelse(test[,3]<= ISS_tab[12,4],14, ifelse(test[,3]<= ISS_tab[13,4],16, ifelse(test[,3]<= ISS_tab[14,4],17, ifelse(test[,3]<= ISS_tab[15,4],18,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                           ifelse(test[,3]<= ISS_tab[16,4],19, ifelse(test[,3]<= ISS_tab[17,4],20, ifelse(test[,3]<= ISS_tab[18,4],21, ifelse(test[,3]<= ISS_tab[19,4],22, ifelse(test[,3]<= ISS_tab[20,4],24,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ifelse(test[,3]<= ISS_tab[21,4],25, ifelse(test[,3]<= ISS_tab[22,4],26, ifelse(test[,3]<= ISS_tab[23,4],27, ifelse(test[,3]<= ISS_tab[24,4],29, ifelse(test[,3]<= ISS_tab[25,4],30,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         ifelse(test[,3]<= ISS_tab[26,4],33, ifelse(test[,3]<= ISS_tab[27,4],34, ifelse(test[,3]<= ISS_tab[28,4],35, ifelse(test[,3]<= ISS_tab[29,4],36, ifelse(test[,3]<= ISS_tab[30,4],38,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                ifelse(test[,3]<= ISS_tab[31,4],41, ifelse(test[,3]<= ISS_tab[32,4],42, ifelse(test[,3]<= ISS_tab[33,4],43, ifelse(test[,3]<= ISS_tab[34,4],45, ifelse(test[,3]<= ISS_tab[35,4],50,54)))))))))))))))))))))))))))))))))))
  #looks good
  
  
  GCS <- test[,4]
  
  test[,4] <- ifelse(test[,4]<= GCS_tab[1,4],3,ifelse(test[,4]<= GCS_tab[2,4],4,ifelse(test[,4]<= GCS_tab[3,4],5,ifelse(test[,4]<= GCS_tab[4,4],6,ifelse(test[,4]<= GCS_tab[5,4],7,ifelse(test[,4]<= GCS_tab[6,4],8,ifelse(test[,4]<= GCS_tab[7,4],9,ifelse(test[,4]<= GCS_tab[8,4],10,ifelse(test[,4]<= GCS_tab[9,4],11,ifelse(test[,4]<= GCS_tab[10,4],12,ifelse(test[,4]<= GCS_tab[11,4],13,ifelse(test[,4]<= GCS_tab[12,4],14,15))))))))))))
  
  
  #Base the blunt trauma / not on the:
  #christensen et al data (jan 1st 2000 - 31st Dec 2005)
  #blunt trauma,  n = 1,365
  #pent trauma, n = 35,564
  
  pent_trauma <- ifelse(runif(length(test[,4]))< 1365/(35564+1365),1,0)
  
  test <- cbind(test, pent_trauma)
  }
  else{
    #Estimate ISS from our data
    ISS <- test[,3]
    
    test[,3] <- ifelse(test[,3]<= ISS_tab[1,4], ISS_tab[1,1], ifelse(test[,3]<= ISS_tab[2,4], ISS_tab[2,1], ifelse(test[,3]<= ISS_tab[3,4], ISS_tab[3,1], ifelse(test[,3]<= ISS_tab[4,4],ISS_tab[4,1], ifelse (test[,3]<= ISS_tab[5,4], ISS_tab[5,1],
                                                                                                                                                                                                               ifelse(test[,3]<= ISS_tab[6,4],ISS_tab[6,1], ifelse (test[,3]<= ISS_tab[7,4],ISS_tab[7,1], ifelse(test[,3]<= ISS_tab[8,4],ISS_tab[8,1], ifelse(test[,3]<= ISS_tab[9,4], ISS_tab[9,1], ifelse(test[,3]<= ISS_tab[10,4],  ISS_tab[10,1],
                                                                                                                                                                                                                                                                                                                                                                                                            ifelse(test[,3]<= ISS_tab[11,4],ISS_tab[11,1], ifelse(test[,3]<= ISS_tab[12,4],ISS_tab[12,1], ifelse(test[,3]<= ISS_tab[13,4],ISS_tab[13,1], ifelse(test[,3]<= ISS_tab[14,4],ISS_tab[14,1], ifelse(test[,3]<= ISS_tab[15,4],ISS_tab[15,1],
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               ifelse(test[,3]<= ISS_tab[16,4],ISS_tab[16,1], ifelse(test[,3]<= ISS_tab[17,4],ISS_tab[17,1], ifelse(test[,3]<= ISS_tab[18,4],ISS_tab[18,1], ifelse(test[,3]<= ISS_tab[19,4],ISS_tab[19,1], ifelse(test[,3]<= ISS_tab[20,4],ISS_tab[20,1],
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ifelse(test[,3]<= ISS_tab[21,4],ISS_tab[21,1], ifelse(test[,3]<= ISS_tab[22,4],ISS_tab[22,1], ifelse(test[,3]<= ISS_tab[23,4],ISS_tab[23,1], ifelse(test[,3]<= ISS_tab[24,4],ISS_tab[24,1], ifelse(test[,3]<= ISS_tab[25,4],ISS_tab[25,1],
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ifelse(test[,3]<= ISS_tab[26,4],ISS_tab[26,1], ifelse(test[,3]<= ISS_tab[27,4],ISS_tab[27,1], ifelse(test[,3]<= ISS_tab[28,4],ISS_tab[28,1], ifelse(test[,3]<= ISS_tab[29,4],ISS_tab[29,1], ifelse(test[,3]<= ISS_tab[30,4],ISS_tab[30,1],
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        ifelse(test[,3]<= ISS_tab[31,4],ISS_tab[31,1], ifelse(test[,3]<= ISS_tab[32,4],ISS_tab[32,1], ifelse(test[,3]<= ISS_tab[33,4],ISS_tab[33,1], ifelse(test[,3]<= ISS_tab[34,4],ISS_tab[34,1], ifelse(test[,3]<= ISS_tab[35,4],ISS_tab[35,1],ifelse(test[,3]<= ISS_tab[36,4],ISS_tab[36,1],ifelse(test[,3]<= ISS_tab[37,4],ISS_tab[37,1],ifelse(test[,3]<= ISS_tab[38,4],ISS_tab[38,1],ifelse(test[,3]<= ISS_tab[39,4],ISS_tab[39,1],ISS_tab[40,1])))))))))))))))))))))))))))))))))))))))
    #GCS
    GCS <- test[,4]
    
    test[,4] <- ifelse(test[,4]<= GCS_tab[1,4],GCS_tab[1,1],ifelse(test[,4]<= GCS_tab[2,4],GCS_tab[2,1],ifelse(test[,4]<= GCS_tab[3,4],GCS_tab[3,1],ifelse(test[,4]<= GCS_tab[4,4],GCS_tab[4,1],ifelse(test[,4]<= GCS_tab[5,4],GCS_tab[5,1],ifelse(test[,4]<= GCS_tab[6,4],GCS_tab[6,1],ifelse(test[,4]<= GCS_tab[7,4],GCS_tab[7,1],ifelse(test[,4]<= GCS_tab[8,4],10,ifelse(test[,4]<= GCS_tab[9,4],GCS_tab[9,1],ifelse(test[,4]<= GCS_tab[10,4],GCS_tab[10,1],ifelse(test[,4]<= GCS_tab[11,4],GCS_tab[11,1],ifelse(test[,4]<= GCS_tab[12,4],GCS_tab[12,1],GCS_tab[13,1]))))))))))))
    
    
    
    
    test[,5] <- ifelse(test[,5]<=blunt_tab[1,4],0,1)
  }
  colnames(test) <- c("Age", "Gender", "ISS", "GCS", "Blunt_trauma")
  #write all data into the test 2 matrix
  #add in pat chars
  test2[,"ID"] <- 1:nrow(test2)
  test2[,"ISS"] <- test[,"ISS"]
  test2[,"Age"] <- test[,"Age"]
  test2[,"Gender"] <- ifelse(test[,"Gender"]==1,0,1)
  test2[,"GCS"] <- test[,"GCS"]
  test2[,"Blunt_trauma"] <- test[,"Blunt_trauma"]
  
  #As we have no data on mCCI, set everyone to have a missing mCCI
  test2[,"CCI"] <- -99
  
  return(test2)
}


#Function to set the value of all parameters in the simulation at the start of the model run.
#Parameter 1 is a switch variable to determine whether or not a PSA is being run
#Second parameter is 
gen_parameters <- function(PSA_switch,PSA_numb, parameters){

  #Set up a matrix to store all parameter values
  param_matrix <- matrix(nrow = ifelse(PSA_switch==1,PSA_numb,1), ncol = nrow(parameters))
  
  colnames(param_matrix) <- c("P_MTC_Tri_pos_ISS_o15","P_MTC_Tri_neg_ISS_o15","P_MTC_Tri_pos_ISS_u16","P_MTC_Tri_neg_ISS_u16","Transfer_nMTC_to_MTC_ISSo15_TP","Transfer_nMTC_to_MTC_ISSo15_TN", "Transfer_nMTC_to_MTC_ISSu16_TP", "Transfer_nMTC_to_MTC_ISSu16_TN", "p_death_hosp_ISSo15_MTC",
                              "RR_p_death_hosp_ISSo15_nMTC", "p_death_hosp_ISSu16", "p_death_y1_ISSo15_MTC", "RR_p_death_y1_nMTC", "p_death_y1_ISSu16", "HR_p_death_lm_ISSo15", "HR_p_death_lm_ISSu15", 
                              "U_ISS_o15_MTC","U_ISS_o15_nMTC","U_ISS_u16_o8","Umult_ISS_u9", "U_genpop_cons", "U_genpop_male", "U_genpop_age", "U_genpop_age_squared", "p_death_hosp_TARN_sqrt_ISS", "p_death_hosp_TARN_ln_ISS",
                              "p_death_hosp_TARN_GCS_3", "p_death_hosp_TARN_GCS_4_5", "p_death_hosp_TARN_GCS_6_8", "p_death_hosp_TARN_GCS_9_12", "p_death_hosp_TARN_GCS_13_14", "p_death_hosp_TARN_GCS_intubated",
                              "p_death_hosp_TARN_CCI_unknown", "p_death_hosp_TARN_CCI_1_5", "p_death_hosp_TARN_CCI_6_10", "p_death_hosp_TARN_CCI_o_10", "p_death_hosp_TARN_age_0_5", "p_death_hosp_TARN_age_6_10",
                              "p_death_hosp_TARN_age_11_15", "p_death_hosp_TARN_age_45_54", "p_death_hosp_TARN_age_55_64", "p_death_hosp_TARN_age_65_74", "p_death_hosp_TARN_gen_f", "p_death_hosp_TARN_age_o_75" ,"p_death_hosp_TARN_age_0_5_gen_f", "p_death_hosp_TARN_age_6_10_gen_f",
                              "p_death_hosp_TARN_age_11_15_gen_f", "p_death_hosp_TARN_age_45_54_gen_f", "p_death_hosp_TARN_age_55_64_gen_f", "p_death_hosp_TARN_age_65_74_gen_f", "p_death_hosp_TARN_age_75_plus_gen_f", 
                              "p_death_hosp_TARN_cons", "p_MTC_ISS_o15_UK", "C_MTC_ISS_o8_u16", "C_MTC_ISS_o15", "C_bluntt_ISS_U10", "C_bluntt_ISS_U17_O_9", "C_bluntt_ISS_U26_O16", "C_bluntt_ISS_O25", "C_pent_ISS_O0_U10",
                              "C_pent_ISS_O9_U16", "C_pent_ISS_O15_U25", "C_pent_ISS_O24_U34", "C_pent_ISS_O34", "C_disch_6m", "C_additional_ambulance", "Increase_lifetime_cost_ISS_o15", "Increase_lifetime_cost_ISS_u15", "TARN_old_Age_0_5", "TARN_old_Age_6_10", "TARN_old_Age_11_15",
                              "TARN_old_Age_45_54", "TARN_old_Age_55_64", "TARN_old_Age_65_75", "TARN_old_Age_over_75", "TARN_old_GCS_9_12", "TARN_old_GCS_6_8", "TARN_old_GCS_4_5", "TARN_old_GCS_3", "TARN_old_GCS_intubated","TARN_old_ISS_SQRT", "TARN_old_ISS_LN", "TARN_old_female",
                              "TARN_old_female_age_0_5", "TARN_old_female_age_6_10", "TARN_old_female_age_11_15", "TARN_old_female_age_45_54", "TARN_old_female_age_55_64", "TARN_old_female_age_65_75", "TARN_old_female_age_75_plus", "TARN_old_constant")
  
  #First parameter, which is the probability of being transfered to an MTC from a non MTC, if the patient's ISS >15 and they have a positive triage rule
  #Step 1, record the name of the parameter in a temproary variable 
  t <- "Transfer_nMTC_to_MTC_ISSo15_TP"
  #Step 2, record the value of the parameter in the simulation
  Transfer_nMTC_to_MTC_ISSo15_TP <-value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- Transfer_nMTC_to_MTC_ISSo15_TP
  
  #Second parameter, which is the probability of being transfered to an MTC from a non MTC, if the patient's ISS >15 and they have a negative triage rule
  t <- "Transfer_nMTC_to_MTC_ISSo15_TN"
  #Step 2, record the value of the parameter in the simulation
  Transfer_nMTC_to_MTC_ISSo15_TN <-value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- Transfer_nMTC_to_MTC_ISSo15_TN
  
  #Thrid parameter, which is he probability of being transfered to an MTC from a non MTC, if the patient's ISS <16 and they have a positive triage rule
  t <- "Transfer_nMTC_to_MTC_ISSu16_TP"
  #Step 2, record the value of the parameter in the simulation
  Transfer_nMTC_to_MTC_ISSu16_TP <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- Transfer_nMTC_to_MTC_ISSu16_TP
  
  #Fourth parameter, which is the probability being transfered to an MTC from a non MTC, if the patient's ISS <16 and they have a negative triage rule
  t <- "Transfer_nMTC_to_MTC_ISSu16_TN"
  #Step 2, record the value of the parameter in the simulation
  Transfer_nMTC_to_MTC_ISSu16_TN <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- Transfer_nMTC_to_MTC_ISSu16_TN
  
  #Probability of within hospital death, if a patient's ISS > 15 and they recieve MTC care
  t <- "p_death_hosp_ISSo15_MTC"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_ISSo15_MTC <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_ISSo15_MTC
  
  
  #RR of within hospital death if ISS > 15 and they do not recieve MTC care
  t <- "RR_p_death_hosp_ISSo15_nMTC"
  #Step 2, record the value of the parameter in the simulation
  RR_p_death_hosp_ISSo15_nMTC <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- RR_p_death_hosp_ISSo15_nMTC
  
  
  #probability of dying in hospital, if ISS < 16
  t <- "p_death_hosp_ISSu16"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_ISSu16 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_ISSu16
  
  #probability of dying between discharge and one year, if the patient's ISS >15 and they recieve treatment at a major trauma centre
  t <- "p_death_y1_ISSo15_MTC"
  #Step 2, record the value of the parameter in the simulation
  p_death_y1_ISSo15_MTC <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_y1_ISSo15_MTC
  
  #Relative risk of death between discharge and one year, if the patient's ISS > 15 and they recieved non-major trauma centre care
  t <- "RR_p_death_y1_nMTC"
  #Step 2, record the value of the parameter in the simulation
  RR_p_death_y1_nMTC <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- RR_p_death_y1_nMTC
  
  #probability of dying between discharge and one year, if the patient's ISS < 16 
  t <- "p_death_y1_ISSu16"
  #Step 2, record the value of the parameter in the simulation
  p_death_y1_ISSu16 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_y1_ISSu16
  
  #Hazard ratio of death for long term modelling, if a patient's ISS > 15
  t<- "HR_p_death_lm_ISSo15"
  #Step 2, record the value of the parameter in the simulation
  HR_p_death_lm_ISSo15 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- HR_p_death_lm_ISSo15
  
  #Hazard ratio of death for long term modelling, if a patient's ISS < 15
  t<- "HR_p_death_lm_ISSu15"
  #Step 2, record the value of the parameter in the simulation
  HR_p_death_lm_ISSu15 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- HR_p_death_lm_ISSu15
  
  #Constant term in the formula for determining age, gender matched utility
  t<- "U_genpop_cons"
  #Step 2, record the value of the parameter in the simulation
  U_genpop_cons <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- U_genpop_cons
  
  
  #Term for the effect of being male for determining age and gender matched utility
  t<- "U_genpop_male"
  #Step 2, record the value of the parameter in the simulation
  U_genpop_male <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- U_genpop_male
  
  #Term for the effect of age for determining age and gender matched utility
  t<- "U_genpop_age"
  #Step 2, record the value of the parameter in the simulation
  U_genpop_age <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- U_genpop_age
  
  #Term for the effect of age squared for determining age and gender matched utility
  t<- "U_genpop_age_squared"
  #Step 2, record the value of the parameter in the simulation
  U_genpop_age_squared <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- U_genpop_age_squared
  
  #Utility score for someone with an ISS over 15, who is sent to a major trama centre
  t<- "U_ISS_o15_MTC"
  #Step 2, record the value of the parameter in the simulation
  U_ISS_o15_MTC <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 4, record the parameter value
  param_matrix[,t] <- U_ISS_o15_MTC
  
  #Utility score for someone with an ISS over 15, who is sent to a non major trama centre
  t<- "U_ISS_o15_nMTC"
  #Step 2, record the value of the parameter in the simulation
  U_ISS_o15_nMTC <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 4, record the parameter value
  param_matrix[,t] <- U_ISS_o15_nMTC
  
  #Utility score for someone with an ISS under 16 but over 8
  t<- "U_ISS_u16_o8"
  #Step 2, record the value of the parameter in the simulation
  U_ISS_u16_o8<- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- U_ISS_u16_o8
  
  #Utility score for someone with an ISS under 16 but over 8
  t<- "Umult_ISS_u9"
  #Step 2, record the value of the parameter in the simulation
  Umult_ISS_u9<- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- Umult_ISS_u9
  
  
  #record TARN parameters
  t<- "p_death_hosp_TARN_sqrt_ISS"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_sqrt_ISS <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_sqrt_ISS
  
  t<- "p_death_hosp_TARN_ln_ISS"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_ln_ISS <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_ln_ISS
  
  t<- "p_death_hosp_TARN_GCS_3"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_GCS_3 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_GCS_3
  
  t<- "p_death_hosp_TARN_GCS_4_5"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_GCS_4_5 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_GCS_4_5
  
  t<- "p_death_hosp_TARN_GCS_6_8"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_GCS_6_8 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_GCS_6_8
  
  t<- "p_death_hosp_TARN_GCS_9_12"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_GCS_9_12 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_GCS_9_12
  
  t<- "p_death_hosp_TARN_GCS_13_14"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_GCS_13_14 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_GCS_13_14
  
  t<- "p_death_hosp_TARN_GCS_intubated"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_GCS_intubated <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_GCS_intubated
  
  t<- "p_death_hosp_TARN_CCI_1_5"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_CCI_1_5 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_CCI_1_5
  
  t<- "p_death_hosp_TARN_CCI_6_10"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_CCI_6_10 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_CCI_6_10
  
  t<- "p_death_hosp_TARN_CCI_unknown"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_CCI_unknown <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_CCI_unknown
  
  t<- "p_death_hosp_TARN_CCI_o_10"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_CCI_o_10 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_CCI_o_10
  
  t<- "p_death_hosp_TARN_age_0_5"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_age_0_5 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_age_0_5
  
  t<- "p_death_hosp_TARN_age_6_10"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_age_6_10 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_age_6_10
  
  t<- "p_death_hosp_TARN_age_11_15"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_age_11_15 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_age_11_15
  
  t<- "p_death_hosp_TARN_age_45_54"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_age_45_54 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_age_45_54
  
  t<- "p_death_hosp_TARN_age_55_64"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_age_55_64 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_age_55_64
  
  t<- "p_death_hosp_TARN_age_65_74"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_age_65_74 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_age_65_74
  
  t<- "p_death_hosp_TARN_age_o_75"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_age_o_75 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_age_o_75
  
  t<- "p_death_hosp_TARN_gen_f"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_gen_f <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_gen_f
  
  t<- "p_death_hosp_TARN_age_0_5_gen_f"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_age_0_5_gen_f <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_age_0_5_gen_f
  
  t<- "p_death_hosp_TARN_age_6_10_gen_f"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_age_6_10_gen_f <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_age_6_10_gen_f
  
  t<- "p_death_hosp_TARN_age_11_15_gen_f"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_age_11_15_gen_f <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_age_11_15_gen_f
  
  t<- "p_death_hosp_TARN_age_45_54_gen_f"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_age_45_54_gen_f <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_age_45_54_gen_f
  
  t<- "p_death_hosp_TARN_age_55_64_gen_f"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_age_55_64_gen_f <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_age_55_64_gen_f
  
  t<- "p_death_hosp_TARN_age_65_74_gen_f"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_age_65_74_gen_f <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_age_65_74_gen_f
  
  t<- "p_death_hosp_TARN_age_75_plus_gen_f"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_age_75_plus_gen_f <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_age_75_plus_gen_f
  
  t<- "p_death_hosp_TARN_cons"
  #Step 2, record the value of the parameter in the simulation
  p_death_hosp_TARN_cons <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_death_hosp_TARN_cons
  
  t<- "p_MTC_ISS_o15_UK"
  #Step 2, record the value of the parameter in the simulation
  p_MTC_ISS_o15_UK <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- p_MTC_ISS_o15_UK
  
  t<- "C_MTC_ISS_o8_u16"
  #Step 2, record the value of the parameter in the simulation
  C_MTC_ISS_o8_u16 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- C_MTC_ISS_o8_u16
  
  t<- "C_MTC_ISS_o15"
  #Step 2, record the value of the parameter in the simulation
  C_MTC_ISS_o15 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- C_MTC_ISS_o15
  
  t<- "C_bluntt_ISS_U10"
  #Step 2, record the value of the parameter in the simulation
  C_bluntt_ISS_U10 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- C_bluntt_ISS_U10
  
  t<- "C_bluntt_ISS_U17_O_9"
  #Step 2, record the value of the parameter in the simulation
  C_bluntt_ISS_U17_O_9 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- C_bluntt_ISS_U17_O_9
  
  t<- "C_bluntt_ISS_U26_O16"
  #Step 2, record the value of the parameter in the simulation
  C_bluntt_ISS_U26_O16 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- C_bluntt_ISS_U26_O16
  
  t<- "C_bluntt_ISS_O25"
  #Step 2, record the value of the parameter in the simulation
  C_bluntt_ISS_O25 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- C_bluntt_ISS_O25
  
  t<- "C_pent_ISS_O0_U10"
  #Step 2, record the value of the parameter in the simulation
  C_pent_ISS_O0_U10 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- C_pent_ISS_O0_U10
  
  t<- "C_pent_ISS_O9_U16"
  #Step 2, record the value of the parameter in the simulation
  C_pent_ISS_O9_U16 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- C_pent_ISS_O9_U16
  
  t<- "C_pent_ISS_O15_U25"
  #Step 2, record the value of the parameter in the simulation
  C_pent_ISS_O15_U25 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- C_pent_ISS_O15_U25
  
  t<- "C_pent_ISS_O24_U34"
  #Step 2, record the value of the parameter in the simulation
  C_pent_ISS_O24_U34 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- C_pent_ISS_O24_U34
  
  t<- "C_pent_ISS_O34"
  #Step 2, record the value of the parameter in the simulation
  C_pent_ISS_O34 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- C_pent_ISS_O34
  
  t<- "C_disch_6m"
  #Step 2, record the value of the parameter in the simulation
  C_disch_6m <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- C_disch_6m
  
  t<- "C_additional_ambulance"
  #Step 2, record the value of the parameter in the simulation
  C_additional_ambulance <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- C_additional_ambulance
  
  
  t<- "P_MTC_Tri_pos_ISS_o15"
  #Step 2, record the value of the parameter in the simulation
  P_MTC_Tri_pos_ISS_o15 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- P_MTC_Tri_pos_ISS_o15
  
  t<- "P_MTC_Tri_neg_ISS_o15"
  #Step 2, record the value of the parameter in the simulation
  P_MTC_Tri_neg_ISS_o15 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- P_MTC_Tri_neg_ISS_o15
  
  t<- "P_MTC_Tri_pos_ISS_u16"
  #Step 2, record the value of the parameter in the simulation
  P_MTC_Tri_pos_ISS_u16 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- P_MTC_Tri_pos_ISS_u16
  
  t<- "P_MTC_Tri_neg_ISS_u16"
  #Step 2, record the value of the parameter in the simulation
  P_MTC_Tri_neg_ISS_u16 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- P_MTC_Tri_neg_ISS_u16
  
  t<- "Increase_lifetime_cost_ISS_o15"
  #Step 2, record the value of the parameter in the simulation
  Increase_lifetime_cost_ISS_o15 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- Increase_lifetime_cost_ISS_o15
  
  t<- "Increase_lifetime_cost_ISS_u15"
  #Step 2, record the value of the parameter in the simulation
  Increase_lifetime_cost_ISS_u15 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- Increase_lifetime_cost_ISS_u15
  
  t<- "TARN_old_Age_0_5"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_Age_0_5 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_Age_0_5
  
  t<- "TARN_old_Age_6_10"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_Age_6_10 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_Age_6_10
  
  t<- "TARN_old_Age_11_15"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_Age_11_15 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_Age_11_15
  
  t<- "TARN_old_Age_45_54"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_Age_45_54 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_Age_45_54
  
  t<- "TARN_old_Age_55_64"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_Age_55_64 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_Age_55_64
  
  t<- "TARN_old_Age_65_75"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_Age_65_75 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_Age_65_75
  
  t<- "TARN_old_Age_over_75"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_Age_over_75 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_Age_over_75
  
  
  t<- "TARN_old_GCS_9_12"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_GCS_9_12 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_GCS_9_12
  
  t<- "TARN_old_GCS_6_8"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_GCS_6_8 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_GCS_6_8
  
  t<- "TARN_old_GCS_4_5"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_GCS_4_5 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_GCS_4_5

  t<- "TARN_old_GCS_3"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_GCS_3 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_GCS_3
  
  t<- "TARN_old_GCS_intubated"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_GCS_intubated <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_GCS_intubated
  
  t<- "TARN_old_ISS_SQRT"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_ISS_SQRT <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_ISS_SQRT
  
  t<- "TARN_old_ISS_LN"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_ISS_LN <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_ISS_LN
  
  t<- "TARN_old_female"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_female <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_female
  
  t<- "TARN_old_female_age_0_5"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_female_age_0_5 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_female_age_0_5
  
  t<- "TARN_old_female_age_6_10"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_female_age_6_10 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_female_age_6_10
  
  t<- "TARN_old_female_age_11_15"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_female_age_11_15 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_female_age_11_15
  
  t<- "TARN_old_female_age_45_54"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_female_age_45_54 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_female_age_45_54
  
  t<- "TARN_old_female_age_55_64"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_female_age_55_64 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_female_age_55_64
  
  t<- "TARN_old_female_age_65_75"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_female_age_65_75 <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_female_age_65_75
  
  t<- "TARN_old_female_age_75_plus"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_female_age_75_plus <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_female_age_75_plus
  
  t<- "TARN_old_constant"
  #Step 2, record the value of the parameter in the simulation
  TARN_old_constant <- value_selector(as.numeric(parameters[t,1]),as.numeric(parameters[t,2]),parameters[t,3],PSA_switch,PSA_numb)
  #Step 3, record the parameter value
  param_matrix[,t] <- TARN_old_constant
  
  return(param_matrix)

}


TARN_mort_pred <- function(pat_chars, parameters, SOUR){
  GCS_3 <- pat_chars[,"GCS"]==3
  GCS_4_to_5 <- pat_chars[,"GCS"]>=4 & pat_chars[,"GCS"]<=5
  GCS_6_to_8 <- pat_chars[,"GCS"]>=6 & pat_chars[,"GCS"]<=8
  GCS_9_to_12 <- pat_chars[,"GCS"]>=9 & pat_chars[,"GCS"]<=12
  GCS_13_to_14 <- pat_chars[,"GCS"]>=13 & pat_chars[,"GCS"]<=14
  GCS_intubated <- pat_chars[,"GCS"]== -99
  
  CCI_1_to_5 <- pat_chars[,"CCI"]>=1 & pat_chars[,"CCI"]<=5
  CCI_6_to_10 <- pat_chars[,"CCI"]>=6 & pat_chars[,"CCI"]<=10
  CCI_over_10 <- pat_chars[,"CCI"]> 10
  CCI_unknown <- pat_chars[,"CCI"]==-99
  
  Age_0_to_5 <- pat_chars[,"Age"]>=0 & pat_chars[,"Age"]<=5
  Age_6_to_10 <- pat_chars[,"Age"]>=6 & pat_chars[,"Age"]<=10
  Age_11_to_15 <- pat_chars[,"Age"]>=11 & pat_chars[,"Age"]<=15
  Age_45_to_54 <- pat_chars[,"Age"]>=45 & pat_chars[,"Age"]<=54
  Age_55_to_64 <- pat_chars[,"Age"]>=55 & pat_chars[,"Age"]<=64
  Age_64_to_74 <- pat_chars[,"Age"]>=65 & pat_chars[,"Age"]<=74
  Age_over_75 <- pat_chars[,"Age"]>= 75
  
  Gender_f <- pat_chars[,"Gender"]==0
  
  ISS_0 <- pat_chars[,"ISS"]==0
  ISS_above_0 <- pat_chars[,"ISS"]>0
  
  #Calculate the fitted value from the TARN mortality prediction equation
  FVs <- parameters[SOUR,"p_death_hosp_TARN_cons"] + parameters[SOUR,"p_death_hosp_TARN_sqrt_ISS"]*(sqrt(10/pat_chars[,"ISS"]) - 0.8686) + 
    parameters[SOUR,"p_death_hosp_TARN_ln_ISS"]*(log(pat_chars[,"ISS"]/10) - 0.2817)+ parameters[SOUR, "p_death_hosp_TARN_GCS_3"]*GCS_3+
    parameters[SOUR,"p_death_hosp_TARN_GCS_4_5"]*GCS_4_to_5 +  parameters[SOUR,"p_death_hosp_TARN_GCS_6_8"]*GCS_6_to_8+
    parameters[SOUR,"p_death_hosp_TARN_GCS_9_12"]*GCS_9_to_12 + parameters[SOUR,"p_death_hosp_TARN_GCS_13_14"]*GCS_13_to_14+ 
    parameters[SOUR,"p_death_hosp_TARN_GCS_intubated"]*GCS_intubated + parameters[SOUR,"p_death_hosp_TARN_CCI_unknown"]*CCI_unknown+
    parameters[SOUR,"p_death_hosp_TARN_CCI_1_5"]*CCI_1_to_5+ parameters[SOUR,"p_death_hosp_TARN_CCI_6_10"]*CCI_6_to_10+
    parameters[SOUR,"p_death_hosp_TARN_CCI_o_10"]*CCI_over_10+
    parameters[SOUR,"p_death_hosp_TARN_age_0_5"]*Age_0_to_5 + parameters[SOUR,"p_death_hosp_TARN_age_6_10"]*Age_6_to_10+
    parameters[SOUR,"p_death_hosp_TARN_age_11_15"]*Age_11_to_15 + parameters[SOUR, "p_death_hosp_TARN_age_45_54"]*Age_45_to_54+
    parameters[SOUR,"p_death_hosp_TARN_age_55_64"]*Age_55_to_64 + parameters[SOUR,"p_death_hosp_TARN_age_65_74"]*Age_64_to_74+
    parameters[SOUR,"p_death_hosp_TARN_age_o_75"]*Age_over_75 + parameters[SOUR,"p_death_hosp_TARN_gen_f"]*Gender_f+ 
    parameters[SOUR,"p_death_hosp_TARN_age_0_5_gen_f"]*Age_0_to_5*Gender_f + parameters[SOUR,"p_death_hosp_TARN_age_6_10_gen_f"]*Age_6_to_10*Gender_f+ 
    parameters[SOUR,"p_death_hosp_TARN_age_11_15_gen_f"]*Age_11_to_15*Gender_f + parameters[SOUR,"p_death_hosp_TARN_age_45_54_gen_f"]*Age_45_to_54*Gender_f+
    parameters[SOUR,"p_death_hosp_TARN_age_55_64_gen_f"]*Age_55_to_64*Gender_f + parameters[SOUR,"p_death_hosp_TARN_age_65_74_gen_f"]*Age_64_to_74*Gender_f+
    parameters[SOUR,"p_death_hosp_TARN_age_75_plus_gen_f"]*Age_over_75*Gender_f 
  
  #Apply the inverse logit function to the fitted values  
  p_surv_temp <- 1/(1+exp(-FVs))
  
  # set the probability of death to 0 if the ISS is 0
  p_death <- ifelse(ISS_0==TRUE, 0, 1- p_surv_temp)
  
  
  
  return(p_death)
  
}

TARN_old_mort_pred <- function(pat_chars, parameters, SOUR){
  
  ### calculate the subgroups for prediction formula
  Age_0_to_5 <- pat_chars[,"Age"]>=0 & pat_chars[,"Age"]<=5
  Age_6_to_10 <- pat_chars[,"Age"]>=6 & pat_chars[,"Age"]<=10
  Age_11_to_15 <- pat_chars[,"Age"]>=11 & pat_chars[,"Age"]<=15
  Age_45_to_54 <- pat_chars[,"Age"]>=45 & pat_chars[,"Age"]<=54
  Age_55_to_64 <- pat_chars[,"Age"]>=55 & pat_chars[,"Age"]<=64
  Age_65_to_75 <- pat_chars[,"Age"]>=65 & pat_chars[,"Age"]<=75
  Age_over_75 <- pat_chars[,"Age"] >= 76
  
  Gender_f <- pat_chars[,"Gender"]==0
  
  GCS_3 <- pat_chars[,"GCS"]==3
  GCS_4_to_5 <- pat_chars[,"GCS"]>=4 & pat_chars[,"GCS"]<=5
  GCS_6_to_8 <- pat_chars[,"GCS"]>=6 & pat_chars[,"GCS"]<=8
  GCS_9_to_12 <- pat_chars[,"GCS"]>=9 & pat_chars[,"GCS"]<=12
  GCS_13_to_14 <- pat_chars[,"GCS"]>=13 & pat_chars[,"GCS"]<=14
  GCS_intubated <- pat_chars[,"GCS"]== -99
  
  
  ISS_0 <- pat_chars[,"ISS"]==0
  ISS_above_0 <- pat_chars[,"ISS"]>0
  
  FVs <- parameters[SOUR,"TARN_old_constant"] + parameters[SOUR,"TARN_old_ISS_SQRT"]*(sqrt(10/pat_chars[,"ISS"]) - 0.953) + 
    parameters[SOUR,"TARN_old_ISS_LN"]*(log(pat_chars[,"ISS"]/10) - 0.0968)+ parameters[SOUR,"TARN_old_GCS_3"]*GCS_3+
    parameters[SOUR,"TARN_old_GCS_4_5"]*GCS_4_to_5 +  parameters[SOUR,"TARN_old_GCS_6_8"]*GCS_6_to_8+ parameters[SOUR,"TARN_old_GCS_9_12"]*GCS_9_to_12 +  
    parameters[SOUR,"TARN_old_GCS_intubated"]*GCS_intubated + parameters[SOUR,"TARN_old_Age_0_5"]*Age_0_to_5 + 
    parameters[SOUR,"TARN_old_Age_6_10"]*Age_6_to_10+ parameters[SOUR,"TARN_old_Age_11_15"]*Age_11_to_15 + parameters[SOUR,"TARN_old_Age_45_54"]*Age_45_to_54+
    parameters[SOUR,"TARN_old_Age_55_64"]*Age_55_to_64 + parameters[SOUR,"TARN_old_Age_65_75"]*Age_65_to_75+
    parameters[SOUR,"TARN_old_Age_over_75"]*Age_over_75 + parameters[SOUR,"TARN_old_female"]*Gender_f+ 
    parameters[SOUR,"TARN_old_female_age_0_5"]*Age_0_to_5*Gender_f + parameters[SOUR,"TARN_old_female_age_6_10"]*Age_6_to_10*Gender_f+ 
    parameters[SOUR,"TARN_old_female_age_11_15"]*Age_11_to_15*Gender_f + parameters[SOUR,"TARN_old_female_age_45_54"]*Age_45_to_54*Gender_f+
    parameters[SOUR,"TARN_old_female_age_55_64"]*Age_55_to_64*Gender_f + parameters[SOUR,"TARN_old_female_age_65_75"]*Age_65_to_75*Gender_f+
    parameters[SOUR,"TARN_old_female_age_75_plus"]*Age_over_75*Gender_f 
  
  #Apply the inverse logit function to the fitted values  
  p_surv <- 1/(1+exp(-FVs))
  
  p_death_temp <- 1 - p_surv
  
  
  # set the probability of death to 0 if the ISS is 0
  p_death <- ifelse(ISS_0==TRUE, 0, p_death_temp)
  
  return(p_death)
  
}

final_dest <- function(strategy, pat_chars, sens, spec, ISS_cutoff_MTC_pos){
  if(strategy=="Manual"){
    
    #set a stream of random numbers, equal to the number of patients
    r <- runif(length(pat_chars[,"ID"]))
    
    #If the patient has an ISS above or equal to the positivity cut-off, assess whether they go to the MTC using 
    pat_chars[,"MTC"] <- ifelse(pat_chars[,"ISS"]>= ISS_cutoff_MTC_pos, ifelse(r[]<sens,1,0), ifelse(r[]<spec,0,1))
  }
}

outcomes <- function(pat_chars, parameters, life_tables, SOUR, strat_name, sensitivity, specificity){
  
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
 
 #estimate the time of death for each patient, as though their ISS is over 15
 #Note this is plus one, because their characteristic is age at baseline and these patients have survived to one year
 Age_at_death_ISS_o15 <-  life_expectancy_ONS(pat_chars[,"Age"]+1,pat_chars[,"Gender"], Life_table_ISS_o_15)
 #estimate the time of death for each patient, as though their ISS is under 16
 Age_at_death_ISS_u16 <-  life_expectancy_ONS(pat_chars[,"Age"]+1,pat_chars[,"Gender"], Life_table_ISS_u_16)
 #Note is age +1 in these calculations, as they must be alive one year after their major trauma to have their long term
 #life expectancy estimated
 
 pat_chars[,"D_1yr_plus"] <- ISS_o15*alive_1yr*Age_at_death_ISS_o15 + ISS_u16*alive_1yr*Age_at_death_ISS_u16
 
 pat_chars[,"Age_death"] <- ifelse(pat_chars[,"D_bl_disch"]==1, pat_chars[,"Age"]+((days_to_discharge*runif(1))/days_in_year), ifelse(pat_chars[,"D_disch_1yr"]==1, pat_chars[,"Age"]+((days_to_discharge+(days_in_year -days_to_discharge)*runif(1))/days_in_year), pat_chars[,"D_1yr_plus"]))
 pat_chars[,"Life_years"] <- pat_chars[,"Age_death"] - pat_chars[,"Age"]
 
 
 return(pat_chars)

}

apply_utils <- function (pat_chars,parameters, SOUR){
  
  #Apply trauma specific utilties
  
  Util_gen_pop_mean_age <- parameters[SOUR,"U_genpop_cons"] + parameters[SOUR,"U_genpop_male"]*mean(pat_chars[,"Gender"]*(91/154)) + 
    parameters[SOUR, "U_genpop_age"]*(61)+parameters[SOUR, "U_genpop_age_squared"]*(61)^2
  #Calculate raw multipliers
  Util_mutl_ISS_o15_MTC <- parameters[SOUR,"U_ISS_o15_MTC"]/Util_gen_pop_mean_age
  Util_mutl_ISS_o15_nMTC <- parameters[SOUR,"U_ISS_o15_nMTC"]/Util_gen_pop_mean_age
  Util_mutl_ISS_u16_o8 <- parameters[SOUR,"U_ISS_u16_o8"]/Util_gen_pop_mean_age
  Util_mult_ISS_u9 <- parameters[SOUR,"Umult_ISS_u9"]
  #constrain the multipliers to be no more than 1 (i.e. you can't be healthier than the general population)
  Util_mutl_ISS_o15_MTC <- ifelse(Util_mutl_ISS_o15_MTC>1,1,Util_mutl_ISS_o15_MTC)
  Util_mutl_ISS_o15_nMTC <- ifelse(Util_mutl_ISS_o15_nMTC>1,1,Util_mutl_ISS_o15_nMTC)
  Util_mutl_ISS_u16 <- ifelse(Util_mutl_ISS_u16_o8>1,1,Util_mutl_ISS_u16_o8)
  
  #create a vector of multpliers relevant to each patient
  #Step 1: create a vector of whether ISS >= 16 or not
  ISS_o15 <- ifelse(pat_chars[,"ISS"]>15,1,0)
  #create a vector of whether ISS >= 9 or not
  ISS_o9 <- ifelse(pat_chars[,"ISS"]>8,1,0)
  #Step 2: create a vector of whether or not the patient has recived major trauma triage or not
  MTC <- ifelse(pat_chars[,"MTC"]==1,1,0)
  #Step 3: determine the appropiate multiplier by the previous vector
  mults <- ifelse(ISS_o15==0,ifelse(ISS_o9==1,Util_mutl_ISS_u16_o8,Util_mult_ISS_u9), ifelse(MTC==0,Util_mutl_ISS_o15_nMTC,Util_mutl_ISS_o15_MTC))

  #Calculate the discounted QALYs
  #The formula for applying continuous discounting to life years is:
  #disc_LY <- exp(-at)
  #where a = ln(1-discount rate), t = time since model entry
  #Ara and brazier formula is:
  # b + c (t+d) + f(t+d)^2 (didn't use e in algebra to ensure I didn't confuse it with exp)
  #Where b = constant + gender coefficent *(1=male, 0=female), c = cofficent for age, t = time since model entry, d = age at model entry, f = coefficent for age ^2
  #therefore, at any given point in time t, a patient's instantaneous discounted utility must be:
  #exp(-at)*(b + c*(t+d) + f*((t+d)^2))
  #Similarly, at any given point in time t, a patient's instantaneous undiscounted utility must be:
  #(b + c*(t+d) + f*((t+d)^2))
  #the indefinate integral (according to wolfram alpha, I was not doing this by hand) for a person's discounted QALYs is:
  #-(exp(-at)*(a^2*b + a^2*c*d+ a^2*c*t + a^2*d^2*f + 2*a^2*d*f*t + a^2*f*t^2 + a*c + 2*a*d*f + 2*a*f*t + 2*f))/a^3 + a constant
  #The indefinate intergral for a person's undiscounted QALYs is:
  #bt+ 0.5*t^2*(c+2*d*f)+d*t*(c+d*f)+(f*t^3)/3 + a constant 
  
  t_start <- 0
  t_end <- ifelse(pat_chars[,"Life_years"]<time_horizon, pat_chars[,"Life_years"], time_horizon)
  a <- log(1+discount_rate_QALYs)
  b <- as.numeric(parameters[SOUR,"U_genpop_cons"] + parameters[SOUR,"U_genpop_male"]*pat_chars[,"Gender"])
  c <- as.numeric(parameters[SOUR,"U_genpop_age"])
  d <- pat_chars[,"Age"]
  f <- as.numeric(parameters[SOUR, "U_genpop_age_squared"])
  
  #Calculate the indefinate integral at the time the patient dies (without the constant)
  t <- t_end
  #Discounted QALYs
  temp_1 <- -(exp(-a*t)*(a^2*b + a^2*c*d + a^2*c*t + a^2*d^2*f + 2*a^2*d*f*t + a^2*f*t^2 + a*c + 2*a*d*f + 2*a*f*t + 2*f))/a^3
  #Undiscounted QALYs
  temp_3 <- b*t + 0.5*(t^2)*(c+2*d*f)+d*t*(c+d*f)+(f*t^3)/3
  #Calculate the indeinate integral at the time the patient enters the model (again without the constant)
  t <- t_start
  #Discounted QALYs
  temp_2 <- -(exp(-a*t)*(a^2*b + a^2*c*d + a^2*c*t + a^2*d^2*f + 2*a^2*d*f*t + a^2*f*t^2 + a*c + 2*a*d*f + 2*a*f*t + 2*f))/a^3
  #Undiscounted QALYs
  temp_4 <- b*t + 0.5*(t^2)*(c+2*d*f)+d*t*(c+d*f)+(f*t^3)/3
  
  #Calculate the definate integral between the time that the patient dies and when they entered the model
  disc_QALYs <- (temp_1 - temp_2)
  undisc_QALYs <- (temp_3 - temp_4)
  
  #Apply multipliers for health status
  disc_QALYs <- disc_QALYs*mults
  undisc_QALYs <- undisc_QALYs*mults

  #Record the discounted QALYs in the pat_chars_matrix
  pat_chars[,"dQALYS"] <- disc_QALYs
  pat_chars[,"QALYS"] <- undisc_QALYs 
  
  return(pat_chars)
}

apply_costs <- function(pat_chars, parameters, SOUR){
  #create a vector of costs of patient's admissions
  #Step 1: create subgroup variables by ISS, trauma type, and whether transfered to MTC
  blunt_trauma <- ifelse(pat_chars[,"Blunt_trauma"]==1,1,0)
  pen_trauma <- ifelse(pat_chars[,"Blunt_trauma"]==0,1,0)
  disch <- ifelse(pat_chars[,"D_bl_disch"]==0,1,0)
  
  pent_ISS_u9 <- pen_trauma*ifelse(pat_chars[,"ISS"]<9,1,0)
  pent_ISS_o8_u16 <- pen_trauma*ifelse(pat_chars[,"ISS"]>8&pat_chars[,"ISS"]<16,1,0)
  pent_ISS_o15_u25 <- pen_trauma*ifelse(pat_chars[,"ISS"]>15&pat_chars[,"ISS"]<25,1,0)
  pent_ISS_o24_u35 <- pen_trauma*ifelse(pat_chars[,"ISS"]>24&pat_chars[,"ISS"]<35,1,0)
  pent_ISS_o34 <- pen_trauma*ifelse(pat_chars[,"ISS"]>34,1,0)
  bluntt_ISS_u10 <- blunt_trauma*ifelse(pat_chars[,"ISS"]<10,1,0)
  bluntt_ISS_u17_o9 <-blunt_trauma*ifelse(pat_chars[,"ISS"]<17&pat_chars[,"ISS"]>9,1,0)
  bluntt_ISS_u26_o16 <-blunt_trauma*ifelse(pat_chars[,"ISS"]<26&pat_chars[,"ISS"]>16,1,0)
  bluntt_ISS_o25 <-blunt_trauma*ifelse(pat_chars[,"ISS"]>25,1,0)
  
  MTC <- ifelse(pat_chars[,"MTC"]==1,1,0)
  MTC_ISS_o8_u16 <- MTC*ifelse(pat_chars[,"ISS"]<16&pat_chars[,"ISS"]>8,1,0)
  MTC_ISS_o15<- MTC*ifelse(pat_chars[,"ISS"]>15,1,0)
  
  MTC_transfer <- ifelse(pat_chars[,"MTC_transfer"]==1,1,0)
  
  #Apply the costs associated with the admission. These include the MTC tariff (if the patient went to an MTC) & treatment costs
  admission_cost <- MTC_ISS_o8_u16*parameters[SOUR,"C_MTC_ISS_o8_u16"]+MTC_ISS_o15*parameters[SOUR,"C_MTC_ISS_o15"]+pent_ISS_u9*parameters[SOUR,"C_pent_ISS_O0_U10"]+
    pent_ISS_o8_u16*parameters[SOUR,"C_pent_ISS_O9_U16"]+pent_ISS_o15_u25*parameters[SOUR,"C_pent_ISS_O15_U25"]+pent_ISS_o24_u35*parameters[SOUR,"C_pent_ISS_O24_U34"]+
    pent_ISS_o34*parameters[SOUR,"C_pent_ISS_O34"]+ bluntt_ISS_u10*parameters[SOUR,"C_bluntt_ISS_U10"]+bluntt_ISS_u17_o9*parameters[SOUR,"C_bluntt_ISS_U17_O_9"]+
    bluntt_ISS_u26_o16*parameters[SOUR,"C_bluntt_ISS_U26_O16"]+bluntt_ISS_o25*parameters[SOUR,"C_bluntt_ISS_O25"]+MTC_transfer*parameters[SOUR,"C_additional_ambulance"]
  
  #Apply the costs associated with being discharged, up until 6 months post-baseline and half of their age-gender matched costs 
  disch_6m_cost <- parameters[SOUR,"C_disch_6m"]*disch
  
  pat_chars[,"Costs"] <- admission_cost + disch_6m_cost 
  
  Frac_future_health_care_costs <- ifelse(pat_chars[,"ISS"]>15, parameters[SOUR,"Increase_lifetime_cost_ISS_o15"],parameters[SOUR,"Increase_lifetime_cost_ISS_u15"])
  
  #add on the costs associated with future health care costs by patient
  for(y in 1:length(pat_chars[,"Life_years"])){
    pat_chars[y,"Costs"] <- pat_chars[y,"Costs"] + ifelse(pat_chars[y,"Life_years"]>0.5,1,0)*Frac_future_health_care_costs[y]*future_costs[ifelse((pat_chars[y,"Age"]+1)<86,pat_chars[y,"Age"]+1,86),ifelse(pat_chars[y,"Gender"]==1,2,3)]
  }
  
  
  
  #Store some temporary variables for discounting
  t0 <- (((days_in_year/2)/days_in_year)-(days_to_discharge/days_in_year))
  t1 <- cont_disc(1,days_to_discharge/days_in_year, ifelse(pat_chars[,"Life_years"]>0.5,0.5, pat_chars[,"Life_years"]), discount_rate_costs)
  t2 <- cont_disc(1,0.5, ifelse(pat_chars[,"Life_years"]<1,pat_chars[,"Life_years"],1), discount_rate_costs)
  #adjust the discounting for the fact that these are not yearly costs, but are instead a cost between discharge
  #and 6 months post-injury
  pat_chars[,"DCosts"] <- admission_cost + disch_6m_cost*(t1/t0) 
  
  #add on the costs associated with future health care costs by patient
  for(y in 1:length(pat_chars[,"Life_years"])){
    pat_chars[y,"DCosts"] <- pat_chars[y,"DCosts"] + ifelse(pat_chars[y,"Life_years"]>0.5,1,0)*Frac_future_health_care_costs[y]*future_costs[ifelse((pat_chars[y,"Age"]+1)<86,pat_chars[y,"Age"]+1,86),ifelse(pat_chars[y,"Gender"]==1,2,3)]*t2[y]
  }
  
  
  
  #do a loop to work out the costs
  for(y in 1:length(pat_chars[,"Life_years"])){
  if(pat_chars[y,"Life_years"]>=0.5){
    for(i in 0:ceiling(pat_chars[y,"Life_years"])){
      #in the very first cycle add on costs for the remaining 6 months that the patient lives
      if(i==0){
        #Work out the time to apply the increase in costs attributable to the trauma incident for. This will be whichever is the least of 1 year or 
        #years lived
        t0 <- min(1,pat_chars[y,"Life_years"])
        #create some temporary varaibles indicating which row (t1) and column (t2) need to be looked up in the future costs spreadsheet
        t1 <- as.numeric(ifelse(pat_chars[y,"Age"]+2+i<86,pat_chars[y,"Age"]+2+i,86))
        t2 <- as.numeric(ifelse(pat_chars[y,"Gender"]==1,2,3))
        
        #If this is the first year, apply the costs for a proportion of the year the individual is alive or for 0.5, whicheve is the smallest
        pat_chars[y,"Costs"] <- pat_chars[y,"Costs"] + (t0 - 0.5)*(Frac_future_health_care_costs[y]-1)*future_costs[t1,t2]
        pat_chars[y,"DCosts"] <- pat_chars[y,"DCosts"] + cont_disc(1,0.5,t0,discount_rate_costs)*(Frac_future_health_care_costs[y]-1)*future_costs[t1,t2]
        
      }
      else if(i==floor(pat_chars[y,"Life_years"])){
        #create some temporary varaibles indicating which row (t1) and column (t2) need to be looked up in the furture costs spreadsheet
        t1 <- as.numeric(ifelse(pat_chars[y,"Age"]+2+i<86,pat_chars[y,"Age"]+2+i,86))
        t2 <- as.numeric(ifelse(pat_chars[y,"Gender"]==1,2,3))
        #If this is in the year that the patient dies only add on a the appropiate proportion of their health care costs
        pat_chars[y,"Costs"] <- pat_chars[y,"Costs"] + (pat_chars[y,"Life_years"]-floor(pat_chars[y,"Life_years"]))*(Frac_future_health_care_costs[y]-1)*future_costs[t1,t2]
        pat_chars[y,"DCosts"] <- pat_chars[y,"DCosts"] + (Frac_future_health_care_costs[y]-1)*future_costs[t1,t2]*cont_disc(1,i,pat_chars[y,"Life_years"],discount_rate_costs)
      }
      else if (i ==ceiling(pat_chars[y,"Life_years"])){
        # if i is greater than the number of life years accrued, end this inner loop
        break
      }
      
      else{
        #otherwise add up the costs within each year
        #create some temporary varaibles indicating which row (t1) and column (t2) need to be looked up in the furture costs spreadsheet
        t1 <- as.numeric(ifelse(pat_chars[y,"Age"]+2+i<86,pat_chars[y,"Age"]+2+i,86))
        t2 <- as.numeric(ifelse(pat_chars[y,"Gender"]==1,2,3))
        #If this is in the year that the patient dies not die add on the full increase in cost for the year
        pat_chars[y,"Costs"] <- pat_chars[y,"Costs"]+(Frac_future_health_care_costs[y]-1)*future_costs[t1,t2]
        pat_chars[y,"DCosts"] <- pat_chars[y,"DCosts"] + (Frac_future_health_care_costs[y]-1)*future_costs[t1,t2]*cont_disc(1,i,i+1,discount_rate_costs)
      }
      
    }
  }
    
    }
  return(pat_chars)
}

### Function to apply triage strategies in the model
triage_strategies <- function(pat_chars, name, sens, spec){
  #Apply the manual strategy, where the sensitivity and specficity of the rule are user
  #defined
  if(name=="manual"){
    major_trauma <- pat_chars[,"ISS"] > 15
    non_mt <- pat_chars[,"ISS"] < 16
    rands <- runif(length(pat_chars[,"ISS"]))
    sens_spec <- ifelse(major_trauma==TRUE, sens, 1-spec)
    temp <- ifelse(rands[]<sens_spec, 1,0)
    pat_chars[,"Triage_rule"] <- temp
    #For the manual rules are based on final destination as the outcome. Therefore, I do
    #not account for compliance
    pat_chars[,"MTC"] <- pat_chars[,"Triage_rule"]
  }
  #further strategies to be added at a later date
  return(pat_chars)
}


model_single_run <- function(pat_chars, parameters, SOUR, life_tables, strat_name, sensitivity, specificity, pop_report){
  
  #add in line of code to generate parameters here
  #estimate the clinical outcomes
  pat_chars <- outcomes(pat_chars, parameters, life_tables, SOUR, strat_name, sensitivity, specificity)
  #apply the utilities
  pat_chars <- apply_utils(pat_chars,parameters, SOUR)
  #apply the costs
  pat_chars <- apply_costs(pat_chars, parameters, SOUR)
  
  #store the key summary statistics in a results table
  #Create a matrix with 1 row
  res_tab <- matrix(nrow = 1, ncol = 12)
 
  #record the sensitivity of the DR
  #First temp vector is a vector indicating that the patient recieved MTC care and their ISS > 15
  t1 <- ifelse(pat_chars[,"MTC"]==1&pat_chars[,"ISS"]>15,1,0)
  #Second temp vector is a vector indicating that the patient had an ISS > 15
  t2 <- ifelse(pat_chars[,"ISS"]>15,1,0)
  #Record the sensitivity of the DR
  res_tab[1,1] <- sum(t1)/sum(t2)
  #Record the specificity of the DR
  #First temp vector is a vector indicating that the patient went to a non-MTC and their ISS < 16
  t1 <- ifelse(pat_chars[,"MTC"]==0&pat_chars[,"ISS"]<16,1,0)
  #Second temp vector is a vector indicating that the patient had an ISS > 15
  t2 <- ifelse(pat_chars[,"ISS"]<16,1,0)
  #Record the sensitivity of the DR
  res_tab[1,2] <- sum(t1)/sum(t2)
  #record the mean number of patients who recieve MTC care
  res_tab[1,3] <- mean(pat_chars[,"MTC"])
  #record the probability of dying between baseline and discharge
  res_tab[1,4] <- mean(pat_chars[,"D_bl_disch"])
  #record the probability of dying between discharge and one year post-accident
  res_tab[1,5] <- mean(pat_chars[,"D_disch_1yr"])
  #record the mean life years lived
  res_tab[1,6] <- mean(pat_chars[,"Life_years"])
  #record the mean undiscounted QALYs
  res_tab[1,7]<- mean(pat_chars[,"QALYS"])
  #record the mean discounted QALYs
  res_tab[1,8] <- mean(pat_chars[,"dQALYS"])
  #record the mean undiscounted costs
  res_tab[1,9] <- mean(pat_chars[,"Costs"])
  #record the mean discounted costs
  res_tab[1,10] <- mean(pat_chars[,"DCosts"])
  #record the percentage of patients with an ISS of 16 or more
  #First temp vector is a vector indicating that the patient had an ISS > 15
  t1 <- ifelse(pat_chars[,"ISS"]>15,1,0)
  #divide the sum of this vector by the number of patients in the population
  t2 <- sum(t1)/length(pat_chars[,"ISS"])
  #record the result
  res_tab[1,11] <- t2
  #record the percentage of patients with an ISS of 9 or more, but less than 16
  t1 <- ifelse(pat_chars[,"ISS"]>8&pat_chars[,"ISS"]<16,1,0)
  #divide the sum of this vector by the number of patients in the population
  t2 <- sum(t1)/length(pat_chars[,"ISS"])
  #record the result
  res_tab[1,12] <- t2
  
  if(pop_report==0){
    return(pat_chars)
  }else{
    return(res_tab)
  }
}


run_simulation <- function(param_inputs, PSA_switch, PSA_numb, pat_numb, strat_name, sensitivity, specificity, pop_report){
  
  #set the random number seed
  set.seed(26090100)
  #Generate pat chars to be 
  pat_chars <- gen_pat_chars(pat_numb, means, covariance, age_tab, gen_tab, ISS_tab, GCS_tab)
  #generate the parameters
  parameters <- gen_parameters(PSA_switch,PSA_numb, param_inputs)
  #As the utility parameters for people with an ISS > 9 are all the same set all the utility samples to be the same 
  parameters[,"U_ISS_o15_nMTC"] <- parameters[,"U_ISS_o15_MTC"]
  parameters[,"U_ISS_u16_o8"] <- parameters[,"U_ISS_o15_MTC"]
  
  #export the parameters, if required (bug checking / SAVI)
  if(Param_export==1){
    write.csv(parameters, file = "parameter_outputs.csv")
  }

  #setup the matrix to store results
  results <- matrix (nrow = ifelse(PSA_switch == 1,PSA_numb,1), ncol=12)
  #create names for the results matrix
  colnames(results) <- c("Sens_DR","Spec_DR", "Number_recieving_MTC_care","proportion_died_before_discharge","proportion_died_between_discharge_and_1_year", "Years_lived",
                              "undiscounted_QALYs", "discounted_QALYs", "undiscounted_Costs", "discounted_Costs", "proportion_ISS_over_16", "proportion_ISS_over_8_under_16")
  
  #set SOUR to 1 to start the simulation
  SOUR <- 1
  #define the strategy
  
  #run the simulation, calling the user defined function to run the model once
  if(pop_report==0){
    test <- model_single_run(pat_chars, parameters, SOUR, life_tabs,strat_name, sensitivity, specificity, pop_report)
    return(test)
  }else if(PSA_switch==0){
    results[SOUR,] <- model_single_run(pat_chars, parameters, SOUR, life_tabs,strat_name, sensitivity, specificity, pop_report)
  }else{
    for(SOUR in 1:PSA_numb){
      results[SOUR,]<- model_single_run(pat_chars, parameters, SOUR, life_tabs,strat_name, sensitivity, specificity, pop_report)
      #make some text appear indicating the current PSA run
      print(SOUR)
      #update the screen so you can see the model hasn't crashed during the PSA
      flush.console()
    }
    return (results)
  } 

}


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

##########################################################
##Check number of patients
if(test_pat_chars=="Yes"){
sens_100_spec_3 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.998, 0.025,0)
sens_95_spec_19 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.948, 0.187,0)
sens_90_spec_58 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.904, 0.584,0)
sens_88_spec_63 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.875, 0.628,0)
sens_75_spec_66 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.746, 0.657,0)
sens_70_spec_70 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.698, 0.701,0)
sens_64_spec_76 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.642, 0.761,0)
sens_57_spec_80 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.57, 0.8,0)
sens_28_spec_89 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.284, 0.886,0)

write.csv(sens_100_spec_3, paste(file_location,"Patient Characteristic Check\\sens_100_spec_3.csv", sep=""))
write.csv(sens_95_spec_19, paste(file_location,"Patient Characteristic Check\\sens_95_spec_19.csv", sep=""))
write.csv(sens_90_spec_58, paste(file_location,"Patient Characteristic Check\\sens_90_spec_58.csv", sep=""))
write.csv(sens_88_spec_63, paste(file_location,"Patient Characteristic Check\\sens_88_spec_63.csv", sep=""))
write.csv(sens_75_spec_66, paste(file_location,"Patient Characteristic Check\\sens_75_spec_66.csv", sep=""))
write.csv(sens_70_spec_70, paste(file_location,"Patient Characteristic Check\\sens_70_spec_70.csv", sep=""))
write.csv(sens_64_spec_76, paste(file_location,"Patient Characteristic Check\\sens_64_spec_76.csv", sep=""))
write.csv(sens_57_spec_80, paste(file_location,"Patient Characteristic Check\\sens_57_spec_80.csv", sep=""))
write.csv(sens_28_spec_89, paste(file_location,"Patient Characteristic Check\\sens_28_spec_89.csv", sep=""))


#run them normally as summaries too for error checking
sens_100_spec_3 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.998, 0.025,1)
sens_95_spec_19 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.948, 0.187,1)
sens_90_spec_58 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.904, 0.584,1)
sens_88_spec_63 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.875, 0.628,1)
sens_75_spec_66 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.746, 0.657,1)
sens_70_spec_70 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.698, 0.701,1)
sens_64_spec_76 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.642, 0.761,1)
sens_57_spec_80 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.57, 0.8,1)
sens_28_spec_89 <- run_simulation(param_data_bc, 0, 1, pat_numb, "manual", 0.284, 0.886,1)

write.csv(sens_100_spec_3, paste(file_location,"Patient Characteristic Check\\sens_100_spec_3_ec.csv", sep=""))
write.csv(sens_95_spec_19, paste(file_location,"Patient Characteristic Check\\sens_95_spec_19_ec.csv", sep=""))
write.csv(sens_90_spec_58, paste(file_location,"Patient Characteristic Check\\sens_90_spec_58_ec.csv", sep=""))
write.csv(sens_88_spec_63, paste(file_location,"Patient Characteristic Check\\sens_88_spec_63_ec.csv", sep=""))
write.csv(sens_75_spec_66, paste(file_location,"Patient Characteristic Check\\sens_75_spec_66_ec.csv", sep=""))
write.csv(sens_70_spec_70, paste(file_location,"Patient Characteristic Check\\sens_70_spec_70_ec.csv", sep=""))
write.csv(sens_64_spec_76, paste(file_location,"Patient Characteristic Check\\sens_64_spec_76_ec.csv", sep=""))
write.csv(sens_57_spec_80, paste(file_location,"Patient Characteristic Check\\sens_57_spec_80_ec.csv", sep=""))
write.csv(sens_28_spec_89, paste(file_location,"Patient Characteristic Check\\sens_28_spec_89_ec.csv", sep=""))

}else{ #if we aren't testing the stability of the results wrt to the number of patients run normal analyses
  
  
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
  write.csv(sens_100_spec_3_PSA, paste(file_location,"PSA results\\sens_100_spec_3_PSA.csv", sep=""))
  use_params_sens_100_spec_3_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_100_spec_3_PSA, "PSA results\\sens_100_spec_3_PSA_params.csv")
  }
  if(PSA_strat == "S95"){
  sens_95_spec_19_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.948, 0.187,1)
  write.csv(sens_95_spec_19_PSA, paste(file_location,"PSA results\\sens_95_spec_19_PSA.csv", sep=""))
  use_params_sens_95_spec_19_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_95_spec_19_PSA, "PSA results\\sens_95_spec_19_PSA_params.csv")
  }
  if(PSA_strat == "S90"){
  sens_90_spec_58_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.904, 0.584,1)
  write.csv(sens_90_spec_58_PSA, paste(file_location,"PSA results\\sens_90_spec_58_PSA.csv", sep=""))
  use_params_sens_90_spec_58_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_90_spec_58_PSA, "PSA results\\sens_90_spec_58_PSA_params.csv")
  }
  if(PSA_strat == "S88"){
  sens_88_spec_63_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.875, 0.628,1)
  write.csv(sens_88_spec_63_PSA, paste(file_location,"PSA results\\sens_88_spec_63_PSA.csv", sep=""))
  use_params_sens_88_spec_63_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_88_spec_63_PSA, "PSA results\\sens_88_spec_63_PSA_params.csv")
  }
  if(PSA_strat == "S75"){
  sens_75_spec_66_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.746, 0.657,1)
  write.csv(sens_75_spec_66_PSA, paste(file_location,"PSA results\\sens_75_spec_66.csv", sep=""))
  use_params_sens_75_spec_66_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_75_spec_66_PSA, "sens_75_spec_66_PSA_params.csv")
  }
  if(PSA_strat == "S70"){
  sens_70_spec_70_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.698, 0.701,1)
  write.csv(sens_70_spec_70_PSA, paste(file_location,"PSA results\\sens_70_spec_70.csv", sep=""))
  use_params_sens_70_spec_70_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_70_spec_70_PSA, "sens_70_spec_70_PSA_params.csv")
  }
  if(PSA_strat == "S64"){
  sens_64_spec_76_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.642, 0.761,1)
  write.csv(sens_64_spec_76_PSA, paste(file_location,"PSA results\\sens_64_spec_76.csv", sep=""))
  use_params_sens_64_spec_76 <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_64_spec_76, "PSA results\\sens_90_spec_58_PSA_params.csv")
  }
  if(PSA_strat == "S57"){
  sens_57_spec_80_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.57, 0.8,1)
  write.csv(sens_57_spec_80_PSA, paste(file_location,"PSA results\\sens_57_spec_80.csv", sep=""))
  use_params_sens_57_spec_80_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_57_spec_80_PSA, "PSA results\\sens_90_spec_58_PSA_params.csv")
  }
  if(PSA_strat == "S28"){
  sens_28_spec_89_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0.284, 0.886,1)
  write.csv(sens_28_spec_89_PSA, paste(file_location,"PSA results\\sens_28_spec_89.csv", sep=""))
  use_params_sens_28_spec_89_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_28_spec_89_PSA, "PSA results\\sens_28_spec_89_PSA_params.csv")
  }
}
#Use the newer TARN mortality equation
TARN_mort_eq <- "New" 
MTCs_in_mort_risk <- "Yes"
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

write.csv(det_analyses,"new TARN equations.csv")
}

#Do a scenario analysis using the old TARN equations, where some patients get MTC equivalent care
param_data_SA2 <- param_data
Proportion_RR_MTC_ISS_o8_u16 <- 0
TARN_mort_eq <- "Old" 
MTCs_in_mort_risk <- "Yes" 
population_source <- "Dutch" 
percent_TARN_cases_reported_ISS_o16 <- 1
percent_TARN_cases_reported_ISS_o9_u16 <- 1

if(PSA_switch==0){
  sens_100_spec_3 <- run_simulation(param_data_SA2, 0, 1, pat_numb, "manual", 0.998, 0.025,1)
  sens_95_spec_19 <- run_simulation(param_data_SA2, 0, 1, pat_numb, "manual", 0.948, 0.187,1)
  sens_90_spec_58 <- run_simulation(param_data_SA2, 0, 1, pat_numb, "manual", 0.904, 0.584,1)
  sens_88_spec_63 <- run_simulation(param_data_SA2, 0, 1, pat_numb, "manual", 0.875, 0.628,1)
  sens_75_spec_66 <- run_simulation(param_data_SA2, 0, 1, pat_numb, "manual", 0.746, 0.657,1)
  sens_70_spec_70 <- run_simulation(param_data_SA2, 0, 1, pat_numb, "manual", 0.698, 0.701,1)
  sens_64_spec_76 <- run_simulation(param_data_SA2, 0, 1, pat_numb, "manual", 0.642, 0.761,1)
  sens_57_spec_80 <- run_simulation(param_data_SA2, 0, 1, pat_numb, "manual", 0.57, 0.8,1)
  sens_28_spec_89 <- run_simulation(param_data_SA2, 0, 1, pat_numb, "manual", 0.284, 0.886,1)
  
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
  
  write.csv(det_analyses,"MTC equivalent care old TARN.csv")
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
  sens_100_spec_10_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 1, 0.1,1)
  write.csv(sens_100_spec_10_PSA, paste(file_location,"PSA results\\sens_100_spec_10_PSA.csv", sep=""))
  use_params_sens_100_spec_10_PSA <- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_100_spec_10_PSA, "PSA results\\sens_100_spec_10_PSA_params.csv")
  }
  if(PSA_strat == "nMTC"){
  sens_0_spec_90_PSA <- run_simulation(param_data_bc, 1, PSA_numb, pat_numb, "manual", 0, 0.9,1)
  write.csv(sens_0_spec_90_PSA, paste(file_location,"PSA results\\sens_0_spec_90_PSA.csv", sep=""))
  use_params_sens_0_spec_90_PSA<- read.csv("parameter_outputs.csv")
  write.csv(use_params_sens_0_spec_90_PSA, "PSA results\\sens_0_spec_90_PSA_params.csv")
}
}
}