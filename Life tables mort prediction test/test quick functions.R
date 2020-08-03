##### call in data that I'll need
pat_chars <- read.csv("D:\\Work from Home\\Local Git\\Life tables mort prediction test\\pat chars.csv")

age <- pat_chars[,"Age"]
gender <- pat_chars[,"Gender"]

life_tables <- read.csv("D:\\Work from Home\\Local Git\\Life tables mort prediction test\\ONSlifetables.csv")

parameters <- read.csv("D:\\Work from Home\\Local Git\\Life tables mort prediction test\\parameters.csv")

##### working , but slow, function
#setting the seed to see if any errors have been introduced
set.seed(1)

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

start_time <- Sys.time()
slow_results <- life_expectancy_ONS(age,gender, life_tables)
end_time <- Sys.time()

elapsed_time <- end_time - start_time
elapsed_time

#####alternatives

set.seed(1)

life_expectancy_ONS_new <- function(age, gender, life_tables){
  #Check that the length of the two vectors (age and gender) is the same
  
  #Generate a temporary vector of the same length to store the results, with all values defaulting to -99.
  #a value of -99 indicates that the age of death has not been determined in the loop
  all_cause_death <- rep(-99, length.out = length(age))
  
  #Loop so that the process is repeated for each patient in the model on the basis of their age
 
    
    #create a vector of 16 to 100
    ages <- seq(from=16, to =101, by =1)

    #create a matrix by gender with the probability of death, and life expectancy by age
    Life_tab_m <- matrix(data = NA, nrow = 86, ncol = 4)
    Life_tab_m[1:86,1] <- ages
    Life_tab_m[1:85,2] <- life_tables[17:101,2]
    
    #set the table so that people aged 101, live the life expectancy of someone aged 100
    Life_tab_m[86,2] <- 1
    #Overwrite 101, with the average life expectancy of the men who survive to age 100
    Life_tab_m[86,1] <- 102.02
    Life_tab_m[1,3] <- 1000
    
    #record the same table for females
    Life_tab_f <- Life_tab_m
    #overwrite the probs of death, with female equivalent probabilities
    Life_tab_f[1:85,2] <- life_tables[17:101,3]
    #change the life expectancy of people who survive to 101 with female data
    Life_tab_f[86,1] <- 102.32
    
    #estimate the risk of death by age for each future age
    for(i in 2:86){
    Life_tab_m[i,3] <- Life_tab_m[i-1,3] - Life_tab_m[i-1,3]*Life_tab_m[i-1,2]
    Life_tab_m[i-1,4] <- (Life_tab_m[i-1,3] - Life_tab_m[i,3])/Life_tab_m[1,3]
    Life_tab_f[i,3] <- Life_tab_f[i-1,3] - Life_tab_f[i-1,3]*Life_tab_f[i-1,2]
    Life_tab_f[i-1,4] <- (Life_tab_f[i-1,3] - Life_tab_f[i,3])/Life_tab_f[1,3]
    
    }
    
    Life_tab_m[86,4] <- Life_tab_m[86,3]/Life_tab_m[1,3]
    Life_tab_f[86,4] <- Life_tab_f[86,3]/Life_tab_f[1,3]
    
    #generate cumulative probability of death column by ainitial starting age
    #start at the starting age of 16
      
      Prob_death_m_16 <-  matrix(data=NA, nrow = (101 -  16+1), ncol = 2)
      Prob_death_m_16[,2] <- Life_tab_m[,1]
      Prob_death_m_16[,1] <- Life_tab_m[,4]
      
      Prob_death_f_16 <-  matrix(data=NA, nrow = (101 -  16+1), ncol = 2)
      Prob_death_f_16[,2] <- Life_tab_f[,1]
      Prob_death_f_16[,1] <- Life_tab_f[,4]
      
      #17
      Prob_death_m_17 <- Prob_death_m_16[(17-15):86,]
      Prob_death_m_17[,1] <- Prob_death_m_17[,1] / sum(Prob_death_m_17[,1])
      Prob_death_m_17[,1] <- cumsum(Prob_death_m_17[,1])
      
      Prob_death_f_17 <- Prob_death_f_16[(17-15):86,]
      Prob_death_f_17[,1] <- Prob_death_f_17[,1] / sum(Prob_death_f_17[,1])
      Prob_death_f_17[,1] <- cumsum(Prob_death_f_17[,1])
      
      #18
      Prob_death_m_18 <- Prob_death_m_16[(18-15):86,]
      Prob_death_m_18[,1] <- Prob_death_m_18[,1] / sum(Prob_death_m_18[,1])
      Prob_death_m_18[,1] <- cumsum(Prob_death_m_18[,1])
      
      Prob_death_f_18 <- Prob_death_f_16[(18-15):86,]
      Prob_death_f_18[,1] <- Prob_death_f_17[,1] / sum(Prob_death_f_18[,1])
      Prob_death_f_18[,1] <- cumsum(Prob_death_f_18[,1])
      
      #19
      Prob_death_m_19 <- Prob_death_m_16[(19-15):86,]
      Prob_death_m_19[,1] <- Prob_death_m_19[,1] / sum(Prob_death_m_19[,1])
      Prob_death_m_19[,1] <- cumsum(Prob_death_m_19[,1])
      
      Prob_death_f_19 <- Prob_death_f_16[(19-15):86,]
      Prob_death_f_19[,1] <- Prob_death_f_19[,1] / sum(Prob_death_f_19[,1])
      Prob_death_f_19[,1] <- cumsum(Prob_death_f_19[,1])
      
      #20
      Prob_death_m_20 <- Prob_death_m_16[(20-15):86,]
      Prob_death_m_20[,1] <- Prob_death_m_20[,1] / sum(Prob_death_m_20[,1])
      Prob_death_m_20[,1] <- cumsum(Prob_death_m_20[,1])
      
      
      Prob_death_f_20 <- Prob_death_f_16[(20-15):86,]
      Prob_death_f_20[,1] <- Prob_death_f_20[,1] / sum(Prob_death_f_20[,1])
      Prob_death_f_20[,1] <- cumsum(Prob_death_f_20[,1])
      
      #21
      Prob_death_m_21 <- Prob_death_m_16[(21-15):86,]
      Prob_death_m_21[,1] <- Prob_death_m_21[,1] / sum(Prob_death_m_21[,1])
      Prob_death_m_21[,1] <- cumsum(Prob_death_m_21[,1])
      
      
      Prob_death_f_21 <- Prob_death_f_16[(21-15):86,]
      Prob_death_f_21[,1] <- Prob_death_f_21[,1] / sum(Prob_death_f_21[,1])
      Prob_death_f_21[,1] <- cumsum(Prob_death_f_21[,1])
      
      #22
      Prob_death_m_22 <- Prob_death_m_16[(22-15):86,]
      Prob_death_m_22[,1] <- Prob_death_m_22[,1] / sum(Prob_death_m_22[,1])
      Prob_death_m_22[,1] <- cumsum(Prob_death_m_22[,1])
      
      
      Prob_death_f_22 <- Prob_death_f_16[(22-15):86,]
      Prob_death_f_22[,1] <- Prob_death_f_22[,1] / sum(Prob_death_f_22[,1])
      Prob_death_f_22[,1] <- cumsum(Prob_death_f_22[,1])
      
      #23
      Prob_death_m_23 <- Prob_death_m_16[(23-15):86,]
      Prob_death_m_23[,1] <- Prob_death_m_23[,1] / sum(Prob_death_m_23[,1])
      Prob_death_m_23[,1] <- cumsum(Prob_death_m_23[,1])
      
      
      Prob_death_f_23 <- Prob_death_f_16[(23-15):86,]
      Prob_death_f_23[,1] <- Prob_death_f_23[,1] / sum(Prob_death_f_23[,1])
      Prob_death_f_23[,1] <- cumsum(Prob_death_f_23[,1])
      
      #24
      Prob_death_m_24 <- Prob_death_m_16[(24-15):86,]
      Prob_death_m_24[,1] <- Prob_death_m_24[,1] / sum(Prob_death_m_24[,1])
      Prob_death_m_24[,1] <- cumsum(Prob_death_m_24[,1])
      
      
      Prob_death_f_24 <- Prob_death_f_16[(24-15):86,]
      Prob_death_f_24[,1] <- Prob_death_f_24[,1] / sum(Prob_death_f_24[,1])
      Prob_death_m_24[,1] <- cumsum(Prob_death_m_24[,1])
      
      
      #25
      Prob_death_m_25 <- Prob_death_m_16[(25-15):86,]
      Prob_death_m_25[,1] <- Prob_death_m_25[,1] / sum(Prob_death_m_25[,1])
      Prob_death_m_25[,1] <- cumsum(Prob_death_m_25[,1])
      
      
      Prob_death_f_25 <- Prob_death_f_16[(25-15):86,]
      Prob_death_f_25[,1] <- Prob_death_f_25[,1] / sum(Prob_death_f_25[,1])      
      Prob_death_f_25[,1] <- cumsum(Prob_death_f_25[,1])

      
      #26
      Prob_death_m_26 <- Prob_death_m_16[(26-15):86,]
      Prob_death_m_26[,1] <- Prob_death_m_26[,1] / sum(Prob_death_m_26[,1])
      Prob_death_m_26[,1] <- cumsum(Prob_death_m_26[,1])
      
      
      Prob_death_f_26 <- Prob_death_f_16[(26-15):86,]
      Prob_death_f_26[,1] <- Prob_death_f_26[,1] / sum(Prob_death_f_26[,1])
      Prob_death_f_26[,1] <- cumsum(Prob_death_f_26[,1])
      
      
      #27
      Prob_death_m_27 <- Prob_death_m_16[(27-15):86,]
      Prob_death_m_27[,1] <- Prob_death_m_27[,1] / sum(Prob_death_m_27[,1])
      Prob_death_m_27[,1] <- cumsum(Prob_death_m_27[,1])
      
      
      Prob_death_f_27 <- Prob_death_f_16[(27-15):86,]
      Prob_death_f_27[,1] <- Prob_death_f_27[,1] / sum(Prob_death_f_27[,1])
      Prob_death_f_27[,1] <- cumsum(Prob_death_f_27[,1])
      
      
      #28
      Prob_death_m_28 <- Prob_death_m_16[(28-15):86,]
      Prob_death_m_28[,1] <- Prob_death_m_28[,1] / sum(Prob_death_m_28[,1])
      Prob_death_m_28[,1] <- cumsum(Prob_death_m_28[,1])
      
      
      Prob_death_f_28 <- Prob_death_f_16[(28-15):86,]
      Prob_death_f_28[,1] <- Prob_death_f_28[,1] / sum(Prob_death_f_28[,1])
      Prob_death_f_28[,1] <- cumsum(Prob_death_f_28[,1])
      
      
      #29
      Prob_death_m_29 <- Prob_death_m_16[(29-15):86,]
      Prob_death_m_29[,1] <- Prob_death_m_29[,1] / sum(Prob_death_m_29[,1])
      Prob_death_m_29[,1] <- cumsum(Prob_death_m_29[,1])
      
      
      Prob_death_f_29 <- Prob_death_f_16[(29-15):86,]
      Prob_death_f_29[,1] <- Prob_death_f_29[,1] / sum(Prob_death_f_29[,1])
      Prob_death_f_29[,1] <- cumsum(Prob_death_f_29[,1])
      
      
      #30
      Prob_death_m_30 <- Prob_death_m_16[(30-15):86,]
      Prob_death_m_30[,1] <- Prob_death_m_30[,1] / sum(Prob_death_m_30[,1])
      Prob_death_m_30[,1] <- cumsum(Prob_death_m_30[,1])
      
      
      Prob_death_f_30 <- Prob_death_f_16[(30-15):86,]
      Prob_death_f_30[,1] <- Prob_death_f_30[,1] / sum(Prob_death_f_30[,1])
      Prob_death_f_30[,1] <- cumsum(Prob_death_f_30[,1])
      
      
      #31
      Prob_death_m_31 <- Prob_death_m_16[(31-15):86,]
      Prob_death_m_31[,1] <- Prob_death_m_31[,1] / sum(Prob_death_m_31[,1])
      Prob_death_m_31[,1] <- cumsum(Prob_death_m_31[,1])
      
      Prob_death_f_31 <- Prob_death_f_16[(31-15):86,]
      Prob_death_f_31[,1] <- Prob_death_f_31[,1] / sum(Prob_death_f_31[,1])
      Prob_death_f_31[,1] <- cumsum(Prob_death_f_31[,1])
      
      #32
      Prob_death_m_32 <- Prob_death_m_16[(32-15):86,]
      Prob_death_m_32[,1] <- Prob_death_m_32[,1] / sum(Prob_death_m_32[,1])
      Prob_death_m_32[,1] <- cumsum(Prob_death_m_32[,1])
      
      Prob_death_f_32 <- Prob_death_f_16[(32-15):86,]
      Prob_death_f_32[,1] <- Prob_death_f_32[,1] / sum(Prob_death_f_32[,1])
      Prob_death_f_32[,1] <- cumsum(Prob_death_f_32[,1])
      
      #33
      Prob_death_m_33 <- Prob_death_m_16[(33-15):86,]
      Prob_death_m_33[,1] <- Prob_death_m_33[,1] / sum(Prob_death_m_33[,1])
      Prob_death_m_33[,1] <- cumsum(Prob_death_m_33[,1])
      
      
      Prob_death_f_33 <- Prob_death_f_16[(33-15):86,]
      Prob_death_f_33[,1] <- Prob_death_f_33[,1] / sum(Prob_death_f_33[,1])
      Prob_death_f_33[,1] <- cumsum(Prob_death_f_33[,1])
      
      
      #34
      Prob_death_m_34 <- Prob_death_m_16[(34-15):86,]
      Prob_death_m_34[,1] <- Prob_death_m_34[,1] / sum(Prob_death_m_34[,1])
      Prob_death_m_34[,1] <- cumsum(Prob_death_m_34[,1])
      
      
      Prob_death_f_34 <- Prob_death_f_16[(34-15):86,]
      Prob_death_f_34[,1] <- Prob_death_f_34[,1] / sum(Prob_death_f_34[,1])
      Prob_death_f_34[,1] <- cumsum(Prob_death_f_34[,1])
      
      
      #35
      Prob_death_m_35 <- Prob_death_m_16[(35-15):86,]
      Prob_death_m_35[,1] <- Prob_death_m_35[,1] / sum(Prob_death_m_35[,1])
      Prob_death_m_35[,1] <- cumsum(Prob_death_m_35[,1])
      
      Prob_death_f_35 <- Prob_death_f_16[(35-15):86,]
      Prob_death_f_35[,1] <- Prob_death_f_35[,1] / sum(Prob_death_f_35[,1])
      Prob_death_f_35[,1] <- cumsum(Prob_death_f_35[,1])
      
      
      #36
      Prob_death_m_36 <- Prob_death_m_16[(36-15):86,]
      Prob_death_m_36[,1] <- Prob_death_m_36[,1] / sum(Prob_death_m_36[,1])
      Prob_death_m_36[,1] <- cumsum(Prob_death_m_36[,1])
      
      
      Prob_death_f_36 <- Prob_death_f_16[(36-15):86,]
      Prob_death_f_36[,1] <- Prob_death_f_36[,1] / sum(Prob_death_f_36[,1])
      Prob_death_f_36[,1] <- cumsum(Prob_death_f_36[,1])
      
      
      #37
      Prob_death_m_37 <- Prob_death_m_16[(37-15):86,]
      Prob_death_m_37[,1] <- Prob_death_m_37[,1] / sum(Prob_death_m_37[,1])
      Prob_death_m_37[,1] <- cumsum(Prob_death_m_37[,1])
      
      
      Prob_death_f_37 <- Prob_death_f_16[(37-15):86,]
      Prob_death_f_37[,1] <- Prob_death_f_37[,1] / sum(Prob_death_f_37[,1])
      Prob_death_f_37[,1] <- cumsum(Prob_death_f_37[,1])
      
      
      #38
      Prob_death_m_38 <- Prob_death_m_16[(38-15):86,]
      Prob_death_m_38[,1] <- Prob_death_m_38[,1] / sum(Prob_death_m_38[,1])
      Prob_death_m_38[,1] <- cumsum(Prob_death_m_38[,1])
      
      Prob_death_f_38 <- Prob_death_f_16[(38-15):86,]
      Prob_death_f_38[,1] <- Prob_death_f_38[,1] / sum(Prob_death_f_38[,1])
      Prob_death_f_38[,1] <- cumsum(Prob_death_f_38[,1])
      
      #39
      Prob_death_m_39 <- Prob_death_m_16[(39-15):86,]
      Prob_death_m_39[,1] <- Prob_death_m_39[,1] / sum(Prob_death_m_39[,1])
      Prob_death_m_39[,1] <- cumsum(Prob_death_m_39[,1])
      
      Prob_death_f_39 <- Prob_death_f_16[(39-15):86,]
      Prob_death_f_39[,1] <- Prob_death_f_39[,1] / sum(Prob_death_f_39[,1])
      Prob_death_f_39[,1] <- cumsum(Prob_death_f_39[,1])
      
      #40
      Prob_death_m_40 <- Prob_death_m_16[(40-15):86,]
      Prob_death_m_40[,1] <- Prob_death_m_40[,1] / sum(Prob_death_m_40[,1])
      Prob_death_m_40[,1] <- cumsum(Prob_death_m_40[,1])
      
      Prob_death_f_40 <- Prob_death_f_16[(40-15):86,]
      Prob_death_f_40[,1] <- Prob_death_f_40[,1] / sum(Prob_death_f_40[,1])
      Prob_death_f_40[,1] <- cumsum(Prob_death_f_40[,1])
      
      #41
      Prob_death_m_41 <- Prob_death_m_16[(41-15):86,]
      Prob_death_m_41[,1] <- Prob_death_m_41[,1] / sum(Prob_death_m_41[,1])
      Prob_death_m_41[,1] <- cumsum(Prob_death_m_41[,1])
      
      Prob_death_f_41 <- Prob_death_f_16[(41-15):86,]
      Prob_death_f_41[,1] <- Prob_death_f_41[,1] / sum(Prob_death_f_41[,1])
      Prob_death_f_41[,1] <- cumsum(Prob_death_f_41[,1])
      
      #42
      Prob_death_m_42 <- Prob_death_m_16[(42-15):86,]
      Prob_death_m_42[,1] <- Prob_death_m_42[,1] / sum(Prob_death_m_42[,1])
      Prob_death_m_42[,1] <- cumsum(Prob_death_m_42[,1])
      
      
      Prob_death_f_42 <- Prob_death_f_16[(42-15):86,]
      Prob_death_f_42[,1] <- Prob_death_f_42[,1] / sum(Prob_death_f_42[,1])
      Prob_death_f_42[,1] <- cumsum(Prob_death_f_42[,1])
      
      
      #43
      Prob_death_m_43 <- Prob_death_m_16[(43-15):86,]
      Prob_death_m_43[,1] <- Prob_death_m_43[,1] / sum(Prob_death_m_43[,1])
      Prob_death_m_43[,1] <- cumsum(Prob_death_m_43[,1])
      
    
      Prob_death_f_43 <- Prob_death_f_16[(43-15):86,]
      Prob_death_f_43[,1] <- Prob_death_f_43[,1] / sum(Prob_death_f_43[,1])
      Prob_death_f_43[,1] <- cumsum(Prob_death_f_43[,1])
      
      
      #44
      Prob_death_m_44 <- Prob_death_m_16[(44-15):86,]
      Prob_death_m_44[,1] <- Prob_death_m_44[,1] / sum(Prob_death_m_44[,1])
      Prob_death_m_44[,1] <- cumsum(Prob_death_m_44[,1])
      
      
      Prob_death_f_44 <- Prob_death_f_16[(44-15):86,]
      Prob_death_f_44[,1] <- Prob_death_f_44[,1] / sum(Prob_death_f_44[,1])
      Prob_death_f_44[,1] <- cumsum(Prob_death_f_44[,1])
      
      
      #45
      Prob_death_m_45 <- Prob_death_m_16[(45-15):86,]
      Prob_death_m_45[,1] <- Prob_death_m_45[,1] / sum(Prob_death_m_45[,1])
      Prob_death_m_45[,1] <- cumsum(Prob_death_m_45[,1])
      
      
      Prob_death_f_45 <- Prob_death_f_16[(45-15):86,]
      Prob_death_f_45[,1] <- Prob_death_f_45[,1] / sum(Prob_death_f_45[,1])
      Prob_death_f_45[,1] <- cumsum(Prob_death_f_45[,1])
      
      
      #46
      Prob_death_m_46 <- Prob_death_m_16[(46-15):86,]
      Prob_death_m_46[,1] <- Prob_death_m_46[,1] / sum(Prob_death_m_46[,1])
      Prob_death_m_46[,1] <- cumsum(Prob_death_m_46[,1])
      
      
      Prob_death_f_46 <- Prob_death_f_16[(46-15):86,]
      Prob_death_f_46[,1] <- Prob_death_f_46[,1] / sum(Prob_death_f_46[,1])
      Prob_death_f_46[,1] <- cumsum(Prob_death_f_46[,1])
      
      
      #47
      Prob_death_m_47 <- Prob_death_m_16[(47-15):86,]
      Prob_death_m_47[,1] <- Prob_death_m_47[,1] / sum(Prob_death_m_47[,1])
      Prob_death_m_47[,1] <- cumsum(Prob_death_m_47[,1])
      
      
      Prob_death_f_47 <- Prob_death_f_16[(47-15):86,]
      Prob_death_f_47[,1] <- Prob_death_f_47[,1] / sum(Prob_death_f_47[,1])
      Prob_death_f_47[,1] <- cumsum(Prob_death_f_47[,1])
      
      
      #48
      Prob_death_m_48 <- Prob_death_m_16[(48-15):86,]
      Prob_death_m_48[,1] <- Prob_death_m_48[,1] / sum(Prob_death_m_48[,1])
      Prob_death_m_48[,1] <- cumsum(Prob_death_m_48[,1])
      
      
      Prob_death_f_48 <- Prob_death_f_16[(48-15):86,]
      Prob_death_f_48[,1] <- Prob_death_f_48[,1] / sum(Prob_death_f_48[,1])
      Prob_death_f_48[,1] <- cumsum(Prob_death_f_48[,1])
      
      
      #49
      Prob_death_m_49 <- Prob_death_m_16[(49-15):86,]
      Prob_death_m_49[,1] <- Prob_death_m_49[,1] / sum(Prob_death_m_49[,1])
      Prob_death_m_49[,1] <- cumsum(Prob_death_m_49[,1])
      
      
      Prob_death_f_49 <- Prob_death_f_16[(49-15):86,]
      Prob_death_f_49[,1] <- Prob_death_f_49[,1] / sum(Prob_death_f_49[,1])
      Prob_death_f_49[,1] <- cumsum(Prob_death_f_49[,1])
      
      
      #50
      Prob_death_m_50 <- Prob_death_m_16[(50-15):86,]
      Prob_death_m_50[,1] <- Prob_death_m_50[,1] / sum(Prob_death_m_50[,1])
      Prob_death_m_50[,1] <- cumsum(Prob_death_m_50[,1])
      
      
      Prob_death_f_50 <- Prob_death_f_16[(50-15):86,]
      Prob_death_f_50[,1] <- Prob_death_f_50[,1] / sum(Prob_death_f_50[,1])
      Prob_death_f_50[,1] <- cumsum(Prob_death_f_50[,1])
      
      
      #51
      Prob_death_m_51 <- Prob_death_m_16[(51-15):86,]
      Prob_death_m_51[,1] <- Prob_death_m_51[,1] / sum(Prob_death_m_51[,1])
      Prob_death_m_51[,1] <- cumsum(Prob_death_m_51[,1])
      
      
      Prob_death_f_51 <- Prob_death_f_16[(51-15):86,]
      Prob_death_f_51[,1] <- Prob_death_f_51[,1] / sum(Prob_death_f_51[,1])
      Prob_death_f_51[,1] <- cumsum(Prob_death_f_51[,1])
      
      
      #52
      Prob_death_m_52 <- Prob_death_m_16[(52-15):86,]
      Prob_death_m_52[,1] <- Prob_death_m_52[,1] / sum(Prob_death_m_52[,1])
      Prob_death_m_52[,1] <- cumsum(Prob_death_m_52[,1])
      
      
      Prob_death_f_52 <- Prob_death_f_16[(52-15):86,]
      Prob_death_f_52[,1] <- Prob_death_f_52[,1] / sum(Prob_death_f_52[,1])
      Prob_death_f_52[,1] <- cumsum(Prob_death_f_52[,1])
      
      
      #53
      Prob_death_m_53 <- Prob_death_m_16[(53-15):86,]
      Prob_death_m_53[,1] <- Prob_death_m_53[,1] / sum(Prob_death_m_53[,1])
      Prob_death_m_53[,1] <- cumsum(Prob_death_m_53[,1])
      
      
      Prob_death_f_53 <- Prob_death_f_16[(53-15):86,]
      Prob_death_f_53[,1] <- Prob_death_f_53[,1] / sum(Prob_death_f_53[,1])
      Prob_death_f_53[,1] <- cumsum(Prob_death_f_53[,1])
      
      
      #54
      Prob_death_m_54 <- Prob_death_m_16[(54-15):86,]
      Prob_death_m_54[,1] <- Prob_death_m_54[,1] / sum(Prob_death_m_54[,1])
      Prob_death_m_54[,1] <- cumsum(Prob_death_m_54[,1])
      
      
      Prob_death_f_54 <- Prob_death_f_16[(54-15):86,]
      Prob_death_f_54[,1] <- Prob_death_f_54[,1] / sum(Prob_death_f_54[,1])
      Prob_death_f_54[,1] <- cumsum(Prob_death_f_54[,1])
      
      
      #55
      Prob_death_m_55 <- Prob_death_m_16[(55-15):86,]
      Prob_death_m_55[,1] <- Prob_death_m_55[,1] / sum(Prob_death_m_55[,1])
      Prob_death_m_54[,1] <- cumsum(Prob_death_m54[,1])
      
      
      Prob_death_f_55 <- Prob_death_f_16[(55-15):86,]
      Prob_death_f_55[,1] <- Prob_death_f_55[,1] / sum(Prob_death_f_55[,1])
      Prob_death_f_55[,1] <- cumsum(Prob_death_f_55[,1])
      
      
      #56
      Prob_death_m_56 <- Prob_death_m_16[(56-15):86,]
      Prob_death_m_56[,1] <- Prob_death_m_56[,1] / sum(Prob_death_m_56[,1])
      
      Prob_death_f_56 <- Prob_death_f_16[(56-15):86,]
      Prob_death_f_56[,1] <- Prob_death_f_56[,1] / sum(Prob_death_f_56[,1])
      
      
      #57
      Prob_death_m_57 <- Prob_death_m_16[(57-15):86,]
      Prob_death_m_57[,1] <- Prob_death_m_57[,1] / sum(Prob_death_m_57[,1])
      
      Prob_death_f_57 <- Prob_death_f_16[(57-15):86,]
      Prob_death_f_57[,1] <- Prob_death_f_57[,1] / sum(Prob_death_f_57[,1])
      
      #58
      Prob_death_m_58 <- Prob_death_m_16[(58-15):86,]
      Prob_death_m_58[,1] <- Prob_death_m_58[,1] / sum(Prob_death_m_58[,1])
      
      Prob_death_f_58 <- Prob_death_f_16[(58-15):86,]
      Prob_death_f_58[,1] <- Prob_death_f_58[,1] / sum(Prob_death_f_58[,1])
      
      #59
      Prob_death_m_59 <- Prob_death_m_16[(59-15):86,]
      Prob_death_m_59[,1] <- Prob_death_m_59[,1] / sum(Prob_death_m_59[,1])
      
      Prob_death_f_59 <- Prob_death_f_16[(59-15):86,]
      Prob_death_f_59[,1] <- Prob_death_f_59[,1] / sum(Prob_death_f_59[,1])
      
      #60
      Prob_death_m_60 <- Prob_death_m_16[(60-15):86,]
      Prob_death_m_60[,1] <- Prob_death_m_60[,1] / sum(Prob_death_m_60[,1])
      
      Prob_death_f_60 <- Prob_death_f_16[(60-15):86,]
      Prob_death_f_60[,1] <- Prob_death_f_60[,1] / sum(Prob_death_f_60[,1])
      
      #61
      Prob_death_m_61 <- Prob_death_m_16[(61-15):86,]
      Prob_death_m_61[,1] <- Prob_death_m_61[,1] / sum(Prob_death_m_61[,1])
      
      Prob_death_f_61 <- Prob_death_f_16[(61-15):86,]
      Prob_death_f_61[,1] <- Prob_death_f_61[,1] / sum(Prob_death_f_61[,1])
      
      
      #62
      Prob_death_m_62 <- Prob_death_m_16[(62-15):86,]
      Prob_death_m_62[,1] <- Prob_death_m_62[,1] / sum(Prob_death_m_62[,1])
      
      Prob_death_f_62 <- Prob_death_f_16[(62-15):86,]
      Prob_death_f_62[,1] <- Prob_death_f_62[,1] / sum(Prob_death_f_62[,1])
      
      #63
      Prob_death_m_63 <- Prob_death_m_16[(63-15):86,]
      Prob_death_m_63[,1] <- Prob_death_m_63[,1] / sum(Prob_death_m_63[,1])
      
      Prob_death_f_63 <- Prob_death_f_16[(63-15):86,]
      Prob_death_f_63[,1] <- Prob_death_f_63[,1] / sum(Prob_death_f_63[,1])
      
      #64
      Prob_death_m_64 <- Prob_death_m_16[(64-15):86,]
      Prob_death_m_64[,1] <- Prob_death_m_64[,1] / sum(Prob_death_m_64[,1])
      
      Prob_death_f_64 <- Prob_death_f_16[(64-15):86,]
      Prob_death_f_64[,1] <- Prob_death_f_64[,1] / sum(Prob_death_f_64[,1])
      
      #65
      Prob_death_m_65 <- Prob_death_m_16[(65-15):86,]
      Prob_death_m_65[,1] <- Prob_death_m_65[,1] / sum(Prob_death_m_65[,1])
      
      Prob_death_f_65 <- Prob_death_f_16[(65-15):86,]
      Prob_death_f_65[,1] <- Prob_death_f_65[,1] / sum(Prob_death_f_65[,1])
      
      #66
      Prob_death_m_66 <- Prob_death_m_16[(66-15):86,]
      Prob_death_m_66[,1] <- Prob_death_m_66[,1] / sum(Prob_death_m_66[,1])
      
      Prob_death_f_66 <- Prob_death_f_16[(66-15):86,]
      Prob_death_f_66[,1] <- Prob_death_f_66[,1] / sum(Prob_death_f_66[,1])
      
      #67
      Prob_death_m_67 <- Prob_death_m_16[(67-15):86,]
      Prob_death_m_67[,1] <- Prob_death_m_67[,1] / sum(Prob_death_m_67[,1])
      
      Prob_death_f_67 <- Prob_death_f_16[(67-15):86,]
      Prob_death_f_67[,1] <- Prob_death_f_67[,1] / sum(Prob_death_f_67[,1])
      
      #68
      Prob_death_m_68 <- Prob_death_m_16[(68-15):86,]
      Prob_death_m_68[,1] <- Prob_death_m_68[,1] / sum(Prob_death_m_68[,1])
      
      Prob_death_f_68 <- Prob_death_f_16[(68-15):86,]
      Prob_death_f_68[,1] <- Prob_death_f_68[,1] / sum(Prob_death_f_68[,1])
      
      
      #69
      Prob_death_m_69 <- Prob_death_m_16[(69-15):86,]
      Prob_death_m_69[,1] <- Prob_death_m_69[,1] / sum(Prob_death_m_69[,1])
      
      Prob_death_f_69 <- Prob_death_f_16[(69-15):86,]
      Prob_death_f_69[,1] <- Prob_death_f_69[,1] / sum(Prob_death_f_69[,1])
      
      #70
      Prob_death_m_70 <- Prob_death_m_16[(70-15):86,]
      Prob_death_m_70[,1] <- Prob_death_m_70[,1] / sum(Prob_death_m_70[,1])
      
      Prob_death_f_70 <- Prob_death_f_16[(70-15):86,]
      Prob_death_f_70[,1] <- Prob_death_f_70[,1] / sum(Prob_death_f_70[,1])
      
      #71
      Prob_death_m_71 <- Prob_death_m_16[(71-15):86,]
      Prob_death_m_71[,1] <- Prob_death_m_71[,1] / sum(Prob_death_m_71[,1])
      
      Prob_death_f_71 <- Prob_death_f_16[(71-15):86,]
      Prob_death_f_71[,1] <- Prob_death_f_71[,1] / sum(Prob_death_f_71[,1])
      
      #72
      Prob_death_m_72 <- Prob_death_m_16[(72-15):86,]
      Prob_death_m_72[,1] <- Prob_death_m_72[,1] / sum(Prob_death_m_72[,1])
      
      Prob_death_f_72 <- Prob_death_f_16[(72-15):86,]
      Prob_death_f_72[,1] <- Prob_death_f_72[,1] / sum(Prob_death_f_72[,1])
      
      #73
      Prob_death_m_73 <- Prob_death_m_16[(73-15):86,]
      Prob_death_m_73[,1] <- Prob_death_m_73[,1] / sum(Prob_death_m_73[,1])
      
      Prob_death_f_73 <- Prob_death_f_16[(73-15):86,]
      Prob_death_f_73[,1] <- Prob_death_f_73[,1] / sum(Prob_death_f_73[,1])
      
      #74
      Prob_death_m_74 <- Prob_death_m_16[(74-15):86,]
      Prob_death_m_74[,1] <- Prob_death_m_74[,1] / sum(Prob_death_m_74[,1])
      
      Prob_death_f_74 <- Prob_death_f_16[(74-15):86,]
      Prob_death_f_74[,1] <- Prob_death_f_74[,1] / sum(Prob_death_f_74[,1])
      
      #75
      Prob_death_m_75 <- Prob_death_m_16[(75-15):86,]
      Prob_death_m_75[,1] <- Prob_death_m_75[,1] / sum(Prob_death_m_75[,1])
      
      Prob_death_f_75 <- Prob_death_f_16[(75-15):86,]
      Prob_death_f_75[,1] <- Prob_death_f_75[,1] / sum(Prob_death_f_75[,1])
      
      #76
      Prob_death_m_76 <- Prob_death_m_16[(76-15):86,]
      Prob_death_m_76[,1] <- Prob_death_m_76[,1] / sum(Prob_death_m_76[,1])
      
      Prob_death_f_76 <- Prob_death_f_16[(76-15):86,]
      Prob_death_f_76[,1] <- Prob_death_f_76[,1] / sum(Prob_death_f_76[,1])
      
      #77
      Prob_death_m_77 <- Prob_death_m_16[(77-15):86,]
      Prob_death_m_77[,1] <- Prob_death_m_77[,1] / sum(Prob_death_m_77[,1])
      
      Prob_death_f_77 <- Prob_death_f_16[(77-15):86,]
      Prob_death_f_77[,1] <- Prob_death_f_77[,1] / sum(Prob_death_f_77[,1])
      
      
      #78
      Prob_death_m_78 <- Prob_death_m_16[(78-15):86,]
      Prob_death_m_78[,1] <- Prob_death_m_78[,1] / sum(Prob_death_m_78[,1])
      
      Prob_death_f_78 <- Prob_death_f_16[(78-15):86,]
      Prob_death_f_78[,1] <- Prob_death_f_78[,1] / sum(Prob_death_f_78[,1])
      
      #79
      Prob_death_m_79 <- Prob_death_m_16[(79-15):86,]
      Prob_death_m_79[,1] <- Prob_death_m_79[,1] / sum(Prob_death_m_79[,1])
      
      Prob_death_f_79 <- Prob_death_f_16[(79-15):86,]
      Prob_death_f_79[,1] <- Prob_death_f_79[,1] / sum(Prob_death_f_79[,1])
      
      #80
      Prob_death_m_80 <- Prob_death_m_16[(80-15):86,]
      Prob_death_m_80[,1] <- Prob_death_m_80[,1] / sum(Prob_death_m_80[,1])
      
      Prob_death_f_80 <- Prob_death_f_16[(80-15):86,]
      Prob_death_f_80[,1] <- Prob_death_f_80[,1] / sum(Prob_death_f_80[,1])
      
      #81
      Prob_death_m_81 <- Prob_death_m_16[(81-15):86,]
      Prob_death_m_81[,1] <- Prob_death_m_81[,1] / sum(Prob_death_m_81[,1])
      
      Prob_death_f_81 <- Prob_death_f_16[(81-15):86,]
      Prob_death_f_81[,1] <- Prob_death_f_81[,1] / sum(Prob_death_f_81[,1])
      
      #82
      Prob_death_m_82 <- Prob_death_m_16[(82-15):86,]
      Prob_death_m_82[,1] <- Prob_death_m_82[,1] / sum(Prob_death_m_82[,1])
      
      Prob_death_f_82 <- Prob_death_f_16[(82-15):86,]
      Prob_death_f_82[,1] <- Prob_death_f_82[,1] / sum(Prob_death_f_82[,1])
      
      #83
      Prob_death_m_83 <- Prob_death_m_16[(83-15):86,]
      Prob_death_m_83[,1] <- Prob_death_m_83[,1] / sum(Prob_death_m_83[,1])
      
      Prob_death_f_83 <- Prob_death_f_16[(83-15):86,]
      Prob_death_f_83[,1] <- Prob_death_f_83[,1] / sum(Prob_death_f_83[,1])
      
      #84
      Prob_death_m_84 <- Prob_death_m_16[(84-15):86,]
      Prob_death_m_84[,1] <- Prob_death_m_84[,1] / sum(Prob_death_m_84[,1])
      
      Prob_death_f_84 <- Prob_death_f_16[(84-15):86,]
      Prob_death_f_84[,1] <- Prob_death_f_84[,1] / sum(Prob_death_f_84[,1])
      
      #85
      Prob_death_m_85 <- Prob_death_m_16[(85-15):86,]
      Prob_death_m_85[,1] <- Prob_death_m_85[,1] / sum(Prob_death_m_85[,1])
      
      Prob_death_f_85 <- Prob_death_f_16[(85-15):86,]
      Prob_death_f_85[,1] <- Prob_death_f_85[,1] / sum(Prob_death_f_85[,1])
      
      #86
      Prob_death_m_86 <- Prob_death_m_16[(86-15):86,]
      Prob_death_m_86[,1] <- Prob_death_m_86[,1] / sum(Prob_death_m_86[,1])
      
      Prob_death_f_86 <- Prob_death_f_16[(86-15):86,]
      Prob_death_f_86[,1] <- Prob_death_f_86[,1] / sum(Prob_death_f_86[,1])
      
      #87
      Prob_death_m_87 <- Prob_death_m_16[(87-15):86,]
      Prob_death_m_87[,1] <- Prob_death_m_87[,1] / sum(Prob_death_m_87[,1])
      
      Prob_death_f_87 <- Prob_death_f_16[(87-15):86,]
      Prob_death_f_87[,1] <- Prob_death_f_87[,1] / sum(Prob_death_f_87[,1])
      
      #88
      Prob_death_m_88 <- Prob_death_m_16[(88-15):86,]
      Prob_death_m_88[,1] <- Prob_death_m_88[,1] / sum(Prob_death_m_88[,1])
      
      Prob_death_f_88 <- Prob_death_f_16[(88-15):86,]
      Prob_death_f_88[,1] <- Prob_death_f_88[,1] / sum(Prob_death_f_88[,1])
      
      #89
      Prob_death_m_89 <- Prob_death_m_16[(89-15):86,]
      Prob_death_m_89[,1] <- Prob_death_m_89[,1] / sum(Prob_death_m_89[,1])
      
      Prob_death_f_89 <- Prob_death_f_16[(89-15):86,]
      Prob_death_f_89[,1] <- Prob_death_f_89[,1] / sum(Prob_death_f_89[,1])
      
      #90
      Prob_death_m_90 <- Prob_death_m_16[(90-15):86,]
      Prob_death_m_90[,1] <- Prob_death_m_90[,1] / sum(Prob_death_m_90[,1])
      
      Prob_death_f_90 <- Prob_death_f_16[(90-15):86,]
      Prob_death_f_90[,1] <- Prob_death_f_90[,1] / sum(Prob_death_f_90[,1])
      
      #91
      Prob_death_m_91 <- Prob_death_m_16[(91-15):86,]
      Prob_death_m_91[,1] <- Prob_death_m_91[,1] / sum(Prob_death_m_91[,1])
      
      Prob_death_f_91 <- Prob_death_f_16[(91-15):86,]
      Prob_death_f_91[,1] <- Prob_death_f_91[,1] / sum(Prob_death_f_91[,1])
      
      #92
      Prob_death_m_92 <- Prob_death_m_16[(92-15):86,]
      Prob_death_m_92[,1] <- Prob_death_m_92[,1] / sum(Prob_death_m_92[,1])
      
      Prob_death_f_92 <- Prob_death_f_16[(92-15):86,]
      Prob_death_f_92[,1] <- Prob_death_f_92[,1] / sum(Prob_death_f_92[,1])
      
      #93
      Prob_death_m_93 <- Prob_death_m_16[(93-15):86,]
      Prob_death_m_93[,1] <- Prob_death_m_93[,1] / sum(Prob_death_m_93[,1])
      
      Prob_death_f_93 <- Prob_death_f_16[(93-15):86,]
      Prob_death_f_93[,1] <- Prob_death_f_93[,1] / sum(Prob_death_f_93[,1])
      
      #94
      Prob_death_m_94 <- Prob_death_m_16[(94-15):86,]
      Prob_death_m_94[,1] <- Prob_death_m_94[,1] / sum(Prob_death_m_94[,1])
      
      Prob_death_f_94 <- Prob_death_f_16[(94-15):86,]
      Prob_death_f_94[,1] <- Prob_death_f_94[,1] / sum(Prob_death_f_94[,1])
      
      #95
      Prob_death_m_95 <- Prob_death_m_16[(95-15):86,]
      Prob_death_m_95[,1] <- Prob_death_m_95[,1] / sum(Prob_death_m_95[,1])
      
      Prob_death_f_95 <- Prob_death_f_16[(95-15):86,]
      Prob_death_f_95[,1] <- Prob_death_f_95[,1] / sum(Prob_death_f_95[,1])
      
      #96
      Prob_death_m_96 <- Prob_death_m_16[(96-15):86,]
      Prob_death_m_96[,1] <- Prob_death_m_96[,1] / sum(Prob_death_m_96[,1])
      
      Prob_death_f_96 <- Prob_death_f_16[(96-15):86,]
      Prob_death_f_96[,1] <- Prob_death_f_96[,1] / sum(Prob_death_f_96[,1])
      
      #97
      Prob_death_m_97 <- Prob_death_m_16[(97-15):86,]
      Prob_death_m_97[,1] <- Prob_death_m_97[,1] / sum(Prob_death_m_97[,1])
      
      Prob_death_f_97 <- Prob_death_f_16[(97-15):86,]
      Prob_death_f_97[,1] <- Prob_death_f_97[,1] / sum(Prob_death_f_97[,1])
      
      #98
      Prob_death_m_98 <- Prob_death_m_16[(98-15):86,]
      Prob_death_m_98[,1] <- Prob_death_m_98[,1] / sum(Prob_death_m_98[,1])
      
      Prob_death_f_98 <- Prob_death_f_16[(98-15):86,]
      Prob_death_f_98[,1] <- Prob_death_f_98[,1] / sum(Prob_death_f_98[,1])
      
      #99
      Prob_death_m_99 <- Prob_death_m_16[(99-15):86,]
      Prob_death_m_99[,1] <- Prob_death_m_99[,1] / sum(Prob_death_m_99[,1])
      
      Prob_death_f_99 <- Prob_death_f_16[(99-15):86,]
      Prob_death_f_99[,1] <- Prob_death_f_99[,1] / sum(Prob_death_f_99[,1])
      
      #100
      Prob_death_m_100 <- Prob_death_m_16[(100-15):86,]
      Prob_death_m_100[,1] <- Prob_death_m_100[,1] / sum(Prob_death_m_100[,1])
      
      Prob_death_f_100 <- Prob_death_f_16[(100-15):86,]
      Prob_death_f_100[,1] <- Prob_death_f_100[,1] / sum(Prob_death_f_100[,1])
      
      #101
      Prob_death_m_101 <- Prob_death_m_16[(101-15):86,]
      Prob_death_m_101[1] <- 1
      
      Prob_death_f_101 <- Prob_death_f_16[(101-15):86,]
      Prob_death_f_101[1] <- 1
      
      
      
    #create a vector of random numbers associated with each patient
      rands <- runif(length(age))
      
      if(age[]==16&gender[]==1){
        age_death <- runif[]
      }
      
      
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
  return(all_cause_death)
}


