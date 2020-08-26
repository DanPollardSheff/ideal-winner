##### call in data that I'll need
pat_chars <- read.csv("D:\\Work from Home\\Local Git\\Life tables mort prediction test\\pat chars.csv")

age <- pat_chars[,"Age"]
gender <- pat_chars[,"Gender"]

life_tables <- read.csv("D:\\Work from Home\\Local Git\\Life tables mort prediction test\\ONSlifetables.csv")

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
 
    #If the patient is 101 or over when this function is run, assume they have the remaining life expectancy of a 100 year old for their gender
    #Life expectancies are 103.02 for men and 103.32 for women (e.e 2.02 nad 2.32 remaining years of life respectively)
    all_cause_death <-    ifelse(age[]>=101, age[] + ifelse(gender[]==1,2.02,2.32),-99 )
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
    Life_tab_f[i,3] <- Life_tab_f[i-1,3] - Life_tab_f[i-1,3]*Life_tab_f[i-1,2]
    }
    
    #generate cumulative probability of death column by age band
    
      
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