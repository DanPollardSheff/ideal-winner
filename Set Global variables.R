#Global variables
PSA_switch <- 1                            #1=run PSA, 0=deterministic
PSA_numb <- 1000                            #number of PSA runs
pat_numb <- 100000                          #number of patients
days_to_discharge <- 30                     #number of days to discharge from hospital
days_in_year <- 365.25                      #number of days in a year
time_horizon <- 100                         #time horizon, years
discount_rate_QALYs <- 0.035                #discount rate, QALYs
discount_rate_costs <- 0.035                #discount rate, costs
Param_export <- 1                           #1=save a copy of PSA parameters

#The proportion of clinical benefit that patients with an ISS of over 8 and under 16 receive
#compared to people with an ISS of 16 or more (0 to 1). Default is no benefit (0)
Proportion_RR_MTC_ISS_o8_u16_hosp <- 0      #Continuous number
Proportion_RR_MTC_ISS_o8_u16_1yr <- 0       #Continuous number

#The proportion of clinical benefit that patients who are initially sent to a non-MTC receive
#compared to people sent straight to an MTC (0 to 1). Default is full benefit (1)
Proportion_RR_MTC_transfer_hosp <- 1        #Continuous number
Proportion_RR_MTC_ISS_transfer_1yr <- 1     #Continuous number


TARN_mort_eq <- "Old"                       #options are new or old. Default is old
MTCs_in_mort_risk <- "No"                   #options are Yes or no. Relates to whether the mort eq is a composite risk score for a 

#population who has / has not been to an MTC or a population who hasn't gone to an MTC. Default is no, as the
#default for the  mortality equation is the Old TARN equation.
percent_TARN_cases_reported_ISS_o16 <- 1    #Continuous number
percent_TARN_cases_reported_ISS_o9_u16 <- 1

population_source <- "Dutch"                #Source of simulated population.Options are UK, 
                                            #Dutch, Dutch_simp. Dutch is the default

population_ISS_over16_only <- "No"          #Option for resampling to produce a population with ISS >= 16. Options are "Yes" or "No". Default is no. 
population_ISS_under16_only <- "No"         #Option for resampling to produce a population with ISS < 16. Can be "Yes" or "No". Default is no.

PSA_rand_no <-  -99                         #random number to determine PSA parameters either -99 (to not reset the seed) or any positive number

scenario <- "_basecase"                         #name to append to saved files 


Eldery_specific_params  <- T                #Takes value T or F. If T model has different parameters for
                                            #elderly (65+) populations

Pead_specific_params    <- F                #Takes value T or F. If T the model has different parameters
                                            #for pediatric (14 and under) population

TARN_22_params          <- T                #Use the TARN 22 parameter estimates? T = TRUE, F = FALSE

#If using TARN22 parameters, set the TARN mortality equation to new and set the option 
#to indicate that MTCs were implemented in the estimates of MTC mortality risk
#This is to avoid manual errors
if(TARN_22_params ==T){
  TARN_mort_eq <- "New" 
  MTCs_in_mort_risk <- "Yes"
}

Util_source             <- "Kruithof"        #Option for the source of the utility values. Either Ahmed or Kruithof. Default is Kruithof