#####Load / install packages
#install.packages("ggplot2")
#install.packages("shiny")
library("ggplot2")
library("shiny")

##Variable controls
#variables
#base case threshold
Threshold <- 20000
#maximum cost-effectiveness threshold to consider
Max_threshold <- 50000
#how many steps you want to take in the cost-effecitveness threshold
Threshold_incr <- 100
#decimal places on your costs, NMBs and ICERS
dp_costs <- 2
#decimal places on your QALYs
dp_QALYS <- 0

#Method switch, determines whether the best strategy on average is compared to the next best strategy on average
#or the best remaining strategy within the PSA run. Default is next best on average. Options are "Average" or "PSA"
#default is Average
Method_Hat <- "Average"

###set the working directory
setwd("C:\\Users\\polla\\Documents\\GitHub\\fantastic-lamp\\Major Trauma Triage v2")

#load my data and organise into a data frame of discounted costs and discounted QALYs
sens_28_spec_89 <- read.csv("PSA\\sens_28_spec_89.csv", header = TRUE, row.names = 1)
sens_57_spec_80 <- read.csv("PSA\\sens_57_spec_80_PSA.csv", header = TRUE, row.names = 1)
sens_64_spec_76 <- read.csv("PSA\\sens_64_spec_76_PSA.csv", header = TRUE, row.names = 1)
sens_70_spec_70 <- read.csv("PSA\\sens_70_spec_70_PSA.csv", header = TRUE, row.names = 1)
sens_75_spec_66 <- read.csv("PSA\\sens_75_spec_66_PSA.csv", header = TRUE, row.names = 1)
sens_88_spec_63 <- read.csv("PSA\\sens_88_spec_63_PSA.csv", header = TRUE, row.names = 1)
sens_90_spec_58 <- read.csv("PSA\\sens_90_spec_58_PSA.csv", header = TRUE, row.names = 1)
sens_95_spec_19 <- read.csv("PSA\\sens_95_spec_19_PSA.csv", header = TRUE, row.names = 1)
sens_100_spec_3 <- read.csv("PSA\\sens_100_spec_3_PSA.csv",header = TRUE, row.names = 1)

#names of columns
column_names <- c("sens_28_spec_89","sens_57_spec_80","sens_64_spec_76","sens_70_spec_70","sens_75_spec_66","sens_88_spec_63",
                 "sens_90_spec_58","sens_95_spec_19","sens_100_spec_3")
#names of rows
row_names

#make a matrix of costs & QALYs
Costs <- matrix(data=NA, nrow = length(sens_28_spec_89$Sens_DR), ncol = length(column_names))
QALYs <- matrix(data=NA, nrow = length(sens_28_spec_89$Sens_DR), ncol = length(column_names))

#name the columns
colnames(Costs) <- column_names
colnames(QALYs) <- column_names

#start writting the discounted costs and QALYs into the data matrices
Costs[,"sens_28_spec_89"] <- sens_28_spec_89$discounted_Costs
Costs[,"sens_57_spec_80"] <- sens_57_spec_80$discounted_Costs
Costs[,"sens_64_spec_76"] <- sens_64_spec_76$discounted_Costs
Costs[,"sens_70_spec_70"] <- sens_70_spec_70$discounted_Costs
Costs[,"sens_75_spec_66"] <- sens_75_spec_66$discounted_Costs
Costs[,"sens_88_spec_63"] <- sens_88_spec_63$discounted_Costs
Costs[,"sens_90_spec_58"] <- sens_90_spec_58$discounted_Costs
Costs[,"sens_95_spec_19"] <- sens_95_spec_19$discounted_Costs
Costs[,"sens_100_spec_3"] <- sens_100_spec_3$discounted_Costs

QALYs[,"sens_28_spec_89"] <- sens_28_spec_89$discounted_QALYs
QALYs[,"sens_57_spec_80"] <- sens_57_spec_80$discounted_QALYs
QALYs[,"sens_64_spec_76"] <- sens_64_spec_76$discounted_QALYs
QALYs[,"sens_70_spec_70"] <- sens_70_spec_70$discounted_QALYs
QALYs[,"sens_75_spec_66"] <- sens_75_spec_66$discounted_QALYs
QALYs[,"sens_88_spec_63"] <- sens_88_spec_63$discounted_QALYs
QALYs[,"sens_90_spec_58"] <- sens_90_spec_58$discounted_QALYs
QALYs[,"sens_95_spec_19"] <- sens_95_spec_19$discounted_QALYs
QALYs[,"sens_100_spec_3"] <- sens_100_spec_3$discounted_QALYs

write.csv(QALYs, "QALY outcomes.csv")
write.csv(Costs, "Cost outcomes.csv")
######Analysis code

### Are the results stable ?#########
#Use the methods in Hatswell et al 2018

#Set up a matrix to store the average NMB for each strategy
#extract the column names from the QALYs matrix
names_columns <- colnames(QALYs)
#populate it as an empty matrix, with 1 row and as many columns as there are in the QALYs matrix
NMB <- matrix(data=NA, nrow = length(QALYs[,1]), ncol = length(QALYs[1,]))
Av_NMB <- matrix(data=NA, nrow = 1, ncol = length(QALYs[1,]))
#Make sure that the column names in the av_NMB matrix 
colnames(NMB) <- names_columns
#fill in the matrix with the mean NMB values
NMB <- (QALYs[]*Threshold - Costs[])

#record the average NMB values for each strategy
Av_NMB <- as.data.frame(colMeans(NMB))



#best strategy is sensitivity sens28_spec_89, the next best strategy is sens 57, spec 80
best_strat <- "sens_28_spec_89"
next_best_strat <- "sens_57_spec_80"

diff_NMB <- NMB[,best_strat]- NMB[,next_best_strat]

mean_diff_NMB <- mean(diff_NMB)
SE_diff_NMB <- sd(diff_NMB)/sqrt(length(diff_NMB))

#create an if statement to print whether more PSAs are needed (or not)
#The rule we are checking is there a less than 2.5% chance that the mean NMB is 0 or less
if((mean_diff_NMB/SE_diff_NMB) < qnorm(0.975,0,1)){
  print("More PSAs are needed")
} else{
  print("Sufficent PSAs have been run")
}

#Give and 95% CI around the mean NMB
NMB_95CI_low <- qnorm(0.025,mean_diff_NMB, SE_diff_NMB)
NMB_95CI_high <- qnorm(0.975, mean_diff_NMB, SE_diff_NMB)

#bring the numbers into a string to print
mean <- paste("£",round(mean_diff_NMB,dp_costs), sep="")
low_95 <-  paste("£",round(NMB_95CI_low,dp_costs), sep="")
high_95 <- paste("£", round(NMB_95CI_high, dp_costs), sep="")


temp <-  paste(mean,low_95, sep = " (")
temp <-  paste(temp,high_95, sep = ", ")
temp <-  paste(temp,"", sep=")")

print(temp)

######### CEAC ####################

#Data objects for manipulation
#set of numbers which contain all cost-effectiveness thresholds that I want to run
thres_seq <- seq(0, Max_threshold, by = Threshold_incr)
#set up a matrix to store which strategy is the most cost-effective
CE <- matrix(data=NA, nrow = length(QALYs[,1]), ncol = length(QALYs[1,]))
NMBs <- matrix(data=NA, nrow = length(QALYs[,1]), ncol = length(QALYs[1,])+1)

#create a matrix to store the results
CEAC_res <- matrix(data=NA, nrow=length(thres_seq), ncol   =  length(QALYs[1,])+1)
#make the first column the threshold sequence
CEAC_res[,1] <- thres_seq



for( i in 1:length(thres_seq)){
  #calculate the NMB
  NMBs[,1:length(QALYs[1,])] <- QALYs[,]*thres_seq[i]-Costs[,]
  #calculate whether each strategy is the most cost-effective
  NMBs[,length(QALYs[1,])+1] <- apply(NMBs[,1:length(QALYs[1,])],1,max)
  for(y in 1:length(NMBs[,1])){
  #determine if the value in each column is the maximum in that row
  CE[y,] <- ifelse(NMBs[y,1:length(QALYs[1,])]==max(NMBs[y,length(QALYs[1,])+1]),1,0)
  }
  CEAC_res[i,2:(length(QALYs[1,])+1)] <- colMeans(CE[])
}
