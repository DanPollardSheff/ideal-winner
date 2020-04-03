#set the working directory
#change this for your file location
file_location <- "~/Work from Home/Matts - Model/pat char check/"



#read in the csv files
sens_100_spec_3 <- read.csv(paste(file_location,"sens_100_spec_3.csv",sep=""))
sens_95_spec_19 <- read.csv(paste(file_location,"sens_95_spec_19.csv",sep=""))
sens_90_spec_58 <- read.csv(paste(file_location,"sens_90_spec_58.csv",sep=""))
sens_88_spec_63 <- read.csv(paste(file_location,"sens_88_spec_63.csv",sep=""))
sens_75_spec_66 <- read.csv(paste(file_location,"sens_75_spec_66.csv",sep=""))
sens_70_spec_70 <- read.csv(paste(file_location,"sens_70_spec_70.csv",sep=""))
sens_64_spec_76 <- read.csv(paste(file_location,"sens_64_spec_76.csv",sep=""))
sens_57_spec_80 <- read.csv(paste(file_location,"sens_57_spec_80.csv",sep=""))
sens_28_spec_89 <- read.csv(paste(file_location,"sens_28_spec_89.csv",sep=""))

#set the cost effectiveness threshold
 threshold <- 20000
 
 #Calculate the net benefit for each strategy 
 NMB_sens_100_spec_3 <- sens_100_spec_3[,"dQALYS"]*threshold - sens_100_spec_3[,"DCosts"]
 NMB_sens_95_spec_19 <- sens_95_spec_19[,"dQALYS"]*threshold - sens_95_spec_19[,"DCosts"]
 NMB_sens_90_spec_58 <- sens_90_spec_58[,"dQALYS"]*threshold - sens_90_spec_58[,"DCosts"]
 NMB_sens_88_spec_63 <- sens_88_spec_63[,"dQALYS"]*threshold - sens_88_spec_63[,"DCosts"]
 NMB_sens_75_spec_66 <- sens_75_spec_66[,"dQALYS"]*threshold - sens_75_spec_66[,"DCosts"]
 NMB_sens_64_spec_76 <- sens_64_spec_76[,"dQALYS"]*threshold - sens_64_spec_76[,"DCosts"]
 NMB_sens_57_spec_80 <- sens_57_spec_80[,"dQALYS"]*threshold - sens_57_spec_80[,"DCosts"]
 NMB_sens_28_spec_89 <- sens_28_spec_89[,"dQALYS"]*threshold - sens_28_spec_89[,"DCosts"] 
 
 #Calculate the mean NMB
 mean_NMB_sens_100_spec_3 <- mean(NMB_sens_100_spec_3[])
 mean_NMB_sens_95_spec_19 <- mean(NMB_sens_95_spec_19[])
 mean_NMB_sens_90_spec_58 <- mean(NMB_sens_90_spec_58[])
 mean_NMB_sens_88_spec_63 <- mean(NMB_sens_88_spec_63[])
 mean_NMB_sens_75_spec_66 <- mean(NMB_sens_75_spec_66[])
 mean_NMB_sens_64_spec_76 <- mean(NMB_sens_64_spec_76[])
 mean_NMB_sens_57_spec_80 <- mean(NMB_sens_57_spec_80[])
 mean_NMB_sens_28_spec_89 <- mean(NMB_sens_28_spec_89[])
 
 #return the maximum NMB
 max(mean_NMB_sens_100_spec_3,mean_NMB_sens_95_spec_19,mean_NMB_sens_90_spec_58,mean_NMB_sens_88_spec_63,mean_NMB_sens_75_spec_66,mean_NMB_sens_64_spec_76,mean_NMB_sens_57_spec_80,mean_NMB_sens_28_spec_89)
 #optimal strategy, sens_28_spec_89
 
 #what is the next best strategy
 max(mean_NMB_sens_100_spec_3,mean_NMB_sens_95_spec_19,mean_NMB_sens_90_spec_58,mean_NMB_sens_88_spec_63,mean_NMB_sens_75_spec_66,mean_NMB_sens_64_spec_76,mean_NMB_sens_57_spec_80)
 #sens 57, spec 80
 
 diff_NMB <- NMB_sens_28_spec_89 - NMB_sens_57_spec_80
 mean(diff_NMB)
 sd(diff_NMB)/sqrt(20000)
 #number of SDs away from no effectt
 mean(diff_NMB)/(sd(diff_NMB)/sqrt(20000))
 #definately stable at 20k pats, 0 is 2.0 SEs away from the mean