###----
#Note for reviewers:
#Below is the pipeline we used to extract light data from the monthly survey data
#A similar pipeline (different input files) was used to extract biweekly light data
#In our paper we only considered total irradiance (auc) and total uv (uv.auc)
#But this pipeline also extracts peak wavelength & maximum spectral energy (data not included in submission)
#Instructions and notes are for future grad students/undergrads in the lab
####----

#Instructions:
#Below is the pipeline I developed to extract  light data from many *.JazIrrad files at one time
#Code for reading & extraction of data is adapted from:
## Aphalo, Pedro J. (2015) The r4photobiology suite. UV4Plants Bulletin, 2015:1, 21-29. DOI:10.19232/uv4pb.2015.1.14

#This pipeline is designed for those with lots of *.JazIrrad files that contain little to no useful metadata
#For example: our lab takes measurements in the field so locations are recorded in a notebook then later entered into a *.csv file
#For those who have useful metadata in your *.JazIrrad files, I encourage you to look at the citation above for other methods to save time & energy

###Versions:
#Use R.3.6.3 and the PhotobiologyInOut library released in 2020
#Newer versions of R or the library may result in erroneous results
#If in doubt double check the contents of a jaz.spct object and compare to the *.JazIrrad File  
#The s.e.irrad column should match the 4th column containing calibrated data (P) of *.JazIrrad Files
#E_irrad values should be similar to integrals calculated in SpectraSuite 


Sys.setenv(TZ = 'UTC')
library(photobiology)
library(photobiologyInOut)
library(pavo)
library(readr)
library(data.table)
library(stringr)

options(tibble.print_max = 5,
        tibble.print_min = 3,
        photobiology.strict.range = NA_integer_)

library(dplyr)
library(tidyr)


###File Organization
## In order to delineate between different sampling times I sorted files into folders 
## Each folder represents a different sampling event, usually containing 3 files to be averaged 
##Whatever naming convention you use, it must be consistent for all folders

#### Metadata
## In a *.csv file record all meta data with the information you used to name folders & files
## Each folder/sampling event should have its own line
## You could decide to have just two columns for file names and folder numbers
## Here, I separated meta data into different columns for ease of later processing
## One column should contain the "xx" designaters for all  files found in the folder
## Any character can be used to separate file names in the File Names Column: eg. 77.78.79 or 77\78\79
## Do not use spaces to separate file names
## Example of entry (spaces indicate different columns): 2021Nov05 Woolsey 4 98.99.00


#read in metadata
temp <- read.csv("~/Research/Community Data Paper/IRR_2017_2021.csv")

colnames(temp)

#Here is where you describe folder names so the code can find the files
#You must place your folders in the photobiologyInOut package folder called extdata!
#In this example, I made a new folder called BGO
#Folders where labeled with method: "Date Site IRR.Point"
#I find putting date first makes it easier to scroll through and find folders later
#(My PC arranges them numerically)
#Example Folder Name: "2021Nov05 Woolsey 4"

temp$folder.names<-paste("extdata/BGO/", temp$Date, " ", temp$Site, " ",temp$IRR.Point, sep="")
temp$folder.names

#split the file category
#Make sure you take into account the max number of files you will have per sampling event
#I overwhelming had 3 files per sampling event, but the max I had was 6
#this will separate everything you have in the file name column into different columns
#Each file name has its own column

temp$File.Name<-as.character(temp$File.Name)
temp<-temp %>% separate(File.Name, c("File1", "File2","File3","File4","File5","File6"))
tempdf<-as.data.frame(temp)

#Function to give files names
#This uses the default naming convention
give.file.name <- function(x) {paste("OUTPUTFILE00",x,".JazIrrad", sep="")}

#name the files
tempdf$File1<- give.file.name(tempdf$File1)
tempdf$File2<- give.file.name(tempdf$File2)
tempdf$File3<- give.file.name(tempdf$File3)
tempdf$File4<- give.file.name(tempdf$File4)
tempdf$File5<- give.file.name(tempdf$File5)
tempdf$File6<- give.file.name(tempdf$File6)

#remove all "OUTPUTFILE00NA.JazIrrad"
tempdf[tempdf == "OUTPUTFILE00NA.JazIrrad"] <- NA

#Get ready to find files in folders
spec.df<-data.frame(tempdf$File1,tempdf$File2,tempdf$File3,tempdf$File4, tempdf$File5, tempdf$File6, tempdf$folder.names)
spec.df<-spec.df%>% distinct()
spec.df <-spec.df %>% rename(File1 = 1)
spec.df <-spec.df %>% rename(File2 = 2)
spec.df <-spec.df %>% rename(File3 = 3)
spec.df <-spec.df %>% rename(File4 = 4)
spec.df <-spec.df %>% rename(File5 = 5)
spec.df <-spec.df %>% rename(File6 = 6)
spec.df <-spec.df %>% rename(folder.names = 7)

#For each file# run the loop below
#jaz.spcttrim <- trim_wl(jaz.spct, range = c(290, 800)) restricts wavelengths to visible & IR spectrum
#spec.df[i,8]<-e_irrad(jaz.spcttrim) finds the area under the curve for visible & IR spectrum
#area under the curve = total irradiance
#Total irradiance (auce) units are: W/m^2/nm
#spec.df[i,9]<-jaz.spcttrim[which.max(jaz.spcttrim$s.e.irrad),1] finds the wavelength with the highest IRR value
#Wavelengths are in nm
# jaz.uv <- trim_wl(jaz.spct, range = c(290, 400)) extracts all wavelengths in the visible uv range
#spec.df[i,10]<-max(jaz.uv$s.e.irrad) finds the maximum value of spectral energy in the UV range
#spec.df[i,11]<-e_irrad(jaz.uv) finds the area under the curve for the uv spectrum
#spec.df[i,12]<-jaz.uv[which.max(jaz.uv$s.e.irrad),1] finds the UV wavelength with the highest IRR value

###Troubleshooting:
#if the loop stops prematurely:
#Double check metadata was entered correctly
#Double check the folder exists and contains the right file numbers
#Make sure folder names do not have a typo

#file1
for (i in seq_len(nrow(spec.df))) {
  if (!is.na(spec.df[i,1])){
    jaz.s.irrad.file <- 
      system.file(spec.df[i,7], spec.df[i,1], 
                  package = "photobiologyInOut", mustWork = TRUE)
    jaz.spct <- read_oo_jazirrad(file = jaz.s.irrad.file)
    jaz.spcttrim <- trim_wl(jaz.spct, range = c(290, 800))
    spec.df[i,8]<-e_irrad(jaz.spcttrim)
    spec.df[i,9]<-jaz.spcttrim[which.max(jaz.spcttrim$s.e.irrad),1]
 #uv range   
    jaz.uv <- trim_wl(jaz.spct, range = c(290, 400))
    spec.df[i,10]<-max(jaz.uv$s.e.irrad)
    spec.df[i,11]<-e_irrad(jaz.uv)
    spec.df[i,12]<-jaz.uv[which.max(jaz.uv$s.e.irrad),1]
   #max uv
  }else
    spec.df[i,1] <- NA
}
# CODE TO RENAME COLUMNS: FILE 1 IS RENAMED
spec.df <-spec.df %>% rename(auc1 = 8)
spec.df <-spec.df %>% rename(peak.wl1 = 9)
spec.df <-spec.df %>% rename(uv.max.se1 = 10)
spec.df <-spec.df %>% rename(uv.auc1 = 11)
spec.df <-spec.df %>% rename(uv.peak1 = 12)

#file2
for (i in seq_len(nrow(spec.df))) {
  if (!is.na(spec.df[i,2])){
    jaz.s.irrad.file <- 
      system.file(spec.df[i,7], spec.df[i,2], 
                  package = "photobiologyInOut", mustWork = TRUE)
    jaz.spct <- read_oo_jazirrad(file = jaz.s.irrad.file)
    jaz.spcttrim <- trim_wl(jaz.spct, range = c(290, 800))
    spec.df[i,13]<-e_irrad(jaz.spcttrim)
    spec.df[i,14]<-jaz.spcttrim[which.max(jaz.spcttrim$s.e.irrad),1]
    #uv
    jaz.uv <- trim_wl(jaz.spct, range = c(290, 400))
    spec.df[i,15]<-max(jaz.uv$s.e.irrad)
    spec.df[i,16]<-e_irrad(jaz.uv)
    spec.df[i,17]<-jaz.uv[which.max(jaz.uv$s.e.irrad),1]
  }else
    spec.df[i,2] <- NA
}
#RENAME COLUMNS
spec.df <-spec.df %>% rename(auc2 = 13)
spec.df <-spec.df %>% rename(peak.wl2 = 14)
spec.df <-spec.df %>% rename(uv.max.se2 = 15)
spec.df <-spec.df %>% rename(uv.auc2 = 16)
spec.df <-spec.df %>% rename(uv.peak2 = 17)

#file3
for (i in seq_len(nrow(spec.df))) {
  if (!is.na(spec.df[i,3])){
    jaz.s.irrad.file <- 
      system.file(spec.df[i,7], spec.df[i,3], 
                  package = "photobiologyInOut", mustWork = TRUE)
    jaz.spct <- read_oo_jazirrad(file = jaz.s.irrad.file)
    jaz.spcttrim <- trim_wl(jaz.spct, range = c(290, 800))
    spec.df[i,18]<-e_irrad(jaz.spcttrim)
    spec.df[i,19]<-jaz.spcttrim[which.max(jaz.spcttrim$s.e.irrad),1]
    #uv
    jaz.uv <- trim_wl(jaz.spct, range = c(290, 400))
    spec.df[i,20]<-max(jaz.uv$s.e.irrad)
    spec.df[i,21]<-e_irrad(jaz.uv)
    spec.df[i,22]<-jaz.uv[which.max(jaz.uv$s.e.irrad),1]
  }else
    spec.df[i,3] <- NA
}
spec.df <-spec.df %>% rename(auc3 = 18)
spec.df <-spec.df %>% rename(peak.wl3 = 19)
spec.df <-spec.df %>% rename(uv.max.se3 = 20)
spec.df <-spec.df %>% rename(uv.auc3 = 21)
spec.df <-spec.df %>% rename(uv.peak3 = 22)

#file4
for (i in seq_len(nrow(spec.df))) {
  if (!is.na(spec.df[i,4])){
    jaz.s.irrad.file <- 
      system.file(spec.df[i,7], spec.df[i,4], 
                  package = "photobiologyInOut", mustWork = TRUE)
    jaz.spct <- read_oo_jazirrad(file = jaz.s.irrad.file)
    jaz.spcttrim <- trim_wl(jaz.spct, range = c(290, 800))
    spec.df[i,23]<-e_irrad(jaz.spcttrim)
    spec.df[i,24]<-jaz.spcttrim[which.max(jaz.spcttrim$s.e.irrad),1]
    #uv
    jaz.uv <- trim_wl(jaz.spct, range = c(290, 400))
    spec.df[i,25]<-max(jaz.uv$s.e.irrad)
    spec.df[i,26]<-e_irrad(jaz.uv)
    spec.df[i,27]<-jaz.uv[which.max(jaz.uv$s.e.irrad),1]
  }else
    spec.df[i,4] <- NA
}
#rename columns
spec.df <-spec.df %>% rename(auc4 = 23)
spec.df <-spec.df %>% rename(peak.wl4 = 24)
spec.df <-spec.df %>% rename(uv.max.se4 = 25)
spec.df <-spec.df %>% rename(uv.auc4 = 26)
spec.df <-spec.df %>% rename(uv.peak4 = 27)

#file5
for (i in seq_len(nrow(spec.df))) {
  if (!is.na(spec.df[i,5])){
    jaz.s.irrad.file <- 
      system.file(spec.df[i,7], spec.df[i,5], 
                  package = "photobiologyInOut", mustWork = TRUE)
    jaz.spct <- read_oo_jazirrad(file = jaz.s.irrad.file)
    jaz.spcttrim <- trim_wl(jaz.spct, range = c(290, 800))
    spec.df[i,28]<-e_irrad(jaz.spcttrim)
    spec.df[i,29]<-jaz.spcttrim[which.max(jaz.spcttrim$s.e.irrad),1]
    #uv
    jaz.uv <- trim_wl(jaz.spct, range = c(290, 400))
    spec.df[i,30]<-max(jaz.uv$s.e.irrad)
    spec.df[i,31]<-e_irrad(jaz.uv)
    spec.df[i,32]<-jaz.uv[which.max(jaz.uv$s.e.irrad),1]
    }else
    spec.df[i,5] <- NA
}
#rename columns
spec.df <-spec.df %>% rename(auc5 = 28)
spec.df <-spec.df %>% rename(peak.wl5 = 29)
spec.df <-spec.df %>% rename(uv.max.se5 = 30)
spec.df <-spec.df %>% rename(uv.auc5 = 31)
spec.df <-spec.df %>% rename(uv.peak5 = 32)

#file6
for (i in seq_len(nrow(spec.df))) {
  if (!is.na(spec.df[i,6])){
    jaz.s.irrad.file <- 
      system.file(spec.df[i,7], spec.df[i,6], 
                  package = "photobiologyInOut", mustWork = TRUE)
    jaz.spct <- read_oo_jazirrad(file = jaz.s.irrad.file)
    jaz.spcttrim <- trim_wl(jaz.spct, range = c(290, 800))
    spec.df[i,33]<-e_irrad(jaz.spcttrim)
    spec.df[i,34]<-jaz.spcttrim[which.max(jaz.spcttrim$s.e.irrad),1]
    #uv
    jaz.uv <- trim_wl(jaz.spct, range = c(290, 400))
    spec.df[i,35]<-max(jaz.uv$s.e.irrad)
    spec.df[i,36]<-e_irrad(jaz.uv)
    spec.df[i,37]<-jaz.uv[which.max(jaz.uv$s.e.irrad),1]
  }else
    spec.df[i,6] <- NA
}
spec.df
###RENAME 
spec.df <-spec.df %>% rename(auc6 = 33)
spec.df <-spec.df %>% rename(peak.wl6 = 34)
spec.df <-spec.df %>% rename(uv.max.se6 = 35)
spec.df <-spec.df %>% rename(uv.auc6 = 36)
spec.df <-spec.df %>% rename(uv.peak6 = 37)

#means of files 1 -5
spec.df$mean.auc <- rowMeans(spec.df[,c(8,13,18,23,28)], na.rm = TRUE)
spec.df$mean.peak.wl <- rowMeans(spec.df[,c(9,14,19,24,29)], na.rm = TRUE)
spec.df$mean.uv.max.se <- rowMeans(spec.df[,c(10,15,20,25,30)], na.rm = TRUE)
spec.df$mean.uv.auc <- rowMeans(spec.df[,c(11,16,21,26,31)], na.rm = TRUE)
spec.df$mean.uv.peak <- rowMeans(spec.df[,c(12,17,22,27,32)], na.rm = TRUE)
#joins site data with means
#add to end
#merge metadata with means
spec.all <- merge(tempdf, spec.df, by=c("File1", "File2","File3","File4","File5","File6","folder.names"))
#write csv
#change file name for you needs
write.csv(spec.all,"~/Research/BGO_spec2.csv")

