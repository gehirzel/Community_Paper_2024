
# Title: "PCA & GLMS for for Community Data Paper"
#author: "Grace Hirzel"
#date: "2023-12-08"
#This script was run with R 4.2.2

#### Libraries ####

library (factoextra)
library(stringr)
library(dplyr)
library(reshape2)
library(ggplot2)

library(ggpmisc)
library(ggpubr) 
library(ggthemes)

#these libraries not strictly necessary for analysis but improve data interpretation
library(corrplot)


#Monthly Surveys #### 


#* Light Data ####

IRM <- read.csv("~/Monthly_Light.csv")
colnames(IRM)
IRM$site.point<-paste(IRM$Site, " ", IRM$IRR.Point, sep="")

#Rename by Site
IRM$site.point<-str_replace(IRM$site.point, "Woolsey 1", "Woolsey")
IRM$site.point<-str_replace(IRM$site.point, "Woolsey 2", "Woolsey")
IRM$site.point<-str_replace(IRM$site.point, "Woolsey 3", "Woolsey")
IRM$site.point<-str_replace(IRM$site.point, "Woolsey 4", "Woolsey")
IRM$site.point<-str_replace(IRM$site.point, "Woolsey 5", "Woolsey")
colnames(IRM)

#calculate relative UV
IRM$relative.uv <- IRM$mean.uv.auc/IRM$mean.auc

#extract month info
IRM$month<- substr(IRM$Date,5,7)
IRM$month <- factor(IRM$month, levels = month.abb)
summary(IRM)

#sort by BGO sunny, shade, & butterfly house sites
IRM$site.simple<-IRM$site.point
IRM$site.simple<-str_replace(IRM$site.simple, "BGO 1", "BGO Sunny")
IRM$site.simple<-str_replace(IRM$site.simple, "BGO 2", "BGO Sunny")
IRM$site.simple<-str_replace(IRM$site.simple, "BGO 3", "Butterfly_House")
IRM$site.simple<-str_replace(IRM$site.simple, "BGO 4", "BGO Shade")
IRM$site.simple<-str_replace(IRM$site.simple, "BGO 5", "BGO Shade")

#Remove Butterfly House Data (not part of study)
IR.Final<-filter(IRM, site.simple != "Butterfly_House")

#match months with abbreviations
IR.Final$abb<-match(IR.Final$month, month.abb)


#*Combine Count & Light Data ######

#create survey ID to match with abundance data
IR.Final$Survey.ID<-paste(IR.Final$Date, IR.Final$Site, sep = " ")

#find mean total irradiance for each survey visit
mo.light <- aggregate(IR.Final$mean.auc, 
                      by = list(IR.Final$Survey.ID), FUN = mean, na.rm=F)
#find mean relative UV for each survey visit
mo.uv <- aggregate(IR.Final$relative.uv, 
                   by = list(IR.Final$Survey.ID), FUN = mean, na.rm=F)

#rename columns
names(mo.light)[1] <- "Survey.ID"
names(mo.light)[2] <- "total_irrad"
names(mo.uv)[1] <- "Survey.ID"
names(mo.uv)[2] <- "rel.uv"

#join average total irradiance & relative uv dataframes
mo.most2<-left_join(mo.light,mo.uv, by=c("Survey.ID"))

#Get abundance & weather data
mo.counts <- read.csv("~/Monthly_Abundance.csv")


#Create Survey ID
mo.counts$Site<-str_replace(mo.counts$Site, "Botanical_Garden", "BGO")
mo.counts$Survey.ID<-paste(mo.counts$Date, mo.counts$Site, sep = " ")
colnames(mo.counts)
#tidy the data
mo.melt = melt(mo.counts, id.vars = c("Survey.ID", "Site", "Date", "Month", "Year", "total_butterflies", "TempF", "Weekly_Precip", "Weekly_Precip_cm"),
               measure.vars = c("nymphalidae", "pieridae", "papilionidae", "hesperiidae", "lycaenidae"))
names(mo.melt)[10] <- "Family"
names(mo.melt)[11] <- "Count"
mo.melt$Family<-str_replace(mo.melt$Family, "nymphalidae", "Nymphalidae")
mo.melt$Family<-str_replace(mo.melt$Family, "pieridae", "Pieridae")
mo.melt$Family<-str_replace(mo.melt$Family, "hesperiidae", "Hesperiidae")
mo.melt$Family<-str_replace(mo.melt$Family, "papilionidae", "Papilionidae")
mo.melt$Family<-str_replace(mo.melt$Family, "lycaenidae", "Lycaenidae")
summary(mo.light)
summary(mo.light$Survey.ID)

#join count, weather & light data
mo<-left_join(mo.melt, mo.most2, by=c("Survey.ID"))

#add numerical months
mo$abb<-match(mo$Month, month.name)
#Convert to Celsius
mo$TempC<-(mo$TempF-32)/ 1.8 


#*Correlation Matrix ####

#convert site(factor) into number to plug into correlation matrix

#convert site(factor) into number to plug into correlation matrix

mo$site.num<-mo$Site 
mo$site.num<-str_replace(mo$site.num, "Woolsey", "2")
mo$site.num<-str_replace(mo$site.num, "BGO", "1")
mo$site.num<-as.numeric(mo$site.num)
colnames(mo)

#include Year, Weekly Precipitation (cm) , Total Irradiance, Relative UV, Abbreviated Month (abb), Temperature (C), & Site
mo.matx <- cor(na.omit(mo[ , c(5,9,12,13,14,15,16)]))

#Nice optional visualization
corrplot(mo.matx, method="pie")

mo.matx.rounded <- round(x= mo.matx, digits=3)
View(mo.matx.rounded)


#*PCA ####

#Runs the Principal Components analysis and shows how much each component explains variance in the data

#Site was left out of PCA (no column 17) for later analysis of effect of site
mo.pca <- princomp(na.omit(mo[ , c(5,9,12,13,14,15)]),
                   scores = TRUE,
                   cor = TRUE
)
#proportion of variance table
summary(mo.pca)


#Details the loadings for each Principal Component ("abb" represents survey month as a numeric14
loadings.mo <- unclass(mo.pca[["loadings"]])
sweep (loadings.mo,
       MARGIN = 2,
       STATS = colSums(loadings.mo),
       FUN = "/")

#Biplot & Scree Plots for monthly data PCA

#visualize Principal Components

library("factoextra")
summary(mo.pca)

screeplot(x = mo.pca,
          npcs = length(mo.pca[["sdev"]]),
          type = "lines")

fviz_pca_biplot(mo.pca, 
                col.var = "red",
                label = "var") +
  theme_classic()

#Preparing the data for generalized linear models:
  
#Bind together PCs of interest
moPC1 <- data.frame(mo.pca[["scores"]])[,1]
moPC2 <- data.frame(mo.pca[["scores"]])[,2]
moPC3<- data.frame(mo.pca[["scores"]])[,3]
moscores1<-cbind(moPC1, moPC2)
moscores<-cbind(moscores1, moPC3)

#Bind abundance data with PC scores
mo1<-na.omit(mo)
mo.PCs <- cbind(mo1, moscores)


#*GLM Visualizations ####

##Prepare the data to create faceted figure
melt.mo <- melt(mo.PCs, id.vars = c("Site", "Date", "Year", "Month", "Survey.ID","total_butterflies"),
                measure.vars = c("moPC1", "moPC2", "moPC3"))

melt.mo.total <- melt(mo.PCs, id.vars = c("Site", "Date", "Year", "Month", "Survey.ID","total_butterflies"),
                      measure.vars = c("moPC1", "moPC2", "moPC3"))
melt.mo.total$Family <- "All Families"
names(melt.mo.total)[6] <- "Count"

melt.mo.fam <- melt(mo.PCs, id.vars = c("Site", "Date", "Year", "Month", "Survey.ID","Family","Count"),
                    measure.vars = c("moPC1", "moPC2", "moPC3"))
mo.facet<-rbind(melt.mo.fam,melt.mo.total)
mo.facet$Family <- factor(mo.facet$Family , levels=c("All Families", "Hesperiidae", "Lycaenidae", "Nymphalidae", "Papilionidae", "Pieridae"))

mo.facet2<-mo.facet %>%
  distinct(.keep_all = TRUE)


##Visualize: Family vs PC linear regressions
# R^2 values appear on graphs
q<-ggplot(data = mo.facet2, aes(x = value, y = Count, color=Site)) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point(size=3) +
  theme_classic() +
  annotate(
    geom = 'segment',
    y = -Inf,
    yend = Inf,
    x = Inf,
    xend = Inf
  ) + 
  scale_color_manual(values = c('#ffc000', '#0072b2')) +
  theme(strip.placement = "outside")
facet(q, facet.by = c("Family","variable"),
      scales = "free",
      panel.labs.background = list(size = 0.1),
      panel.labs.font = list(color = "black"),
      panel.labs.font.x = list( color = "black"))

#Visualize: Family vs Site boxplot

ggplot(mo.facet2,aes(Family, Count, fill=Site)) +
  stat_boxplot(size=2,outlier.size=5) + 
  theme_classic() +  
  theme(legend.title = element_blank(),legend.text=element_text(size=25), 
        axis.text=element_text(size=10, angle = 45, hjust=1, face="bold"),
        axis.text.x=element_text(margin = margin(t=0, b=0)),
        axis.title.y = element_text(size=28,face="bold",
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(size=28,face="bold",
                                    margin = margin(t = 20, r = 0, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size=32, face="bold")) +
  scale_fill_manual(values = c('#ffc000', '#0072b2')) +
  labs(x="Family", y= "# of Butterflies")
#* GLMs ####

#filter for each site
mo.BGO<-filter(mo.PCs,Site=="BGO")
mo.Woolsey<-filter(mo.PCs,Site=="Woolsey")
summary(mo.BGO)
summary(mo.Woolsey)

#filter for each family
mo.hes<-filter(mo.PCs,Family=="Hesperiidae")
mo.pap<-filter(mo.PCs,Family=="Papilionidae")
mo.nym<-filter(mo.PCs,Family=="Nymphalidae")
mo.pie<-filter(mo.PCs,Family=="Pieridae")
mo.lyc<-filter(mo.PCs,Family=="Lycaenidae")

#filter for each family BGO
BGO.hes<-filter(mo.BGO,Family=="Hesperiidae")
BGO.pap<-filter(mo.BGO,Family=="Papilionidae")
BGO.nym<-filter(mo.BGO,Family=="Nymphalidae")
BGO.pie<-filter(mo.BGO,Family=="Pieridae")
BGO.lyc<-filter(mo.BGO,Family=="Lycaenidae")

#filter for each family Woolsey
Wol.hes<-filter(mo.Woolsey,Family=="Hesperiidae")
Wol.pap<-filter(mo.Woolsey,Family=="Papilionidae")
Wol.nym<-filter(mo.Woolsey,Family=="Nymphalidae")
Wol.pie<-filter(mo.Woolsey,Family=="Pieridae")
Wol.lyc<-filter(mo.Woolsey,Family=="Lycaenidae")


mo.Pcsd<-mo.PCs %>%
  distinct(Survey.ID, .keep_all = TRUE)
BGO.Pcsd<-mo.BGO %>%
  distinct(Survey.ID, .keep_all = TRUE)
Wol.Pcsd<-mo.Woolsey %>%
  distinct(Survey.ID, .keep_all = TRUE)
summary(BGO.Pcsd)
summary(Wol.Pcsd)
#GLM for all butterflies  
all.glm <- glm(total_butterflies ~ moPC1+moPC2+moPC3 + Site, data=mo.Pcsd, family=quasipoisson)
summary(all.glm)
anova(all.glm, test = "F")

mo.Pcsd$YearF<-as.factor(mo.Pcsd$Year)
Year.aov<-aov(total_butterflies ~ YearF, data = mo.Pcsd)
summary(Year.aov)

BGO.glm <- glm(total_butterflies ~ moPC1+moPC2+moPC3, data=BGO.Pcsd, family=quasipoisson)
summary(BGO.glm)
anova(BGO.glm, test = "F")


BGOh.glm <- glm(Count ~ moPC1+moPC2+moPC3, data=BGO.hes, family=quasipoisson)
summary(BGOh.glm)
anova(BGOh.glm, test = "F")

BGOl.glm <- glm(Count ~ moPC1+moPC2+moPC3, data=BGO.lyc, family=quasipoisson)
summary(BGOl.glm)
anova(BGOl.glm, test = "F")

BGOn.glm <- glm(Count ~ moPC1+moPC2+moPC3, data=BGO.nym, family=quasipoisson)
anova(BGOn.glm, test = "F")

BGOp.glm <- glm(Count ~ moPC1+moPC2+moPC3, data=BGO.pap, family=quasipoisson)
anova(BGOp.glm, test = "F")

BGOpi.glm <- glm(Count ~ moPC1+moPC2+moPC3, data=BGO.pie, family=quasipoisson)
anova(BGOpi.glm, test = "F")

Wol.glm <- glm(total_butterflies ~ moPC1+moPC2+moPC3, data=Wol.Pcsd, family=quasipoisson)
summary(Wol.glm)
anova(Wol.glm, test = "F")

Wolh.glm <- glm(Count ~ moPC1+moPC2+moPC3, data=Wol.hes, family=quasipoisson)
summary(Wolh.glm)
anova(Wolh.glm, test = "F")

Woll.glm <- glm(Count ~ moPC1+moPC2+moPC3, data=Wol.lyc, family=quasipoisson)
summary(Woll.glm)
anova(Woll.glm, test = "F")

Woln.glm <- glm(Count ~ moPC1+moPC2+moPC3, data=Wol.nym, family=quasipoisson)
anova(Woln.glm, test = "F")

Wolp.glm <- glm(Count ~ moPC1+moPC2+moPC3, data=Wol.pap, family=quasipoisson)
anova(Wolp.glm, test = "F")

Wolpi.glm <- glm(Count ~ moPC1+moPC2+moPC3, data=Wol.pie, family=quasipoisson)
anova(Wolpi.glm, test = "F")



#Biweekly Surveys ####
#*Correlation Matrix ####

#####Biweekly Irrad Measurements#######
IR1.2 <- read.csv("Biweekly_LightFinal.csv")
IR1.2$Survey.ID <-str_sub(IR1.2$folder.names,19,-7)
IR1.2$Date <-str_sub(IR1.2$folder.names,19,27)
IR1.2$Irr.ID <-str_sub(IR1.2$folder.names,-1,-1)

#agg by site point
auc.cp <- aggregate(IR1.2$mean.auc,
                    by = list(IR1.2$Survey.ID,IR1.2$Date), FUN = mean, na.rm=FALSE)

IR1.2$relative.uv <- IR1.2$mean.uv.auc/IR1.2$mean.auc
rel.uv.cp <- aggregate(IR1.2$relative.uv,
                       by = list(IR1.2$Survey.ID,IR1.2$Date) , FUN = mean, na.rm=FALSE)

IR.cp<-left_join(auc.cp,rel.uv.cp, by=c("Group.1","Group.2"))
colnames(IR.cp)
IR.cp


names(IR.cp)[1] <- "Survey.ID"
names(IR.cp)[2] <- "Date"
names(IR.cp)[3] <- "auc"
names(IR.cp)[4] <- "rel.uv"


summary(IR.cp)
IR.cp[!duplicated(IR.cp), ]

###########Count & Weather Data######

bi <- read.csv("~/Biweekly_Count_Weather_Final.csv")

bi$Survey.ID<-paste(bi$Date, bi$Site, sep = " ")

bi.melt = melt(bi, id.vars = c("Survey.ID", "Site", "Date", "Time.of.Day", "Year", "Julian_Week", "Total", "Temp", "This_week.s_total_precip_in", "Transect"),
               measure.vars = c("Nymphalidae", "Pieridae", "Papilionidae", "Hesperiidae", "Lycaenidae"))

######Temp Data######

bi$TempC<-(bi$Temp-32)/ 1.8
bi.temp <- aggregate(bi$TempC, 
                     by = list(bi$Survey.ID, bi$Julian_Week) , FUN = mean, na.rm=FALSE)

names(bi.temp)[1] <- "Survey.ID"
names(bi.temp)[2] <- "Week"
names(bi.temp)[3] <- "TempC"

bi.temp$Week.factor<-as.factor(bi.temp$Week)



#####Precip Data######
bi.precip <- aggregate(bi$This_week.s_total_precip_in, 
                       by = list(bi$Year, bi$Julian_Week) , FUN = mean, na.rm=FALSE)


names(bi.precip)[1] <- "Year"
names(bi.precip)[2] <- "Week"
names(bi.precip)[3] <- "Precip"

bi.precip$PrecipCM<-bi.precip$Precip*2.54


#bi.precip$Julian_Week<-as.factor(bi.precip$Julian_Week)
bi.precip$Weekf<-as.factor(bi.precip$Week)

########## Combine Count. Light & Weather######
bi.light <- aggregate(bi.melt$value, 
                      by = list(bi.melt$Year, bi.melt$Julian_Week, bi.melt$Site, bi.melt$Date, bi.melt$variable), FUN = sum, na.rm=FALSE)
bi.light2 <- aggregate(bi.melt$Total, 
                       by = list(bi.melt$Year, bi.melt$Julian_Week, bi.melt$Site, bi.melt$Date, bi.melt$variable), FUN = sum, na.rm=FALSE)

names(bi.light)[1] <- "Year"
names(bi.light)[2] <- "Week"
names(bi.light)[3] <- "Site"
names(bi.light)[4] <- "Date"
names(bi.light)[5] <- "Family"
names(bi.light)[6] <- "Fam.Count"

names(bi.light2)[1] <- "Year"
names(bi.light2)[2] <- "Week"
names(bi.light2)[3] <- "Site"
names(bi.light2)[4] <- "Date"
names(bi.light2)[5] <- "Family"
names(bi.light2)[6] <- "Total.Count"

bi.light$Survey.ID<-paste(bi.light$Date, bi.light$Site, sep = " ")
bi.light2$Survey.ID<-paste(bi.light2$Date, bi.light2$Site, sep = " ")
bi.most<-left_join(bi.light,bi.light2, by=c("Survey.ID","Week","Year", "Site", "Date", "Family"))

bi.temp$Week<-as.numeric(bi.temp$Week)
bi.precip$Week<-as.numeric(bi.precip$Week)
bi.most2<-left_join(bi.most,bi.temp, by=c("Survey.ID","Week"))
bi.most3<-left_join(bi.most2,bi.precip, by=c("Year","Week"))


bi.all<-left_join(bi.most3,IR.cp, by=c("Date"))

summary(bi.all)

#convert site(factor) into number to plug into correlation matrix

bi.all$site.num<-bi.all$Site
bi.all$site.num<-str_replace(bi.all$site.num, "Woolsey", "1")
bi.all$site.num<-str_replace(bi.all$site.num, "Stump", "2")
bi.all$site.num<-str_replace(bi.all$site.num, "Chesney", "3")
bi.all$site.num<-as.numeric(bi.all$site.num)

#include Year, Weekly Precipitation (cm) , Total Irradiance (auc), Relative UV, Survey Week, Temperature (C), & Site
bi.matx <- cor(na.omit(bi.all[ , c(1,2,9,12,15,16,17)]))

#Nice optional visualization
corrplot(bi.matx, method="pie")

bi.matx.rounded <- round(x= bi.matx, digits=3)
View(bi.matx.rounded)
#write.csv(bi.matx.rounded, "week_corr_MayTest.csv", row.names=TRUE)

#*PCA ####
#Site was left out of PCA (no column 21) for later analysis of effect of site
bi.pca <- princomp(na.omit(bi.all[ , c(1,2,9,12,15,16)]),
                   scores = TRUE,
                   cor = TRUE
)
#proportion of variance table
summary(bi.pca)


loadings.bi <- unclass(bi.pca[["loadings"]])
sweep (loadings.bi,
       MARGIN = 2,
       STATS = colSums(loadings.bi),
       FUN = "/")

#Biplot & Scree Plots for biweekly data PCA
#visualize Principal Components

screeplot(x = bi.pca,
          npcs = length(bi.pca[["sdev"]]),
          type = "lines")

fviz_pca_biplot(bi.pca, 
                col.var = "red",
                label = "var") +
  theme_classic()




#Preparing the data for generalized linear models:

#Bind together PCs of interest
biPC1 <- data.frame(bi.pca[["scores"]])[,1]
biPC2 <- data.frame(bi.pca[["scores"]])[,2]
biPC3 <- data.frame(bi.pca[["scores"]])[,3]
biPC4 <- data.frame(bi.pca[["scores"]])[,4]
biscores<-cbind(biPC1, biPC2, biPC3, biPC4)

#Bind abundance data with PC scores
bi1<-na.omit(bi.all)
bi.PCs <- cbind(bi1, biscores)

#*GLM Visualizations ####
##Prepare the data to create faceted figure
melt.bi.total <- melt(bi.PCs, id.vars = c("Site", "Date", "Year", "Week", "Survey.ID.x","Total.Count"),
                      measure.vars = c("biPC1", "biPC2", "biPC3", "biPC4"))
melt.bi.total$Family <- "All Families"
names(melt.bi.total)[6] <- "Fam.Count"

melt.bi.fam <- melt(bi.PCs, id.vars = c("Site", "Date", "Year", "Week", "Survey.ID.x","Family","Fam.Count"),
                    measure.vars = c("biPC1", "biPC2", "biPC3","biPC4"))
bi.facet<-rbind(melt.bi.fam,melt.bi.total)
bi.facet$Family <- factor(bi.facet$Family , levels=c("All Families","Hesperiidae", "Lycaenidae", "Nymphalidae", "Papilionidae", "Pieridae"))

bi.facet2<-bi.facet %>%
  distinct(.keep_all = TRUE)

#Run all sites together

q<-ggplot(data = bi.facet2, aes(x = value, y = Fam.Count)) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point(size=3) +
  theme_base()+
  theme(strip.placement = "outside")

facet(q, facet.by = c("Family","variable"),
      scales = "free",
      panel.labs.background = list(size = 0.5),
      panel.labs.font = list(color = "black"),
      panel.labs.font.x = list( color = "black"))

q<-ggplot(data = bi.facet2, aes(x = value, y = Fam.Count, color = Site)) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point(size=3) +
  theme_base()+
  theme(strip.placement = "outside")
facet(q, facet.by = c("Family","variable"),
      scales = "free",
      panel.labs.background = list(size = 0.5),
      panel.labs.font = list(color = "black"),
      panel.labs.font.x = list( color = "black"))

##Visualize: Family vs PC linear regressions
# R^2 values appear on graphs
#Sites Separated 

q<-ggplot(data = bi.facet2, aes(x = value, y = Fam.Count, color=Site)) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point(size=3) +
  theme_base()+
  scale_color_manual(values = c('#d55e00', '#f0e442', '#0072b2')) +
  theme(strip.placement = "outside")

facet(q, facet.by = c("Family","variable"),
      scales = "free",
      panel.labs.background = list(size = 0.5),
      panel.labs.font = list(color = "black"),
      panel.labs.font.x = list( color = "black"))

#Visualize: Family vs Site boxplot
ggplot(bi.facet2,aes(x=factor(Family, level=c("All Families", "Hesperiidae", "Lycaenidae", "Nymphalidae", "Papilionidae", "Pieridae")),Fam.Count, fill=Site)) +
  stat_boxplot(size=2,outlier.size=5) + 
  theme_classic() +  
  theme(legend.title = element_blank(),legend.text=element_text(size=25), 
        axis.text=element_text(size=28, face="bold"),
        axis.title.y = element_text(size=28,face="bold",
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(size=28,face="bold",
                                    margin = margin(t = 20, r = 0, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size=32, face="bold")) +
  scale_fill_manual(values = c('#d55e00', '#f0e442', '#0072b2')) +
  labs(x="Family", y= "# of Butterflies")

#* GLMs ####
colnames(bi.PCs)[colnames(bi.PCs ) == 'Survey.ID.x'] <- 'Survey.ID'

bi1.hes<-filter(bi.PCs,Family=="Hesperiidae")
bi1.pap<-filter(bi.PCs,Family=="Papilionidae")
bi1.nym<-filter(bi.PCs,Family=="Nymphalidae")
bi1.pie<-filter(bi.PCs,Family=="Pieridae")
bi1.lyc<-filter(bi.PCs,Family=="Lycaenidae")

bi.Pcsd<-bi.PCs %>%
  distinct(Survey.ID, .keep_all = TRUE)
summary(bi.Pcsd)

#GLM for all butterflies  
bi.tot <- glm(Total.Count ~ biPC1 + biPC2 + biPC3 + biPC4 + Site, data=bi.Pcsd, family=quasipoisson)
summary(bi.tot)
anova(bi.tot, test = "F")

#GLM for Hesperiidae
bi.hes.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4 + Site, data=bi1.hes, family=quasipoisson)
summary(bi.hes.glm)
anova(bi.hes.glm, test = "F")

#GLM for Lycaenidae
bi.lyc.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4 + Site, data=bi1.lyc, family=quasipoisson)
summary(bi.lyc.glm)
anova(bi.lyc.glm, test = "F")

#GLM for Nymphalidae 
bi.nym.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4 + Site, data=bi1.nym, family=quasipoisson)
summary(bi.nym.glm)
anova(bi.nym.glm, test = "F")

#GLM Papilionidae
bi.pap.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4 + Site, data=bi1.pap, family=quasipoisson)
summary(bi.pap.glm)
anova(bi.pap.glm, test = "F")

#GLM for Pieridae
bi.pie.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4 + Site, data=bi1.pie, family=quasipoisson)
summary(bi.pie.glm)
anova(bi.pie.glm, test = "F")

#Behavior####

bi.beh <- read.csv("C:/Users/lenny/Cloud-Drive/Community/Data/processed_beh_transects_2024Feb.csv")
#bi.beh <- read.csv(
bi.beh$Survey.ID<-paste(bi.beh$Date, bi.beh$Site, sep = " ")
bi.beh$Year<-as.numeric(bi.beh$Year)

bi.combo<-left_join(bi.beh,bi.PCs, by=c("Year","Site","Family","Week"))
summary(bi.combo$Behavior.Processed)
bi.combo$Behavior.Processed<-as.factor(bi.combo$Behavior.Processed)

#File for JMP: "beh_for_JMP_May15.csv"
melt.beh <- melt(bi.combo, id.vars = c("Site", "Year", "Week", "Survey.ID","Behavior.Processed","Family"),
                 measure.vars = c("biPC1", "biPC2", "biPC3", "biPC4"))

melt.beh$Family <- factor(melt.beh$Family , levels=c("Nymphalidae", "Hesperiidae", "Lycaenidae", "Pieridae", "Papilionidae"))
write.csv(melt.beh, "beh_for_JMP_Feb25.csv")


beh.tot <- glm(Behavior.Processed ~ biPC1 + biPC2 + biPC3 + biPC4 + Site, data=bi.combo, family=binomial)
summary(beh.tot)
anova(beh.tot, test = "Chisq")


beh.hes<-filter(bi.combo,Family=="Hesperiidae")
beh.pap<-filter(bi.combo,Family=="Papilionidae")
beh.nym<-filter(bi.combo,Family=="Nymphalidae")
beh.pie<-filter(bi.combo,Family=="Pieridae")
beh.lyc<-filter(bi.combo,Family=="Lycaenidae")
#Hesperiidae
beh.glm.hes <- glm(Behavior.Processed ~ biPC1 + biPC2 + biPC3 + biPC4 + Site, data=beh.hes, family=binomial)
summary(beh.glm.hes)
anova(beh.glm.hes, test = "Chisq")


#Nymphalidae
beh.glm.nym <- glm(Behavior.Processed ~ biPC1 + biPC2 + biPC3 + biPC4 + Site, data=beh.nym, family=binomial)
summary(beh.glm.nym)
anova(beh.glm.nym, test = "Chisq")


#Lycaenidae
beh.glm.lyc <- glm(Behavior.Processed ~ biPC1 + biPC2 + biPC3 + biPC4 + Site, data=beh.lyc, family=binomial)
summary(beh.glm.lyc)
anova(beh.glm.lyc, test = "Chisq")


###Not enough data
#Pieridae
beh.glm.pie <- glm(Behavior.Processed ~ biPC1 + biPC2 + biPC3 + biPC4 + Site, data=beh.pie, family=binomial)
summary(beh.glm.pie)
anova(beh.glm.pie, test = "Chisq")

#Papilionidae
beh.glm.pap <- glm(Behavior.Processed ~ biPC1 + biPC2 + biPC3 + biPC4 + Site, data=beh.pap, family=binomial)
summary(beh.glm.pap)
anova(beh.glm.pap, test = "Chisq")

###Supporting Analyses####

bi.Pcsd<-bi.PCs %>%
  distinct(Survey.ID, .keep_all = TRUE)


Wol.bi<-filter(bi.PCs,Site=="Woolsey")
Che.bi<-filter(bi.PCs,Site=="Chesney")
Stu.bi<-filter(bi.PCs,Site=="Stump")


Wolb.hes<-filter(Wol.bi,Family=="Hesperiidae")
Wolb.pap<-filter(Wol.bi,Family=="Papilionidae")
Wolb.nym<-filter(Wol.bi,Family=="Nymphalidae")
Wolb.pie<-filter(Wol.bi,Family=="Pieridae")
Wolb.lyc<-filter(Wol.bi,Family=="Lycaenidae")

Cheb.hes<-filter(Che.bi,Family=="Hesperiidae")
Cheb.pap<-filter(Che.bi,Family=="Papilionidae")
Cheb.nym<-filter(Che.bi,Family=="Nymphalidae")
Cheb.pie<-filter(Che.bi,Family=="Pieridae")
Cheb.lyc<-filter(Che.bi,Family=="Lycaenidae")

Stub.hes<-filter(Stu.bi,Family=="Hesperiidae")
Stub.pap<-filter(Stu.bi,Family=="Papilionidae")
Stub.nym<-filter(Stu.bi,Family=="Nymphalidae")
Stub.pie<-filter(Stu.bi,Family=="Pieridae")
Stub.lyc<-filter(Stu.bi,Family=="Lycaenidae")


#GLM for all butterflies  
bi.tot <- glm(Total.Count ~ biPC1 + biPC2 + biPC3 + biPC4, data=Wol.bi, family=quasipoisson)
summary(bi.tot)
anova(bi.tot, test = "F")

bi.tot <- glm(Total.Count ~ biPC1 + biPC2 + biPC3 + biPC4, data=Che.bi, family=quasipoisson)
summary(bi.tot)
anova(bi.tot, test = "F")

bi.tot <- glm(Total.Count ~ biPC1 + biPC2 + biPC3 + biPC4, data=Stu.bi, family=quasipoisson)
summary(bi.tot)
anova(bi.tot, test = "F")

#GLM for Hesperiidae
bi.hes.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4, data=Wolb.hes, family=quasipoisson)
summary(bi.hes.glm)
anova(bi.hes.glm, test = "F")

bi.hes.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4, data=Cheb.hes, family=quasipoisson)
summary(bi.hes.glm)
anova(bi.hes.glm, test = "F")

bi.hes.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4, data=Stub.hes, family=quasipoisson)
summary(bi.hes.glm)
anova(bi.hes.glm, test = "F")

#GLM for Lycaenidae
bi.lyc.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4, data=Wolb.lyc, family=quasipoisson)
summary(bi.lyc.glm)
anova(bi.lyc.glm, test = "F")

bi.lyc.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4, data=Cheb.lyc, family=quasipoisson)
summary(bi.lyc.glm)
anova(bi.lyc.glm, test = "F")

bi.lyc.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4, data=Stub.lyc, family=quasipoisson)
summary(bi.lyc.glm)
anova(bi.lyc.glm, test = "F")


#GLM for Nymphalidae 
bi.nym.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4, data=Wolb.nym, family=quasipoisson)
summary(bi.nym.glm)
anova(bi.nym.glm, test = "F")

bi.nym.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4, data=Cheb.nym, family=quasipoisson)
summary(bi.nym.glm)
anova(bi.nym.glm, test = "F")

bi.nym.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4, data=Stub.nym, family=quasipoisson)
summary(bi.nym.glm)
anova(bi.nym.glm, test = "F")

#GLM Papilionidae
bi.pap.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4, data=Wolb.pap, family=quasipoisson)
summary(bi.pap.glm)
anova(bi.pap.glm, test = "F")

bi.pap.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4, data=Cheb.pap, family=quasipoisson)
summary(bi.pap.glm)
anova(bi.pap.glm, test = "F")

bi.pap.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4, data=Stub.pap, family=quasipoisson)
summary(bi.pap.glm)
anova(bi.pap.glm, test = "F")

#GLM for Pieridae
bi.pie.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4, data=Wolb.pie, family=quasipoisson)
summary(bi.pie.glm)
anova(bi.pie.glm, test = "F")

bi.pie.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4, data=Cheb.pie, family=quasipoisson)
summary(bi.pie.glm)
anova(bi.pie.glm, test = "F")

bi.pie.glm<-glm(Fam.Count~biPC1 + biPC2 + biPC3 + biPC4, data=Stub.pie, family=quasipoisson)
summary(bi.pie.glm)
anova(bi.pie.glm, test = "F")

                                              