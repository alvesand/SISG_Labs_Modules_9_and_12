#This script illustrates how to estimate heritability with the ANOVA
# A small example is given with the Beetle Data (SAS, Saxton)
#install.packages("lme4)

rm(list=ls())
library(lme4)

setwd("C:\\Users\\Badger\\Desktop\\Labs SISG modules 9 and 12")

beetle <- read.table("beetle data.txt", header=TRUE)

str(beetle)
head(beetle)

beetle$sire <- factor(beetle$sire)
beetle$dam <- factor(beetle$dam)

fecund<-beetle[which(beetle$trait=="fec"),]

# ANOVA 
m2<-aov(value ~ 1 + sire + sire:dam, data=fecund)
summary(m2)

attach(fecund)
mytable <- table(sire,dam)

# Approximately n=5 progeny by mating
# Average of q=4.58 dams per sire

# E[MS]
# Sire:      Ve + nVd + nqVs
# Dam/Sire:  Ve + nVd
# Residual:  Ve


###############
# REML estimates

m1 <-lmer(value ~ 1 + (1|sire) + (1|sire:dam), data=fecund)
summary(m1)


# Mass

mass<-beetle[which(beetle$trait=="mass"),]

# ANOVA 
m4<-aov(value ~ 1 + sire + sire:dam, data=mass)
summary(m4)

attach(mass)
mytable <- table(sire,dam)

# Approximately n=5 progeny by mating
# Average of q=4.58 dams per sire

# E[MS]
# Sire:      Ve + nVd + nqVs
# Dam/Sire:  Ve + nVd
# Residual:  Ve


###############
# REML estimates

m3 <-lmer(value ~ 1 + (1|sire) + (1|sire:dam), data=mass)
summary(m3)

