#This scripts gives some examples on how to perform a mixed model analysis in R 
#Download the packages required for this Lab (unblock the comments by deleting #):
#install.packages("rrBLUP")
#install.packages("EMMREML")
#install.packages("MCMCglmm")
#install.packages("pedigree")
#install.packages("nadiv")

#Import the data
getwd()
setwd("C:/Users/Badger/Desktop/Labs SISG modules 9 and 12")
ysim = read.table("ysim.txt", h = T)
#Build a pedigree file for our data
ped = ysim[,1:3]
head(ped)


#---------------------------------------------------------------------------------#
# option 1: using our function programmed from scratch
#---------------------------------------------------------------------------------#
#import the function
source("C:\\Users\\Badger\\Desktop\\Labs SISG modules 9 and 12\\EM_algorithm.R")

#For using this function we will have compute the A (additive relationship matrix) inverse 
#besides the incidence matrices for the fixed effects (X) and animal effect (Z)

#A simple function for creating incidence matrices
incidence = function(x){
  temp = data.frame(table(x))
  temp = subset(temp, Freq >0)
  n_class = dim(temp)[1]
  inc = matrix(0, length(x), n_class)
  colnames(inc) = temp[,1]
  for (i in 1:length(x)){
    j = which(x[i]==temp[,1])
    inc[i,j] = 1
  }
  return(inc)
}

#---------------------------------------------------------------------------#
#-------------- Create the A inverse ---------------------------------------#
#---------------------------------------------------------------------------#
library(pedigree)
Ai = makeAinv(ped=ysim[,1:3])
Ai <- read.table("Ainv.txt")
nInd <- nrow(ysim[,1:3])
Ainv <- matrix(0,nrow = nInd,ncol = nInd)
Ainv[as.matrix(Ai[,1:2])] <- Ai[,3]
dd <- diag(Ainv)
Ainv <- Ainv + t(Ainv)
diag(Ainv) <- dd
dim(Ainv)
colnames(Ainv) = rownames(Ainv) <- ysim$animal

ysim2 = subset(ysim, y>0) #Maintain only animals with phenotype

#Building the X matrix
X=cbind(incidence(ysim2$sex),ysim2$feed)
#See how it looks like
head(X)

#Building the Z matrix
Z=incidence(ysim$animal)
Z = Z[-which(ysim$y==0),] #delete the animals without information for Z
head(Z[1:10,1:10])
head(Z[1:10,121:131])
dim(Z)
#check if A and Z are conformable
dim(Z) == dim(Ainv)

#Getting an idea for starting points within the parameter space
summary(lm(y~as.factor(sex)+feed,data=ysim[which(ysim$y>0),]))
guess1 = (28.22^2)/2
guess2 = guess1

#Perform the analysis
#G and R will be the genetic and residual variances for the model after convergence
mme = emreml610(maxrounds = 1000,guess = c(guess1,guess1),criter=10^-5,Ainv=Ainv,X=X,Z=Z,y=ysim2$y)

#Extract the solutions
blue = mme$blue
blup = mme$blup
mme$G
mme$R
mme$G/(mme$G+mme$R) #Compute the estimated heritability

#Check the correlation with true breeding values
par(mfrow =c(1,1))
plot(blup$ebv, ysim$tbv, xlab = "EBV", ylab= "TBV", main = "", pch = 16, col = rgb(0.5,0.5,0.9,0.6))
cor(blup$ebv, ysim$tbv)
mean(blup$acc)

#---------------------------------------------------------------------------------#
# Option 2: Using EMMREML for Variance Component Estimation (VCE) using the REML method 
# and the EMMA algorithm
#---------------------------------------------------------------------------------#
library(EMMREML)
#This package actually needs A not its inverse, so we are back solving Ainv
mme2 = emmreml(y=ysim2$y, X=X, Z=Z, K=solve(Ainv),varbetahat=T,varuhat=T, PEVuhat=T, test=FALSE)
mme2$Vu
mme2$Ve
mme2$Vu/(mme2$Vu+mme2$Ve) #Compute the estimated heritability
blue2 = mme2$betahat
blup2 = mme2$uhat

#---------------------------------------------------------------------------------
# Option 3: Using rrBlup for VCE with the Maximum Likelihood (ML) method
#---------------------------------------------------------------------------------
library(rrBLUP)
mme3 = mixed.solve(y=ysim2$y, Z=Z, K=solve(Ainv), X=X, method="ML",
                   bounds=c(1e-09, 1e+09), SE=T, return.Hinv=FALSE)

mme3$Vu
mme3$Ve
mme3$Vu/(mme3$Vu+mme3$Ve) #Compute the estimated heritability

blup3 = mme3$u
blue3 = mme3$beta


#---------------------------------------------------------------------------------
# Option 4: Using MCMCglmm for VCE via Bayesian Inference
#---------------------------------------------------------------------------------

#Define a flat prior for the variance components
flat_prior <- list(R = list(V = 1, nu = 0), G = list(G1 = list(V = 1, nu = 0)))

library(MCMCglmm)
library(nadiv)

head(ysim)
ysim = subset(ysim, y!=0) #Delete the missing data 
head(ysim)

ped1 = prepPed(ped)#Prepare the pedigree data
ysim$sex1 =  ifelse(ysim$sex==1,"M","F")#create factor levels for the sex 
mme4 = MCMCglmm(y~sex1+feed, random =~ animal, data = ysim, pedigree = ped1, verbose=T,
                nitt=20000, thin=5, burnin=10000, pr=T, family = "gaussian", prior =flat_prior)

#Estimate heritability
h2 = mme4$VCV[,1]/(mme4$VCV[,1]+mme4$VCV[,2])

#Visualize the heritability a posteriori estimates
par(mfrow=c(1,2))
plot(as.numeric(h2), type = "l", ylab = "h2 values", main = "Trace plot")
plot(density(h2), main = "h2 marginal distribution")
polygon(density(h2), col = rgb(0.2,0.4,0.8,0.6))
abline(v = mean(h2), col = "red", lty = 2)

#Visualize the solutions
head(mme4$Sol[,1:10])

#Create a matrix for storing the blups
blup4 = data.frame("animal" = rep(0,nrow(ped1)),
                  "ebv" = rep(0,nrow(ped1)), "sd" = rep(0,nrow(ped1)))

#Average the solutions for the blup
for(i in 1:(ncol(mme4$Sol)-3)){
  blup4[i,1] = ped$animal[i]
  blup4[i,2] = mean(mme4$Sol[,i+3])
  blup4[i,3] = sd(mme4$Sol[,i+3])
}

head(blup4)
head(blup3)
#---------------------------------------------------------------------------------#
# Comparing all solutions
#---------------------------------------------------------------------------------#
results = data.frame(blup1 = blup$ebv, blup2 = blup2, blup3 = blup3, blup4 = blup4$ebv)
cor(results)

