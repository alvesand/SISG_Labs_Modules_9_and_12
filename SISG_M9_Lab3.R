#This script creates a function for simulating the effects of response to selection based on the 
#on a trait influenced by the fixed effects of sex, percentage of protein in the diet and 
#the genetic additive as random effect 
#The simulation model is y = mu + prot% + sex + u + e, with u ~ N(0, A*varu) and e ~ N(0, I*vare)

#Install the necessary packages for this simulation
#install.packages("pedigree")

library(pedigree)

#Run all code from line 10 to line 164 together
#Let's do the simulation within a function to avoid unnecessary data at end of the process
sim_fun = function(h2, i_m, i_f, ngen, mu = mu, vare, seed){
  #-------------- Set the population parameters ----------------------------------#
  set.seed(seed) # Set a seed for reproducibility sake
  mu = mu #Global mean
  sol_sex = c(25, -25) #Sex effects
  nsire = 20 #number of sires on the base population
  slope = .75 #slope for the protein percentage on the diet
  ndam = 100 #number of dams on the base population
  vare = vare #The desired residual variance per generation
  h2 = h2 #Trait heritability
  varu = h2*vare/(1-h2) #Expected Additive variance per generation
  ngen = ngen #number of generations simulated
  i_m = i_m #Selection intensity for males
  i_f = i_f #Selection intensity for females
  protein = c(10.0,12.5,15,17.5,20,22.5,25) #Indicates the levels of protein on the animal's diet
  
  #--- Create empty arrays to store the information for the base population ----#
  sex = matrix(0,nsire+ndam,1) 
  u= matrix(0,nsire+ndam,1)
  ped = matrix(0,nsire+ndam,4)
  y = matrix(0,nsire+ndam,1)
  feed = matrix(0,nsire+ndam) 
  #------------------------------------------------------------------------------#

  #Loop the data for each animal in the base population
  for (j in 1:(nsire+ndam)){
    ped[j,1] = j #Animal ID on the ped file
    u[j] = rnorm(n = 1,mean = 0,sd = sqrt(varu)) #Sample the breeding value from the global additive variance
    if(j<=nsire){sex[j]=1}else{sex[j]=2} #Set the sex of the animal
  }

  #--- Create a data.frame indicating the animals available for mating in the expanded population --#
  #n_offs will count the number of offspring per sire or dam 
  #y is the phenotype (0 or unavailable for the base population)
  #f is the inbreeding coefficient (assumed to be 0 for the base population)

  #Split the base population data by sex
  dam_data = data.frame(dam = ped[which(sex==2),1], n_offs = 0, gen = 0, y = 0, f = 0)
  rownames(dam_data) = dam_data$dam

  sire_data = data.frame(sire = ped[which(sex==1),1], n_offs = 0, gen = 0, y = 0, f = 0)
  rownames(sire_data) = sire_data$sire

  #All animals will be allowed to mate for the Generation 1
  dams_ava = dam_data[,c(1,5)] 
  sires_ava = sire_data[,c(1,5)]

  #Prepare the ped file with all necessary information
  ped = data.frame(ped)
  ped$f = 0
  colnames(ped) = c("anim", "sire","dam","gen","f")
  #-----------------------------------------------------------------------------------------------#

  #Looping for each generation in the expanded population
  for (j in 1:ngen){
    gen = j
    nanim = nrow(dams_ava) #Set the number of animals per generation (Equal to the number of available dams)
  
    #This data.frame will store the offspring data for the current generation
    sel_cand = data.frame("anim" = rep(0,nanim),"n_offs" = rep(0,nanim),
                        "gen" = rep(0,nanim),"y" = rep(0,nanim),"f" = rep(0,nanim), sex = rep(0, nanim))
  
    #Nested loop for creating the offspring data 
    for (i in 1:nanim){
      ind = sample(1:nrow(dams_ava),1) #Sample the animal's dam 
      dam = dams_ava[ind,1] #Store the dam id
      f_dam = dams_ava[ind,2] #Store the dam inbreeding coefficient
      dams_ava = dams_ava[-ind,] #Delete the sampled dam (dams are allowed to have only 1 offspring per generation)
      sire = sires_ava[sample(1:nrow(sires_ava),1),1] #Same process for the animal's sire
      f_sire = sires_ava[sample(1:nrow(sires_ava),1),2] 
    
      anim = nrow(ped) + 1 #Store the id of the born animal
      dam_data[as.character(dam), 2] =  dam_data[as.character(dam), 2] + 1 #Add 1 offspring for the sampled dam
      sire_data[as.character(sire), 2] =  sire_data[as.character(sire), 2] + 1 #Add 1 offspring for the sampled sire
      f = calcInbreeding(ped[,1:3])[nrow(ped)] #Computes the inbreeding coefficient for the new animal (Requires the package pedigree)
      ped = rbind(ped, c(anim, sire, dam, gen,f)) #Allocate the offspring information on the ped data
      sex_offs = rbinom(n=1, size=1, prob=0.5)+1 #Sample the offspring sex from a binomial distribution (1 = M, 2 = F)
      sex = rbind(sex, c(sex_offs)) #Allocate the offspring sex on the sex array 
      n = nrow(sex) #Store the current row number for the simulated data
    
      #Compute the offspring TBV (True Breeding Value), with bvi = 1/2*(EBV_Si + EBV_Di) + Mi
      #In which Mi is the Mendelian deviation for the animal i, with Mi ~ N(0, 0.5*(1-0.5*(f_sire+f_dam))*varu)
      #Without inbreed parents this reduces to Mi ~ N(0, 0.5*varu)
      u = rbind(u, c(0.5*u[sire] + 0.5*u[dam] + rnorm(n = 1,mean = 0, sd = sqrt(0.5*(1-0.5*(f_sire+f_dam))*varu))))
    
      #Sample and store the protein percentage that the animal will receive on the diet
      feed = rbind(feed, protein[sample(1:7,1)])
    
      #Simulate the animal phenotype 
      phen = as.numeric(mu + sol_sex[sex_offs] + (feed[n]*slope) + u[n] + rnorm(n = 1,mean = 0,sd = sqrt(vare)))
      y = rbind(y, c(phen)) #Store the phenotype
    
      #Store all relevant data for the animal on the selection candidates data
      sel_cand$anim[i] = anim; sel_cand$gen[i] = gen; sel_cand$y[i] = phen; 
      sel_cand$f[i] = f; sel_cand$sex[i] = sex_offs
    
    }
  
    #Split the selection candidates according to sex
    m_data = subset(sel_cand, sex == 1)
    f_data = subset(sel_cand, sex == 2)
    
    #Order the males data according to the observed performance (Phenotypic Selection)
    ord = order(m_data$y, decreasing = T)
    m_data = m_data[ord,]
    m_data = m_data[1:round(i_m*nrow(m_data)),-6] #Apply the selection intensity for the male selection candidates
  
    #Same process for the females
    ord = order(f_data$y, decreasing = T) 
    f_data = f_data[ord,]
    f_data = f_data[1:round(i_f*nrow(f_data)),-6] #-6 will delete the sex (to make sure that the rbind function will work)
  
    #rename the rows and columns for the data of selected individuals
    colnames(f_data) = c("dam", "n_offs", "gen", "y", "f" )
    colnames(m_data) = c("sire", "n_offs", "gen", "y", "f" )
  
    rownames(f_data) = f_data$dam
    rownames(m_data) = m_data$sire
  
    #Sires will be culled after 3 generations
    sire_data = sire_data[which(gen-sire_data$gen <= 3),]
    sire_data = rbind(sire_data, m_data) #Add the selected individuals as new young sires
  
    #Add the selected females to the dams data
    dam_data = rbind(dam_data, f_data)
    #Compute the probability of the dams being culled given their number of offspring (roughly their age)
    p_cull= 1-(0.95*exp(-0.15*dam_data$n_offs)) #The higher the age the higher the risk of being culled
  
    #Create an variable indicating if the dam was culled 
    culling = NULL
    for (k in 1:nrow(dam_data)){
      culling[k] = rbinom(n=1,size = 1,prob = p_cull[k]) #Sample it from a binomial distribution 
    }
  
  #Delete the culled dams 
  dam_data = dam_data[which(culling==0),]
  
  #Update the information on the available sires and dams
  dams_ava = dam_data[,c(1,5)]
  sires_ava = sire_data[,c(1,5)]
  rownames(dam_data) = dam_data$dam
  rownames(sire_data) = sire_data$sire
  #Begin next generation 
  
  }

  #Packing the simulated data
  ysim = data.frame(ped,sex, feed, u, y)
  colnames(ysim) = c("animal","sire","dam", "gen", "f", "sex", "feed", "tbv","y")

  return(ysim)
}

#Lets use the created function to simulate a trait with heritability of 0.3
#Selection intensity of 0.1 for the males and 0.5 for the females
#The overall mean will be 396.4 and the residual variance 550
#The trait will be selected for 10 generations

ysim = sim_fun(h2=0.3, i_m = 0.1, i_f = 0.5, ngen = 10, mu = 396, vare = 550, seed = 182)

#-----------------------------------------
#Some sanity checks
#----------------------------------------#
head(ysim)
tail(ysim)

#Plot the histograms for the Breeding values and Phenotypes
par(mfrow=c(2,1))
hist(ysim$tbv, breaks = 30, xlab = "True Breeding Values", main = "", col = rgb(0.8,0.3,0.5,0.7)) #EBV mean is no longer 0 on the current population
hist(ysim$y[which(ysim$y>0)], breaks=30, main = "", xlab = "Phenotypes", col = rgb(0.1,0.6,0.5,0.7))

var(ysim$tbv[which(ysim$gen==10)]) #Check if the additive variance per generation approaches the desired
summary(ysim$f) #Check the inbreeding coefficient statistics 


#Plot the Selection response per generation on the simulated data
ngen = 10
mu_ebv = NULL
for (i in 0:ngen){
  mu_ebv[i] = mean(ysim$tbv[which(ysim$gen==i)])
}

mu_phen = NULL
for (i in 1:ngen){
  mu_phen[i] = mean(ysim$y[which(ysim$gen==i)])
}

f = NULL
for (i in 0:ngen){
  f[i] = mean(ysim$f[which(ysim$gen==i)])
}

varu = NULL
for (i in 0:ngen){
  varu[i] = var(ysim$tbv[which(ysim$gen==i)])
}

#Look at the following plots:
par(mfrow=c(2,2))
plot(mu_phen, type = "l", lty = 2, col ="violet", xlab = "Generations", ylab = "Average Performance")
plot(mu_ebv, type = "l", lty = 2, col ="deepskyblue1", xlab = "Generations", ylab = "Average TBV")
plot(f, type = "l", lty = 2, col ="darkgray", xlab = "Generations", ylab = "Average Population Inbreeding")
plot(varu, type = "l", lty = 2, col ="cyan3", xlab = "Generations", ylab = "Additive genetic variance")

#Now answer the questions
#1) What happens with the EBV (breeding value) through the generations? why? 
#2) What is happening with the population average inbreeding? What is the reason for that?
#3) What is happening with the Additive genetic variance? Do a search in the internet about the Bulmer effect to learn more about it.
#4) Try to simulate a new trait with different heritability and selection intensity for males and females. Look what happens



