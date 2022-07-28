#This script illustrates how to compute the additive relationship matrix (A) using different R packages
#Install the packages if necessary by unblocking the next few comments
#install.packages("pedigree")
#install.packages ("kinship2")
#install.packages("pedigreemm")
#install.packages("reshape")
#intall.packages("gdata")

#Change the directory where you want to work in
setwd("C:\\Users\\Badger\\Desktop\\Labs SISG modules 9 and 12")

# Example of numerator relationship matrix
# Pedigree
# ID Sire Dam
#  1   -   - 
#  2   -   -
#  3   1   2
#  4   1   -
#  5   4   3
#  6   5   2

id <- c(1,2,3,4,5,6) # individual ID
dam <- c(NA,NA,2,NA,3,2) # Dams; NA for unknown
sire <- c(NA,NA,1,1,4,5) # Sires; NA for unknown
ped <- data.frame(id,dam,sire)

library("pedigree")
library("reshape2")
library("gdata")
makeA(ped, which = c(rep(FALSE,0), rep(TRUE,6)))
A1 <- read.table("A.txt")
A1 <- data.matrix(dcast(A1, V1 ~ V2, value.var="V3")[,-1])
upperTriangle(A1, byrow=TRUE)=lowerTriangle(A1)

library("kinship2")
A2 <- kinship2::kinship(id=id, dadid=sire, momid=dam)*2

library("pedigreemm")
A3 <- as.matrix(crossprod(relfactor(pedigreemm::pedigree(sire=sire, dam=dam, label=id))))
A4 <- as.matrix(pedigreemm::getA(pedigreemm::pedigree(sire=sire, dam=dam, label=id)))



