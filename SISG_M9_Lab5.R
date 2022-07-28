#This script gives a short example on how to perform a QTL mapping analysis
#install.packages("qtl")
library(qtl)

#Read the data
data<-read.cross("csv", "\\Users\\Badger\\Desktop\\Labs SISG modules 9 and 12", "IntMap.csv", estimate.map = FALSE)

# summary of the data
summary(data)

# description of the data
boxplot(phenotype ~ sex, data = data$pheno, horizontal = TRUE, xlab = "phenotype", ylab = "sex", col = c("blue", "green"))

anova(aov(phenotype ~ sex, data = data$pheno))

# interval mapping with sex as an additive covariate
data<-calc.genoprob(data, step = 1, error.prob = 0.001)

qtl_woc<-scanone(data, method = "hk")
a_qtl<-scanone(data, addcovar = data$pheno[,2], method = "hk")

plot(qtl_woc, a_qtl, col = c("blue", "red"), lty = c(1,2), ylab = "LOD score")
plot(a_qtl - qtl_woc, col = c("blue"), ylab = "LOD with cov. - LOD without cov.")

# permutation test with sex as an additive covariate
perm_woc<-scanone(data, method = "hk", n.perm = 300, verbose = FALSE)

a_perm<-scanone(data, addcovar = data$pheno[,2], method = "hk", n.perm = 300, verbose = FALSE)

par(mfrow = c(1,2), las = 1)

plot(perm_woc, breaks = seq(0,6,len = 50))

plot(a_perm, breaks = seq(0,6,len = 50))

both_a_qtl<-c(qtl_woc, a_qtl, labels = c("nocov", "cov"))

both_a_perm<-cbind(perm_woc, a_perm, labels = c("nocov", "cov"))

# significance thresholds with sex as an additive covariate

summary(both_a_perm, alpha = c(0.1, 0.05, 0.01))

summary(both_a_qtl, perms = both_a_perm, format = "allpeaks", alpha = 0.20, pvalues = TRUE)

####################################################

####################################################

# interval mapping with sex as an interactive covariate

data<-calc.genoprob(data, step = 1, error.prob = 0.001)

i_qtl<-scanone(data, addcovar = data$pheno[,2], intcovar = data$pheno[,2], method = "hk")

plot(a_qtl, i_qtl, col = c("blue", "red"), lty = c(1,2), ylab = "LOD score")

plot(i_qtl - a_qtl, col = c("blue"), ylab = "LODi score")

# permutation test with sex as an interactive covariate

set.seed(54955149)

a_perm<-scanone(data, addcovar = data$pheno[,2], method = "hk", n.perm = 1000, verbose = FALSE)

set.seed(54955149)

i_perm<-scanone(data, addcovar = data$pheno[,2], intcovar = data$pheno[,2], method = "hk", n.perm = 1000, verbose = FALSE)

both_i_qtl<-c(i_qtl, i_qtl - a_qtl, labels = c("f", "i"))

both_i_perm<-cbind(i_perm, i_perm - a_perm, labels = c("f", "i"))

# significance thresholds with sex as an additive covariate

summary(both_i_perm, alpha = c(0.1, 0.05, 0.01))

summary(both_i_qtl, perms = both_i_perm, format = "allpeaks", alpha = 0.20, pvalues = TRUE)

####################################################

####################################################

# separate analysis for each sex - interval mapping

data<-calc.genoprob(data, step = 1, error.prob = 0.001)

m_qtl<-scanone(subset(data, ind = data$pheno[,2] == 1), method = "hk")

f_qtl<- scanone(subset(data, ind = data$pheno[,2] == 0), method = "hk")

plot(m_qtl, f_qtl, col = c("blue", "red"), lty = c(1,2), ylab = "LOD score")

# separate analysis for each sex - permutation test

m_perm<-scanone(subset(data, ind = data$pheno[,2] == 1), method = "hk", n.perm = 1000, verbose = FALSE)

f_perm<-scanone(subset(data, ind = data$pheno[,2] == 0), method = "hk", n.perm = 1000, verbose = FALSE)

both_qtl<-c(m_qtl, f_qtl, labels = c("male", "female"))

both_perm<-cbind(m_perm, f_perm, labels = c("male", "female"))

# separate analysis for each sex - significance thresholds

summary(both_perm, alpha = c(0.1, 0.05, 0.01))

summary(both_qtl, perms = both_perm, format = "allpeaks", alpha = 0.20, pvalues = TRUE)
