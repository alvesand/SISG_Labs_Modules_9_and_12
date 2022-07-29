#############################################################################
#This script illustrates how the pictures in pages 43 and 44 from IntroQG-seattle-2022-Lecture01_0
#can be obtained. Those plots depict the additive and dominance variance in function of the
#allele frequencies for a single locus
#According to Falconer & Mackay (1996) the additive and dominance variances
#for a single locus a given by VA = 2pq(a + d(q-p))^2 and (2pqd)^2
#in which p and q are the allele frequencies, a and d are the additive and dominance 
#effects

#Get a range of frequencies for the first allele
p = seq(0,1,0.01)
#Assuming Hardy-Weinberg Equilibrium p+q must be 1
q = 1-p

#First lets assume no dominance
a = 10
d = 0
VA = 2*p*q*(a + d*(q-p))^2
plot(VA~q, type = "l", lty = 2, ylab = "Variance") 
legend(legend = "VA", lty = 2, "topleft")

#Assuming complete dominance (d = a)
a = 10
d = 10
VA = 2*p*q*(a + d*(q-p))^2
VD = (2*p*q*d)^2

plot(VA~q, type = "l", lty = 2, ylab = "Variance", ylim = c(min(VA, VD), max(VA,VD))) 
points(VD~q, type = "l", lty = 2, col = "red") 
legend(legend = c("VA", "VD"), lty = c(2,2), col = c("black", "red"), "topleft")

#Questions
#1) What happens with the Additive variance under no dominance and complete dominance?
#2) Change the size of the additive effect for the first plot and see what happens with the VA scale
#3) Plot the graph for a situation of over-dominance (d > a)
#4) Plot the graph considering a = 0 and d = 5, what happens with the additive variance in function of the
#allelic frequencies?

