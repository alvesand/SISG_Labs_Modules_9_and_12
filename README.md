# **Instructions for the R Labs (Modules 9 and 12 of the SISG)**

This repository stores the R scripts for the modeules 9 (Quantitative Genetics) and 12 (Mixed Models in Quantitative Genetics) of the 27th Summer Institute in Statistical Genetics (SISG).

# **About the R labs**

The R scripts are divided according to the Module (9 or 12), feel free to follow along all scripts and try to understand every piece of code that you are running. 
**Discousure:** The scripts come with **no** warranty and are not free of potential bugs and typos. Feel free to communicate us any potential issue.   

## Module 9 - Quantitative Genetics   

* **SISG_M9_Lab1:** This script illustrates how the pictures in pages 43 and 44 from IntroQG-seattle-2022-Lecture01_0 can be obtained. Those plots depict the additive and dominance variance in function of the allele frequencies for a single locus.   

* **SISG_M9_Lab2:** This script illustrates how to compute the additive relationship matrix (A) using different R packages.  

* **SISG_M9_Lab3:** This script creates a function for simulating the effects of response to selection based on a trait influenced by the fixed effects of sex, percentage of protein in the diet and the genetic additive as random effect. Make sure to run the code from lines 13 to 164 together. This will create the simulation function. After created, the function accepts 7 arguments: h2 -> the trait heritability; i_m -> percentage of the best males being selected; i_f -> percentage of the best females being selected; ngen -> the number of generations simulated; mu -> the trait mean; vare -> the residual variance for the trait; seed -> the seed if you want to replicate always the same results. The result is a dataframe with 9 colummns: Animal, Sire, Dam, Gen, f, sex, feed, tbv, y. In this dataframe 0 means no observation, f represents the animal's inbreeding coefficient, tbv is the simulated true breeding value; feed is a covariable indicating the percentage of protein in the diet for that specific animal. *The main goalof this script is to think about the genetic consequences of selection in the population.*  

* **SISG_M9_Lab4:** This script illustrates how to estimate heritability with the simplest method (ANOVA). A small example is given with the Beetle Data (SAS, Saxton).  

* **SISG_M9_Lab5:**  This script gives a short example on how to perform a QTL mapping analysis.


