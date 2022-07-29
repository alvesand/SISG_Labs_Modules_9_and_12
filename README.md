# **Instructions for the R Labs (Modules 9 and 12 of the SISG)**

This repository stores the R scripts for the modeules 9 (Quantitative Genetics) and 12 (Mixed Models in Quantitative Genetics) of the 27th Summer Institute in Statistical Genetics (SISG).

# **About the R labs**

The R scripts are divided according to the Module (9 or 12), feel free to follow along all scripts and try to understand every piece of code that you are running. Make sure to download every package needed for the labs, also, please remember to change the working directory when necessary. To download the files go to the green tab Code -> Download ZIP
**Discousure:** The R scripts and this README document come with **no** warranty and are not free of potential bugs and typos. Feel free to communicate us any potential issue.   

## Module 9 - QUANTITATIVE GENETICS   

* **SISG_M9_Lab1:** This script illustrates how the pictures in pages 43 and 44 from IntroQG-seattle-2022-Lecture01_0 can be obtained. Those plots depict the additive and dominance variance in function of the allele frequencies for a single locus.   

* **SISG_M9_Lab2:** This script illustrates how to compute the additive relationship matrix (A) using different R packages.  

* **SISG_M9_Lab3:** This script creates a function for simulating a trait under selection. The simulated trait is influenced by the fixed effects of sex, percentage of protein in the diet and the genetic additive as random effect. Make sure to run the code from lines 13 to 164 together, this will create the simulation function. After created, the function accepts 7 arguments: h2 -> the trait heritability; i_m -> percentage of the best males being selected; i_f -> percentage of the best females being selected; ngen -> the number of generations simulated; mu -> the trait mean; vare -> the residual variance for the trait; seed -> the seed if you want to replicate always the same results. The result is a dataframe with 9 colummns: Animal, Sire, Dam, Gen, f, sex, feed, tbv, y. In this dataframe 0 means no observation, f represents the animal's inbreeding coefficient, tbv is the simulated true breeding value; feed is a covariable indicating the percentage of protein in the diet for that specific animal. *The main goal of this script is make the students think about the genetic consequences of selection in the population.*  

* **SISG_M9_Lab4:** This script illustrates how to estimate heritability with the simplest method (ANOVA). A small example is given with the Beetle Data (SAS, Saxton).  

* **SISG_M9_Lab5:**  This script gives a short example on how to perform a QTL mapping analysis.

## Module 12 - MIXED MODELS IN QUANTITATIVE GENETICS  

* **SISG_M12_Lab1:** Some basic commands to perform different matrix operations in R.

* **SISG_M12_Lab2:** This script gives a simple and small example on how the mixed model equations are set and solved.  

* **SISG_M12_Lab3:** This scripts gives some examples on how to perform a complete mixed model analysis in R. The following animal model is assumed: y = Xb + Zu + e.

* **SISG_M12_Lab4:** This script gives some examples on how to perform a Genomic Prediction Analysis in R. The main package used is the [BGLR](https://cran.r-project.org/web/packages/BGLR/BGLR.pdf).

* **SISG_M12_Lab5:** This script gives some examples on how to perform a Genome Wide Association Study (GWAS) in R. Different strategies of data correction are used to compare the influence on the results. The main package for this Lab is the [gaston](https://cran.r-project.org/web/packages/gaston/gaston.pdf).

# **About the data used**

Some Labs use external dataset that are provided in this repository. The Necessary datasets are:

* **beetle data.txt**: A txt data containing the pedigree information (Sire, Dam, Progeny), the trait observation (value) and the trait identification. The dataset contains two traits: fec and mass. 

* **IntMap.csv**: This is a csv document containing the data generated for a simulated backcross scheme with N = 400 animals. For more information please refer to [Broman and Sen (2009)](https://link.springer.com/book/10.1007/978-0-387-92125-9). The first row contains the information about the phenotype, sex and the markers identification. The second row indicates the chromosome ID. Third row indicates the marker position (cM). The genotypes are coded as A (homozygous animals) and H (Heterozygous animals), sex is coded as 0 (Male) and 1 (female).  

* **ysim.txt**: A dataset simulated by the function created with SISG_M9_Lab3. The following arguments were used: h2=0.3, i_m = 0.1, i_f = 0.5, ngen = 10, mu = 396, vare = 550, seed = 182.  

* **EM_algorithm.R**: This extra R script illustrate how to program an algorithm for variance component estimation from scratch. The programmed algorithm is the EM (Expectation Maximization) that uses the REML (Restricted Maximum Likelihehood) as estimation method. 

* **pheno.txt:** A simulated dataset containing the animal ID, Sex, and observed value (N = 4000). 

* **PED.txt:** The pedigree data for pheno.txt (11,000 animals in the pedigree)

* **GEN.txt:** The genotype data for pheno.txt. Includes the genotype information for the animals 110001 to 111000 (the last 1000 animals in PED.txt). The animals are genotyped for 5000 SNP markers. The AA, AB, and BB genotypes are coded as 0, 1 and 2, respectively. 

* **phenotype.txt:** This dataset contains the animal ID and 10 simulated traits. All traits are uncorrelated and have different genomic regions affecting them. 

* **genotype.txt:** The genotypes for animals in *phenotype.txt* each column gives the animals' full genotype. The animals' ID are given in the columns of the first row. There are 14,500 genotyped loci for 1,000 animals. 

* **snp_map.txt:** the map information for the *genotype.txt* file. This archive contains three columns: SNP (the SNP ID), Chromosome (The Chromosome ID), Position (THe position in base pairs). 








