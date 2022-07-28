#######################################################################################
###########################                        ####################################
###########################       GWAS in R        ####################################
###########################                        ####################################
#######################################################################################

#--> Set the directory:
setwd("C:/Users/ligia/Downloads/GWAS")


#--> Packages:
#install.packages("data.table")
#install.packages("qqman")
#install.packages("gaston")
#install.packages("FactoMineR")
#install.packages("factoextra")
#install.packages("rgl")
library("data.table") # Fast and memory efficient file reader and writer
library("qqman")      # Q-Q and Manhattan Plots for GWAS Data
library("gaston")     # Genetic Data Handling (QC, GRM, LD, PCA) & Linear Mixed Models
library("FactoMineR") # Multivariate Exploratory Data Analysis and Data Mining
library("factoextra") # Extract and Visualize the Results of Multivariate Data Analyses
library("rgl")        # 3D Visualization Using OpenGL
  
  
#--> Read files:
GEN <- fread("Genotype.txt",header=T)
MAP <- fread("SNP_Map.txt",                      
             header=T,                           
             col.names=c("SNP","Chromosome","Position"), # Rename all 4 columns
             stringsAsFactors=FALSE)                     # Keep character column as character
PHE <- fread("Phenotype.txt", header = T)        

#--> Some adjustments in the files:
#GEN <- subset(GEN, GEN$SNP %in% MAP$SNP)         # Keep only those SNPs that are in the MAP file
#GEN <- GEN[order(match(GEN$SNP,MAP$SNP)),]       # Sort GEN data so its SNP column match with the SNP column on MAP data
all(GEN$SNP==MAP$SNP)                            # Check if both columns are in the same order
PHE <- subset(PHE, Animal %in% colnames(GEN))    # Keep only genotyped animals
Keep_Birds <- as.character(PHE$Animal)           # Save phenotyped birds IDs that were genotyped
GEN <- GEN[, ..Keep_Birds]                       # Keep genotypes of phenotyped birds
PHE <- PHE[order(match(PHE$Animal,colnames(GEN))),] # Sort phenotype data
all(PHE$Animal==colnames(GEN))                   # Check if both columns are in the same order
GEN[GEN=="AA"] <- 0                              # Recode AA to 0
GEN[GEN=="AB"] <- 1                              # Recode AB to 1
GEN[GEN=="BA"] <- 1                              # Recode BA to 1
GEN[GEN=="BB"] <- 2                              # Recode BB to 2
GEN[GEN=="--"] <- NA                             # Recode -- to NA (missing)
GEN <- t(GEN); dim(GEN)                          # Convert from data.table to a tranposed matrix
mode(GEN) <- "numeric"                           # Convert from character to numeric
colnames(GEN) <- MAP$SNP                         # Rename the columns with SNP names


#--> Files format:
FAM=data.frame(famid=PHE$Animal,                 # PLINK .fam format file
               id=PHE$Animal,
               father=0,
               mother=0,
               sex=0,
               pheno=0)
BIM=data.frame(chr=as.numeric(MAP$Chromosome),   # PLINK .bim format file
               id=MAP$SNP,
               dist=0,
               pos=as.numeric(MAP$Position),
               A1="A",
               A2="B")

# Up to here, GEN is just a matrix with animals in rows and SNP in columns
# In the next step, the function as.bed.matrix will use the GEN matrix, FAM and BIM files
# to create a new GEN data, which will no longer be just a matrix, it will be a bed matrix.
# A bed.matrix is a S4 class for SNP genotype matrices used for gaston package.

options(gaston.autosomes = 1:29)                # Set the number of chr
GEN <- as.bed.matrix(x=GEN,                     # x is the numeric matrix of SNPs
                     fam=FAM,                   # PLINK .fam format file
                     bim=BIM)                   # PLINK .bim format file


#--> Quality controls:

# Hardy-Weinberg Proportion exact test
GEN <- set.hwe(GEN, method="exact", verbose=FALSE)

GEN <- select.snps(x = GEN,            # x is a bed.matrix data 
                   maf>0.01 &          # maf: minor allele frequency
                     callrate>0.95 &   # callrate: SNP callrat
                     hwe>1e-5)         # hwe: p-values for each SNP


#--> Save Clean Genotypes:
write.bed.matrix(GEN, "GEN_Clean")


#--> Genomic relatioship matrix:
standardize(GEN) <- 'p' 
G <- GRM(GEN)           


#--> Principal Component Analysis (PCA)
res.pca <- PCA(G, ncp = 3, graph = FALSE)
res.hcpc <- HCPC(res.pca, graph=FALSE, min=2)
Groups <- res.hcpc$data.clust$clust
loadings <- sweep(res.pca$var$coord,2,sqrt(res.pca$eig[1:3,1]),FUN="/")
#plot3d(loadings[,1:3], col=Groups)
plot(loadings[,1:2], ylab="PC 2", xlab="PC 1", col=Groups)
PHE$PC1 <- as.numeric(loadings[,1])             # Add the 1st PC to PHE data
PHE$PC2 <- as.numeric(loadings[,2])             # Add the 2nd PC to PHE data
PHE$PC3 <- as.numeric(loadings[,3])             # Add the 3rd PC to PHE data

head(PHE)
#--> Fixed Effects
PHE$Sex <- sample(c("M","F"), nrow(PHE), replace=TRUE) # Fake Sex
PHE$CG <- sample(100:110, nrow(PHE), replace=TRUE)     # Fake CG
PHE$Sex <- as.factor(PHE$Sex)                    # Set the class of Sex to factor
PHE$CG <- as.factor(PHE$CG)                      # Set the class of CG to factor
X <- model.matrix(Phen1 ~ CG + Sex, data=PHE)    # Creates a design (or model) matrix
PHE$Phen_1 <- as.numeric(PHE$Phen_1)               # Set the class of Phen1 to numeric


#--> Bonferroni threshold
Thr <- -log10(0.05/mean(table(GEN@snps$chr)))   


#--------------------------------------------------------------------------------------------------#
#--> Association test using Linear Model - no Fixed Effect
lm.no.FE  <- association.test(x=GEN,          # A bed.matrix
                                Y=PHE$Phen_1,   # The phenotype vector
                                method="lm")    # Method to use: "lm" for (generalized) linear model
# Manhattan plot
qqman::manhattan(lm.no.FE, chr="chr", bp="pos", p="p", snp="id", main="lm - no Fixed Effects",
          cex=0.5, cex.axis=0.7, col=c("gray45","red4"), 
          suggestiveline=Thr, genomewideline=F, logp=T)
# QQ plot
qqman::qq(lm.no.FE$p)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
#--> Association test using Linear Model - with Fixed Effect
lm.FE  <- association.test(x=GEN,          # A bed.matrix
                             Y=PHE$Phen_1,    # The phenotype vector
                             X=X,            # A covariable (design) matrix. "X <- model.matrix(Phen1 ~ CG + Sex, data=PHE)"
                             method="lm")    # Method to use: "lm" for (generalized) linear model,
                                             # and "lmm" for (generalized) linear mixed model
# Manhattan plot
qqman::manhattan(lm.FE, chr="chr", bp="pos", p="p", snp="id", main="lm - with Fixed Effects",
          cex=0.5, cex.axis=0.7, col=c("gray45","red4"), 
          suggestiveline=Thr, genomewideline=F, logp=T)
# QQ plot
qqman::qq(lm.FE$p)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
#--> Association test using Linear Model with Fixed Effects and PC
X.PC <- model.matrix(Phen_1 ~ CG + Sex + PC1 + PC2 + PC3, data=PHE)
FE.PC  <- association.test(GEN, Y=PHE$Phen_1, X=X.PC, method="lm", test="wald")
# Manhattan plot
qqman::manhattan(FE.PC, chr="chr", bp="pos", p="p", snp="id", main="lm - with Fixed Effects and PC",
          cex=0.5, cex.axis=0.7, col=c("gray45","red4"), 
          suggestiveline=Thr, genomewideline=F, logp=T)
# QQ plot
qqman::qq(FE.PC$p)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
#--> Association test using a Linear Mixed Model - no Fixed Effect
MM.no.FE <- association.test(GEN, Y=PHE$Phen1, K=G,        method="lmm")
# Manhattan plot
qqman::manhattan(MM.no.FE, chr="chr", bp="pos", p="p", snp="id",  main="LMM - no Fixed Effects",
                 cex=0.5, cex.axis=0.7, col=c("gray45","red4"), 
                 suggestiveline=Thr, genomewideline=F, logp=T)
# QQ plot
qqman::qq(MM.no.FE$p)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
#--> Association test using a Linear Mixed Model - with Fixed Effect
MM.FE <- association.test(x=GEN,                # A bed.matrix
                                Y=PHE$Phen9,    # The phenotype vector
                                X=X,            # A covariable (design) matrix
                                K=G,            # A Genetic Relationship Matrix
                                method="lmm")   # Method to use: "lm" for (generalized) linear model
                                                # and "lmm" for (generalized) linear mixed model
# Manhattan plot
qqman::manhattan(MM.FE, chr="chr", bp="pos", p="p", snp="id",  main="LMM - with Fixed Effects",
                 cex=0.5, cex.axis=0.7, col=c("gray45","red4"), 
                 suggestiveline=Thr, genomewideline=F, logp=T)
# QQ plot
qqman::qq(MM.FE$p)
#--------------------------------------------------------------------------------------------------#


# gene mapping - Exemple
gm <- MM.FE[,c("id","chr","pos","p")]
gm$start_pos<-gm$pos - 50
gm$end_pos<-gm$pos + 50
gm$p_bonferroni<-p.adjust(gm$p, method = "bonferroni", n = length(gm$p))
top<-gm[order(gm$p_bonferroni),]
print(head(top,10),row.names = FALSE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("biomaRt")
library("biomaRt")

ensembl <- useMart("ENSEMBL_MART_ENSEMBL")             # loading genes database
#searchDatasets(mart = ensembl, pattern = "gallus")    # searching the data base
ensembl <- useDataset("ggallus_gene_ensembl", ensembl) # loading database
#attributes = listAttributes(ensembl)                  # to see all attributes
attributes <- c("ensembl_gene_id","chromosome_name","start_position","end_position",
                "external_gene_name","description")
filter<-"chromosomal_region"

### pic in chr 13
value1<-c("13:8895:8995")
gene.list1<-getBM(attributes=attributes, filters=filter, values=value1, mart=ensembl)
gene.list1


