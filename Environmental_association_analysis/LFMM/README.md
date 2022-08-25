## Detecting SNPs under selection using the LFMM analysis


LFMM tests for association between loci and environmental variables. LFMM can take into account the effect of population strcuture by introducing hidden latent factors in the model and identifying nonrandom associations between SNPs and environmental variables.
The optimal number of genetic clusters in the dataset can determine the number of latent factors in the model. Using the _sNMF_ function of the R package _LEA_, we can explore population structure and choose the optimal K.


We will run this analysis in R:

```
library(LEA)    
library(qvalue)   
library(lfmm) 

##### import input files #####

#read population allele frequency data
gen <- read.table("pop_alleleFreq.txt")[,-1]

#create lfmm format from vcf (the vcf is in our working directory) for sNMF analysis

vcf2lfmm(input.file = "./Qff.vcf")

#import env data
All_env <- read.delim("Qff_allBioclim.txt")

#choose the non correlated env variables:
env_subset <- All_env[,c("bio_3","bio_5","bio_8","bio_9","bio_12","bio_19")]

##### assessing population structure using snmf to choose the optimal K  ######

project = NULL
project = snmf("Qff.lfmm", K = 1:10, entropy = TRUE, repetitions = 20, project = "new")
plot(project, col = "blue", pch = 19, cex = 1.2)  ##I will use k=3 

##### running lfmm #####

Qff.lfmm <- lfmm_ridge(Y = gen, X = env_subset, K = 3) #The ridge estimates are based on minimimizing a regularized least-squares problem with an L2 penalty.

# performs association testing using the fitted model:
pv <- lfmm_test(Y = gen, X = env_subset, lfmm = Qff.lfmm, calibrate="gif")

#The histogram of test significance values is expected to be flat, with a peak near zero.
#A QQ-plot is displayed as follows.

pvalues <- pv$calibrated.pvalue 
qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)

#check genomic inflation factor which gives us a sense for how well the model has accounted for confounding factors in the data.
pv$gif #An appropriately calibrated set of tests will have a GIF of around 1.

#check how applying the GIF to the pvalues can change the pvalue distribution for each env variables

hist(pv$pvalue[,1], main="Unadjusted p-values bio_3")        
hist(pv$calibrated.pvalue[,1], main="GIF-adjusted p-values bio_3") #we should see a flat histogram with a peak for smaller pvalues (near zero)
                                                                   
hist(pv$pvalue[,2], main="Unadjusted p-values bio_5")        
hist(pv$calibrated.pvalue[,2], main="GIF-adjusted p-values bio_5") 

hist(pv$pvalue[,3], main="Unadjusted p-values bio_8")        
hist(pv$calibrated.pvalue[,3], main="GIF-adjusted p-values bio_8") 

hist(pv$pvalue[,4], main="Unadjusted p-values bio_9")        
hist(pv$calibrated.pvalue[,4], main="GIF-adjusted p-values bio_9") 

hist(pv$pvalue[,5], main="Unadjusted p-values bio_12")        
hist(pv$calibrated.pvalue[,5], main="GIF-adjusted p-values bio_12") 

hist(pv$pvalue[,6], main="Unadjusted p-values bio_19")        
hist(pv$calibrated.pvalue[,6], main="GIF-adjusted p-values bio_19") 

#### Estimate adjusted p-values ####
pvalues <- pv$pvalue 
zs <- pv$score  

summary(zs) #shows the env variables that were used in the analysis and their associated z scores


scaffolds <- row.names(qv)
calibratedpvalues <- as.data.frame(pv$calibrated.pvalue)


qv_bio_3 <- as.data.frame(qvalue(calibratedpvalues$bio_3)$qvalues)
candidtaes_bio_3 <- scaffolds[which(qv_bio_3 < 0.1)]
candidtaes_bio_3

qv_bio_5 <- as.data.frame(qvalue(calibratedpvalues$bio_5)$qvalues)
candidtaes_bio_5 <- scaffolds[which(qv_bio_5 < 0.1)]
candidtaes_bio_5

qv_bio_8 <- as.data.frame(qvalue(calibratedpvalues$bio_8)$qvalues)
candidtaes_bio_8 <- scaffolds[which(qv_bio_8 < 0.1)]
candidtaes_bio_8

qv_bio_9 <- as.data.frame(qvalue(calibratedpvalues$bio_9)$qvalues)
candidtaes_bio_9 <- scaffolds[which(qv_bio_9 < 0.1)]
candidtaes_bio_9

qv_bio_12 <- as.data.frame(qvalue(calibratedpvalues$bio_12)$qvalues)
candidtaes_bio_12 <- scaffolds[which(qv_bio_12 < 0.1)]
candidtaes_bio_12

qv_bio_19 <- as.data.frame(qvalue(calibratedpvalues$bio_19)$qvalues)
candidtaes_bio_19 <- scaffolds[which(qv_bio_19 < 0.1)]
candidtaes_bio_19


allsnps = list(bio_3 = candidtaes_bio_3, bio_5 = candidtaes_bio_5, 
         bio_8 = candidtaes_bio_8, bio_9 = candidtaes_bio_9, 
         bio_12 = candidtaes_bio_12, bio_19 = candidtaes_bio_19)

capture.output(allsnps, file = "LFMM_adaptive_SNPs.txt")    #save the list of outlier snps

#Using K=3, the default GIF correction, and an FDR threshold of 0.10, we detect 42 candidate SNPs 
#under selection in response to 6 bioclimatic variables. 33 of these SNPs are unique and will
#be considered as the final LFMM candidate SNP list based on all 6 env variables.
```

Using the R package 'ggplot2' we can create Manhattan plot to show the distribution of candidate SNPs for each bioclimatic variable:

```
library(ggplot2)
library(dplyr)

scaffolds <- read.table("LFMM_allSNPs_CHRinfo6526.txt", header = TRUE) #import chromosome and snp position information for the 6526 SNPs

#we will plot candidates under bio_3, by using the same code you can plot SNPs for other bioclim variables too.
pvals_df <- as.data.frame(pvalues)
pval_bio_3 <- as.data.frame(-log10(pvals_df$bio_3))
pval_bio_3 <- cbind(pval_bio_3,scaffolds)
pval_bio_3 <- rename(pval_bio_3, log10pval=`-log10(pvals_df$bio_3)`)
outliers_bio_3 <- read.table("bio_3_outliers.txt")
outliers_bio_3_vec <- as.vector(outliers_bio_3$V2)

selected_snps3 <- pval_bio_3[pval_bio_3$SNP %in% outliers_bio_3_vec, ] #selecting candidates SNPs under bio_3

ggplot(pval_bio_3, aes(x=MRK, y = abs(log10pval), color = CHR)) + ##abs <- absolute value of plog
  geom_point(show.legend = FALSE, alpha = 0.8, size = 3) +
  scale_color_gradient(low="#D19C1D", high="#472C1B") +
  geom_point(selected_snps3, mapping = aes(x=MRK, y = abs(log10pval)), color = "#EB4511", size = 4) +
  theme_minimal() +
  ggtitle("bio_3")
```

Here is the Manhattan plot:

![bio3](https://user-images.githubusercontent.com/13001264/186780672-7fd7fee2-4d11-4f98-8dd7-fa5db4bef68e.png)



#### References
https://bookdown.org/hhwagner1/LandGenCourse_book/WE_11.html

https://cran.r-project.org/web/packages/lfmm/vignettes/lfmm.html




