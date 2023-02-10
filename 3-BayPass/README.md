## Identifying loci associated with invasive status of populations
We identified outlier SNPs in each main invasive cluster (clusters were identified using spatial genomic analyses) using the contrast statistic (C2) in BayPass.
BayPass accounts for the impact of shared demographic history and population structure when identifying variants that have allele frequencies that significantly deviate from the expected multivariate normal distribution.
We ran BayPass with default parameters, using three different ecotype files in three separate runs to independently compare the Queensland native populations to each main invasive cluster. 

```
module load Python
module load ifort

python ./reshaper_baypass.py bialSNP_MAF_geno_LD.vcf Qfly_popmap.txt Qfly.geno  #converting vcf to BayPass geno using reshaper_baypass.py (https://gitlab.com/YDorant/Toolbox/-/blob/master/reshaper_baypass.py)
./i_baypass -gfile Qff.geno -contrastfile Qff_nativevsAlice.ecotype -efile Qff_nativevsAlice.ecotype -outprefix Qff_nativevsAlice -nthreads 6
./i_baypass -gfile Qff.geno -contrastfile Qff_nativevsExpanded.ecotype -efile Qff_nativevsExpanded.ecotype -outprefix Qff_nativevsExpanded -nthreads 6
./i_baypass -gfile Qff.geno -contrastfile Qff_nativevsIslands.ecotype -efile Qff_nativevsIslands.ecotype -outprefix Qff_nativevsIslands -nthreads 6
```

##### Sanity check
Before selecting the candidate SNPs, it's crucial to examine the p-value behavior to prevent our outcomes from being impacted by issues related to multiple testing. We will examine the distribution of p-values by creating a histogram in R. The ideal p-value histogram should approximate a uniform distribution for higher p-values.
```
Alice.C2=read.table("Qff_nativevsAlice_summary_contrast.out",h=T) ##plot the histogram for all three runs
hist(10**(-1*Alice.C2$log10.1.pval.),freq=F,breaks=40)
abline(h=1)
```

## Choosing candidate SNPs and creating Manhattan plots
`Qff_nativevsAlice_summary_contrast.out` contains posterior mean of the C2 contrast statistics (M_C2), standard deviation of C2 contrast statistics (SD_C2), calibrated estimator of C2 statistics (C2_std) and its corrected p value (log10(1/pval)). We will use this file to calculate false discovery rates (FDR 5%) in R and choose candidate SNPs. The procedure bellow was conducted for all three BayPass runs.
```
library(qvalue)
library(dplyr)
library(ggplot2)

##import the BayPass output file and scaffold list
Alice.C2=read.table("Qff_nativevsAlice_summary_contrast.out",h=T)
scaffolds = read.table("scaffold_list.txt") ##I previuosly extracted scaffold names from VCF using bash commands: cat bialSNP_MAF_geno_LD.vcf | grep -v "#" | cut -f3 > scaffold_list.txt
Alice.C2 = as.data.frame(cbind(Alice.C2, scaffolds))

##make a simple Manhattan plot to check the distribution of outlier SNPs
plot(Alice.C2$log10.1.pval.)
abline(h=3,lty=2) _#0.001 p--value theshold_

##calculate q-value for each SNP and choose SNPs with q-values < 0.05
pvalues <- as.data.frame(10^-Alice.C2[,6]) ##convert -log p values to normal p values
Alice.C2 <- cbind(pvalues,Alice.C2)  ##append it to the df
pvalues <- as.vector(Alice.C2$`10^-Alice.C2[, 6]`)
qval <- qvalue(p = pvalues)
plot(qval$qvalues)
abline(h=0.05)  ##qvalues < 0.05 are highly sig.

qvalues = as.data.frame(qval$qvalues)

Alice.C2 <- cbind(Alice.C2,qvalues)
Alice.C2 = rename(Alice.C2, qval=`qval$qvalues`, LocusName = V1)

##extract those SNPs with q value < 0.05
selected_SNPs = Alice.C2[Alice.C2$qval < 0.05, ]
write.table(selected_SNPs,"BayPassOutliers_nativevsAlice_5prctFDR.txt", sep = "\t")

##Manhattan plot
scaffolds2 = read.table("../../Qff_chromosomes-chrNumbers.txt", header = TRUE) #Here I changed scaffold names to chromosome numbers (CHR 1-5 and number 6 for unplaced scaffolds)
Alice.C2_2 <- cbind(Alice.C2, scaffolds2)

ggplot(Alice.C2_2, aes(x=MRK, y = log10.1.pval., color = CHR)) + 
  geom_point(show.legend = FALSE, alpha = 1, size = 2) +
  scale_color_gradient(low="#D19C1D", high="#472C1B") +
  geom_point(selected_SNPs, mapping = aes(x=MRK, y = log10.1.pval.), color = "#EB4511", size = 2) +
  theme_classic()
  
 ## we will also make a Venn diagram to check the common SNPs between the three comparisons
library(ggvenn)
 
Alice = read.table("./native vs Alice/Qff_BayPassOutliers_Alice_5prctFDR_103list.txt")
Southern = read.table("./native vs expanded/Qff_BayPassOutliers_Expanded_5prctFDR_2SNPs_list.txt")
Islands = read.table("./native vs Islands/Qff_BayPassOutliers_Islands_5prctFDR_17list.txt")

x= list(Alice=Alice$V1, Southern=Southern$V1, Islands=Islands$V1)
ggvenn(x, fill_color = c("#F8D210", "#593F1F", "#DEE3CA"), fill_alpha = 0.7, stroke_size = 0.2, set_name_size = 4, stroke_color = "black", show_percentage = FALSE)
```
