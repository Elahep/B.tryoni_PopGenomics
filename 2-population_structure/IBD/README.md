To identify patterns of gene flow, isolation, and discrete population structure we used EEMS and conStruct.

## EEMS (Estimating Effective Migration Surfaces)

EEMS implements the concept of isolation by resistance, which characterises variation in migration rates between adjacent locations in a stepping-stone model, to identify barriers to gene flow across the landscape using geographically indexed genetic data. To run EEMS we need three input files:
_.coord_ file including the coordinates of each individual sample; _.outer_ file containing an enclosed habitat polygon; _.diffs_ containing the matrix of genetic dissimilarities.  

To prepare a habitat polygon based on the geographic distribution of our sampling sites, we used the Google Maps API v3 tool (http://www.birdtheme.org/useful/v3tool.html). We excluded the Tahiti Island population to obtain an enclosed polygon around Australia. 

To prepare the genetic dissimilarity matrix, we used the bed2diffs_v1 R function (distributed with the EEMS software):

```
# First converting VCF to str:
library(vcfR)
library(poppr)
library(radiator)

snp_vcf2 = read.vcfR("D:/fastq R1/EP4/M1n2_FINAL/depth_miss39prcnt_noCamp4.vcf")
snp_vcf2
pop.data2 = read.table("D:/fastq R1/EP4/M1n2_FINAL/Onpopmap_noCamp4.txt", header = F)
gl.snp2 <- vcfR2genlight(snp_vcf2)
pop(gl.snp2) <- rep(pop.data2$V2)

genomic_converter(data = gl.snp2, output = "structure")

# Now using bed2diffs_v1
library(adegenet)
library(stringr)
bed2diffs_v1 <- function(Geno) {
  nIndiv <- nrow(Geno)
  nSites <- ncol(Geno)
  Diffs <- matrix(0, nIndiv, nIndiv)
  
  for (i in seq(nIndiv - 1)) {
    for (j in seq(i + 1, nIndiv)) {
      x <- Geno[i, ]
      y <- Geno[j, ]
      Diffs[i, j] <- mean((x - y)^2, na.rm = TRUE)
      Diffs[j, i] <- Diffs[i, j]
    }
  }
  Diffs
}

data2 <- read.structure("Qff.str",
                       onerowperind = FALSE, 
                       n.ind = 293, n.loc = 6707, col.lab = 1, col.pop = 2, ask = FALSE)
data2

Geno2 <- data2@tab
## Keep SNPs that are observed to be bi-allelic.
multi.loci <- names(which(data2@loc.n.all != 2))
## Explanation: 
## Suppose we want to remove locus, `L100` with alleles, `L100.00`, `L100.01`, `L100.02`, 
## then detect columns whose names match the regular expression `^L100\\.\\d+$`
multi.cols <- which(grepl(paste0("^", multi.loci, "\\.\\d+$", collapse = "|"), colnames(Geno2)))
if (length(multi.cols)) Geno2 <- Geno2[, - multi.cols]
dim(Geno2)

stopifnot(identical(data2@type, 'codom'))
Geno2 <- Geno2[, str_detect(colnames(Geno2), "\\.01$")]

Diffs_v1 <- bed2diffs_v1(Geno2)

##checking eigenvalues. Check that the dissimilarity matrix has 
##one positive eigenvalue and nIndiv-1 negative eigenvalues, as required by 
##a full-rank Euclidean distance matrix. Choose the diff matrix that meets this criterion.

eigneV1 = eigen(Diffs_v1)
sort(eigneV1$values)

Diffs_v1 <- round(Diffs_v1, digits = 293)

sort(round(eigen(Diffs_v1)$values, digits = 2))

write.table(Diffs_v1, "Qff.diffs", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
```


EEMS can be run on a terminal. We ran two independent MCMC runs for 5 million generations using different seeds. We prepared the _params_Qff_chain1.ini_ input file which contains the following parameters:

```
the datapath = ../data/Qff
mcmcpath = ../data/Qff-chain1
nIndiv = 293
nSites = 6707
nDemes = 200
seed = 1112004 
numMCMCIter = 5000000
numBurnIter = 2000000
numThinIter = 9999
qEffctProposalS2 = 0.03
qSeedsProposalS2 = 0.01
mEffctProposalS2 = 3
mrateMuProposalS2 = 0.1
mSeedsProposalS2 = 0.01
```

Then, we ran the program using the following command:

```
module load Eigen
module load Boost
./runeems_snps --params params_Qff_chain1.ini
```

## ConStruct 

To identify whether genetic divergence is represented as discrete clusters or continuous clines, we used the R package conStruct v. 1.0.5. We need two input files for running conStruct: a matrix of allele frequency for each individual and a matrix of pairwise geographic distances between all samples.

A matrix of allele frequency was created using structure2conStruct function of the R package conStruct:

```
library(conStruct)
conStruct.data <- structure2conStruct(infile = "./Qff.str",
                                      onerowperind = FALSE,
                                      start.loci = 3,
                                      start.samples = 1,
                                      missing.datum = -9,
                                      outfile = "./myConStructData")
is.matrix(conStruct.data)

pop.data.matrix <- matrix(NA,nrow=28,ncol=ncol(conStruct.data))
pop.index = c(1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,18,18,18,18,18,18,18,18,18,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,21,21,21,21,21,21,21,21,21,22,22,22,22,22,22,22,22,22,22,22,22,22,23,23,23,23,23,23,23,23,24,24,24,24,24,24,24,24,24,24,25,25,25,25,25,25,25,25,25,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,27,27,27,27,27,27,27,27,27,27,28,28,28,28,28,28,28,28)
for(i in 1:nrow(pop.data.matrix)){
  pop.data.matrix[i,] <- colMeans(
    conStruct.data[
      which(pop.index==i),,
      drop=FALSE
    ],na.rm=TRUE
  )
}

is.matrix(pop.data.matrix)
write.table(pop.data.matrix, "construct", sep = "\t")
```

A matrix of pairwise geographic distances was created by calculating pairwise great-circle distance between sampling coordinates:

```
# Reading matrix of coordinates for each individual (i.e., 301 rows of lat. and long. format)
coords <- as.matrix(read.table("Qff_coordinates_latlon.txt", header = FALSE)) 

# Calculate pairwise distances using rdist.earth
distances <- rdist.earth(coords, miles = FALSE)

write.table(distances, file = "distances.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```

We then ran spatial and non-spatial models for K=1-6, comparing them using cross-validation:

```
library(conStruct)
library(parallel)
library(foreach)
library(doParallel)

setwd("/nesi/nobackup/XXXX/Eli/conStruct_2")
allele_freq <- as.matrix(read.table("./construct"))
coords <- as.matrix(read.table("./Qff_coordinates_latlon.txt"))
distances <- as.matrix(read.table("./distances.txt", header = TRUE))

ncpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))

cluster <- makeCluster(ncpus)

registerDoParallel(cluster)

# Running spatial vs non-spatial model comparison using a cross-validation analysis:

my.xvals <- x.validation(train.prop = 0.9,
                         n.reps = 3,
                         K = 1:6,
                         freqs = allele_freq,
                         data.partitions = NULL,
                         geoDist = distances,
                         coords = coords,
                         prefix = "Qf_test",
                         n.iter = 5000,
                         make.figs = TRUE,
                         save.files = TRUE,
                         parallel = TRUE,
                         n.nodes = ncpus)

stopCluster(cluster)

## Visualization of results:

sp.results <- as.matrix(
                read.table("Qf_test_sp_xval_results.txt",
                           header = TRUE,
                           stringsAsFactors = FALSE)
               )
nsp.results <- as.matrix(
                read.table("Qf_test_nsp_xval_results.txt",
                           header = TRUE,
                           stringsAsFactors = FALSE)
               )


write.table(sp.results, "./sp.results.txt", sep = "\t")
write.table(nsp.results, "./nsp.results.txt", sep = "\t")

## Plotting the output:

# First, get the 95% confidence intervals for the spatial and nonspatial
#  models over values of K (mean +/- 1.96 the standard error)

sp.CIs <- apply(sp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})

# Then, plot cross-validation results for K=1:7 with 5 replicates
pdf(file="./cross_validation_plots.pdf")
par(mfrow=c(1,2))
plot(rowMeans(sp.results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(sp.results,nsp.results),
     main="cross-validation results")
    points(rowMeans(nsp.results),col="green",pch=19)

#   Finally, visualize results for the spatial model
#   separately with its confidence interval bars.

#   note that you could do the same with the spatial model, 
#   but the confidence intervals don't really show up 
#   because the differences between predictive accuracies
#   across values of K are so large.

plot(rowMeans(sp.results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(sp.CIs),
     main="spatial cross-validation results")
segments(x0 = 1:nrow(sp.results),
         y0 = sp.CIs[1,],
         x1 = 1:nrow(sp.results),
         y1 = sp.CIs[2,],
         col = "blue",lwd=2)

dev.off()
```

This gave us the following plot, indicating that the spatial model is favored over the non-spatial model, so IBD is probably a feature of our data:


![cross_validation_plots_page-0001](https://github.com/Elahep/B.tryoni_PopGenomics/assets/13001264/88d1b463-7be7-4003-ad56-45f652af7309)


We now focused on the spatial model, by running it separately for 10,000 iterations to calculate layer contribution. We set a threshold criterion of 0.02 for the covariate contribution, selecting the k value that fell above this cutoff as the most suitable.

```
library(conStruct)
setwd("/nesi/nobackup/XXXX/Eli/sp_models")
allele_freq <- as.matrix(read.table("./construct"))
coords <- as.matrix(read.table("./Qff_coordinates_latlon.txt"))
distances <- as.matrix(read.table("./distances.txt", header = TRUE))

# Running spatial model for K=1-6 and more iterations
for (k in 1:6) {
    sp_k <- conStruct(spatial = TRUE, 
                      K = k, 
                      freqs = allele_freq,
                      geoDist = distances, 
                      coords = coords,
                      prefix = paste0("spK", k), 
                      n.iter = 10000, 
                      n.chains = 1, 
                      make.figs = TRUE, 
                      save.files = TRUE)
}

## Calculating layer contribution:

# Loop through output files generated by conStruct 
#   runs with K=1 through 6 and calculating the 
#   layer contributions for each layer in each run  

layer.contributions <- matrix(NA,nrow=6,ncol=6)

# Load the conStruct.results.Robj and data.block.Robj
#   files saved at the end of a conStruct run
load("spK1_conStruct.results.Robj")
load("spK1_data.block.Robj")

# Calculate layer contributions
layer.contributions[,1] <- c(calculate.layer.contribution(conStruct.results[[1]],data.block),rep(0,5)) #rep = (0,max_k-1)
tmp <- conStruct.results[[1]]$MAP$admix.proportions

for(i in 2:6){
  # load the conStruct.results.Robj and data.block.Robj
  #   files saved at the end of a conStruct run
  load(sprintf("spK%d_conStruct.results.Robj",i))
  load(sprintf("spK%d_data.block.Robj",i))
  
  # match layers up across runs to keep plotting colors consistent
  # for the same layers in different runs
  tmp.order <- match.layers.x.runs(tmp,conStruct.results[[1]]$MAP$admix.proportions)  
  
  # Calculate layer contributions
  layer.contributions[,i] <- c(calculate.layer.contribution(conStruct.results=conStruct.results[[1]],
                                                            data.block=data.block,
                                                            layer.order=tmp.order),
                               rep(0,6-i))
  tmp <- conStruct.results[[1]]$MAP$admix.proportions[,tmp.order]
}

write.table(layer.contributions, "layer_contribution_sp.txt", sep = "\t")

## Plotting layer contributions across values of k
pdf(file = "layer_contribution_sp.pdf")
barplot(layer.contributions,
        col=c("orchid","goldenrod", "deeppink", "royalblue", "firebrick", "darkgreen", "steelblue"),
        xlab="",
        ylab="layer contributions",
        names.arg=paste0("K=",1:6))
dev.off()
```

#### References
EEMS and conStruct manuals and GitHub pages.

















