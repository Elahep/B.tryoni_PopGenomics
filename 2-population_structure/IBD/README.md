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
write.table(pop.data.matrix, "indv_matrix.txt", sep = "\t")
```

A matrix of pairwise geographic distances was created by calculating pairwise great-circle distance between sampling coordinates:

```
# Reading matrix of coordinates for each individual (i.e., 301 rows of lat. and long. format)
coords <- as.matrix(read.table("Qff.coord", header = FALSE)) 

# Calculate pairwise distances using rdist.earth
distances <- rdist.earth(coords, miles = FALSE)

write.table(distances, file = "distances.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```



























