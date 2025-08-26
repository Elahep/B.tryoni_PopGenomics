## Preparing input files for environmental association analysis
To determine the significance of predictor variables, including climatic and spatial variables, in explaining the alteration of allele frequency in the native and introduced ranges of the Queensland fruit fly, we conducted a gradient forest (GF) analysis. The GF analysis employed random forest machine learning algorithms to incrementally assess the combined impact of environmental variables on the variance of allele frequency. To do this, we first need to prepare input files including environmental variables for our sampling sites and population allele frequencies.

For the gradient forest check instructions in Amy's GitHub at https://github.com/AlaMetTyr/Queensland-fruit-fly-gradient-forest-analysis

### Extract environmental variables from Bioclimatic raster files

Bioclimatic variables can be downloaded from the <a href="https://www.worldclim.org/data/bioclim.html" title="WorldClim database" >WorldClim database</a>. We have downloaded raster files at a resolution of 30 s degrees for all the 19 variables. We used R to extract climatic data for each sampling site of the Queensland fruit fly from these raster files:

```
module load R-Geo/4.1.0-gimkl-2020a
R
library(raster)
library(rgdal)
bio_1 <- raster("wc2.1_30s_bio_1.tif")
bio_2 <- raster("wc2.1_30s_bio_2.tif")
bio_3 <- raster("wc2.1_30s_bio_3.tif")
bio_4 <- raster("wc2.1_30s_bio_4.tif")
bio_5 <- raster("wc2.1_30s_bio_5.tif")
bio_6 <- raster("wc2.1_30s_bio_6.tif")
bio_7 <- raster("wc2.1_30s_bio_7.tif")
bio_8 <- raster("wc2.1_30s_bio_8.tif")
bio_9 <- raster("wc2.1_30s_bio_9.tif")
bio_10 <- raster("wc2.1_30s_bio_10.tif")
bio_11 <- raster("wc2.1_30s_bio_11.tif")
bio_12 <- raster("wc2.1_30s_bio_12.tif")
bio_13 <- raster("wc2.1_30s_bio_13.tif")
bio_14 <- raster("wc2.1_30s_bio_14.tif")
bio_15 <- raster("wc2.1_30s_bio_15.tif")
bio_16 <- raster("wc2.1_30s_bio_16.tif")
bio_17 <- raster("wc2.1_30s_bio_17.tif")
bio_18 <- raster("wc2.1_30s_bio_18.tif")
bio_19 <- raster("wc2.1_30s_bio_19.tif")
places <- read.delim("Qff_coordinates.txt", header = TRUE)
coords<-data.frame(lon=places[,2], lat=places[,3])
coords$lon <- as.numeric(coords$lon)
coords$lat <- as.numeric(coords$lat)
coordinates(coords) <- c("lon","lat")
#double check the distribution of points on the map:
library(maps)
map()
points(coords, pch=16)
dev.off()
#use "extract" function to get values for our coordinates:
val_bio_1 <- extract(x=bio_1, y=coords) # this will give a vector of length = number of populations of the value of bio_1 for each sampling location (or population).
val_bio_2 <- extract(x=bio_2, y=coords)
val_bio_3 <- extract(x=bio_3, y=coords)
val_bio_4 <- extract(x=bio_4, y=coords)
val_bio_5 <- extract(x=bio_5, y=coords)
val_bio_6 <- extract(x=bio_6, y=coords)
val_bio_7 <- extract(x=bio_7, y=coords)
val_bio_8 <- extract(x=bio_8, y=coords)
val_bio_9 <- extract(x=bio_9, y=coords)
val_bio_10 <- extract(x=bio_10, y=coords)
val_bio_11 <- extract(x=bio_11, y=coords)
val_bio_12 <- extract(x=bio_12, y=coords)
val_bio_13 <- extract(x=bio_13, y=coords)
val_bio_14 <- extract(x=bio_14, y=coords)
val_bio_15 <- extract(x=bio_15, y=coords)
val_bio_16 <- extract(x=bio_16, y=coords)
val_bio_17 <- extract(x=bio_17, y=coords)
val_bio_18 <- extract(x=bio_18, y=coords)
val_bio_19 <- extract(x=bio_19, y=coords)
df <- data.frame(val_bio_1,val_bio_2,val_bio_3,val_bio_4,val_bio_5,val_bio_6,val_bio_7,val_bio_8,val_bio_9,val_bio_10,val_bio_11,val_bio_12,val_bio_13,val_bio_14,val_bio_15,val_bio_16,val_bio_17,val_bio_18,val_bio_19)
write.table(df, "Qff_allBioclim.txt", sep = "\t")
## We will remove highly correlated (|r| > 0.75) bioclimatic variables:
pairs.panels(df, scale=T)  # by looking at the correlation coefficients, we reomved highly correlated varaibles and retained 6 variables for further analysis
##remaining variables:
#bio_3
#bio_5
#bio_8
#bio_9
#bio_12
#bio_19
```


### Create population allele frequencies from vcf file

We will use R to read vcf, convert it to a genotype matrix (individuals in rows and SNPs in columns) and use the genotype matrix to create population allele frequencies:


```
library(vcfR)
library(adegenet)
data <- read.vcfR("./Qff_allSNPs.vcf")
geno <- extract.gt(data)
dim(geno)
G <- geno
G[geno %in% c("0/0")] <- 0  # 0 when the individual is homozygous for the major allele
G[geno  %in% c("0/1")] <- 1  # 1 when the individual is heterozygous
G[geno %in% c("1/1")] <- 2  # 2 when the individual is homozygous for the second (or alternative) allele
gen <- t(G)   # this is our genotype matrix
dim(gen)
sum(is.na(gen))  # we can see there are missing data here. We will impute them after calculating population allele frequencies
Genotypes <- as.data.frame(gen)
popmap <- read.table("./Qff_popmap_bysite.txt")

Genotypes <- Genotypes[match(popmap$V1, row.names(Genotypes), nomatch = 0),] ##doing this to make sure all the individual in the popmap file are present in our Genotype file.
write.table(Genotypes,"./vegan RDA/Genotypes.txt", sep = "\t") ##the aggregate
#function below doesn't like the : in my loci names(in columns). so I'll 
#save the genotype table and in notepad replace : with _.

Genotypes2 <- read.delim("./vegan RDA/Genotypes.txt") #read the table with : replaced with _


AllFreq <- aggregate(Genotypes2, by = list(popmap$V2), function(x) mean(x, na.rm = T)/2)  # Estimating population allele frequencies
row.names(AllFreq) <- as.character(AllFreq$Group.1)
## imputing missing genotypes with the median of the locus allele frequencies across all populations
for (i in 1:ncol(AllFreq))
{
  AllFreq[which(is.na(AllFreq[,i])),i] <- median(AllFreq[-which(is.na(AllFreq[,i])),i], na.rm=TRUE)
}
## filtering on MAF again (although we had done it during SNP calling, we have now done missing data imputation and that could have resulted in alleles with small MAF)
freq_mean <- colMeans(AllFreq[,-1])
AllFreq <- AllFreq[,-which(freq_mean>=0.95 | freq_mean<=0.05)] # this final genetic matrix includes
                                                               # allele frequency data for 28 populations and 6526 SNPs
write.table(AllFreq, "pop_alleleFreq.txt", sep = "\t")  # save population allele frequency dataset. This will be used in the Gradient Forest, RDA and LFMM analyses.
```


****************
###### References
<a href="https://popgen.nescent.org/2018-03-27_RDA_GEA.html" title="Detecting multilocus adaptation using Redundancy Analysis (RDA)" >Detecting multilocus adaptation using Redundancy Analysis (RDA)</a>


<a href="https://github.com/Capblancq/RDA-landscape-genomics" title="RDA applications for landscape genomics" >RDA applications for landscape genomics</a>
