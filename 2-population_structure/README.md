## Population structure and diversity

We used the _populations_ module of Stacks to calculate population genomic indices including overall nucleotide diversity (π), average observed and expected heterozygosity, and inbreeding coefficient (Fis) for each local population.
For this purpose, we need a population map file which contains the list of each individual sample and its corresponding local population:

```
head Qfly_popmap.txt
Weipa_1	Weipa
Weipa_2	Weipa
Weipa_3	Weipa
Weipa_4	Weipa
Weipa_5	Weipa
Mapoon_1	Mapoon
Mapoon_2	Mapoon
Mapoon_3	Mapoon
Mapoon_4	Mapoon
...
```

Now let's use Stacks to calculate the indices:
```
module load Stacks
populations -V bialSNP_MAF_geno_LD.vcf -O ./ -M Qfly_popmap.txt 
```

Next, we will use R to conduct a PCA and calculate between-population Fst. Here I have also included the code we used to create the PCA and Fst heatmap plots.

```
library(vcfR)
library(adegenet)
library(StAMPP)
library(ggplot2)
library(corrplot)

##running PCA
snp_vcf2 = read.vcfR("./bialSNP_MAF_geno_LD.vcf")
pop.data2 = read.table("./Qfly_popmap.txt", header = F)
gl.snp2 <- vcfR2genlight(snp_vcf2)
pop(gl.snp2) <- rep(pop.data2$V2)
snp.pca2 <- glPca(gl.snp2, nf = 10)

##write PCA scores
snp.pca.scores2 <- as.data.frame(snp.pca2$scores)
snp.pca.scores2$pop <- pop(gl.snp2)
write.table(snp.pca.scores2, "./Qfly_adegenetPCA.txt", sep = "\t")

##export list of eigen values and percentage variances for each PC
#eigen values for PCs
eig.val<-snp.pca2$eig
eig.val
#percentages of variance for each PC
eig.perc <- 100*snp.pca2$eig/sum(snp.pca2$eig)
eig.perc
eigen<-data.frame(eig.val,eig.perc)
eigen
#writing file with both
write.csv(eigen,file="./Qfly_adegenetPCA_eigen-summary.csv",row.names=TRUE,quote=FALSE)

##PCA plotting
data2 = read.delim("Qfly_adegenetPCA.txt") #I have manually added population information to this file prior to loading it
mycol = c("#f1c039","#f37d21", "#51692d", "#56ba32")
ggplot(data2, aes(x=PC1, y=PC2, color=pop2)) +
  geom_point(size = 2) + 
  scale_color_manual(values=mycol) +
  theme_classic()+
  xlab("PC1 (1.99%)") +
  ylab("PC2 (1.34%)")


##calculating Fst between populations
Qfly_Fst <- stamppFst(gl.snp2, nboots = 100, percent = 95, nclusters = 6)
Fst <- Qfly_Fst$Fsts
pFst <- Qfly_Fst$Pvalues
write.table(Fst, "Fst.txt", sep="\t")
write.table(pFst, "Fst_pvalue.txt", sep="\t")

##creating heatmap
# Melt the correlation matrix
Q_Fst <- as.matrix(read.table("Fst.txt"))
QflyFs <- melt(Q_Fst, na.rm = TRUE)
summary(QflyFs$value)
ggplot(data = QflyFs, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#ffd60a", high = "#001d3d", mid = "#4e9de6", 
                       midpoint = 0.056, limit = c(0.005,0.11), space = "Lab", 
                       name="Pairwise Fst") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
  ```
  
  
We have also used R to anaylze population structure using a non-negative matrix factorisation (sNMF) analysis  and create admixture plots:
```
library(LEA)
library(pophelper)

##creating input files
vcf2geno(input.file = "bialSNP_MAF_geno_LD.vcf", output.file = "Qff.geno")

##snmf clustering
projectalpha = NULL
projectalpha = snmf("Qff.geno", K = 1:10, repetitions = 50, entropy = T, CPU = 8, project = "new")

# plot cross-entropy criterion for all runs in the snmf project
pdf(file = "./cross_ent_alphadefualt.pdf")
plot(projectalpha, col = "maroon4", pch = 19, cex = 1.2)
dev.off()

best2 = which.min(cross.entropy(projectalpha, K = 2))
best2
best3 = which.min(cross.entropy(projectalpha, K = 3))
best3
best4 = which.min(cross.entropy(projectalpha, K = 4))
best4
best5 = which.min(cross.entropy(projectalpha, K = 5))
best5
best6 = which.min(cross.entropy(projectalpha, K = 6))
best6
best7 = which.min(cross.entropy(projectalpha, K = 7))
best7
best8 = which.min(cross.entropy(projectalpha, K = 8))
best8
best9 = which.min(cross.entropy(projectalpha, K = 9))
best9
best10 = which.min(cross.entropy(projectalpha, K = 10))
best10

##creating admixture plots. For this, you need to first create a new folder (All_Qfiles) and move the Q files with "best" entropies from the LEA runs into it. 
sfiles <- list.files(path=("./All_Qfiles"), full.names=T)
slist <- readQ(files=sfiles)
plotQ(qlist=slist[2],imgtype = "pdf",
       height = 1.5, clustercol = c("#51692d", "#f1c039"), dpi = 1200, exportpath = "./")
 plotQ(qlist=slist[3],imgtype = "pdf",
       height = 1.5, clustercol = c("#51692d","#56ba32","#f1c039"), dpi = 1200, exportpath = "./")
 plotQ(qlist=slist[4],imgtype = "pdf",
       height = 1.5, clustercol = c("#51692d","#56ba32","#f1c039","#f37d21"), dpi = 1200, exportpath = "./")
 plotQ(qlist=slist[5],imgtype = "pdf",
       height = 1.5, clustercol = c("#f37d21","#51692d","#f1c039","#56ba32","#a63838"), dpi = 1200, exportpath = "./")
 plotQ(qlist=slist[6],imgtype = "pdf",
       height = 1.5, clustercol = c("#a63838","#f1c039","#ecbcab","#56ba32","#51692d","#f37d21"), dpi = 1200, exportpath = "./")
 plotQ(qlist=slist[7],imgtype = "pdf",
       height = 1.5, clustercol = c("#51692d","#a63838","#f1c039","#ecbcab","#caf291","#56ba32","#f37d21"), dpi = 1200, exportpath = "./")
 plotQ(qlist=slist[8],imgtype = "pdf",
       height = 1.5, clustercol = c("#3d87db","#a63838","#51692d","#f1c039","#ecbcab","#caf291","#56ba32","#f37d21"), dpi = 1200, exportpath = "./")
 plotQ(qlist=slist[9],imgtype = "pdf",
       height = 1.5, clustercol = c("#f1c039","#a63838","#3d87db","#ecbcab","#caf291","#56ba32","#0000FF","#51692d","#f37d21"), dpi = 1200, exportpath = "./")
 plotQ(qlist=slist[1],imgtype = "pdf",
       height = 1.5, clustercol = c("#dd00ff","#51692d","#caf291","#3d87db","#ecbcab","#a63838","#56ba32","#0000FF","#f1c039","#f37d21"), dpi = 1200, exportpath = "./") 
```

We used DivMigrate in R to assess the directionality and magnitude of gene flow between populations using Nei’s Gst (Nei, 1973) as a measure of population genetic differentiation. Here I have also provided the code to create the gene flow heatmap.
```
library(diveRsity)
library(corrplot)

divMigrate(infile="./H2_bialSNP_MAF05_geno_LD50502_388_rename_reorder.gen", outfile="Gst_5000boot", stat="gst", para=TRUE, plot_network=TRUE)
dev.off()

##before plotting the heatmap, I have manually edited the Gst output matrix of divMigrate, making sure the matrix is tab delimited and population names are correct.
data <- as.matrix(read.table("divMigrate_Gst.txt"))
corrplot(data, method="color", order = "alphabet", is.corr = FALSE)
```
