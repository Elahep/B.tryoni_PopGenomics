## Functional annotation of candidate genes
We investigated genes associated with BayPass candidate SNPs using BEDtools and the UCSC Table Browser (https://genome.ucsc.edu/).

We first need to extract the outlier SNPs from our full-SNP VCF file using VCFtools. For that, we need a .txt file containing the list of these outliers with the scaffold and SNP position information for each SNP. 

```
module load VCFtools

#this procedure will be done on the candidates from all three BayPass runs on the three main invasive clusters
vcftools --vcf bialSNP_MAF_geno_LD.vcf --snps nativevsAlice_candidates.txt --recode --recode-INFO-all
mv out.recode.vcf nativevsAlice_outliers.vcf
```

We now need to convert .vcf to .bed file. "bed" format is a text file format that represents genomic data (scaffold or chromosome names and SNP coordinates) in the form of columns separated by spaces or tabs.

```
module load BEDOPS
vcf2bed < nativevsAlice_outliers.vcf > nativevsAlice_outliers.bed
```

The .bed file can be utilized to obtain flanking sequence data using the reference genome. BEDTools can perform this task, but it requires a designated "genome file." To create the BEDTools genome file from the reference genome fasta index file, I employed the following bash script.
```
awk -v OFS='\t' {'print $1,$2'} Btry.fa.fai > Btry.txt   #Btry.txt is the special bedtools genome format (to get fai.fa file use: samtools faidx Btry.fa)
```

Now all the necessary input files are ready to use BEDTools to extract flanking regions (e.g. 10,000 base pairs) around the outlier SNPs:
```
module load BEDTools
bedtools flank -i nativevsAlice_outliers.bed -g Btry.txt -b 10000 | cut -f1-3 > nativevsAlice_outliers_10kb_flank.bed
```

The resulting .bed file contains sequences at flanking 10 kb of candidate (or outlier) SNPs. I then uploaded these sequences in the UCSC Table Browser (https://genome.ucsc.edu/) to retrieve transcript IDs associated with these sequences.


![image](https://user-images.githubusercontent.com/13001264/217693489-964b35aa-6460-447c-9161-b2d80260b7b9.png)

The output of the Table Browser tool contains information about the assembly ID used for extracting exons, the transcript ID and range (scaffold and coordinates) for each exon. We only need the transcript ID for each exon. I used simple bash scripts to extract only the lines that contain the transcript IDs and use notepad+ to manually edit the final file and remove any extra character/information from the list of transcript IDs.


## GO terms analysis
We performed a Gene Ontology (GO) enrichment test to assess the gene functions related to putatively adaptive loci. We first need to get the genomic GO annotation for all loci (candidate and non-candidate) and create a complete set of GO terms as a reference for the enrichment analysis. This task can be done using Interproscan, but before that we need to get the transcript ID for all loci (candidate and non-candidate) and we will do that using snpEff:

```
module load Miniconda3
snpEff -c snpEff.config mygenome ./bialSNP_MAF_geno_LD.vcf > Qfly_FullData.anno.vcf
cat snpEff_genes.txt | grep -v "#" | cut -f3 | uniq > Qfly_FullData_transIDfromsnpeff.txt
```
This list of all transcript IDs should be uploaded to NCBI Batch Entrez tool to search against the nucleotide database and obtain the DNA sequences associated with the transcript IDs. The DNA sequences (in fasta format) can now be used in InterProScan to obtain the associated GO terms:

```
module load InterProScan
interproscan.sh -i Qfly_FullData_batchentrez.fasta -t n --goterms
```

Using the above command, InterProScan will output multiple files (such as .gff3, .json and .tsv). We only need the .tsv file to extract GO terms for different transcript IDs. When we look at the .tsv output file, we can see that InterProScan has not found GO terms for some transcript IDs. Also, there can be multiple GO terms suggested for one transcript ID. Lets first remove the lines that do not contain any GO terms. In bash:

```
grep -w "GO" Qfly_FullData_batchentrez.fasta.tsv | cut -f1,14 > Qfly_FullData_GOlist.txt 
```

The resulting file contains two columns: in the first column we have transcript IDs of all of our genes (candidate and non-candidate) and in the second column we have GO terms associated with each transcript ID. Note that multiple GO terms in each row are separated by "|".

Now all the necessary files are ready to be used in GO terms analysis in R. Note that we ran this code for candidate SNPs from native vs Alice Springs and native vs Pacific Islands. Candidates from native vs Expanded ranges were not associated with any exons so we couldn't use them in the GO terms analysis.

```
library(topGO)


##import tab delimited file of all genes and their GO terms
geneID2GO <- readMappings("./InterProScan/Qfly_FullData_GOlist.txt")  
geneID2GO$XM_014420108.2  #check the GO terms for some of the transcript IDs
str(head(geneID2GO))

geneNames <- names(geneID2GO)
head(geneNames)

##import transcript IDs of outlier SNPs:
interesting_genes = read.table("./nativevsAlice_outliers_transcriptIDs.txt")
myInterestingGenes <- as.vector(interesting_genes$V1) ##import the vector of your genes of interest. This should be just the ID of the genes.
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)

##create a topGo object for the "Biological Processes" terms (BP):
GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO) 

##get the list of significant genes
sig_genes = sigGenes(GOdata_BP) 

## return the GO graph
graph(GOdata_BP)

##GO enrichment test
resultFisher <- runTest(GOdata_BP, algorithm="weight01", statistic="fisher")  ##resultFisher is a TopGoresult objetc. You can see p values in this object using the "score" function: 
resultFisher  #this shows how many GO terms are significant
allRes <- GenTable(GOdata_BP, raw.p.value = resultFisher, classicFisher = resultFisher, ranksOf = "classicFisher", Fis = resultFisher, topNodes = length(resultFisher@score)) 
allRes
write.table(allRes, "nativevsAlice_GOresults_BP.txt", sep = "\t")  ##save the results for the BP terms.


##now we can do GO enrichment term analysis using different ontology. For example, we will analyze for the "Molecular Function" terms here (MF):
GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_MF
resultFisher_MF <- runTest(GOdata_MF, algorithm="weight01", statistic="fisher")
resultFisher_MF  #this shows how many GO terms are significant
allRes_MF <- GenTable(GOdata_MF, raw.p.value = resultFisher_MF, classicFisher = resultFisher_MF, ranksOf = "classicFisher", Fis = resultFisher_MF, topNodes = length(resultFisher_MF@score)) 
allRes_MF
write.table(allRes_MF, "nativevsAlice_GOresults_MF.txt", sep = "\t")

##Similarly, we will do the analysis for the "Cellular Component" terms as well (CC):
GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_CC
resultFisher_CC <- runTest(GOdata_CC, algorithm="weight01", statistic="fisher")
resultFisher_CC  #this shows how many GO terms are significant
allRes_CC <- GenTable(GOdata_CC, raw.p.value = resultFisher_CC, classicFisher = resultFisher_CC, ranksOf = "classicFisher", Fis = resultFisher_CC, topNodes = length(resultFisher_CC@score)) 
allRes_CC
write.table(allRes_CC, "nativevsAlice_GOresults_CC.txt", sep = "\t")
```
