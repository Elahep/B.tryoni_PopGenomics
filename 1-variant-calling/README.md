## Calling and filtering SNPs

As the raw DarTseq data provided by <a href="https://www.nature.com/articles/s41598-020-67397-5" title="Popa‑Báez et al (2020)" >Popa‑Báez et al (2020)</a> was already trimmed and quality-controlled, we have used them directly for SNP calling.


We will start by downloading the latest whole genome assembly of _D. tryoni_ from NCBI (BioProject number: PRJNA560467) and using ```tar -xf``` to extarct .tar file and get the genome assembly fasta file.

## Aligning raw reads against reference genome
We use BWA to do reference genome alignmnet. First, we need to index the reference genome to help the aligner work more efficiently by finding all exact matches to a sequenced read using a single lookup in the given data structure (see <a href="https://www.frontiersin.org/articles/10.3389/fpls.2021.657240/full#:~:text=Indexing%20the%20reference%20genome%20allows,to%20finish%20for%20each%20read." title="Musich et al (2021)" >Musich et al (2021)</a> for more details):
```
module load BWA
bwa index Btry.fa
```

Now the actual alignment can be done using BWA-MEM default parameters and the M and R flags to mark secondary reads and add read groups, respectively:
```
for i in {4..301};
do
bwa mem -t 10 ./ref_genome/Btry.fa -M -R "@RG\tID:Btry_${i}\tLB:Btry_${i}\tPL:ILLUMINA\tPM:HISEQ\tSM:Btry_${i}" Btry_${i}.fastq.gz | gzip > Btry_${i}.sam.gz 
done
```

## SNP calling
To call variants, we will use mpileup and call modules of BCFtools with the default calling method, outputting only variant positions, and including variants with minimum mapping quality and minimum base quality scores of 20.
```
module load BCFtools
bcftools mpileup -Ou -f ./ref_genome/Btry.fa --bam-list bam_list_all.txt --min-MQ 20 --min-BQ 20 | bcftools call -mv -Ob -o calls_all_samples.bcf
```

And further filtering SNPs and retain only high-quality, uncorrelated SNPs using BCFtools and PLINK2:
```
module load BCFtools
module load PLINK
bcftools view -m2 -M2 -v snps -O b -o bialSNP.bcf calls_all_samples.bcf  ##removing all indels and non-biallelic SNPs
bcftools filter -i 'MAF > 0.05' -O v -o bialSNP_MAF.vcf bialSNP.bcf ##applying MAF threshold of 0.05
plink2 --vcf bialSNP_MAF.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --geno 0.1 --recode vcf --out bialSNP_MAF_geno
plink2 --vcf bialSNP_MAF_geno.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --indep-pairwise 50 5 0.2
plink2 --vcf bialSNP_MAF_geno.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --exclude plink.prune.out --recode vcf --out bialSNP_MAF_geno_LD
```

