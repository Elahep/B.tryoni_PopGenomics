
# Identifying genomic signal of local adaptation in the invasive Queensland fruit fly, _Bactrocera tryoni_

Here we will use the DarTseq data generated by <a href="https://www.nature.com/articles/s41598-020-67397-5" title="Popa‑Báez et al (2020)" >Popa‑Báez et al (2020)</a> from the native and expanded ranges of _B. tryoni_ to understand population genomic underpinnigs of invasion success. We use genome-wide scan to identify highly differentiated genomic variants (SNPs) associated with the invasive status of _B. tryoni_ populations. Additionally, we determined which environmental variables have a key role in shaping allelic trends across the distribution range of _B. tryoni_.


Codes in this repository have been used in:

Parvizi, E., Vaughan, A.L., Dhami, M.K. McGaughran, A. Genomic signals of local adaptation across climatically heterogenous habitats in an invasive tropical fruit fly (Bactrocera tryoni). Heredity (2023). https://doi.org/10.1038/s41437-023-00657-y


***************

## Table of contents
### 1-variant_calling:   
Contains scripts for running BWA MEM, BCFtools and PLINK

### 2-population_structure:
Contains scripts for PCA, sNMF, StAMPP, pophelper, and tests of IBD (EEMS and ConStruct)

### 3-BayPass
Contains script for running the C2-statistic analysis in BayPass

### 4-GOterms
Contains scripts for running Interproscan and topGO 

### 5-environmental_association
Contains scripts for processing bioclimatic variables + gradient forest analysis
