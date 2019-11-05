# pptc-pdx-ethnicity-inference

**This repository contains the data and methods used to infer approximate ethnic backgrounds for 252 patient-derived xenograft (PDX) models from the Pediatric Preclinical Testing Consortium (PPTC) using SNP array genotyping data. Published in Cell Reports: https://doi.org/10.1016/j.celrep.2019.09.071 (see Figure S2)**

**Authors:** Laura Egolf, Zalman Vaksman, Jo Lynne Rokita (2018)

Note that genotyping data derived from tumors are not ideal for inferring ethnicity, so these methods and results are meant to serve only as an approximation.

## Datasets used:
- PDX SNP array data (Illumina Final Report files from GenomeStudio): Deposited in https://figshare.com/account/home#/projects/38147
- HapMap 3 (release 2): Downloaded from ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/latest_phaseIII_ncbi_b36/plink_format/

## Scripts:
- *pdx_pca.sh:* Bash script used to convert the Illumina Final Report genotype files to PLINK format, merge these with HapMap 3, and run the PCA
- *plot_pdx_pca.R:* R script used to plot the first two PCs and assign an approximate ethnicity to each sample

## Supporting files:
- *filelist_InfiniumOmniExpress-24v1-2_A1.txt* and *filelist_humanomniexpress-24-v1-1-a.txt:* List of SNP array data files for the 254 samples (separated by chip type)
- *snps_to_exclude.txt:* List of 414 SNPs with problematic allele coding that caused errors in PLINK
- *ethnicity_coordinates_40kSNPs.txt:* Coordinates used for assigning samples to general ethnicity groups
- *2019-02-09-all-hist-colors.txt:* Hexadecimal color codes used for plotting samples according to tumor histotype
- *2019-07-25-clin.txt:* Clinical annotation for PDX samples

## Software used:
- R version 3.4.3
- PLINK 1.07
- PLINK 1.9

## Output files:
- PCA.plink.eigenval
- PCA.plink.eigenvec
- PCA.plink.eigenvec.var
- inferred_ethnicities_40kSNPs.txt
