#!/usr/bin/bash


plink10="~/plink1.07/plink-1.07-x86_64/plink"      # Directory path for PLINK 1.07
plink19="~/plink1.9/plink"     # Directory path for PLINK 1.9
hapmap="hapmap3_r2_b36_fwd.consensus.qc.poly" 


### Convert Illumina Final Report to PLINK format 

## For HumanOmni chip

# Create map file (use one sample from each chip type)
awk -F "\t" 'NR>11 {print $3 "\t" $2 "\t0\t" $4}' PPTC-IC-9320GCT-D-1.SNPFinalReport.txt > humanomniexpress-24-v1-1-a.map

# Create lgen and fam files
while read filename; do
    # lgen file - one line per SNP
    awk -F "\t" 'NR>11 {print "0\t" $1 "\t" $2 "\t" $9 "\t" $10}' $filename >> humanomniexpress-24-v1-1-a.lgen
        # Use Allele1 - Plus and Allele2 - Plus
    # fam file - one line per file
    awk -F "\t" 'NR==12 {print "0\t" $1 "\t0\t0\t0\t0"}' $filename >> humanomniexpress-24-v1-1-a.fam
done < filelist_humanomniexpress-24-v1-1-a.txt

# Generate ped/map from map/lgen/fam
$plink10 --noweb --missing-genotype - --output-missing-genotype 0 --lfile humanomniexpress-24-v1-1-a --recode --out humanomniexpress-24-v1-1-a.plink 

## For InfiniumOmni chip

# Create map file (use one sample from each chip type)
awk -F "\t" 'NR>11 {print $3 "\t" $2 "\t0\t" $4}' PPTC-AF0D-XTP1-B-1-0-D-1.SNPFinalReports.txt > InfiniumOmniExpress-24v1-2_A1.map

# Create lgen and fam files
while read filename; do
    # lgen file - one line per SNP
    awk -F "\t" 'NR>11 {print "0\t" $1 "\t" $2 "\t" $9 "\t" $10}' $filename >> InfiniumOmniExpress-24v1-2_A1.lgen
        # Use Allele1 - Plus and Allele2 - Plus
    # fam file - one line per file
    awk -F "\t" 'NR==12 {print "0\t" $1 "\t0\t0\t0\t0"}' $filename >> InfiniumOmniExpress-24v1-2_A1.fam
done < filelist_InfiniumOmniExpress-24v1-2_A1.txt

# Generate ped/map from map/lgen/fam
$plink10 --noweb --missing-genotype - --output-missing-genotype 0 --lfile InfiniumOmniExpress-24v1-2_A1 --recode --out InfiniumOmniExpress-24v1-2_A1.plink
    # Had to switch to PLINK 1.07 due to errors caused by "-" genotypes in 1.9

## Merge chip types

# Remove 414 SNPs with problematic allele coding; convert to bed/bim/fam
$plink19 --noweb --file humanomniexpress-24-v1-1-a.plink --exclude snps_to_exclude.txt --make-bed --out humanomniexpress-24-v1-1-a.excluded-snps.plink
$plink19 --noweb --file InfiniumOmniExpress-24v1-2_A1.plink --exclude snps_to_exclude.txt --make-bed --out InfiniumOmniExpress-24v1-2_A1.excluded-snps.plink

# Run merge command
$plink19 --noweb --bfile humanomniexpress-24-v1-1-a.excluded-snps.plink --bmerge InfiniumOmniExpress-24v1-2_A1.excluded-snps.plink --out pdx.plink


### SNP Quality Control

# Remove sex chromosomes, SNPs with MAF <1%, call rate <90% or Hardy-Weinberg equilibrium surpassing 0.00005
$plink19 --noweb --bfile pdx.plink --chr 1-22 --maf 0.1 --geno 0.1 --hwe 0.00005 --make-bed --out pdx.qc.plink


### Merge with HapMap and prune non-informative SNPs

# Extract the intersecting snps
$plink19 --noweb --bfile pdx.qc.plink --extract $hapmap.bim --make-bed --out pdx.extract.plink
$plink19 --noweb --bfile $hapmap --extract pdx.qc.plink.bim --make-bed --out hapmap.extract.plink

# Merge HapMap and PDX files
$plink19 --noweb --bfile pdx.extract.plink --bmerge hapmap.extract.plink --make-bed --out pdx_hapmap.plink

# Pruning step; gets rid of non-informative or low coverage snps 
$plink19 --noweb --bfile pdx_hapmap.plink --indep-pairwise 50 5 0.1 --out pdx_hapmap_ind-pair.plink

## Extract the the pruned set of snps
$plink19 --noweb --bfile pdx_hapmap.plink -extract pdx_hapmap_ind-pair.plink.prune.in --make-bed --out pdx_hapmap_prune.plink


### Run PCA 
$plink19 --bfile pdx_hapmap_prune.plink --pca 20 header var-wts --threads 8 --out PCA.plink


