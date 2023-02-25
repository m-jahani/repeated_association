# input files:

## annuus_maf_0.03_Rvalue_NFFD

> **Spearsman correlation R coeficients for Annuus SNPs and NFFD climate variable** 
> - column1: species ID 
> - column2: variable ID.
> - column3: SNP ID in format of chromosomeID:Position.
> - column4: Corresponding window.
> - column5: Crrelation coeficient value.

## argophyllus_maf_0.03_Rvalue_NFFD

> **Spearsman correlation R coeficients for Argophyllus SNPs and NFFD climate variable**
> - column1: species ID
> - column2: variable ID
> - column3: SNP ID in format of chromosomeID:Position
> - column4: Corresponding window
> - column5: Crrelation coeficient value

## 5kbwindow_recombination

> **recombination rate file**
>  - column1: window ID
>  - column2: recombination rate

## Annuus_threshold.90

> **significant LD threshld for different window distance, base on null distribution cinstructed on 10000 random window withthe corresponding distance**
>  - column1: window distance (bp)
>  - column2: LD threshold

## Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt

> **Genetic Map**
> - column1: Chromosome
> - column2: Physical Position (bp)
> - column3: Genetic Position (cM)

## Annuus.tranche90.snp.env.90.bi.remappedHa412HO.beagle.vcf

> **standard VCF"

# Scripts:

## TOPCAN.R

> **identifies windows of the genome within each species that showed strong signatures of association**
>
> Runs as:
> 
> Rscript TOPCAN.R \\
> annuus_maf_0.03_Rvalue_NFFD /\
> argophyllus_maf_0.03_Rvalue_NFFD /\
> R /\
> 5kbwindow_recombination /\
> 0.99 /\
> 0.9999 /\
> ~ /\
> 150


