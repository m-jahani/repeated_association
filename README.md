# input files:

### annuus_maf_0.03_Rvalue_NFFD

> **Spearsman correlation R coeficients for Annuus SNPs and NFFD climate variable** 
> - column1: species ID 
> - column2: variable ID.
> - column3: SNP ID in format of chromosomeID:Position.
> - column4: Corresponding window.
> - column5: Crrelation coeficient value.

### argophyllus_maf_0.03_Rvalue_NFFD

> **Spearsman correlation R coeficients for Argophyllus SNPs and NFFD climate variable**
> - column1: species ID
> - column2: variable ID
> - column3: SNP ID in format of chromosomeID:Position
> - column4: Corresponding window
> - column5: Crrelation coeficient value

### 5kbwindow_recombination

> **recombination rate file**
>  - column1: window ID
>  - column2: recombination rate

### Annuus_threshold.90

> **significant LD threshld for different window distance, base on null distribution cinstructed on 10000 random window withthe corresponding distance**
>  - column1: window distance (bp)
>  - column2: LD threshold

### Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt

> **Genetic Map**
> - column1: Chromosome
> - column2: Physical Position (bp)
> - column3: Genetic Position (cM)

### Annuus.tranche90.snp.env.90.bi.remappedHa412HO.beagle.vcf

> **standard VCF**

# Scripts:

### TOPCAN.R

> **identifies windows of the genome within each species that showed strong signatures of association**
>
> Runs as:
> 
> Rscript TOPCAN.R \\ # input file\
> annuus_maf_0.03_Rvalue_NFFD \\ # input file\
> argophyllus_maf_0.03_Rvalue_NFFD \\ # input file\
> R \\ # R: correlation BF: BayesFactor P: Pvalue\
> 5kbwindow_recombination \\ # input file\
> 0.99 \\ # cut of value for assocition test, 0.01 for top percentile or 10 of BayesFactor > 10\
> 0.9999 \\ # quantile threshold of the binomial expectation \
> ~ \\ # save directory \
> 150 # number of threads
>
> output:
>
>> annuus_maf_0.03_Rvalue_NFFD_argophyllus_maf_0.03_Rvalue_NFFD_topcandidate
>> **information on number of total snps and outlier snps for each window and whether the window is a top candidate**


### NULLW.R

> **identifies windows of Repeated Association (WRAs) between pairs of species**
**
>
> Runs as:
>
> Rscript NULLW.R \
> annuus_maf_0.03_Rvalue_NFFD \\ # input file\
> argophyllus_maf_0.03_Rvalue_NFFD \\ # input file\
> R \\ # R: correlation BF: BayesFactor P: Pvalue\
> 5kbwindow_recombination \\ # input file\
> annuus_maf_0.03_Rvalue_NFFD_argophyllus_maf_0.03_Rvalue_NFFD_topcandidate \\ #TOPCAN.R output
> 0.05 \\ # threshold of q-value in FDR test \ 
> ~  \\ # save directory \
> 150 # number of threads
>
> output:
>
>> Annuus_Argophyllus_NFFD_nullw_result_FDR_0.05





