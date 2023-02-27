# Input files:

### annuus_maf_0.03_Rvalue_NFFD
```
Spearsman correlation R Coefficients for Annuus SNPs and NFFD climate variable
 - column1: species ID 
 - column2: variable ID.
 - column3: SNP ID in format of chromosomeID:Position.
 - column4: Corresponding window.
 - column5: Crrelation Coefficient value.
```

### argophyllus_maf_0.03_Rvalue_NFFD
```
Spearsman correlation R Coefficients for Argophyllus SNPs and NFFD climate variable
 - column1: species ID
 - column2: variable ID
 - column3: SNP ID in format of chromosomeID:Position
 - column4: Corresponding window
 - column5: Crrelation Coefficient value
 ```

### 5kbwindow_recombination
```
Recombination rate file
  - column1: window ID
  - column2: recombination rate
```
### Annuus_threshold.90
```
significant LD threshld for different window distance, base on null distribution cinstructed on 10000 random window withthe corresponding distance
 - column1: window distance (bp)
 - column2: LD threshold
```

### Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt
```
Genetic Map
- column1: Chromosome
- column2: Physical Position (bp)
- column3: Genetic Position (cM)
```
### Annuus.tranche90.snp.env.90.bi.remappedHa412HO.beagle.vcf
```
Standard VCF
```
# Scripts:

### TOPCAN.R
```
Identifies windows of the genome within each species that showed strong signatures of association

 Runs as:
 
 Rscript TOPCAN.R \\\
 annuus_maf_0.03_Rvalue_NFFD \\ # input file\
 argophyllus_maf_0.03_Rvalue_NFFD \\ # input file\
 R \\ # R: correlation BF: BayesFactor P: Pvalue\
 5kbwindow_recombination \\ # input file\
 0.99 \\ # cut of value for assocition test, 0.01 for top percentile or 10 of BayesFactor > 10\
 0.9999 \\ # quantile threshold of the binomial expectation \
 ~ \\ # save directory \
 150 # number of threads
```

### Output:
```
 annuus_maf_0.03_Rvalue_NFFD_argophyllus_maf_0.03_Rvalue_NFFD_topcandidate
 annuus_maf_0.03_Rvalue_NFFD_topcandidate
 argophyllus_maf_0.03_Rvalue_NFFD_topcandidate
```


### NULLW.R
```
Identifies windows of Repeated Association (WRAs) between pairs of species

Runs as:

Rscript NULLW.R \\\
annuus_maf_0.03_Rvalue_NFFD \\ # input file\
argophyllus_maf_0.03_Rvalue_NFFD \\ # input file\
R \\ # R: correlation BF: BayesFactor P: Pvalue\
5kbwindow_recombination \\ # input file\
annuus_maf_0.03_Rvalue_NFFD_argophyllus_maf_0.03_Rvalue_NFFD_topcandidate \\ #TOPCAN.R output \
0.05 \\ # threshold of q-value in FDR test \ 
~  \\ # save directory \
150 # number of threads
```

### Output:
```
Annuus_Argophyllus_NFFD_nullw_result_FDR_0.05
Argophyllus_Annuus_NFFD_nullw_result_FDR_0.05
Argophyllus_Annuus_NFFD_nullw_result
Annuus_Argophyllus_NFFD_nullw_result
```

### LDCLUSTER.R
```
Forms Clusters of Repeated Association (CRAs) base on LD and genetic distance between them

This script needs plink installed to run
 
Runs as:

Rscript LDCLUSTER.R \\\
Annuus_Argophyllus_NFFD_nullw_result_FDR_0.05 \\ #NULLW.R output \
annuus_maf_0.03_Rvalue_NFFD \\ # input file\
Annuus.tranche90.snp.env.90.bi.remappedHa412HO.beagle.vcf \\ # input file\
Annuus_threshold.90 \\ # input file\
Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt \\ # input file\
5 \\ # genetic distance threshold \
~  \\ # save directory \
25 # number of threads
```

### Output:
```
Annuus_Argophyllus_NFFD_FDR_0.05.clustering_result
```
