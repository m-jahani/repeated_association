# Repeatability of Adaptation in Sunflowers




## Input files:

### annuus_maf_0.03_Rvalue_NFFD
```
Spearsman correlation R Coefficients for Annuus SNPs and NFFD climate variable

 - column1: species ID 
 - column2: variable ID.
 - column3: SNP ID in format of chromosomeID:Position.
 - column4: Corresponding window.
 - column5: Correlation Coefficient value.
```

### argophyllus_maf_0.03_Rvalue_NFFD
```
Spearsman correlation R Coefficients for Argophyllus SNPs and NFFD climate variable

 - column1: species ID
 - column2: variable ID
 - column3: SNP ID in format of chromosomeID:Position
 - column4: Corresponding window
 - column5: Correlation Coefficient value
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

### HAPLOBLOCKS
```
The position of haploblocks in the sunflower genome

- column1: Chromosome
- column2: Start Position 
- column3: End Position 
- column4: Species
- column5: Haploblock ID
```

### Ha412.genome

```
length of each chromosome in Ha412.V2 reference genome

- column1: Chromosome
- column2: length (bp)
```

## Scripts:

### TOPCAN.R

Detects genomic windows within each species that exhibit strong association signals.

### Runs as:
 ```
 Rscript TOPCAN.R \
 annuus_maf_0.03_Rvalue_NFFD \ # input file
 argophyllus_maf_0.03_Rvalue_NFFD \ # input file
 R \ # R: correlation BF: BayesFactor P: Pvalue
 5kbwindow_recombination \ # input file
 0.95 \ # cut of value for assocition test, 0.01 for top percentile or 10 of BayesFactor > 10
 0.9999 \ # quantile threshold of the binomial expectation 
 ~ \ # save directory 
 150 # number of threads
```

### Output:
```
 annuus_maf_0.03_Rvalue_NFFD_argophyllus_maf_0.03_Rvalue_NFFD_topcandidate
 annuus_maf_0.03_Rvalue_NFFD_topcandidate
 argophyllus_maf_0.03_Rvalue_NFFD_topcandidate
```


### NULLW.R
Detects windows of Repeated Association (WRAs) between pairs of species.

#### Runs as:
```
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
Forms Clusters of Repeated Association (CRAs) base on LD and genetic distance between them

#### Runs as:
needs plink installed to run
```
Rscript LDCLUSTER.R \
Annuus_Argophyllus_NFFD_nullw_result_FDR_0.05 \ #NULLW.R output \
annuus_maf_0.03_Rvalue_NFFD \ # input file
Annuus.tranche90.snp.env.90.bi.remappedHa412HO.beagle.vcf \ # input file
Annuus_threshold.90 \ # input file
Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt \ # input file
1 \ # genetic distance threshold 
~  \ # save directory 
25 # number of threads
```
Rscript LDCLUSTER.R \
Argophyllus_Annuus_NFFD_nullw_result_FDR_0.05 \ #NULLW.R output \
argophyllus_maf_0.03_Rvalue_NFFD \ # input file
Annuus.tranche90.snp.env.90.bi.remappedHa412HO.beagle.vcf \ # input file
Annuus_threshold.90 \ # input file
Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt \ # input file
1 \ # genetic distance threshold 
~  \ # save directory 
25 # number of threads
```

### Output:
```
Annuus_Argophyllus_NFFD_FDR_0.05.clustering_result

Argophyllus_Annuus_NFFD_FDR_0.05.clustering_result
```

### Merge LDCLUSTER.R results

```
cat Annuus_Argophyllus_NFFD_FDR_0.05.clustering_result Argophyllus_Annuus_NFFD_FDR_0.05.clustering_result > Annuus_Argophyllus_NFFD_FDR_0.05.clustering_result_merged
```

### OVERLAP.R
Identifies haploblocks that overlap with clusters of repeated association.

#### Runs as:
```
Rscript OVERLAP.R \
HAPLOBLOCKS \  # input file
Annuus_Argophyllus_NFFD_FDR_0.05.clustering_result_merged \  # merged repeated association clusters
Annuus_Argophyllus \  # IDs of two species being compared (seprated by underscore)
~ \  # save directory
100 # number of threads
```

### Output:
```
Annuus_Argophyllus_NFFD_FDR_0.05.clustering_result_merged_overalp_haploblocks
```

### OVERLAP.R
Evaluates whether the length of overlap and the number of overlapped repeated clusters with haploblocks significantly deviate from the null expectation.

### Runs as:
```
Rscript OVERLAP_PVALUE.R \
HAPLOBLOCKS \  # input file
Annuus_Argophyllus_NFFD_FDR_0.05.clustering_result_merged_overalp_haploblocks \ #CRA and haploblock ovrlap result from OVERLAP.R script
Annuus_Argophyllus \ # IDs of two species being compared (seprated by underscore)
Ha412.genome \ # input file
repeated_association \
100 # number of threads
```

### Output:
```
Annuus_Argophyllus_NFFD_FDR_0.05.clustering_result_merged_overalp_haploblocks_Pvalue
```
