# Hawk

## Overview
**Hawk** is a lightweight tool designed to facilitate the selection of elite plants and animals based on known Quantitative Trait Loci (QTLs). It uses our [custom PanMap file](https://www.biorxiv.org/content/10.1101/2025.03.28.645917v1.full) and QTL BED file to identify lines carrying favorable alleles at key loci, enabling marker-assisted selection for precision breeding.

## Requirements
### [KhufuEnv](https://github.com/w-korani/KhufuEnv)

## Inputs
### 1. PanMap file
Hawk is designed to be used with our [custom PanMap file](https://www.biorxiv.org/content/10.1101/2025.03.28.645917v1.full). This file format is similar in concept to a HapMap, but it displays variants against a pangenome graph rather than a linear reference genome. 

Alleles in the PanMap are encoded using integers and the actual allele sequences are stored in a corresponding FASTA file. 

To convert other input formats to our PanMap format and to generate the corresponding FASTA file, we recommend using our [KhufuEnv](https://github.com/w-korani/KhufuEnv) tools, `vcf2panmap` or `hapmap2panmap` and `panmapGetFasta`.

**Khufu PanMap Example**
![F2 large](https://github.com/user-attachments/assets/aa84fc5a-816d-43fa-937a-2bb500dc90a8)
*Image Source: [Wright, Hallie C, Catherine E. M. Davis, Josh Clevenger, and Walid Korani. “KhufuEnv, an Auxiliary Toolkit for Building Computational Pipelines for Plant and Animal Breeding.” bioRxiv, January 1, 2025, 2025.03.28.645917.](https://doi.org/10.1101/2025.03.28.645917)*

### 2. QTL BED file
This is an extended BED file format. It follows the standard 3-column BED structure, with an additional fourth column which includes the donor parent and associated trait. 

| Column # | Description                                      | Example           |
|----------|--------------------------------------------------|-------------------|
|1         | Chromosome name                                  | `TRv2Chr.12`      |
|2         | Start position                                   | `1500000`         |
|3         | End position                                     | `2000000`         |
|4         | Comma-delimited string: donor_parent,trait_name  | `Ascasubi,smut`   |


**Note:** For traits with multiple QTL regions, there must be a unique string in the fourth column as in the example below.
```bash
TRv2Chr.12  1500000    2000000    Ascasubi,smut
TRv2Chr.15  40000000   150000000  Lariat,SclerBlight
TRv2Chr.19  152000000  155000000  Lariat,HO1
TRv2Chr.19  200000000  205000000  Lariat,HO2
```
## Outputs
### 1. Similarity Summary Table

Tab-delimited file containing similarity scores for each sample and QTL information. The first line is a header starting with `##` listing the traits used, `trait1,trait2,trait3`.

| Column  | Description                                                               | Example           |
|----------|---------------------------------------------------------------------------|-------------------|
|1         | **Line ID**                                                                   | `x167.196.02`      |
|2         | **Average similarity:** Calculated by averaging of all the trait similarity scores listed in column 3. Reflects overall similarity to the donor parent across all traits listed in the header of the file. | `0.71`      |
|3         | **Similarity for each trait:** Comma-delimited list of similarity scores for each trait (in the order shown in the header line). A value of `NA` may appear if the similarity score could not be calculated due to the number of variants failing to meet the minimum depth set.                                                              | `0.86,0.57,0.71`      |
|4         | **Number of variants per trait:** Comma-delimited list of the number of variants found per trait region (in the order shown in the header line).                                                        | `14,453,7`      |
|5         | **Range of the trait region:** Start and stop positions of the trait regions (in the order shown in the header line).                                                            | `1594957-1964360,40138292-149634268,152278500-153931210`      |
|6         | **Length of each region:** Length of each trait region in bp (in the order shown in the header line).                                                          | `369403,109495976,1652710`      |
|7         | **Background parent similarity and variant count:** Similarity score to the background parent and the number of variants after dropping all QTL.                                                          | `0.65;2602`      |

**Example**
```bash
##           smut,SclerBlight,HO                                                                                                                      
ID           AverageSimilarity    Similarity      NumberOfVariants  Range                                                   len                       Lariat
x167.196.02  0.71                 0.86,0.57,0.71  14,453,7          1594957-1964360,40138292-149634268,152278500-153931210  369403,109495976,1652710  0.65;2602
x132.2.02    0.69                 1.00,0.37,NA    10,210,4          1823422-1957028,40138292-148149473,152278500-154784694  133606,108011181,2506194  0.60;1224
x168.211.01  0.66                 0.94,0.55,0.50  16,311,6          1673886-1964171,40138292-148716487,152278500-154784694  290285,108578195,2506194  0.62;1850
x135.60.01   0.60                 0.92,0.49,0.40  13,380,5          1634190-1964171,40132422-149581379,152366038-154784694  329981,109448957,2418656  0.62;1783
x137.83.01   0.58                 0.57,0.56,0.60  14,593,5          1594957-1964360,40132422-149634268,152278500-154784694  369403,109501846,2506194  0.69;3027
x126.8       0.56                 1.00,0.48,0.20  10,296,5          1767392-1957028,46505718-148716487,152366038-154784694  189636,102210769,2418656  0.62;1442
x169.214.01  0.56                 0.50,0.62,0.57  6,342,7           1594957-1956545,40132422-148747441,152278500-154784694  361588,108615019,2506194  0.65;2152
x169.215.01  0.56                 0.64,0.65,0.50  11,510,6          1634190-1964360,40934005-149634268,152366038-154784694  330170,108700263,2418656  0.67;2927
x132.14.01   0.52                 0.53,0.51,NA    15,279,3          1673886-1964171,42380806-149581379,152366038-154784694  290285,107200573,2418656  0.65;1632
x160.133.01  0.36                 0.20,0.52,NA    10,365,4          1822625-1964171,40223802-148800184,152366038-153931210  141546,108576382,1565172  0.59;1774
```

### 2. Graph (optional)

A dot plot graph may be generated for visualization by setting the `--graph` option to `1`. 


## Usage

### Example Usage
```bash
hawk.sh -panmap test.panmap -bed test.bed -graph 1 -MinDep 5 -dominance 1 -o test.hawk -background Lariat
```

### Options

| Options/Flags |   Description                                                                |
| -----------------| -----------------------------------------------------------------------------|
| `-h`           | Display the help menu.                                                       |
| `-t`           | Number of threads to use.                                                    |
| `-o`           | Output file name.                                                            |
| `-u`          | Set to `1` to require donor allele uniqueness or `0` to allow consideration of all donor alleles. Default is `1`.                                             |
| `--bed`        | Input QTL BED file.                                                          |
| `--panmap`     | Input PanMap file.                                                           |
| `--background`     | Comma-delimited list of background parents. ex: `parent1,parent2`. Default is no background parents. Set to `All` to use all parents.          |
| `--dominance`  | Dominance of the allele (float). Set to `NA` to exclude heterozygotes. Default dominance is `1`.                       |
| `--MinDep`     | Minimum number of variants required to calculate similarity. Default used is `5`.                  |
| `--graph`      | Set to `1` to generate graph, or `0` to skip. Default is `0`.                 |
| `--clean`      | Set to `1` to clean intermediate files, or `0` to keep. Default is `1`.       |

## Citation
