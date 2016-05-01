gtCNV
=====
Genotyping Copy Nuber Variation with Machine Learning. 

*A resource for human whole genome next-generation sequencing libraries.* 

## Getting Started

clone the repository from github and run :grinning:

``` 
$ git clone git@github.com:dantaki/gtCNV.git
```
### Prerequisites 

gtCNV requires **python 2.7**i :snake:

gtCNV has been tested on Linux and MacOS with [bioconda](https://bioconda.github.io/)

**Required python libraries**

* numpy
* pandas
* pybedtools
* pysam 
* scikit-learn

### Installation

No installation required. Ensure prerequisites are met and run :grinning:

```
$ python gtCNV -i sample_info.in -b CNV.bed [ OPTIONAL -o gtCNV_genotypes.vcf -g hg19 --pre gtCNV_preprocessing/ ]
```
### Inputs

#### 1. Sample information < -i >

##### Must be tab-delimited

ID | BAM PATH | Gender [M/F]
--- | --- | --- 
NA12878 | /home/usr/bam/NA12878_BWAMEM.bam | F
HG00096 | /home/usr/bam/HG00096_BWAMEM.bam | M

##### :heavy_exclamation_mark: **BAM files must be BWA-MEM aligned** :heavy_exclamation_mark:


#### 2. BED file < -b > 

##### Must be tab-delimited

CHROM | START | END | TYPE [DEL/DUP]
--- | --- | --- | --- 
chr1 | 1000 | 2000 | DEL 
chr2 | 3500 | 4500 | DUP
chr2 | 5000 | 5300 | DEL_ALU
chr3 | 1000 | 2000 | DUP_mCNV

#####:heavy_exclamation_mark: **CNV type must contain either 'DEL' or 'DUP'** :heavy_exclamation_mark:

## Options

Display options

```
$ python gtCNV --help
```

Flag | Description
--- | ------------
-i | Sample information input
-b | BED file of CNVs
-c | Number of samples to run in parallel. Limited by available CPUs
-g | Reference Genome Build [ hg19, hg38 ]. Default is hg19
-s | random seed for genomic shuffling. Used in preprocessing
-o | VCF output 
--pre | Preprocessing output directory. Skips preprocessing if completed
--feats | Feature output directory. Skips feature extraction if completed

## Useage 

gtCNV is designed for human whole genome next-generation sequencing libraries. Given a list of CNV positions, gtCNV returns an annotated VCF with predicted copy number states.


## Credits

####Author:

* Danny Antaki
    * dantaki@ucsd.edu
* William Brandler

## History

[gtCNV version 0.1](https://github.com/dantaki/gtCNV/tree/Version-0.1) used in Brander et al. *AJHG* 2016 ([DOI](http://dx.doi.org/10.1016/j.ajhg.2016.02.018) PMID:    27018473)

## Credits

#### Acknowlegements:

* Prateek Tandon 
* Jonathan Sebat
    * Sebat Lab http://sebatlab.ucsd.edu/index.php/software-data

## More information

gtCNV was trained on 27 high coverage genomes with validated genotypes from the phase 3 intgrated structural variation release (doi:10.1038/nature15394;PMID:     26432246). 
## License

gtCNV
    Copyright (C) 2015  Danny Antaki

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
