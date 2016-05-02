gtCNV
=====
Genotyping Copy Nuber Variation with Machine Learning. 

*A resource for whole genome next-generation sequencing libraries.* 

## Getting Started

Download from github
``` 
git clone git@github.com:dantaki/gtCNV.git
```
Download resource files. Available compressed files include .tar.gz and .zip 

:floppy_disk: .tar.gz
```
cd gtCNV/
wget https://www.dropbox.com/s/holuoc1v35ecmuo/gtCNV_v1.0_resources.tar.gz
tar zxvf gtCNV_v1.0_resources.tar.gz
```
:floppy_disk: .zip
```
cd gtCNV/ 
wget https://www.dropbox.com/s/k6r3c8rynur2hjo/gtCNV_v1.0_resources.zip
unzip gtCNV_v1.0_resources.zip
```
### Prerequisites 

gtCNV requires **python 2.7** :snake:

gtCNV has been tested on Linux and MacOS with [bioconda](https://bioconda.github.io/)

**Required python libraries**

* numpy
* pandas
* pybedtools 
* pysam 0.9.0 
* scikit-learn


:wrench: pybedtools requires bedtools. **gtCNV supports bedtools version [2.25.0](https://github.com/arq5x/bedtools2/releases) or later**

### Installation

Download resource files and the correct versions of bedtools and pysam. 

Now run! :grinning:

```
python gtCNV -i sample_info.in -b CNV.bed [ OPTIONAL -o gtCNV_genotypes.vcf -g hg19 --pre gtCNV_preprocessing/ ]
```
### Inputs

#### 1. Sample information < -i >

##### Must be tab-delimited. No headers

##### :heavy_exclamation_mark: **BAM files must be BWA-MEM aligned** :heavy_exclamation_mark:

ID | BAM PATH | Gender [M/F]
--- | --- | --- 
NA12878 | /home/usr/bam/NA12878_BWAMEM.bam | F
HG00096 | /home/usr/bam/HG00096_BWAMEM.bam | M

#### 2. BED file < -b > 

##### Must be tab-delimited. No headers

#####:heavy_exclamation_mark: **CNV type must contain either 'DEL' or 'DUP'** :heavy_exclamation_mark:

CHROM | START | END | TYPE [DEL/DUP]
--- | --- | --- | --- 
chr1 | 1000 | 2000 | DEL 
chr2 | 3500 | 4500 | DUP
chr2 | 5000 | 5300 | DEL_ALU
chr3 | 1000 | 2000 | DUP_mCNV

## Options

Display options

```
python gtCNV --help
```

Flag | Description
--- | ------------
-i | Sample information input
-b | BED file of CNVs
-c | Number of samples to run in parallel. Limited by available CPUs
-g | Reference Genome Build [ hg19, hg38 ]. Default is hg19
-s | Random seed for genome shuffling. Used in preprocessing
-o | VCF output 
--pre | Preprocessing output directory. Skips preprocessing if completed
--feats | Feature output directory. Skips feature extraction if completed

## Tutorial

Refer to README.md in tutorial/ directory for help

```
#After downloading and decompressing the resource files
cd gtCNV/
python gtCNV -i tutorial/tutorial.in -b tutorial/tutorial.bed -o tutorial_genotypes.vcf
```

## Usage 

* gtCNV is designed for human whole genome next-generation sequencing libraries. Given a list of CNV positions, gtCNV returns an annotated VCF with predicted copy number states.

* The training set included 27 high coverage genomes for deletions and 2,503 low coverage genomes for duplications from the [1000 Genomes Project](http://www.1000genomes.org/).Validated genotypes were obtained from the [phase 3 integrated structural variation call set](ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/integrated_sv_map/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz).([DOI:10.1038/nature15394](http://dx.doi.org/10.1038%2Fnature15394); PMID:    26432246). 

   * A copy of the raw training is located here ` gtCNV/resources/training_sets/gtCNV_raw_training-sets.zip`

* gtCNV performs a preprocessing step before genotyping; preprocessing output is located in ` gtCNV/gtCNV_preprocessing/` directory. 
   * You may wish to run gtCNV on new positions using the same samples. 
   * Pass the preprocessing directory in the command to skip this step 
      * `python gtCNV -i sample_info.in -b new_cnv.bed --pre gtCNV_preprocessing/`

* Features for genotyping include coverage, discordant paired-ends, and split reads. 
   * BAM files must be BWA-MEM aligned to annotate split-reads. 
   * Raw features are located in the ` gtCNV/gtCNV_genotypes/` directory. 

* CNVs with high coverages (normalized coverage >5 /estimated autosome copy number >10) are omitted. Such loci genotype poorly and interfere with the SVM model. 

* Output is in VCF format. 
   * Median Phred-adjusted non-reference likelihoods are reported in the QUAL column
   * Positions are annotated based on their overlap to genes, repeats, and 1000 Genomes phase 3 CNV

* Suggested filters in the VCF were determined using [svtoolkit intensity rank sum annotator](http://gatkforums.broadinstitute.org/gatk/discussion/2715/documentation-for-intensityranksum-annotator)
   * deletions: PASS at >= 12 median non-reference likelihood 
      * FDR ~ 1% at 5% allele frequency
   * duplications: PASS at >= 10 median non-reference likelihood
      * FDR ~ 3% at 5% allele frequency

* **NOTE:** non-allelic homologous recombination derived duplications are problematic for gtCNV
   * We recommend a more conservative median non-reference cutoff
      * NAHR duplications: PASS at >= 6 median non-reference likelihood
   * gtCNV Version 2.0 will address this issue. 

## Credits

####Author:

* Danny Antaki
    * dantaki@ucsd.edu
* William Brandler

## History

[gtCNV version 0.1](https://github.com/dantaki/gtCNV/tree/Version-0.1) used in Brander et al. *AJHG* 2016 ([DOI](http://dx.doi.org/10.1016/j.ajhg.2016.02.018) PMID:    27018473)

#### Acknowledgements: 

* Prateek Tandon 
* Jonathan Sebat
    * Sebat Lab http://sebatlab.ucsd.edu/index.php/software-data

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

Contact
------
dantaki@ucsd.edu 
