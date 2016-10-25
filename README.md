gtCNV
=====
Genotype Copy Number Variation with Machine Learning.

## Getting Started
#### 1: Download Source Code :floppy_disk:
```
wget http://downloads.sourceforge.net/project/gtcnv/gtcnv-2.3.zip # gtcnv-2.3.tar.gz also available
unzip gtcnv-v2.3.zip
```
#### 2: Configure Environment
Run `configure.pl` to define install location and paths to FASTA assemblies
```
cd gtcnv-v2.3/
perl configure.pl # follow the instructions
```
#### 3: Compile from Source
```
python setup.py install # ignore numpy warnings
```
## Options
`gtCNV --help`

Flag | Description
--- | -------------
-i \| -in | Sample information input
-r \| -cnv | CNVs to genotype. BED or VCF
-c \| -cpu | Number of samples to run in parallel. Limited by available CPUs
-g \| -genome | Reference Genome Build [ hg19, hg38 ]. Default: hg19
-pcrfree | GC content normalization for PCR-free libraries
-s \| -seed | Random seed for genome shuffling. Used in preprocessing
-o \| -out | output
-pre | Preprocessing output directory. *Skips preprocessing*
-feats | Feature output directory. *Skips feature extraction*
#### Example
```
# genotype cnv
gtCNV -i CEU.in -r cnv.vcf -g hg38 -o CEU_cnv_genotypes.vcf

# genotype cnv2 skipping preprocessing
gtCNV -i CEU.in -r cnv2.vcf -g hg38 -o CEU_cnv2_genotypes.vcf -pre gtCNV_preprocessing/

# produce a VCF of one individual skipping feature extraction
head -n 1 CEU.in >sub.in
gtCNV -i sub.in -r cnv.vcf -g hg38 -o sub_CEU_cnv_genotypes.vcf -pre gtCNV_preprocessing/ -feats gtCNV_features/

# output is in gtCNV_genotypes/
ls gtCNV_gentypes/*
    CEU_cnv_genotypes.txt # Tab-delimited genotypes
    CEU_cnv_genotypes.vcf # VCF formatted genotypes
    ...
```
*Output VCF comes with gene annotations and other useful statistics*
## Inputs
#### Sample information < -i >
Tab-delimited file containing sample information. Gender can also be encoded as 1 for M and 2 for F

ID | BAM PATH |  VCF PATH | Gender [M/F]
--- | --- | --- | ---
NA12878 | /bam/NA12878_BWAMEM.bam | /vcf/NA12878_GATK_HC.vcf | F
HG00096 | /bam/HG00096_BWAMEM.bam | /vcf/HG00096_GATK_HC.vcf | M
#### CNVs to genotype < -r >
* BED format
  * Tab-delimited: first four columns
    1. Chromosome
    2. Start
    3. End
    4. Type: DEL | DUP
* VCF format
  * SVTYPE= DEL | DUP
  * Must have END=

## Usage
* gtCNV is designed for paired-end, short-read whole genomes
* Whole genome alignments from [1000 Genomes Project](http://www.1000genomes.org/) were used for training with validated genotypes from the phase 3 integrated structural variation call set [DOI:10.1038/nature15394](http://dx.doi.org/10.1038%2Fnature15394)
* Features for genotyping include coverage, discordant paired-ends, split-reads, and heterozygous allelic depth ratio.
   * BAM files must have supplementary alignment tags (SA).
   * SNV VCF must contain Allelic Depth. gtCNV can accomodate [GATK Haplotype Caller](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) and [FreeBayes](https://github.com/ekg/freebayes) VCFs.`
* gtCNV implements a bi-allelic model with a copy number range of 0-4
* Output is in VCF format
   * Median Phred-adjusted ALT likelihoods are reported in the QUAL column
   * Annotation for variants include genes, 1000 Genomes CNV,segmental duplications, RepeatMasker, and more
* Filters are set according to validated [*de novo* mutations](http://dx.doi.org/10.1016/j.ajhg.2016.02.018). Hence, they may be strict
   * `FT` format field is a sample filter. **1 indicates PASS**
   * `FILTER` column is PASS if at least 90% of the population passed sample-wise filters
* CNVs with outlier coverage features with an estimated autosome copy number >10 cannot be genotyped.

---

### Requirements:
* [python 2.7](https://www.python.org/)
  * [cython](https://github.com/cython/cython)
  * [numpy](http://www.numpy.org/)
  * [pandas](http://pandas.pydata.org/)
  * [pybedtools](https://daler.github.io/pybedtools/)
  * [pysam 0.9+](https://github.com/pysam-developers/pysam)

* [bedtools 2.25.0](https://github.com/arq5x/bedtools2/releases) or later

*gtCNV requires python 2.7*

*gtCNV has been tested on Linux and MacOS with [bioconda](https://bioconda.github.io/)*

---

## Credits

####Author:

* Danny Antaki
    * dantaki@ucsd.edu
    
#### Acknowledgements:
* William Brandler
* Jonathan Sebat
    * Sebat Lab http://sebatlab.ucsd.edu/index.php/software-data

-------    

### History
[gtCNV version 0.1](https://github.com/dantaki/gtCNV/tree/Version-0.1) used in Brander et al. *AJHG* 2016 ([DOI](http://dx.doi.org/10.1016/j.ajhg.2016.02.018) PMID:    27018473)

------

### License
MIT License

Copyright (c) 2016 Danny Antaki

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

##### Contact
dantaki@ucsd.edu
