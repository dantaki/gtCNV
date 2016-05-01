gtCNV
=====
Genotyping Copy Nuber Variation by Machine Learning

## Installation

clone the repository from github and run :grinning:

gtCNV implements these python libraries

* numpy
* pandas
* pybedtools
* pysam <version 0.8.4>
* scikit-learn

## Usage

gtCNV is designed for human whole genome next-generation sequencing libraries. Given a list of CNV positions, gtCNV will return an annotated VCF with predicted copy number states.

## Inputs

#### Sample information < -i >

Must be tab-delimited

ID | BAM PATH | Gender [M/F]
--- | --- | --- 
NA12878 | /home/usr/bam/NA12878_BWAMEM.bam | F
HG00096 | /home/usr/bam/HG00096_BWAMEM.bam | M

**BAM files must be BWA-MEM aligned**

#### BED file 

Must be tab-delimited

CHROM | START | END | TYPE [DEL/DUP]
--- | --- | --- | --- 
chr1 | 1000 | 2000 | DEL 
chr2 | 3500 | 4500 | DUP
chr2 | 5000 | 5300 | DEL_ALU
chr3 | 1000 | 2000 | DUP_mCNV

## Options

Flag | Description
--- | ------------
-i | Sample information input
-b | BED file of CNVs
-c | Number of samples to run in parallel. Limited by available CPUs
-g | Reference Genome Build [ hg19, hg38 ]. Default is hg19
-s | random seed for genomic shuffling. Used in preprocessing
-o | VCF output 
--pre | Preprocessing output directory. If preprocessing has already been completed
--feats | Feature output directory. If feature extraction has already been completed

## Credits

####Author:

* Danny Antaki
    * dantaki@ucsd.edu
* William Brandler

gtCNV integrates coverage, discordant paired-ends, and split reads in genotype prediction. 

For a performance benchmark of gtCNV on 27 high coverage whole genomes (1000 Genomes Project) please refer to this poster presented at ASHG and WCPG (2015)

http://www.dropbox.com/s/09abkh9jsihpe0l/wcpg_antaki_poster.pdf

gtCNV was trained on 27 high coverage genomes with validated genotypes from the phase 3 intgrated structural variation release (doi:10.1038/nature15394;PMID:     26432246). 

Training sets are included under 'resources'; the CNV callset can be found here: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz

Parallelization sample-wise is possible with the optional -c flag following an integer of the number of available CPUs

requires python 2.7 and following python libraries: pysam, pybedtools, pandas, numpy, scikit-learn

## Quick Start
For a quick example read the HOWTO in the 'tutorial' directory

         ############################# WARNING ############################

         gtCNV currently only supports short-read libraries aligned to hg19

            We are working on adding options for different genome builds

            In practice, gtCNV has trouble genotyping CNVs <1kb in size

                  We will resolve this issue in future versions

## Inputs

gtCNV requires the following:
	
        1. list of BAM files with full path
	

        2. BED file of CNV positions. Tab delimited. Formated as chr    start    end    type
                CNV type must contain "DEL" or "DUP" to label losses and gains respectively 
	
## Useage

To view options and help

        $ gtCNV --help

gtCNV is a two step process. You MUST run preprocessing before genotyping step

1. Preprocessing:
       	Estimate coverage, insert size, and read length distriubtions. 

        $ gtCNV --preprocess -b bam.list [ --cpu INT, --out preprocessing.out, --seed INT ] 

2. Genotyping: 
        Copy number prediction with support vector machines
	
        $ gtCNV --genotype -b bam.list -i cnvs.bed --pre gtCNV_preprocessing_out/preprocessing.out [ --cpu INT, --out genotypes.out ] 


## History

gtCNV version 0.1 

## Credits

Author: 
        
        Danny Antaki
	
        dantaki@ucsd.edu

Acknowlegements:
        
        William Brandler
	
        Jonathan Sebat
	
        Sebat Lab http://sebatlab.ucsd.edu/index.php/software-data

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
