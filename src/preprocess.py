import argparse
from argparse import RawTextHelpFormatter
import sys
import os
import glob
import pybedtools
import pysam
import numpy as np
from multiprocessing import Pool
#######################
def bamList (fh):
	bam = []
        with open(fh) as f: 
        	bam = [(line.rstrip('\n')) for line in f ]
	return bam
def randomChr (seed):
	bed_master = pybedtools.BedTool('chr1 1 2',from_string=True)
	bed = pybedtools.BedTool()
	for chr_ref in glob.glob('resources/hg19_chr_lengths/*chr*'):
		bed_master = bed_master.cat(bed.random(l=100000,n=100,seed=seed,g=chr_ref),force_truncate=True,postmerge=False)
	maskbed = pybedtools.BedTool("resources/flags_segDup_unmappable.bed")
	masked = bed_master.subtract(maskbed)
	return masked.sort()
def gtCNV_stats(bed,bam,bamfh):
	read_stats = {}
	chr_size=0
	bed = bed.merge()
	rc={}
	for i in bed:
		(c,s,e) = i
		c = c.replace("chr","")
		region = str(c+":"+s+"-"+e)
		chr_size += (int(e)-int(s)+1)
		for read in bam.fetch(region=region):
			if (read.is_proper_pair == False or read.is_qcfail== True or read.is_duplicate == True): continue
			if(read.mapq < 40): continue
			if read_stats.get(str(read.qname)+str(read.is_read1)) != None: continue 
			read_stats[str(read.qname)+str(read.is_read1)] = ((read.qlen,abs(read.isize)))
	return read_stats,chr_size
def MAD (a,c=0.6745):
	if len(a) == 1:
        	d = np.median(a)
        	m = np.median(np.fabs(a - d) / c)
	else:
		d = np.median(a)
		m = np.median(np.fabs(a - d) / c)
	return m
def gtCNV (bamfh,bed,out):
	bam = pysam.AlignmentFile(bamfh,"rb")
	fhs = bamfh.split("/")
	id = fhs[-1].replace(".bam","")
	ofh = open(out,'a')
	genome_cov=[]
	genome_read_length=[]
	genome_insert_size=[]
	genome_mad=[]
	genome_size=0
	for chr in ( 'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'):
		chr_bed = bed.filter(lambda x: x.chrom==chr)
		(chr_stats,chr_size) = gtCNV_stats(chr_bed,bam,bamfh)
		read_count = len(chr_stats)
		read_length = []
		insert_size = []
		for read in chr_stats:
			(rl,isz) = chr_stats[read]
			read_length.append(rl)
			insert_size.append(isz)
		if read_count == 0: out = (bamfh,chr,'0','0','0','0',str(chr_size))
		else : 
			norm = (float(read_count)/float(chr_size))*np.median(read_length)
			out = (bamfh,chr,str(norm),str(np.median(read_length)),str(np.median(insert_size)),str(MAD(insert_size)),str(chr_size))
			genome_cov.append(norm)
			genome_read_length.append(np.median(read_length))
			genome_insert_size.append(np.median(insert_size))
			genome_mad.append(MAD(insert_size))
			genome_size += chr_size
		ofh.write('\t'.join(out))
		ofh.write('\n')
	out = (bamfh,"GENOME",str(np.median(genome_cov)),str(np.median(genome_read_length)),str(np.median(genome_insert_size)),str(np.median(genome_mad)),str(genome_size))
	ofh.write('\t'.join(out))
	ofh.write('\n')
	ofh.close()
	bam.close()
#######################
if __name__ == '__main__':
	# options 
	splash ='        __________________   ____    __\n_______ __  /__  ____/__  | / /_ |  / /\n__  __ `/  __/  /    __   |/ /__ | / / \n_  /_/ // /_ / /___  _  /|  / __ |/ /  \n_\__, / \__/ \____/  /_/ |_/  _____/   \n/____/ PREPROCESSING\n\n\n-----------------------\nREQUIRED FOR GENOTYPING\n-----------------------\n\nReturns average read length, insert size, and chromosome coverage\n'
	parser = argparse.ArgumentParser(description=splash,formatter_class=RawTextHelpFormatter)
	parser.add_argument('-b','--bam', help='list of bam files with full path',required=True)
	parser.add_argument('-c','--cpu', help='parallelize sample wise. 1 per cpu. DEFAULT=1',required=False,default=1,type=int)
	parser.add_argument('-o','--out', help='outfile path',required=False,default="gtCNV_preprocessing.out",type=str)
	parser.add_argument('-s','--seed', help='integer seed for genomic suffling',required=False,default=42,type=int)

	args = parser.parse_args()
	bam = args.bam
	cores = args.cpu
	ofh = args.out
	seed = args.seed
	outdir = 'gtCNV_preprocessing_out/'
	bamfiles = bamList(bam)
	bed = randomChr(seed)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	ofh = outdir + ofh
	head = ("bam","chr","cov","read_length_median","insert_size_median","insert_size_MAD","chr_bp_parsed")
        outfh = open(ofh,'w')
        outfh.write('\t'.join(head))
        outfh.write('\n')
        outfh.close()
	if cores > 1 :
		pool = Pool(processes=cores)
		for bamfh in bamfiles:
			pool.apply_async(gtCNV, args=(bamfh,bed,ofh) )
		pool.close()
		pool.join()
	else: 
		for bamfh in bamfiles:
			gtCNV(bamfh,bed,ofh)

