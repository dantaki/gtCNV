#!/usr/env python
import sys
import os
import glob
import pybedtools
import pysam
import numpy as np
def randomChr (seed,gen):
	bed_master = pybedtools.BedTool('chr1 1 2',from_string=True)
	bed = pybedtools.BedTool()
	for chr_ref in glob.glob('resources/'+gen+'_chr_lengths/*chr*'):
		bed_master = bed_master.cat(bed.random(l=100000,n=100,seed=seed,g=chr_ref),force_truncate=True,postmerge=False)
	maskbed = pybedtools.BedTool('resources/'+gen+'_unmapped.bed')
	masked = bed_master.subtract(maskbed)
	return masked.sort()
def bamHead(bam):
        chrFlag=False
        bamhead = bam.header['SQ']
        for i in  bamhead:
                if str(i['SN']).find("chr") != -1: chrFlag=True
        return chrFlag
def gtCNV_stats(bed,bam,chrFlag):
	read_stats = {}
	chr_size=0
	bed = bed.merge()
	rc={}
	for i in bed:
		(c,s,e) = i
		if chrFlag == False: c = c.replace("chr","")
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
def gtCNV_preprocess (iid,bamfh,bed,out):
	bam = pysam.AlignmentFile(bamfh,"rb")
	chrFlag = bamHead(bam)
	ofh = open(out,'a')
	genome_cov=[]
	genome_read_length=[]
	genome_insert_size=[]
	genome_mad=[]
	genome_size=0
	chroms = [ 'chr1','chr2','chr3','chr4','chr5','chr6',
		   'chr7','chr8','chr9','chr10','chr11','chr12',
		   'chr13','chr14','chr15','chr16','chr17','chr18',
		   'chr19','chr20','chr21','chr22','chrX','chrY'
		 ]
	for chr in chroms:
		chr_bed = bed.filter(lambda x: x.chrom==chr)
		(chr_stats,chr_size) = gtCNV_stats(chr_bed,bam,chrFlag)
		read_count = len(chr_stats)
		read_length = []
		insert_size = []
		for read in chr_stats:
			(rl,isz) = chr_stats[read]
			read_length.append(rl)
			insert_size.append(isz)
		if read_count == 0: 
			ofh.write('\t'.join(map(str,(iid,chr,'0','0','0','0',chr_size)))+'\n') 
		else : 
			norm = (float(read_count)/float(chr_size))*np.median(read_length)
			out = ( iid,chr,str(norm),
			        str(np.median(read_length)),
				str(np.median(insert_size)),
				str(MAD(insert_size)),
				str(chr_size)
			      )
			genome_cov.append(norm)
			genome_read_length.append(np.median(read_length))
			genome_insert_size.append(np.median(insert_size))
			genome_mad.append(MAD(insert_size))
			genome_size += chr_size
			ofh.write('\t'.join(out)+'\n')
	ofh.write('\t'.join( (  iid,"GENOME",
			        str(np.median(genome_cov)),
				str(np.median(genome_read_length)),
				str(np.median(genome_insert_size)),
				str(np.median(genome_mad)),
				str(genome_size)
	      		     ))+'\n')	
	ofh.write('\t'.join(out))
	ofh.write('\n')
	ofh.close()
	bam.close()
def preprocess(iid,bamfh,ofh,gen,seed):
	outdir = os.getcwd()+'/gtCNV_preprocessing/'
       	bed = randomChr(seed,gen)
        if not os.path.exists(outdir): os.makedirs(outdir)
        ofh = outdir + ofh
        outfh = open(ofh,'w')
        outfh.write('\t'.join(("id","chr","cov","read_length_median","insert_size_median","insert_size_MAD","chr_bp_parsed"))+'\n')
        outfh.close()
	gtCNV_preprocess(iid,bamfh,bed,ofh)
