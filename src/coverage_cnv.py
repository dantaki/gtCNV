#!/usr/env python
import pysam
import numpy as np
def normalize_coverage(read_count,span,chr_cov,read_length): return ((read_count/span)*read_length)/chr_cov
def count_reads(cnv_list,bam,ci,chrFlag):
	"""count reads within CNV span for coverage estimation"""
	read_count=0
	bp_span=0
	for (c,s,e) in cnv_list:
		"""if the chromosome format lacks chr"""
		if chrFlag == False: c = c.replace('chr','')
		region= str(c+':'+s+'-'+e)
		bp_span+=int(e)-int(s)+1
		"""count each read within the span"""
		for read in bam.fetch(region=region):
			"""skip noninformative reads
			-- courtsey of svtyper - https://github.com/hall-lab/svtyper --"""
			if (read.is_reverse == read.mate_is_reverse
			   or read.is_proper_pair == False
			   or read.is_qcfail == True
			   or read.mapping_quality < 10
			   or read.is_secondary
			   or read.is_unmapped
			   or read.mate_is_unmapped
			   or read.is_duplicate
			   or abs(read.tlen) >= ci
			   or read.tid != read.rnext):
				continue
			read_count+=1
	return(read_count,bp_span)
def depth_of_coverage(cnv_list,bamfh,chrFlag,read_length):
	"""return median depth of coverage for CNV <= 1kb"""
	pos_doc={}
	for (c,s,e) in cnv_list:
		if chrFlag == False: c = c.replace('chr','')
		region= str(c+':'+s+'-'+e)
		depth_result = pysam.depth("-a", "-Q" "40", "-r", region, "-l", str(read_length-10), bamfh)
		str_flag=0
		if isinstance(depth_result,str):
			depth_result = depth_result.split('\n')
			str_flag=1
		for x in depth_result:
			r = x.rstrip('\n').split('\t')
			if str_flag == 1:
				if len(r)!=3: continue
				pos_doc[int(r[1])]=int(r[2])
			else: pos_doc[int(r[1])]=int(r[2])
	if len(pos_doc) !=0:
		temp=[]
		for x in pos_doc: temp.append(float(pos_doc[x]))
		return np.median(temp)
	else: return 0
