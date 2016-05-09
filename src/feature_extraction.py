#!/usr/env python
import os
from  cnv_filter import  filter_calls
import coverage_cnv as coverage
from discordant_split_cnv import discordant_split_cnv
from preprocess import bamHead
import pysam
import bed
class Prepro:
	"""preprocessing results"""
	def __init__(self,pre):
		chr_cov = {}
		insert_size={}
		mad = {}
		read_length = {}
		with open(pre) as f:
			next(f)
			for l in f:
				r = l.rstrip('\n').split('\t')
				chr_cov[(str(r[0]),str(r[1]))]=float(r[2])
				if r[1] == "GENOME":
			       		mad[r[0]]=float(r[5])
					insert_size[r[0]]=float(r[4])
					read_length[r[0]]=float(r[3])
		self.chr=chr_cov
		self.insert_size=insert_size
		self.insert_mad=mad
		self.read_len=read_length
def extract_feats(iid,bamfh,bedfh,prefh,gender,out,gen):
	cnv = bed.Bed(bedfh).cnv
	(sorted_cnv,master_cnv) = filter_calls(cnv,gen)
	bam = pysam.AlignmentFile(bamfh,"rb")
	chrFlag = bamHead(bam)
	pre = Prepro(prefh)
	insert_size = pre.insert_size[iid]
	insert_mad = pre.insert_mad[iid]
	outdir = os.getcwd()+'/gtCNV_genotypes/'
	if not os.path.exists(outdir): os.makedirs(outdir)
	ofh = open(outdir+out,'w')
	out_head = '\t'.join(('chr','start','end','type','size','id','coverage','discordant_ratio','split_ratio'))
	ofh.write(out_head+'\n')
	"""iterate for each CNV"""
	for call in sorted_cnv:
		(c,s,e,cl) = call
		(cnv_span,flank_span,windows) = master_cnv[call]
		size = int(e)-int(s)+1
		ci = insert_size+(5*insert_mad)
		(discordant,split,concordant) = discordant_split_cnv(flank_span,bam,size,ci,windows,chrFlag)
		if size > 1000:
			"""read count coverage estimation"""
			(read_count,bp_span) = coverage.count_reads(cnv_span,bam,ci,chrFlag)
			cov = coverage.normalize_coverage(float(read_count),bp_span,pre.chr[(iid,c)],pre.read_len[iid])
		else:
			"""median depth of coverage estimation"""
			median_doc = coverage.depth_of_coverage(cnv_span,bamfh,chrFlag,pre.read_len[iid])
			cov = median_doc/pre.chr[(iid,c)]
		if float(concordant) == 0.0:
			discordant_ratio = str(round(float(discordant)/1.0,3))
			split_ratio = str(round(float(split)/1.0,3))
		else:
			discordant_ratio = str(round(float(discordant)/float(concordant),3))
			split_ratio = str(round(float(split)/float(concordant),3))
		ofh.write('\t'.join((   str(c),str(s),str(e),str(cl),str(size),str(iid),
			     	    str(round(float(cov),3)),
				    discordant_ratio, split_ratio
    			 	))+'\n')
	bam.close()
	ofh.close()
