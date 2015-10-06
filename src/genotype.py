import os
import svm as gtSVM
import argparse
from argparse import RawTextHelpFormatter
import sys
import pybedtools
import pysam
import numpy as np
import pandas as pd
from multiprocessing import Pool
#######################
def bamList (fh):
	bam = []
        with open(fh) as f: 
		bam = [(line.rstrip('\n')) for line in f ]
        return bam
def bedcnv (fh):
	cnv = []
	with open(fh) as f:
		for l in f:
			r = l.rstrip('\n').split('\t')
			cnv.append((r[0],r[1],r[2],r[3]))
	return cnv 
def annotatecnv (cnv):
	tag = 1
	annot=[]
	master={}
	for r in cnv:
		(c,s,e,cl) = r
		annot.append((c,s,e,cl,tag))
		master[tag]=(c,s,e,cl)
		tag+=1
	return (annot,master)
def expandcnv (cnv):
	lens={}
	with open("resources/hg19.chrom.sizes") as f:
		for l in f:
			r = l.rstrip('\n').split('\t')
			lens[r[0]]=r[1]
	dpesr=[]
	dpesrRef={}
	for r in cnv:
		(c,s,e,cl,t) = r
		t = int(t)
		if(lens.get(c) == None): continue
		chrlen = lens[c]
		s = int(s) 
		e = int(e)
		s1 = s-500
		e1 = s+500
		s2 = e-500
		e2 = e+500
		if(s1 < 0): s1=0 
		if(e1 < 0): e1=0
		if(s2 > chrlen): s2=chrlen
		if(e2 > chrlen): e2=chrlen
		dpesr.append((c,s1,e1,cl,t))
		dpesr.append((c,s2,e2,cl,t))
		dpesrRef[t] = (s1,e1,s2,e2) 
	return (dpesr,dpesrRef)
def bedconvert (cnv):
	cnvbed = pybedtools.BedTool(cnv)
	return cnvbed.sort()
def isPAR(bed):
                parbed = pybedtools.BedTool("resources/par_hg19.bed")
		if len(parbed.intersect(bed,f=0.5,wb=True)) > 0: return True
		else: return False 
def mergePos (cnv):
	bed = pybedtools.BedTool(cnv)
	bed = bed.merge()
	list=[]
	for (c,s,e) in bed: list.append((c,s,e))
	return list	
def maskBed (bed):
	maskbed = pybedtools.BedTool("resources/flags_segDup_unmappable.bed")
	masked = bed.subtract(maskbed)
	return masked.sort()
def hash_bed (bed):
	hash={}
	for r in bed:
		(c,s,e,cl,t) = r
    		t = int(t)
		if  hash.get(t) == None:
			temp_list = []
			temp_list.append((c,s,e,cl))
			hash[t] = temp_list
    		else:
        		hash[t].append((c,s,e,cl))
	return hash	
def hash_union(hash1,hash2):
	union=[]
	for k in hash1:
		if hash2.get(k) == None: continue 
		union.append(k)
	return union
def sortBed (union,master):
	tempbed=[]
	for i in union:
		(c,s,e,cl) = master[i]
		tempbed.append((c,s,e,cl,i))
	bed = pybedtools.BedTool(tempbed)
	sorted = bed.sort()
	ordered = []
	for r in sorted:
		(c,s,e,cl,i) = r
		ordered.append(int(i))
	return ordered
def read_count(cov_list,size,bamfh,bam):
	read_count=0
	bp_spanned=0
        read_mapped = {}
	sizeFlag=False
	for pos in cov_list:
		(c,s,e) = pos
		c = c.replace("chr","")
		region = str(c+":"+s+"-"+e)
		bp_spanned += int(e)-int(s)+1
		#if size < 1000:
		if size < 0:
			sizeFlag=True
			depth_result = pysam.depth("-Q" "40", "-r", region, bamfh) 
			for x in depth_result:
				r = x.rstrip('\n').split('\t')
				read_mapped[int(r[1])]=int(r[2])
		else:
			for read in bam.fetch(region=region):
				if (read.is_reverse == read.mate_is_reverse  or read.is_proper_pair == False or read.is_qcfail== True or read.is_duplicate == True or read.mapping_quality < 40): continue
				read_mapped[str(read.qname)+str(read.is_read1)]=1
	if(sizeFlag==True):
        	if(len(read_mapped) != 0):
			temp=[]
			for i in read_mapped: temp.append(float(read_mapped[i]))
			read_count =  np.median(temp)
	else:
		read_count = len(read_mapped)
	return (float(read_count),bp_spanned)
def dpesr_count(dpesr_list,ci,dpesr_window_pos,bamfh):
	(s1,e1,s2,e2) = dpesr_window_pos
	discordant_count=0
	split_count=0
	discordant_reads={}
	split_reads={}
	for pos in dpesr_list:
		(c,s,e,cl) = pos
		c = c.replace("chr","")
		region = str(c+":"+s+"-"+e)
		for read in bamfh.fetch(region=region,until_eof=True):
			if read.is_qcfail == True or read.is_duplicate == True or read.mapping_quality < 20 or read.is_reverse == read.mate_is_reverse: continue
			read_name = str(read.qname)+str(read.is_read1)
			if abs(read.tlen) >= ci:
                                if read.get_overlap(s1-1,e1-1) > 0:
                                	if( s2 <= read.mpos+1 <= e2): discordant_reads[read_name] = 1
                               	if read.get_overlap(s2-1,e2-1) > 0:
                      			if( s1 <= read.mpos+1 <= e1): discordant_reads[read_name] = 1
			if read.is_secondary == True:
                        	second_align = read.get_tag("SA").split(',')
                                if second_align[0] != c: continue
                                second_align[1] = int(second_align[1])
                                if (s1 <= read.pos+1 <= e1) and (s2 <= second_align[1] <= e2): split_reads[read_name] = 1
                               	if (s2 <= read.pos+1 <= e2) and (s1 <= second_align[1] <= e1): split_reads[read_name] = 1
	discordant_count = len(discordant_reads)
        split_count = len(split_reads)
	return(discordant_count,split_count)
def preprocess(pre):
	chr_cov = {}
	auto_cov = []
	insert_size={}
	mad = {}
	read_length = {}
	gender = {}
	with open(pre) as f:
		next(f)
		for l in f:
			r = l.rstrip('\n').split('\t')
			chr_cov_key = str(r[0])+str(r[1])
			chr_cov[chr_cov_key]=float(r[2])
			if r[1] != "GENOME" and r[1] != "chrX" and r[2] != "chrY": auto_cov.append(float(r[2]))
			if r[1] == "chrY": 
				if np.mean(auto_cov)/float(r[2]) > 0.3: gender[r[0]]="Male"
				else: gender[r[0]]="Female"
			if r[1] == "GENOME":
				mad[r[0]]=float(r[5])
				insert_size[r[0]]=float(r[4])
				read_length[r[0]]=float(r[3])
	return(chr_cov,insert_size,mad,read_length,gender)
def normalize_coverage( (read_count,span),size,chr_cov, read_length):
	if size < 1000:
		return read_count/chr_cov
	else:
		return ((read_count/span)*read_length)/chr_cov
def gtCNV (bamfh,union_cnv,hash_cov,hash_dpesr,dpesr_window,master_cnv,out,pre):
	bam = pysam.AlignmentFile(bamfh,"rb")
	(chr_cov,insert_size,mad,read_length,gender) = pre
	fhs = bamfh.split("/")
	id = fhs[-1].replace(".bam","")
	if os.path.isfile(out):
		ofh = open(out,'a')
	else:
		ofh = open(out,'w')
		head = ("chr","start","end","size","type","id","coverage","discordant_PE","split_reads","copy_number","REFERENCE LIKELIHOOD","NONREFERENCE LIKELIHOOD","PHRED REFERENCE","PHRED NONREFERENCE")
		ofh.write('\t'.join(head))
		ofh.write('\n')
	is_mad = float(mad[bamfh])
	autodel = [] 
	autodup = []
	sexdel = []
	sexdup = []
	for i in union_cnv:
		(c,s,e,cl) = master_cnv[i]
		autosome = True
		if gender[bamfh] == "Male" and (c == "chrX"  or c == "chrY"): autosome = isPAR(bedconvert([(c,s,e)]))
		size = int(e)-int(s)+1
		ci = insert_size[bamfh]+(5*is_mad)
		cov_list = hash_cov[i]
		cov_list = mergePos(cov_list)
		dpesr_list = hash_dpesr[i]
		dpesr_window_pos = dpesr_window[i]	
		cov  = normalize_coverage(read_count(cov_list,size,bamfh,bam),size,chr_cov[str(bamfh)+str(c)],read_length[bamfh])
		(dpe,sr) = dpesr_count(dpesr_list,ci,dpesr_window_pos,bam)	
		out = (c,str(s),str(e),str(size),cl,id,str(cov),str(dpe),str(sr))
		if autosome == True and cl.find('DEL') != -1: autodel.append(out)
		if autosome == True and cl.find('DUP') != -1: autodup.append(out)
		if autosome == False and cl.find('DEL') != -1: sexdel.append(out)
		if autosome == False and cl.find('DUP') != -1: sexdup.append(out)
	outdf = pd.DataFrame()
	if len(autodel) > 0: 
		autodel = gtSVM.autosome_del_svm(pd.DataFrame(autodel))
		outdf = outdf.append(autodel)
	if len(sexdel) > 0: 
		sexdel = gtSVM.sexchr_del_svm(pd.DataFrame(sexdel))
		outdf = outdf.append(sexdel)
	if len(autodup) > 0:
		autodup = gtSVM.autosome_dup_svm(pd.DataFrame(autodup))
		outdf = outdf.append(autodup)
	if len(sexdup) > 0: 
		sexdup = gtSVM.sexchr_dup_svm(pd.DataFrame(sexdup))
		outdf = outdf.append(sexdup)
	outdf.to_csv(ofh,sep='\t',header=False,index=False,mode='a')
	ofh.close()
	bam.close()
#######################
if __name__ == '__main__':
	splash ='        __________________   ____    __\n_______ __  /__  ____/__  | / /_ |  / /\n__  __ `/  __/  /    __   |/ /__ | / / \n_  /_/ // /_ / /___  _  /|  / __ |/ /  \n_\__, / \__/ \____/  /_/ |_/  _____/   \n/____/ GENOTYPING\n\n\n\nReturns CNV genotypes\n'
	parser = argparse.ArgumentParser(description=splash,formatter_class=RawTextHelpFormatter)
        parser.add_argument('-b','--bam', help='list of bam files with full path',required=True)
	parser.add_argument('-i','--bed', help='BEDfile of CNVs. Tab delimited. CHROM    START    END    TYPE',type=str,required=True)
	parser.add_argument('-p','--pre', help='gtCNV preprocessing output. DEFAULT=gtCNV_proprocessing_out/gtCNV_preprocessing.out',type=str,required=False,default="gtCNV_proprocessing_out/gtCNV_preprocessing.out")
        parser.add_argument('-c','--cpu', help='parallelize sample wise. 1 per cpu. DEFAULT=1',required=False,default=1,type=int)
        parser.add_argument('-o','--out', help='outfile path',required=False,default="gtCNV_genotypes.out",type=str)
	args = parser.parse_args()
        bam = args.bam
        bed = args.bed
	cores = args.cpu
        ofh = args.out
	prefh=args.pre
        outdir = 'gtCNV_genotype_out/'
	ofh = outdir+ofh
	if not os.path.exists(outdir):
                os.makedirs(outdir)
	bamfiles = bamList(bam)
	pre = preprocess(prefh)
	cnv=[]
	master_cnv={}
	cnv = bedcnv(bed)	
	(cnv,master_cnv) = annotatecnv(cnv)
	(cnv_dpesr,dpesr_window) = expandcnv(cnv)
	cnvbed = bedconvert(cnv)
	mask_cov = maskBed(cnvbed)
	cnvbed = bedconvert(cnv_dpesr)
	mask_dpesr = maskBed(cnvbed)
	hash_cov = hash_bed(mask_cov)
	hash_dpesr = hash_bed(mask_dpesr)
	union_cnv = hash_union(hash_cov,hash_dpesr)
	union_cnv = sortBed(union_cnv,master_cnv)
	if cores > 1 :
		pool = Pool(processes=cores)
		for bamfh in bamfiles:
			pool.apply_async(gtCNV, args=(bamfh,union_cnv,hash_cov,hash_dpesr,dpesr_window,master_cnv,ofh,pre) )
		pool.close()
		pool.join()
	else: 
		for bamfh in bamfiles:
			gtCNV(bamfh,union_cnv,hash_cov,hash_dpesr,dpesr_window,master_cnv,ofh,pre)
