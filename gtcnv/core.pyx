#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
import os,sys,pysam,datetime
from math import ceil
from glob import glob
import pybedtools as pbed
from time import time
from datetime import timedelta	
import pandas as pd
import numpy as np
cimport numpy as np
try:from configparser import ConfigParser
except ImportError:from ConfigParser import ConfigParser
try: import cPickle as pickle
except ImportError: import pickle
np.import_array()
def errFH(fh):
	sys.stderr.write('ERROR {} NOT FOUND\n'.format(fh))
	sys.exit(1)
cdef warnFH(fh): sys.stderr.write('WARNING {} NOT FOUND\n'.format(fh))
cdef checkCall(cl):
	cdef bint c
	c=0
	if 'DEL' in cl or 'DUP' in cl: c=1
	return c
def Bed(fh):
	cdef unsigned int s
	cdef unsigned int e
	cdef bint v
	v=0
	c=''
	cl=''
	cnv=[]
	if not os.path.isfile(fh): errFH(fh)
	header=open(fh).readline().rstrip()
	if 'VCF' in header: v=1
	if v==True:
		with open(fh) as f:
			for l in f:
				if l.startswith('#'): continue
				s=0
				e=0
				r= l.rstrip().split('\t')
				c=str(r[0])
				if not c.startswith('chr'): c='chr'+c
				s=int(r[1])
				for i in r[7].split(';'):
					if 'SVTYPE=' in i: cl=i.replace('SVTYPE=','')
					if 'END=' in i and 'CI' not in i and i.startswith('END='): e=int(i.replace('END=',''))
				if s<=0 or s>=e or e<=0: continue
				if checkCall(cl)==True:cnv.append((c,s,e,cl))
	else:
		with open(fh) as f:
			for l in f:
				if l.startswith('#'): continue
				s=0
				e=0
				r=l.rstrip().split('\t')
				c=str(r[0])
				s=int(r[1])
				e=int(r[2])
				if s<=0 or s>=e or e<=0: continue
				cl=r[3]
				if not c.startswith('chr'): c='chr'+c	
				if checkCall(cl)==True: cnv.append((c,s,e,cl))
	scnv=[]
	for x in pbed.BedTool(cnv).sort():
		(c,ss,ee,cl) = x
		s=int(ss)
		e=int(ee)
		scnv.append((c,s,e,cl))	
	return scnv
def get_path():
	try:
		root = __file__
		if os.path.islink(root): root = os.path.realpath(root)
		return os.path.dirname(os.path.abspath(root))
	except:
		sys.stderr.write('ERROR {} NOT FOUND\n'.format(__file__))
		sys.exit(1)
def check_in(fh):
	cdef bint skip
	bams={}
	vcfs={}
	gender={}
	if not os.path.isfile(fh): errFH(fh)
	genders = ['M','F']
	with open(fh) as f:
		for l in f:
			r = l.rstrip('\n').split('\t')
			if len(r) != 4:
				print 'WARNING: {} not formatted correctly! Please refer to the documentation at www.github.com/dantaki/gtCNV\n'.format(l)
			else:
				(iid,bamfh,vcffh,sx) = r
				skip=0
				if not os.path.isfile(bamfh):
					skip=1
					warnFH(bamfh)
				if not os.path.isfile(vcffh):
					skip=1
					warnFH(vcffh)
				if sx == '1': sx = 'M'
				if sx == '2': sx = 'F'
				if sx not in genders:
					print 'WARNING: {} not an accepted gender. Accepted genders: 1|M 2|F>\n'
					skip=1
				if skip==True: continue
				bams[iid]=bamfh
				vcfs[iid]=vcffh
				gender[iid]=sx
	return bams,vcfs,gender
cdef get_chroma():
	chroms=[]
	for x in range(1,23): chroms.append('chr'+str(x))
	chroms.append('chrX')
	chroms.append('chrY')
	return chroms
def check_cnv(cnv,gen):
	cdef unsigned int s
	sz={}
	with open(get_path()+'/resources/{}.bedtools.genome'.format(gen)) as f:
		for l in f:
			r = l.rstrip('\n').split('\t')
			s=int(r[1])
			sz[r[0]]=s
	chroms=get_chroma()
	raw=[]
	ok=[]
	for (c,s,e,cl) in cnv:
		raw.append((c,s,e,cl))
		if c not in chroms: print 'WARNING: {} not accepted chromosome format!\nAccepted chromosomes:\n{}'.format(c,'\n'.join(chroms))
		else:
			if e > sz[c]: print 'WARNING: {}:{}-{} {} has position greater than {} length!'.format(c,s,e,cl,c)
			else: ok.append((c,s,e,cl))
	return raw,ok
def reportTime(init_time,s):
	f=[]
	s = ' '+s+' <time elapsed: {}>'.format(timedelta(seconds=int(time())-init_time))
	for x in range(len(s)+2): f.append('-')
	print ''.join(f)+'\n'+s+'\n'+''.join(f)
cdef fasta_config(gen):
	cfg = get_path()+'/config/gtcnv.ini'
	if not os.path.isfile(cfg):
		sys.stderr.write('ERROR {} NOT FOUND\n'.format(cfg))
		sys.exit(1)
	conf = ConfigParser()
	conf.read(cfg)
	return conf.get('FASTA_PATHS',gen)
cdef FASTA_head(fh):
	chr_flag=False
	with open (fh,'r') as f:
		l = f.next()
		if '>' not in l:
			sys.stderr.write('ERROR {} MUST CONTAIN ">CHROM" in the first line\n'.format(fh))
			sys.exit(1)
		if '>chr' in l: chr_flag=True
	return chr_flag
cdef append_list2dict(a):
	d={}
	for i in a:
		(c,s,e,cl,tag) = i
		if d.get(tag) == None: d[tag]=[(c,s,e)]
		else: d[tag].append((c,s,e))
	return d
cdef annot_toList(cnv):
	cnv_list = []
	cnv.sort()
	for i in cnv:
		(c,s,e,cl,tag) = i
		cnv_list.append((c,s,e,cl,tag))
	return cnv_list
cdef filter_calls(cnv,gen):
	cdef unsigned int tag
	"""add unique annotation to each CNV"""
	tag = 1
	annot_cnv=[]
	cnv_dict={}
	start_flank=[]
	end_flank=[]
	for r in cnv:
		cnv_dict[tag]=r
		(c,s,e,cl) = r
		annot_cnv.append((c,s,e,cl,tag))
		start_flank.append((c,s,s,cl,tag))
		end_flank.append((c,e,e,cl,tag))
		tag+=1
	"""add 500bp to each start and end position"""
	start_flank500bp=pbed.BedTool(start_flank).slop(b=500,g=get_path()+'/resources/'+gen+'.bedtools.genome')
	end_flank500bp=pbed.BedTool(end_flank).slop(b=500,g=get_path()+'/resources/'+gen+'.bedtools.genome')
	start_flank_dict=dict((i[4],(i[1],i[2])) for i in start_flank500bp)
	end_flank_dict=dict((i[4],(i[1],i[2])) for i in end_flank500bp)
	"""mask regions to segmental duplications and unmappable regions"""
	masked_cnv=annot_toList(pbed.BedTool(annot_cnv).subtract(pbed.BedTool(get_path()+'/resources/'+gen+'_unmapped.bed.gz')))
	masked_start_flank=annot_toList(start_flank500bp.subtract(pbed.BedTool(get_path()+'/resources/'+gen+'_unmapped.bed.gz')))
	masked_end_flank=annot_toList(end_flank500bp.subtract(pbed.BedTool(get_path()+'/resources/'+gen+'_unmapped.bed.gz')))
	"""take the union of the two masked call sets"""
	if len(masked_cnv) == 0:
		sys.stderr.write('ERROR: CNVs completed masked by excluded regions\n')
		sys.exit(1)
	cnv_tag = list(zip(*masked_cnv)[4])
	passed_tag = list(set(cnv_tag))
	masked_cnv_dict=append_list2dict(masked_cnv)
	masked_start_flank_dict=append_list2dict(masked_start_flank)
	masked_end_flank_dict=append_list2dict(masked_end_flank)
	master_cnv={}
	for i in passed_tag:
		windows = (start_flank_dict[i],end_flank_dict[i])
		cnv_span = masked_cnv_dict[i]
		flank_span=None
		if masked_start_flank_dict.get(i) != None and masked_end_flank_dict.get(i) != None:
			flank_span = masked_start_flank_dict[i] + masked_end_flank_dict[i]
		master_cnv[cnv_dict[int(i)]]=(cnv_span,flank_span,windows)
	return master_cnv
cdef gc_norm(pcrfree):
	pcr_dict={}
	with open(get_path()+'/resources/GC_content_reference.txt','r') as f:
		for l in f:
			(tag,GC,norm) = l.rstrip('\n').split('\t')
			key = 'DOC'
			if 'READCOUNT' in tag: key ='RC'
			if pcrfree == True:
				if 'PCRFREE' in tag: pcr_dict[(key,int(GC))]=float(norm)
			else:
				if 'PCRFREE' not in tag: pcr_dict[(key,int(GC))]=float(norm)
	return pcr_dict
cdef gc_cnv(cnv,master_cnv,gen):
	cnv_gc={}
	FASTA = fasta_config(gen)
	chrFlag = FASTA_head(FASTA)
	for call in cnv:
		if master_cnv.get(call)==None: continue
		cnvs = master_cnv[call][0]
		if chrFlag == False: cnvs = [tuple(map(lambda y: str.replace(y,'chr',''),x)) for x in cnvs]
		GC_content = int(5 * round(float(GC(pbed.BedTool(cnvs).nucleotide_content(fi=FASTA)))/5))
		cnv_gc[call]=GC_content
	return cnv_gc
cdef GC(nuc):
	cdef int A=0
	cdef int C=0
	cdef int G=0
	cdef int T=0
	cdef int N=0
	for x in nuc:
		A+=int(x[5])
		C+=int(x[6])
		G+=int(x[7])
		T+=int(x[8])
		N+=int(x[9])
	return ((G+C)*100)/(A+C+G+T+N)
cdef MAD(a, double c=0.6745): return np.around(np.median(np.fabs(a-np.median(a))/c),decimals=3)
cdef normalize_chr_cov(double read_count, double chr_size, read_length): return np.around( (read_count/chr_size)*np.median(read_length),decimals=2)
cdef normalize_coverage(double read_count, double span, double chr_cov, double read_length): return (((read_count/span)*read_length)/chr_cov)
cdef count_reads(cnv_list,bam,ci,chrFlag):
	"""count reads within CNV span for coverage estimation"""
	cdef unsigned int read_count
	cdef unsigned int bp_span
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
cdef depth_of_coverage(cnv_list,bamfh,chrFlag,read_length):
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
cdef is_discordant(read,windows,ci,matepos):
	"""returns True if read is discordant"""
	((s1,e1),(s2,e2)) = windows
	if abs(read.tlen) >= ci:
		"""insert size is greater than 5MAD from median"""
		if (int(s1) <= read.pos+1 <= int(e1) and int(s2) <= int(matepos)+1 <= int(e2)) or (int(s2) <= read.pos+1 <= int(e2) and int(s1) <= int(matepos)+1 <= int(e1)):
			return True
		else:
			return False
cdef is_split(read,windows,c):
	"""returns True if read is split"""
	((s1,e1),(s2,e2)) = windows
	if read.is_secondary == True:
		second_align = read.get_tag("SA").split(',')
		if second_align[0] != c:
			"""return False if maps to other chromosome"""
			return False
		else:
			second_align[1] = int(second_align[1])
			if ( ((int(s1) <= read.pos+1 <= int(e1)) and (int(s2) <= second_align[1] <= int(e2)))
			     or ((int(s2) <= read.pos+1 <= int(e2)) and (int(s1) <= second_align[1] <= int(e1)))):
				"""secondary alignment must be near the opposite breakpoint of the primary alignment"""
				return True
			else:
				return False
cdef discordant_split_cnv(flank_list,bam,size,ci,windows,chrFlag):
	"""count discordant paired ends and split reads"""
	cdef unsigned int discordant_count
	cdef unsigned int split_count
	cdef unsigned int concordant_count	
	discordant_count=0
	split_count=0
	concordant_count=0
	if flank_list==None: return (0,0,0)
	else:
		for (c,s,e) in flank_list:
			if chrFlag == False : c = c.replace("chr","")
			region = str(c+":"+s+"-"+e)
			for read in bam.fetch(region=region,until_eof=True):
				"""skip noninformative reads
				- courtsey of svtyper - https://github.com/hall-lab/svtyper --"""
				if (read.is_qcfail
				    or read.is_unmapped
				    or read.mate_is_unmapped
				    or read.tid != read.rnext
				    or read.is_reverse == read.mate_is_reverse
				    or read.is_duplicate): continue
				else:
					mate_pos = read.pnext
					if is_discordant(read,windows,ci,mate_pos) == True: discordant_count+=1
					if is_split(read,windows,c)  == True: split_count+=1
					if (read.is_proper_pair
					    and read.mapping_quality >= 10
					    and not read.is_secondary
					    and abs(read.tlen) < ci):
						concordant_count+=1
					else: continue
		
		return(discordant_count,split_count,concordant_count)
cdef bamHead(bam):
	cdef bint chrFlag
	chrFlag=0
	bamhead = bam.header['SQ']
	for i in  bamhead:
		if str(i['SN']).find("chr") != -1: chrFlag=1
	return chrFlag
cdef randomChr (seed,gen):
	bed_master = pbed.BedTool('chr1 1 2',from_string=True)
	bed = pbed.BedTool()
	for chr_ref in glob(get_path()+'/resources/'+gen+'_chr_lengths/*chr*'):
		bed_master = bed_master.cat(bed.random(l=100000,n=100,seed=seed,g=chr_ref),force_truncate=True,postmerge=False)
	maskbed = pbed.BedTool(get_path()+'/resources/'+gen+'_unmapped.bed.gz')
	return bed_master.subtract(maskbed).sort().merge()
cdef gtCNV_stats(bed,bam,chrFlag):
	cdef unsigned int chr_size
	read_stats = {}
	chr_size=0
	for (c,s,e) in bed:
		if chrFlag == False: c = c.replace("chr","")
		region = str(c+":"+s+"-"+e)
		chr_size += (int(e)-int(s)+1)
		for read in bam.fetch(region=region):
			if (read.is_proper_pair == False or read.is_qcfail== True or read.is_duplicate == True or read.mapq < 40 or read_stats.get(str(read.qname)+str(read.is_read1)) != None): continue
			read_stats[str(read.qname)+str(read.is_read1)] = ((read.qlen,abs(read.isize)))
	return read_stats,chr_size
cdef gtCNV_preprocess (iid,bamfh,vcffh,bed,out,gen):
	cdef double read_count
	cdef double genome_snp_cov
	cdef unsigned int genome_SNP
	cdef unsigned int genome_size
	cdef unsigned int chr_snp
	cdef unsigned int chr_size
	cdef double chr_snp_cov
	bam = pysam.AlignmentFile(bamfh,"rb")
	chrFlag = bamHead(bam)
	ofh = open(out,'a')
	genome_cov=[]
	genome_read_length=[]
	genome_insert_size=[]
	genome_mad=[]
	genome_size=0
	chroms=get_chroma()
	snp_cov = gtVCF(vcffh,iid,get_path()+'/resources/'+gen+'_unmapped.bed.gz').preprocessVCF()
	for chromo in chroms:
		filt_bed = bed.filter(lambda x: x.chrom==chromo)
		chr_bed=[tuple(x) for x in filt_bed]
		(chr_stats,chr_size) = gtCNV_stats(chr_bed,bam,chrFlag)
		chr_snp=0
		chr_snp_cov=0
		if snp_cov.get(chromo) != None: chr_snp_cov, chr_snp = snp_cov[chromo]
		read_count = len(chr_stats)
		read_length = []
		insert_size = []
		for read in chr_stats:
			(rl,isz) = chr_stats[read]
			read_length.append(rl)
			insert_size.append(isz)
		if read_count == 0:
			ofh.write('\t'.join(map(str,(iid,chromo,'0','0','0','0',chr_size,chr_snp_cov,chr_snp)))+'\n')
		else :
			norm = normalize_chr_cov(read_count,chr_size,np.array(read_length))
			out = ( iid,chromo,str(norm),
				str(np.median(read_length)),
				str(np.median(insert_size)),
				str(MAD(np.array(insert_size))),
				str(chr_size),str(chr_snp_cov),str(chr_snp)
			      )
			genome_cov.append(norm)
			genome_read_length.append(np.median(read_length))
			genome_insert_size.append(np.median(insert_size))
			genome_mad.append(MAD(np.array(insert_size)))
			genome_size += chr_size
			ofh.write('\t'.join(out)+'\n')
	genome_snp_cov=0
	genome_SNP=0
	if snp_cov.get('GENOME') != None: genome_snp_cov, genome_SNP = snp_cov['GENOME']
	ofh.write('\t'.join( (  iid,"GENOME",
				str(np.median(genome_cov)),
				str(np.median(genome_read_length)),
				str(np.median(genome_insert_size)),
				str(np.median(genome_mad)),
				str(genome_size),
				str(genome_snp_cov),
				str(genome_SNP)
			     ))+'\n')
	ofh.close()
	bam.close()
class gtVCF():
	def __init__(self,FH=None,iid=None,mask=None):
		cdef unsigned int depth=0
		cdef double a=0.0
		cdef double b=0.0
		cdef double ratio=0.0
		cdef int IND=-9
		cdef int DP_IND=-9
		cdef int GT_IND=-9
		cdef int AD_IND=-9
		cdef int s=0
		snps=[]
		hets=[]
		SNP={}
		HET={}
		CHROMS=get_chroma()
		if not os.path.isfile(FH): errFH(FH)
		lines=[]
		with open(FH,'r') as f:
			for l in f:
				if not l.startswith('#'): break
				else:
					lines.append(l.rstrip('\n'))
		self.iid=iid
		self.fh=FH
		self.isHC = False
		self.isFB = False
		if [x for x in lines if 'HaplotypeCaller' in x]: self.isHC=True
		if [x for x in lines if 'freeBayes' in x]: self.isFB=True
		if self.isHC == self.isFB:
			sys.stderr.write('ERROR: {} IS AMBIGUOUS. VCF HEADER MUST CONTAIN EITHER <HaplotypeCaller> or <freeBayes>\n'.format(FH))
			sys.exit(1)
		self.formatTag=None
		if self.isFB == True: self.formatTag = 'DPR'
		if self.isHC == True: self.formatTag = 'AD'
		with open(FH) as f:
			for l in f:
				if l.startswith('##'): continue
				elif l.startswith('#CHROM'):
					r=l.rstrip('\n').split('\t')
					for x in xrange(len(r)):
						if r[x] in iid: IND=x
					if IND == -9:
						sys.stderr.write('ERROR {} NOT FOUND IN VCF HEADER: {}\n'.format(iid,l.rstrip('\n')))
						sys.exit(1)
				else:
					r=l.rstrip('\n').split('\t')
					try: r[IND]
					except IndexError:
						sys.stderr.write('ERROR {} MALFORMED. skipping ...\n'.format(l.rstrip('\n')))
						continue
					if len(r) < IND: continue
					c = r[0]
					s = int(r[1])
					if 'chr' not in c: c='chr'+str(c)
					c = c.replace('chr23','chrX').replace('chr24','chrY')
					if c not in CHROMS: continue
					FORMAT = r[8].split(':')
					DP_IND=GT_IND=AD_IND=-9
					for x in xrange(len(FORMAT)):
						if FORMAT[x] == 'DP': DP_IND=x
						if FORMAT[x] == 'GT': GT_IND=x
						if FORMAT[x] == self.formatTag: AD_IND=x
					if DP_IND == -9 or GT_IND == -9: continue
					sample_info = r[IND].split(':')
					killer=False
					for x in xrange(9,len(r)):
						temp_info = r[x].split(':')
						if temp_info[GT_IND] == './.' or temp_info[GT_IND] == '.' or temp_info[GT_IND] == '.|.':killer=True
					if killer == True: continue
					if sample_info[DP_IND] == '.' or sample_info[AD_IND] == '.': continue
					depth=int(sample_info[DP_IND])
					sample_info[GT_IND]=sample_info[GT_IND].replace('|','/')
					gt = [int(x) for x in sample_info[GT_IND].split('/')]
					snps.append((c,s,s,depth))
					if len(gt) != 2: continue
					if gt[0] != gt[1]:
						if AD_IND == -9: continue
						allelic_depth = [float(x) for x in sample_info[AD_IND].split(',')]
						a = allelic_depth[0]
						b = allelic_depth[1]
						if a == 0.0 or b == 0.0: continue
						ratio = b/a
						if b > a: ratio=a/b
						hets.append((c,s,s,ratio))
		LociMask=pbed.BedTool(mask).sort()
		snps = [tuple(x) for x in pbed.BedTool(snps).sort().intersect(LociMask,wa=True,v=True,sorted=True).sort()]
		hets = [tuple(x) for x in pbed.BedTool(hets).sort().intersect(LociMask,wa=True,v=True,sorted=True).sort()]
		for chrom in CHROMS:
			snpsub = [(int(x[1]),int(x[3])) for x in snps if x[0]==chrom]
			hetsub = [(int(x[1]),float(x[3])) for x in hets if x[0]==chrom]
			SNP[chrom]=snpsub
			HET[chrom]=hetsub
		self.SNP=SNP
		self.HET=HET
	def preprocessVCF(self):
		snp_cov={}
		CHROMS = get_chroma()
		genome=[]
		for chrom in CHROMS:
			if self.SNP.get(chrom) == None: snp_cov[chrom]=(0,0)
			else:
				if len(self.SNP[chrom])==0: snp_cov[chrom]=(0,0)
				else: snp_cov[chrom]=(np.median([x[1] for x in self.SNP[chrom]]),len(self.SNP[chrom]))
				genome+=[x[1] for x in self.SNP[chrom]]
		snp_cov['GENOME']=(np.median(genome),len(genome))
		return snp_cov	
	def extract_SNP_features(self,sorted_cnv,master_cnv,snp_cov):
		cdef int s=0
		cdef int e=0
		SNP_feats={}
		HET_feats={}
		for call in sorted_cnv:
			if master_cnv.get(call)==None: continue
			(cnv_span,flank_span,windows) = master_cnv[call]
			snp_cov_locus=[]
			het_ratio_locus=[]
			for x in cnv_span:
				c = x[0]
				norm_cov = snp_cov[(self.iid,c)]
				s,e = map(int,(x[1],x[2]))
				if self.SNP.get(c) == None: continue
				SNPsub = self.SNP[c]
				for x in SNPsub:
					if x[0] >= s and x[0] <= e: snp_cov_locus.append(x[1])
					elif x[0] > e: break	
				if self.HET.get(c) == None: continue
				HETsub = self.HET[c]
				for x in HETsub:
					if x[0] >= s and x[0] <= e: het_ratio_locus.append(x[1])
					elif x[0] > e: break
			if len(snp_cov_locus)==0: SNP_feats[call]=(float('nan'),0)
			else: SNP_feats[call]=(np.nanmedian(snp_cov_locus)/norm_cov,len(snp_cov_locus))
			if len(het_ratio_locus)==0: HET_feats[call]=(float('nan'),0)
			else: HET_feats[call]=(np.nanmedian(het_ratio_locus),len(het_ratio_locus))
		return SNP_feats,HET_feats
def preprocess(iid,bamfh,vcffh,ofh,gen,seed):
	outdir = os.getcwd()+'/gtCNV_preprocessing/'
	bed = randomChr(seed,gen)
	if not os.path.exists(outdir): os.makedirs(outdir)
	ofh = outdir + ofh
	outfh = open(ofh,'w')
	outfh.write('\t'.join(("#ID","CHROM","COVERAGE_MEDIAN","READ_LENGTH_MEDIAN","INSERT_SIZE_MEDIAN","INSERT_SIZE_MAD","BAM_BP_PARSED","SNP_COVERAGE_MEDIAN","SNP_PARSED"))+'\n')
	outfh.close()
	gtCNV_preprocess(iid,bamfh,vcffh,bed,ofh,gen)
class Prepro(object):
	def __init__(self,pre):
		chr_cov={}
		insert_size={}
		mad={}
		read_length={}
		snp_cov={}
		with open(pre) as f:
			next(f)
			for l in f:
				r = l.rstrip('\n').split('\t')
				chr_cov[(str(r[0]),str(r[1]))]=float(r[2])
				snp_cov[(r[0],r[1])]=float(r[7])
				if r[1] == "GENOME":
					mad[r[0]]=float(r[5])
					insert_size[r[0]]=float(r[4])
					read_length[r[0]]=float(r[3])
		self.chr_cov=chr_cov
		self.insert_size=insert_size
		self.insert_mad=mad
		self.read_len=read_length
		self.snp_cov=snp_cov
def extract_feats(iid,bamfh,vcffh,cnv,prefh,gender,out,gen,pcrfree):
	master_cnv = filter_calls(cnv,gen)
	gc_dict = gc_norm(pcrfree)
	cnv_gc=gc_cnv(cnv,master_cnv,gen)
	bam = pysam.AlignmentFile(bamfh,"rb")
	chrFlag = bamHead(bam)
	pre = Prepro(prefh)
	insert_size = pre.insert_size[iid]
	insert_mad = pre.insert_mad[iid]
	snp_cov = pre.snp_cov
	outdir = os.getcwd()+'/gtCNV_features/'
	if not os.path.exists(outdir): os.makedirs(outdir)
	ofh = open(outdir+out,'w')
	out_head = '\t'.join(('chr','start','end','type','size','id','coverage','coverage_GCcorrected','discordant_ratio','split_ratio','SNP_coverage','HET_ratio','SNPs','HET_SNPs'))
	ofh.write(out_head+'\n')
	SNPs, HETs = gtVCF(vcffh,iid,get_path()+'/resources/'+gen+'_unmapped.bed.gz').extract_SNP_features(cnv,master_cnv,snp_cov)
	for call in cnv:
		(c,s,e,cl) = call
		if SNPs.get(call)==None or HETs.get(call)==None or cnv_gc.get(call)==None or master_cnv.get(call)==None: continue
		snp_coverage, snps = SNPs[call]
		het_ratio, hets = HETs[call]
		GC_content = cnv_gc[call]
		(cnv_span,flank_span,windows) = master_cnv[call]
		size = int(e)-int(s)+1
		ci = insert_size+(5*insert_mad)
		discordant,split,concordant = discordant_split_cnv(flank_span,bam,size,ci,windows,chrFlag)
		if size > 1000:
			"""read count coverage estimation"""
			(read_count,bp_span) = count_reads(cnv_span,bam,ci,chrFlag)
			cov=float('nan')
			if pre.chr_cov[(iid,c)] != 0: cov = normalize_coverage(read_count,bp_span,pre.chr_cov[(iid,c)],pre.read_len[iid])
			GCnorm_ratio = 1.0
			if gc_dict.get(('RC',GC_content))==None:
				if gc_dict.get(('DOC',GC_content)) !=None: GCnorm_ratio=gc_dict[('DOC',GC_content)]
			else: GCnorm_ratio = gc_dict[('RC',GC_content)]
			GCcov = cov * GCnorm_ratio
		else:
			"""median depth of coverage estimation"""
			median_doc = depth_of_coverage(cnv_span,bamfh,chrFlag,pre.read_len[iid])
			cov=float('nan')
			if pre.chr_cov[(iid,c)]!= 0: cov = median_doc/pre.chr_cov[(iid,c)]
			GCnorm_ratio=1.0
			if gc_dict.get(('DOC',GC_content)) != None: GCnorm_ratio=gc_dict[('DOC',GC_content)]
			GCcov = cov * GCnorm_ratio
		if float(concordant) == 0.0:
			discordant_ratio = str(round(float(discordant)/1.0,3))
			split_ratio = str(round(float(split)/1.0,3))
		else:
			discordant_ratio = str(round(float(discordant)/float(concordant),3))
			split_ratio = str(round(float(split)/float(concordant),3))
		ofh.write('\t'.join((   str(c),str(s),str(e),str(cl),str(size),str(iid),
				    str(round(float(cov),3)), str(round(float(GCcov),3)),
				    discordant_ratio, split_ratio,
				    str(np.around(snp_coverage,decimals=3)),str(snps),str(np.around(het_ratio,decimals=3)),str(hets) ))+'\n')
	bam.close()
	ofh.close()
cdef returnPAR(cnv,gen):
	out=[]
	results = pbed.BedTool(cnv).intersect(pbed.BedTool(get_path()+'/resources/par_'+gen+'.bed'),f=0.5,wa=True,u=True)
	if len(results) > 0:
		for x in results: out.append(tuple(x))
	return out
cdef removePAR(cnv,gen):
	out=[]
	results= pbed.BedTool(cnv).subtract(pbed.BedTool(get_path()+'/resources/par_'+gen+'.bed'),f=0.5,A=True)
	if len(results) > 0:
		for x in results.sort(): out.append(tuple(x))
	return out
def likFilter(x,lik,f):
	'''
	de novo filter recommendations
	'''
	cdef unsigned int sz
	cdef double dpe
	cdef double sr
	cdef unsigned short passflag
	clf=x[len(x)-1]
	sz=int(x[3])
	dpe=float(x[7])
	sr=float(x[8])
	passflag=1
	if 'DEL' in clf and 'HEMIZYGOUS' not in clf:
		if (dpe!=0 or sr!=0) and lik<12: passflag=0
		elif dpe==0 and sr==0:
			if f=='r':
				if sz<=1000: passflag=0
				elif (1000<sz<=5000) and lik<20: passflag=0
				elif sz>5000 and lik<18: passflag=0
			elif f=='n':
				if sz<=1000: passflag=0
				elif (1000<sz<=3000) and lik<24: passflag=0
				elif sz>3000 and lik<20: passflag=0
	elif 'DUP_PairedEnd' in clf and lik <11: passflag=0
	elif 'DUP_HetAllelicDepth' in clf:
		if (dpe!=0 or sr!=0) and lik<10: passflag=0
		elif dpe==0 and sr==0:
			if f=='r':
				if sz<=3000: passflag=0
				elif (3000<sz<=100000) and lik<15: passflag=0
				elif sz>100000 and lik<13: passflag=0
			elif f=='n':
				if sz<=3000: passflag=0
				elif sz>3000 and lik<10: passflag=0
	elif 'HEMIZYGOUS_DEL' in clf and lik<8: passflag=0
	elif 'HEMIZYGOUS_DUP' in clf:
		if sz<=5000: passflag=0
		elif sz>5000 and lik<10: passflag=0
	return passflag
def genotype(feats,sex,gen,out):
	REF={}
	NON={}
	GQ={}
	HEMI={}
	FILT={}
	males = [k for k in sex if sex[k] == 'M']
	sex_chrom = ['chrX','chrY']
	dels = [ k for k in feats if 'DEL' in k[3]]
	dups = [ k for k in feats if 'DUP' in k[3]]
	autosome_dels = [k for k in dels if k[0] not in sex_chrom]
	autosome_dups = [k for k in dups if k[0] not in sex_chrom]
	sexchr_dels = [k for k in dels if k[0] in sex_chrom]
	sexchr_dups = [k for k in dups if k[0] in sex_chrom]
	del_par=returnPAR(sexchr_dels,gen)
	dup_par=returnPAR(sexchr_dups,gen)
	if len(del_par) > 0:
		for k in del_par: autosome_dels.append(k)
	if len(dup_par) > 0:
		for k in dup_par: autosome_dups.append(k)
	male_sexchr_del=[]
	male_sexchr_dup=[]
	for k in removePAR(sexchr_dels,gen):
		if k[5] not in males: autosome_dels.append(k)
		else: male_sexchr_del.append(k)
	for k in removePAR(sexchr_dups,gen):
		if k[5] not in males: autosome_dups.append(k)
		else: male_sexchr_dup.append(k)
	head = ['chr','start','end','type','size','id','old_covr','covr_GC','dpe','sr','SNP_coverage','SNPs','HET_ratio','HET_SNPs']
	pd.options.mode.chained_assignment = None
	genos=[]
	outdir = os.getcwd()+'/gtCNV_genotypes/'
	if not os.path.exists(outdir): os.makedirs(outdir)
	ofh = open(outdir+out,'w')
	ofh.write('\t'.join(('#CHROM','START','END','LENGTH','TYPE','ID','COVERAGE','DISCORDANT_PAIRED-END','SPLIT_READS','SNP_COVERAGE','SNPs','HETEROZYGOUS_ALLELIC_DEPTH_RATIO','HETEROZYGOUS_SNPs','COPY_NUMBER','GTL_HOM-ALT','GTL_HET','GTL_HOM-REF','GTL_ALT','PHRED_REF','PHRED_ALT','CLASSIFIER'))+'\n')
	if len(autosome_dels) > 0:
		autosome_del_df = pd.DataFrame(autosome_dels)
		autosome_del_df.columns=head
		for x in autosome_del_svm(autosome_del_df).values:
			ofh.write('\t'.join(map(str,x))+'\n')
			k = (x[0],x[1],x[2],x[4],x[5])
			kk = (x[0],x[1],x[2],x[4])
			x[6] = format(float(x[6])*2,'.2f')
			GQ[k]= ','.join((format(float(x[14]),'.2f'),format(float(x[15]),'.2f'),format(float(x[16]),'.2f')))
			if int(x[13]) == 2:
				x[13]='0/0'
				if x[5] not in males and x[0] == 'chrY': x[13]='.'
				else:
					lik=ceil(float(x[18]))
					FILT[k]=likFilter(x,lik,'r')
					if REF.get(kk)==None: REF[kk]=[lik]
					else: REF[kk].append(lik)
			else:
				if int(x[13])==1: x[13]='0/1'
				else: x[13]='1/1'
				if x[5] not in males and x[0] == 'chrY': x[13]='.'
				else:
					lik=ceil(float(x[19]))
					FILT[k]=likFilter(x,lik,'n')
					if NON.get(kk)==None: NON[kk]=[lik]
					else: NON[kk].append(lik)
			genos.append(tuple(x))
	if len(autosome_dups) > 0:
		autosome_dup_df = pd.DataFrame(autosome_dups)
		autosome_dup_df.columns=head
		for x in autosome_dup_svm(autosome_dup_df).values:
			ofh.write('\t'.join(map(str,x))+'\n')
			x[6] = format(float(x[6])*2,'.2f')
			GQ[k]= ','.join((format(float(x[14]),'.2f'),format(float(x[15]),'.2f'),format(float(x[16]),'.2f')))
			if int(x[13]) == 2:
				x[13]='0/0'
				if x[5] not in males and x[0] == 'chrY': x[13]='.'
				else:
					lik=ceil(float(x[18]))
					FILT[k]=likFilter(x,lik,'r')
					if REF.get(kk)==None: REF[kk]=[lik]
					else: REF[kk].append(lik)
			else:
				if int(x[13])==3: x[13]='0/1'
				else: x[13]='1/1'
				if x[5] not in males and x[0] == 'chrY': x[13]='.'
				else:
					lik=ceil(float(x[19]))
					FILT[k]=likFilter(x,lik,'n')
					if NON.get(kk)==None: NON[kk]=[lik]
					else: NON[kk].append(lik)
			genos.append(tuple(x))
	if len(male_sexchr_del) > 0:
		sexchr_del_df = pd.DataFrame(male_sexchr_del)
		sexchr_del_df.columns=head
		for x in sexchr_del_svm(sexchr_del_df).values:
			ofh.write('\t'.join(map(str,x))+'\n')
			x[6] = format(float(x[6]),'.2f')
			GQ[k]= ','.join((format(float(x[15]),'.2f'),format(float(x[16]),'.2f')))
			HEMI[kk]=1
			if int(x[13]) == 1:
				x[13]='0'
				lik=ceil(float(x[18]))
				FILT[k]=likFilter(x,lik,'r')
				if REF.get(kk)==None: REF[kk]=[lik]
				else: REF[kk].append(lik)
			else:
				x[13]='1'
				lik=ceil(float(x[19]))
				FILT[k]=likFilter(x,lik,'n')
				if NON.get(kk)==None: NON[kk]=[lik]
				else: NON[kk].append(lik)
			genos.append(tuple(x))
	if len(male_sexchr_dup) > 0:
		sexchr_dup_df = pd.DataFrame(male_sexchr_dup)
		sexchr_dup_df.columns=head
		for x in sexchr_dup_svm(sexchr_dup_df).values:
			ofh.write('\t'.join(map(str,x))+'\n')
			x[6] = format(float(x[6])*2,'.2f')
			GQ[k]= ','.join((format(float(x[15]),'.2f'),format(float(x[16]),'.2f')))
			HEMI[kk]=1
			if int(x[13]) == 1:
				x[13]='0'
				lik=ceil(float(x[18]))
				if x[5] not in males and x[0] == 'chrY':
					x[13]='.'
					FILT[k]=0
				if x[5] in males:
					FILT[k]=likFilter(x,lik,'r')	
					if REF.get(kk)==None: REF[kk]=[lik]
					else: REF[kk].append(lik)
			else:
				x[13]='1'
				lik=ceil(float(x[19]))
				if x[5] not in males and x[0] == 'chrY':
						x[13]='.'
						FILT[k]=0	
				if x[5] in males:
					FILT[k]=likFilter(x,lik,'n')	
					if NON.get(kk)==None: NON[kk]=[lik]
					else: NON[kk].append(lik)
			genos.append(tuple(x))
	refmed={}
	nonmed={}
	ofh.close()
	for x in REF: refmed[x]=int(np.median(REF[x]))
	for x in NON: nonmed[x]=int(np.median(NON[x]))
	return genos,refmed,nonmed,GQ,HEMI,FILT
cdef del_autosome_tab(df,lik,clfname):
	df['classif']=clfname
	cn=[]
	for x in lik:
		if x[0] > x[1] and x[0] > x[2]: cn.append(0)
		elif x[1] > x[0] and x[1] > x[2]: cn.append(1)
		elif x[2] > x[0] and x[2] > x[1]: cn.append(2)
		else: cn.append(2)
	df['copy_number']=cn	
	df['likHOM'] = lik[:,0]
	df['likHET'] = lik[:,1]
	df['likREF'] = lik[:,2]
	df['NONREF'] = 1 - df['likREF']
	df['PHRED_REF'] = -10.0 * np.log10(1.0-df['likREF'])
	df['PHRED_NONREF'] = -10.0 * np.log10(1.0-df['NONREF'])
	df= df[['chr','start','end','size','type','id','covr','dpe','sr','SNP_coverage','SNPs','HET_ratio','HET_SNPs','copy_number','likHOM','likHET','likREF','NONREF','PHRED_REF','PHRED_NONREF','classif']]	
	return df
cdef del_sexchr_tab (df,lik,clfname):
	df['classif']=clfname
	cn=[]
	for x in lik:
		if x[0] > x[1]: cn.append(0)
		elif x[1] > x[0]: cn.append(1)
		else: cn.append(1)
	df['copy_number']=cn
	df['likVAR'] = lik[:,0]
	df['likREF'] = lik[:,1]
	df['NONREF'] = 1 - df['likREF']
	df['PHRED_REF'] = -10.0 * np.log10(1.0-df['likREF'])
	df['PHRED_NONREF'] = -10.0 * np.log10(1.0-df['NONREF'])	
	df= df[['chr','start','end','size','type','id','covr','dpe','sr','SNP_coverage','SNPs','HET_ratio','HET_SNPs','copy_number','likVAR','likVAR','likREF','NONREF','PHRED_REF','PHRED_NONREF','classif']]
	return df
cdef dup_autosome_tab(df,lik,clfname):
	df['classif']=clfname
	cn=[]
	for x in lik:
		if x[0] > x[1] and x[0] > x[2]: cn.append(4)
		elif x[1] > x[0] and x[1] > x[2]: cn.append(3)
		elif x[2] > x[0] and x[2] > x[1]: cn.append(2)
		else: cn.append(2)	
	df['copy_number']=cn
	df['likCN4'] = lik[:,0]
	df['likCN3'] = lik[:,1]
	df['likREF'] = lik[:,2]
	df['NONREF'] = 1 - df['likREF']
	df['PHRED_REF'] = -10.0 * np.log10(1.0-df['likREF'])
	df['PHRED_NONREF'] = -10.0 * np.log10(1.0-df['NONREF'])
	df= df[['chr','start','end','size','type','id','covr','dpe','sr','SNP_coverage','SNPs','HET_ratio','HET_SNPs','copy_number','likCN4','likCN3','likREF','NONREF','PHRED_REF','PHRED_NONREF','classif']]
	return df
cdef dup_sexchr_tab(df,lik,clfname):
	df['classif']=clfname
	cn=[]
	for x in lik:
		if x[0] > x[1]: cn.append(2)
		elif x[1] > x[0]: cn.append(1)
		else: cn.append(2)
	df['copy_number']=cn
	df['likVAR'] = lik[:,0]
	df['likREF'] = lik[:,1]
	df['NONREF'] = 1 - df['likREF']
	df['PHRED_REF'] = -10.0 * np.log10(1.0-df['likREF'])
	df['PHRED_NONREF'] = -10.0 * np.log10(1.0-df['NONREF'])
	df= df[['chr','start','end','size','type','id','covr','dpe','sr','SNP_coverage','SNPs','HET_ratio','HET_SNPs','copy_number','likVAR','likVAR','likREF','NONREF','PHRED_REF','PHRED_NONREF','classif',]]
	return df
cdef merge3d (df):
	return np.vstack(np.asarray( (df['covr'],df['dpe'],df['sr']), order = 'C', dtype='float' )).T
cdef prepHET(df):
	X = np.vstack(np.asarray((df['covr'],df['HET_ratio']),order='C',dtype='float')).T
	return df,X
cdef prepPE(df):
	X = merge3d(df)
	return df,X
cdef prep3d(df):
	df[['covr_GC','dpe','sr']]=df[['covr_GC','dpe','sr']].astype(float)
	df[['size']]=df[['size']].astype(int)
	df=df[df['covr_GC'] < 5.0]
	df['covr']=df['covr_GC']
	X = merge3d(df)
	return df,X
cdef prep4d(df):
	df[['covr_GC','dpe','sr']]=df[['covr_GC','dpe','sr']].astype(float)
	df[['size']]=df[['size']].astype(int)
	df = df[df['covr_GC'] < 5.0]
	df['covr']=df['covr_GC']
	big = df[df['size'] > 1000]
	sma = df[df['size'] <= 1000]
	X1 = merge3d(big)
	X2 = merge3d(sma)
	return big,sma,X1,X2
cdef autosome_del_svm(df):
	(bigdf,smadf,BIG,SMA) = prep4d(df)
	final=pd.DataFrame()
	if len(BIG) > 0:
		clf=''
		with open(get_path()+'/resources/training_sets/1000Genomes_HighCov_autosome_deletion_LRG-CLF.pkl','rb') as f: clf=pickle.load(f)
		lik = clf.predict_proba(BIG)
		final=final.append(del_autosome_tab(bigdf,lik,'DEL_gt1KB'))
	if len(SMA) > 0:
		clf=''
		with open(get_path()+'/resources/training_sets/1000Genomes_HighCov_autosome_deletion_SMA-CLF.pkl','rb') as f: clf=pickle.load(f)
		lik = clf.predict_proba(SMA)
		final=final.append(del_autosome_tab(smadf,lik,'DEL_lt1KB'))
	return final
cdef sexchr_del_svm(df):
	(df,X) = prep3d(df)
	clf=''
	with open(get_path()+'/resources/training_sets/1000Genomes_HighCov_sexchr_deletion_CLF.pkl','rb') as f: clf=pickle.load(f)
	lik = clf.predict_proba(X)
	return del_sexchr_tab(df,lik,'HEMIZYGOUS_DEL')
cdef autosome_dup_svm(df):
	snpCLF=pd.DataFrame()
	peCLF=pd.DataFrame()
	final=pd.DataFrame()
	df[['SNPs','HET_SNPs']]=df[['SNPs','HET_SNPs']].astype(int)
	df[['covr_GC','dpe','sr']]=df[['covr_GC','dpe','sr']].astype(float)
	mean_covr=[]
	for ind,ele in df.iterrows():
		if ele['SNPs'] > 0: mean_covr.append(np.mean((ele['covr_GC'],float(ele['SNP_coverage']))))
		else: mean_covr.append(ele['covr_GC'])
	df['covr']=mean_covr
	df=df[df['covr'] < 5.0]
	snpList=[]
	peList=[]
	for ind,ele in df.iterrows():
		if ele['dpe']==0 and ele['sr']==0 and ele['HET_SNPs'] > 0: snpList.append(ind)
		else: peList.append(ind)
	snpCLF=df.loc[snpList]
	snpCLF[['HET_ratio']]=snpCLF[['HET_ratio']].astype(float)
	peCLF=df.loc[peList]
	X=None
	clf=None
	temp=pd.DataFrame()
	if len(snpCLF) >0:
		temp,X = prepHET(snpCLF)
		with open(get_path()+'/resources/training_sets/1000Genomes_HighCov_autosome_duplication_SNP-CLF.pkl','rb') as f: clf=pickle.load(f)
		lik = clf.predict_proba(X)
		final=final.append(dup_autosome_tab(temp,lik,'DUP_HetAllelicDepth'))
	clf=None
	X=None
	temp=pd.DataFrame()
	if len(peCLF) >0:
		temp,X= prepPE(peCLF)
		with open(get_path()+'/resources/training_sets/1000Genomes_LowCov_autosome_duplication_PE-CLF.pkl','rb') as f: clf=pickle.load(f)
		lik = clf.predict_proba(X)
		final=final.append(dup_autosome_tab(temp,lik,'DUP_PairedEnd'))
	return final
cdef sexchr_dup_svm(df):
	(df,X) = prep3d(df)
	clf=''
	with open(get_path()+'/resources/training_sets/1000Genomes_LowCov_sexchr_duplication_CLF.pkl','rb') as f: clf=pickle.load(f)
	lik = clf.predict_proba(X)
	return dup_sexchr_tab(df,lik,'HEMIZYGOUS_DUP')
cdef recipOver(x):
	cdef unsigned int s1
	cdef unsigned int s2
	cdef unsigned int e1
	cdef unsigned int e2
	(s1,s2,e1,e2) = x
	cdef double ovr
	cdef double ovr1
	cdef double ovr2
	sz1 = e1-s1+1
	sz2 = e2-s2+1
	smax=max(s1,s2)
	emin=min(e1,e2)
	ovr = emin-smax+1
	return min((ovr/sz1),(ovr/sz2))
cdef cytobandOverlap(cnv,gen):
	cyto={}
	try:
		for (c,s,e,c2,s2,e2,band,t,ovr) in cnv.intersect(get_path()+'/resources/annotation_files/{}_cytoband.bed'.format(gen),wao=True):
			if cyto.get((c,s,e))==None: cyto[(c,s,e)]=c2.replace('chr','')+band
			else: cyto[(c,s,e)]=cyto[(c,s,e)]+','+c2.replace('chr','')+band
	except pbed.cbedtools.MalformedBedLineError: cyto[(c,s,e)]='NA' 	
	return cyto
cdef flagsOverlap(cnv,gen):
	flags={}
	cdef int ovr
	for x in cnv.intersect(get_path()+'/resources/annotation_files/{}_flags.bed.gz'.format(gen),wao=True):
		x=list(x)
		x[len(x)-1]=int(x[len(x)-1])
		if(x[(len(x)-1)]==0): continue
		(c,s,e,c2,s2,e2,flag,ovr) = x
		if flags.get((c,s,e,flag))==None: flags[(c,s,e,flag)]=int(ovr)
		else: flags[(c,s,e,flag)]+=int(ovr)
	for (c,s,e,f) in flags: flags[(c,s,e,f)]= float(flags[(c,s,e,f)])/(int(e)-int(s)+1.0)
	return flags
cdef meiOverlap(cnv,gen):
	cdef double rovr
	mei={}
	maxovr={}
	for x in cnv.intersect(get_path()+'/resources/annotation_files/{}_repeatmasker.bed.gz'.format(gen), f=0.8, F=0.8, wa=True, wb=True):
		rovr=0
		rovr=recipOver(map(int,(x[1],x[4],x[2],x[5])))
		k = (x[0],x[1],x[2])
		name = '{}:{}:{}'.format(x[7],x[8],x[6])	
		if maxovr.get(k)==None:
			maxovr[k]=rovr
			mei[k]=name
		elif maxovr.get(k) != None and rovr > maxovr[k]:
			maxovr[k]=rovr
			mei[k]=name
		else: continue
	return mei,maxovr
cdef thouGenOverlap(genos,gen):
	thouGen={}
	dels = pbed.BedTool([(x[0],x[1],x[2],x[3]) for x in genos if 'DEL' in x[3]]).sort()
	dups = pbed.BedTool([(x[0],x[1],x[2],x[3]) for x in genos if 'DUP' in x[3]]).sort()
	for x in dels.intersect(get_path()+'/resources/annotation_files/{}_1000Genomes_DEL.bed'.format(gen), f=0.5, F=0.5, wao=True):
		x=list(x)
		if int(x[len(x)-1])==0: continue
		else: thouGen[(x[0],x[1],x[2],x[3])]=(x[len(x)-2],format(float(x[len(x)-1])/(int(x[2])-int(x[1])+1.0),'.2f'))
	for x in dups.intersect(get_path()+'/resources/annotation_files/{}_1000Genomes_DUP.bed'.format(gen), f=0.5, F=0.5, wao=True):
		x=list(x)
		if int(x[len(x)-1])==0: continue
		else: thouGen[(x[0],x[1],x[2],x[3])]=(x[len(x)-2],format(float(x[len(x)-1])/(int(x[2])-int(x[1])+1.0),'.2f'))
	return thouGen
cdef geneOverlap(cnv,gen):
	cdef unsigned int exonTotal
	cdef unsigned int exonNum
	cdef unsigned int intronTotal
	cdef unsigned int intronNum
	cdef unsigned int exoncnt
	cdef unsigned int exontot
	cdef unsigned int introncnt
	cdef unsigned int introntot	
	GENEFH=get_path()+'/resources/annotation_files/{}_genes.bed.gz'.format(gen)
	genes={}
	geneTEMP={}
	for (c,s,e,c2,s2,e2,gene) in cnv.intersect(GENEFH, wa=True,wb=True):
		if geneTEMP.get((c,s,e))==None: geneTEMP[(c,s,e)]=[gene]
		elif gene not in geneTEMP[(c,s,e)]: geneTEMP[(c,s,e)].append(gene)
		else: continue
	for x in geneTEMP:
		ori={}
		exons={}
		introns={}
		trx={}
		genlist=[]
		for y in geneTEMP[x]:
			gene = y.split(',')
			ori[gene[0]]=gene[1]
			if 'exon' in y:
				trx[gene[0]]=gene[2]
				exonTotal = int(gene[len(gene)-1].split('/').pop())
				exonNum = int(gene[len(gene)-1].split('/').pop(0).replace('exon_',''))
				exons[(gene[0],'TOTAL')]=exonTotal
				if exons.get(gene[0])==None:
					exons[gene[0]]=1
					exons[(gene[0],exonNum)]=1
				elif exons.get((gene[0],exonNum)) == None:
					exons[gene[0]]+=1
					exons[(gene[0],exonNum)]=1
				else: continue
			elif 'UTR3' in y:
				trx[gene[0]]=gene[2]
				exons[(gene[0],'UTR3')]=1
			elif 'UTR5' in y:
				trx[gene[0]]=gene[2]
				exons[(gene[0],'UTR5')]=1
			elif 'intron' in y:
				trx[gene[0]]=gene[2]
				intronTotal = int(gene[len(gene)-1].split('/').pop())
				intronNum = int(gene[len(gene)-1].split('/').pop(0).replace('intron_',''))
				introns[(gene[0],'TOTAL')]=intronTotal
				if introns.get(gene[0])==None:
					introns[gene[0]]=1
					introns[(gene[0],intronNum)]=1
				elif introns.get((gene[0],intronNum)) == None:
					introns[gene[0]]+=1
					introns[(gene[0],intronNum)]=1
				else: continue
			elif 'stream' in y:
				trx[gene[0]]=gene[2]
				exons[(gene[0],gene[3])]=1
			else: continue
		for y in trx:
			if ori.get(y) == None: continue
			orient=ori[y]
			exoncnt=0
			exontot=0
			if exons.get(y) != None: exoncnt=exons[y]
			if exons.get((y,'TOTAL')) != None: exontot=exons[(y,'TOTAL')]
			if orient == '+':
				if exons.get((y,'upstream_1kb')) == 1 and exons.get((y,'downstream_1kb')) != 1:
					genlist.append(','.join(map(str,(y,trx[y],'upstream_1kb'))))
				if exons.get((y,'UTR3'))!=1 and exons.get((y,'UTR5')) == 1:
					genlist.append(','.join(map(str,(y,trx[y],'UTR5'))))
				if exoncnt != 0: genlist.append(','.join(map(str,(y,trx[y],'exon_{}/{}'.format(exoncnt,exontot)))))
				if exons.get(y) == None or exoncnt != exontot:
					introncnt=0
					introntot=0
					if introns.get(y) != None: introncnt=introns[y]
					if introns.get((y,'TOTAL'))!=None: introntot=introns[(y,'TOTAL')]
					if introncnt !=0: genlist.append(','.join(map(str,(y,trx[y],'intron_{}/{}'.format(introncnt,introntot)))))
				if exons.get((y,'UTR3'))==1 and exons.get((y,'UTR5')) != 1:
					genlist.append(','.join(map(str,(y,trx[y],'UTR3'))))
				if exons.get((y,'upstream_1kb')) != 1 and exons.get((y,'downstream_1kb')) == 1:
					genlist.append(','.join(map(str,(y,trx[y],'downstream_1kb'))))
			else:
				if exons.get((y,'upstream_1kb')) != 1 and exons.get((y,'downstream_1kb')) == 1:
					genlist.append(','.join(map(str,(y,trx[y],'downstream_1kb'))))
				if exons.get((y,'UTR3'))==1 and exons.get((y,'UTR5')) != 1:
					genlist.append(','.join(map(str,(y,trx[y],'UTR3'))))
				if exoncnt != 0: genlist.append(','.join(map(str,(y,trx[y],'exon_{}/{}'.format(exoncnt,exontot)))))
				if exons.get(y) == None or exoncnt != exontot:
					introncnt=0
					introntot=0
					if introns.get(y) != None: introncnt=introns[y]
					if introns.get((y,'TOTAL'))!=None: introntot=introns[(y,'TOTAL')]
					if introncnt !=0: genlist.append(','.join(map(str,(y,trx[y],'intron_{}/{}'.format(introncnt,introntot)))))
				if exons.get((y,'UTR3'))!=1 and exons.get((y,'UTR5')) == 1:
					genlist.append(','.join(map(str,(y,trx[y],'UTR5'))))
				if exons.get((y,'upstream_1kb')) == 1 and exons.get((y,'downstream_1kb')) != 1:
					genlist.append(','.join(map(str,(y,trx[y],'upstream_1kb'))))
		if len(genlist)>=1 : genes[x]='|'.join(genlist)
	return genes
def annotate(raw,genos,gen,REF,NON,GQ,OFH,sex,hemi,filt):
	cdef unsigned int CNV_ID
	cdef unsigned int SZ
	genes={}
	cyto={}
	mei={}
	flags={}
	thouGen={}
	FILT={}
	AF={}
	GT={}
	refcnt={}
	males = [k for k in sex if sex[k] == 'M']
	IID = list(set([x[5] for x in genos]))
	IID.sort(key=str.lower)
	cnv=pbed.BedTool(list(set([(x[0],x[1],x[2]) for x in raw]))).sort()
	cyto=cytobandOverlap(cnv,gen)
	flags=flagsOverlap(cnv,gen)
	mei,meiRO=meiOverlap(cnv,gen)
	thouGen=thouGenOverlap(raw,gen)
	genes=geneOverlap(cnv,gen)
	for x in genos:
		k = (x[0],x[1],x[2],x[4])
		filty=0
		if filt.get((x[0],x[1],x[2],x[4],x[5]))!=None: filty=filt[(x[0],x[1],x[2],x[4],x[5])]
		if GQ.get((x[0],x[1],x[2],x[4],x[5]))!=None:
			lik=format(float(x[len(x)-2]),'.2f')
			if '1' not in x[13]:
				lik=format(float(x[len(x)-3]),'.2f')
				if refcnt.get(k)==None: refcnt[k]=1
				else: refcnt[k]+=1
			if '1' in x[13]:
				if AF.get((k,x[5]))==None:
					if AF.get(k)==None: AF[k]=1
					else: AF[k]+=1
					AF[(k,x[5])]=1
				lik=format(float(x[len(x)-2]),'.2f')
			v = ':'.join(map(str,(x[13],
					filty,
					format(float(x[7]),'.2f'), #DPE
					format(float(x[8]),'.2f'), #SR
					format(float(x[6]),'.2f'), #COVR
					format(float(x[9]),'.2f'),  #SNP COVR
					int(x[10]), #SNPs
					format(float(x[11]),'.2f'), #HET ratio
					int(x[12]), #HETs
					lik,
					GQ[(x[0],x[1],x[2],x[4],x[5])])))
			if FILT.get(k)==None: FILT[k]=[filty]
			else : FILT[k].append(filty)
			GT[(k,x[5])]=v
	for x in genos:
		k = (x[0],x[1],x[2],x[4])
		for iid in IID:
			if GT.get((k,iid))==None:
				if x[0] == 'chrY' and iid not in males: GT[(k,iid)]='.'
				elif x[0] == 'chrY' and iid in males: GT[(k,iid)]='0'
				elif x[0] == 'chrX' and iid in males and hemi.get(k) != None: GT[(k,iid)]='0'
				else: GT[(k,iid)]='0/0'
	VCF=[]
	CNV_ID=1
	for (c,s,e,cl) in pbed.BedTool(list(set([(x[0],x[1],x[2],x[3]) for x in raw]))).sort():
		SZ=int(e)-int(s)+1
		TYPE=cl
		THOUGEN_ID='NA'
		THOUGEN_OVR='NA'
		SEGD=0.00
		ABPTS=0.00
		UNMAP=0.00
		CENTMER=0.00
		STR=0.00
		meiName="NA"
		REPEATMASKER=0.00
		GENE='intergenic'
		CYTOB='NA'
		MEDREF='NA'
		QUAL='NA'
		ALLELE='NA'
		GTS=[]
		PASS='PASS'
		fail=[]
		cnvref=0
		if refcnt.get((c,s,e,cl))!=None: cnvref=refcnt[(c,s,e,cl)]
		if mei.get((c,s,e))!=None: 
			meiName=mei[(c,s,e)]
			TYPE=TYPE+':'+meiName
		DESX='{}:{}-{}_{}'.format(c,s,e,TYPE.replace('DEL','deletion').replace('DUP','duplication'))
		if cyto.get((c,s,e))!=None: CYTOB=cyto[(c,s,e)]
		if flags.get((c,s,e,'abparts'))!=None: ABPTS=format(flags[(c,s,e,'abparts')],'.2f')
		if flags.get((c,s,e,'centromere'))!=None: CENTMER=format(flags[(c,s,e,'centromere')],'.2f')
		if flags.get((c,s,e,'segDup'))!=None: SEGD=format(flags[(c,s,e,'segDup')],'.2f')
		if flags.get((c,s,e,'STR'))!=None: STR=format(flags[(c,s,e,'STR')],'.2f')
		if flags.get((c,s,e,'unmapable'))!=None: UNMAP=format(flags[(c,s,e,'unmapable')],'.2f')
		if thouGen.get((c,s,e,cl))!=None: THOUGEN_ID,THOUGEN_OVR = thouGen[(c,s,e,cl)]
		if meiRO.get((c,s,e))!=None: REPEATMASKER=format(meiRO[(c,s,e)],'.2f')
		if genes.get((c,s,e))!=None: GENE=genes[(c,s,e)]
		if REF.get((c,s,e,cl))!=None: MEDREF=REF[(c,s,e,cl)]
		if NON.get((c,s,e,cl))!=None: QUAL=NON[(c,s,e,cl)]
		if AF.get((c,s,e,cl))!=None: ALLELE=format(AF[(c,s,e,cl)]/float(len(IID)),'.2f')
		if float(ABPTS) >= 0.5: fail.append('ABPARTS')
		if float(CENTMER) >= 0.5: fail.append('CENTROMERE')
		if float(SEGD) >= 0.5: fail.append('SEGDUP')
		if float(STR) >= 0.5: fail.append('STR')
		if float(UNMAP) >= 0.5: fail.append('UNMAPABLE')
		if cnvref == len(IID): fail.append('ALLREF')
		if FILT.get((c,s,e,cl))!= None: 
			gtPass = len([x for x in FILT[(c,s,e,cl)] if x==1])
			if gtPass/len(FILT[(c,s,e,cl)]) < 0.9: fail.append('GENOTYPE-FAIL')
		else: fail.append('FAIL')
		if len(fail) > 0: PASS = ','.join(fail)
		INFO = 'END={};SVTYPE={};SVLEN={};AF={};CYTOBAND={};REPEATMASKER={},{};1000G_ID={};1000G_OVERLAP={};DESCRIPTION={};GENES={};ABPARTS={};CENTROMERE={};SEGDUP={};STR={};UNMAPABLE={};MEDREFGL={}'.format(e,TYPE,SZ,ALLELE,CYTOB,meiName,REPEATMASKER,THOUGEN_ID,THOUGEN_OVR,DESX,GENE,ABPTS,CENTMER,SEGD,STR,UNMAP,MEDREF)
		for x in IID:
			if GT.get(((c,s,e,cl),x))!=None: GTS.append(GT[((c,s,e,cl),x)])
			else: GTS.append('./.')
		out= '\t'.join(map(str,(c,s,CNV_ID,'.','<'+TYPE+'>',QUAL,PASS,INFO,'GT:FT:PE:SR:CN:SC:SP:AR:HT:SQ:GL','\t'.join(GTS))))
		VCF.append(out)
		CNV_ID+=1
	date=[]
	date.append(datetime.date.today())
	VCFHEAD=[	'##fileformat=VCFv4.1',
			'##fileDate={}'.format(date[0]),
			'##reference={}'.format(fasta_config(gen)),
			'##GTCNV_CMD="{}"'.format(' '.join(map(str,sys.argv[:]))),
			'##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">',
			'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
			'##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
			'##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency,in the range (0,1)">',
			'##INFO=<ID=CYTOBAND,Number=.,Type=String,Description="cytoBand(s) of the variant">',
			'##INFO=<ID=REPEATMASKER,Number=2,Type=String,Description="Name and reciprocal overlap of RepeatMasker vairant">',
			'##INFO=<ID=1000G_ID,Number=1,Type=String,Description="1000 Genomes Phase 3 integrated SV callset variant name">',
			'##INFO=<ID=1000G_ID,Number=1,Type=Float,Description="Overlap to 1000 Genomes Phase 3 variant, in the range (0,1)">',
			'##INFO=<ID=DESCRIPTION,Number=1,Type=String,Description="Verbose description of SV"',
			'##INFO=<ID=GENES,Number=1,Type=String,Description="Genes within this SV locus, pipe-separated by transcripts with comma-separated details including the refGene name>"',
			'##INFO=<ID=ABPARTS,Number=1,Type=Float,Description="Parts of antibodies overlap, in the range (0,1)">',
			'##INFO=<ID=CENTROMERE,Number=1,Type=Float,Description="Centromere overlap, in the range (0,1)">',
			'##INFO=<ID=SEGDUP,Number=1,Type=Float,Description="Segmental duplication overlap, in the range (0,1)">',
			'##INFO=<ID=STR,Number=1,Type=Float,Description="Short Tandem Repeat overlap, in the range (0,1)">',
			'##INFO=<ID=UNMAPABLE,Number=1,Type=Float,Description="Overlap with regions unmapable using 100bp reads with respect to the hg19 assembly, in the range (0,1)">',
			'##INFO=<ID=MEDREFGL,Number=1,Type=Float,Description="Median Phred-scaled genotype likelihood for individuals genotyped homozygous reference">',
			'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
			'##FORMAT=<ID=FT,Number=1,Type=Integer,Description="Sample filter, 1 indicates PASS">',
			'##FORMAT=<ID=PE,Number=1,Type=Float,Description="Ratio of Discordant Paired-Ends to Concordant Paired-Ends">',
			'##FOTMAT=<ID=SR,Number=1,Type=Float,Description="Ratio of Split reads to Concordant Paired-Ends">',
			'##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy Number">',
			'##FORMAT=<ID=SC,Number=1,Type=Float,Description="SNP normalized coverage">',
			'##FORMAT=<ID=SP,Number=1,Type=Float,Description="Number of SNPs within locus">',
			'##FORMAT=<ID=AR,Number=1,Type=Float,Description="Heterozygous Allelic Depth ratio">',
			'##FORMAT=<ID=HT,Number=1,Type=Float,Description="Number of Heterozygous SNPs">',
			'##FORMAT=<ID=SQ,Number=1,Type=Float,Description="Phred-scaled genotype likelihood">',
			'##FORMAT=<ID=GL,Number=2|3,Type=Float,Description="Phred-scaled genotype likelihood; homozygous alt, heterogygous alt, homozygous ref">',
			'##ALT=<ID=DEL,Description="Deletion, if 80% reciprocal overlap with RepeatMasker element, the class, name, and family are given separated by colons">'
			'##ALT=<ID=DUP,Description="Duplication, if 80% reciprocal overlap with RepeatMasker element, the class, name, and family are given separated by colons">',
		]
	with open(get_path()+'/resources/{}.chrom.sizes'.format(gen)) as f:
		for l in f:
			(chrom,size) = l.rstrip('\n').split('\t')
			VCFHEAD.append('##contig=<ID={},length={}>'.format(chrom.replace('chr',''),size))
	outdir = os.getcwd()+'/gtCNV_genotypes/'
	if not os.path.exists(outdir): os.makedirs(outdir)
	outfh = open(outdir+OFH,'w')
	outfh.write('\n'.join(VCFHEAD)+'\n')
	outfh.write('\t'.join(('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','\t'.join(IID)))+'\n')
	outfh.write('\n'.join(VCF)+'\n')
	outfh.close()
