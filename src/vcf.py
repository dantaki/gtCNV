#!/usr/env python
from pybedtools import BedTool as Bed
from numpy import median
import datetime
def cytobandOverlap(cnv,gen):
	cyto={}
	for (c,s,e,c2,s2,e2,band,t,ovr) in cnv.intersect('resources/annotation_files/{}_cytoband.bed'.format(gen),wao=True):
		if cyto.get((c,s,e))==None: cyto[(c,s,e)]=c2.replace('chr','')+band
               	else: cyto[(c,s,e)]=cyto[(c,s,e)]+'|'+c2.replace('chr','')+band
	return cyto
def flagsOverlap(cnv,gen):
	flags={}
	for x  in cnv.intersect('resources/annotation_files/{}_flags.bed'.format(gen),wao=True):
                if(int(x[-1])==0): continue
                (c,s,e,c2,s2,e2,flag,ovr) = x
                if flags.get((c,s,e,flag))==None: flags[(c,s,e,flag)]=int(ovr)
                else: flags[(c,s,e,flag)]+=int(ovr)
	for (c,s,e,f) in flags: flags[(c,s,e,f)]= float(flags[(c,s,e,f)])/(int(e)-int(s)+1.0)
	return flags
def meiOverlap(cnv,gen):
	mei={}
	maxovr={}
        for x in cnv.intersect('resources/annotation_files/{}_MEI.bed'.format(gen), f=0.8, F=0.8, wao=True,):
                if int(x[-1])==0: continue
                elif maxovr.get((x[0],x[1],x[2]))==None:
                        maxovr[(x[0],x[1],x[2])]=int(x[-1])
                        mei[(x[0],x[1],x[2])]=x[-2]
                elif maxovr.get((x[0],x[1],x[2])) != None and int(x[-1]) > maxovr[(x[0],x[1],x[2])]:
                        maxovr[(x[0],x[1],x[2])]=int(x[-1])
                        mei[(x[0],x[1],x[2])]=x[-2]
                else: continue
	return mei
def thouGenOverlap(genos,gen):
	thouGen={}
	dels = Bed([(x[0],x[1],x[2],x[4]) for x in genos if 'DEL' in x[4]]).sort()
        dups = Bed([(x[0],x[1],x[2],x[4]) for x in genos if 'DUP' in x[4]]).sort()
	for x in dels.intersect('resources/annotation_files/{}_1000Genomes_DEL.bed'.format(gen), f=0.5, F=0.5, wao=True):
                if int(x[-1])==0: continue
                else: thouGen[(x[0],x[1],x[2],x[3])]=(x[-2],format(float(x[-1])/(int(x[2])-int(x[1])+1.0),'.2f'))
	for x in dups.intersect('resources/annotation_files/{}_1000Genomes_DUP.bed'.format(gen), f=0.5, F=0.5, wao=True):
                if int(x[-1])==0: continue
                else: thouGen[(x[0],x[1],x[2],x[3])]=(x[-2],format(float(x[-1])/(int(x[2])-int(x[1])+1.0),'.2f'))
	return thouGen
def geneOverlap(cnv,genos,gen):
	genes={}
	starts=Bed([(x[0],x[1],x[1]) for x in genos]).sort()
        ends=Bed([(x[0],x[2],x[2]) for x in genos]).sort()
        GENEFH='resources/annotation_files/{}_genes.bed'.format(gen)
        for (c,s,e,c2,s2,e2,gene) in starts.intersect(GENEFH,wa=True,wb=True):
		x = gene.split(',')
                gene = ','.join(map(str,(x[0],x[2],x[3])))
                if genes.get((c,s))==None: genes[(c,s)]=gene
                elif gene not in genes[(c,s)]: genes[(c,s)] = genes[(c,s)]+'|'+gene
                else: continue
        for (c,s,e,c2,s2,e2,gene) in ends.intersect(GENEFH,wa=True,wb=True):
                x = gene.split(',')
                gene = ','.join(map(str,(x[0],x[2],x[3])))
                if genes.get((c,e))==None: genes[(c,e)]=gene
                elif gene not in genes[(c,e)]: genes[(c,e)] = genes[(c,e)]+'|'+gene
                else: continue
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
                                exonTotal = int(gene[-1].split('/').pop())
                                exonNum = int(gene[-1].split('/').pop(0).replace('exon_',''))
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
                                intronTotal = int(gene[-1].split('/').pop())
                                intronNum = int(gene[-1].split('/').pop(0).replace('intron_',''))
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
def annotate(genos,gen,REF,NON,GQ,OFH,sex,hemi):
	genes={}
	cyto={}
	mei={}
	flags={}
	thouGen={}
	AF={}
	GT={}
	refcnt={}
	males = [k for k in sex if sex[k] == 'M']
	IID = list(set([x[5] for x in genos]))
	IID.sort(key=str.lower)
	cnv=Bed(list(set([(x[0],x[1],x[2]) for x in genos]))).sort()
	cyto=cytobandOverlap(cnv,gen)
	flags=flagsOverlap(cnv,gen)
	mei=meiOverlap(cnv,gen)
	thouGen=thouGenOverlap(genos,gen)
	genes=geneOverlap(cnv,genos,gen)
	for x in genos:
		k = (x[0],x[1],x[2],x[4])
		if GQ.get((x[0],x[1],x[2],x[4],x[5]))!=None:
			lik=format(float(x[-2]),'.2f')
			if '1' not in x[9]: 
				if refcnt.get(k)==None: refcnt[k]=1
				else: refcnt[k]+=1
			if '1' in x[9]:
				if AF.get((k,x[5]))==None:
					if AF.get(k)==None: AF[k]=1
					else: AF[k]+=1
					AF[(k,x[5])]=1
				lik=format(float(x[-1]),'.2f')
			v = ':'.join((x[9],format(float(x[7]),'.2f'),format(float(x[8]),'.2f'),format(float(x[6]),'.2f'),lik,GQ[(x[0],x[1],x[2],x[4],x[5])]))	
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
	for (c,s,e,cl) in Bed(list(set([(x[0],x[1],x[2],x[4]) for x in genos]))).sort():
		SZ=int(e)-int(s)+1
		TYPE=cl
		THOUGEN_ID='NA'
		THOUGEN_OVR='NA'
		SEGD=0.00
		ABPTS=0.00
		UNMAP=0.00
		CENTMER=0.00
		STR=0.00
		BP1='intergenic'
		BP2='intergenic'
		GENE='intergenic'
		CYTOB='NA'
		MEDREF='NA'
		QUAL='NA'
		ALLELE='NA'
		STRANDS='+-'
		GTS=[]
		PASS='PASS'
		fail=[]
		cnvref=0
		if refcnt.get((c,s,e,cl))!=None: cnvref=refcnt[(c,s,e,cl)]
		if 'DUP' in cl: STRANDS='-+'
		if mei.get((c,s,e))!=None: TYPE=TYPE+':'+mei[(c,s,e)]
		DESX='{}:{}-{}_{}'.format(c,s,e,TYPE.replace('DEL','deletion').replace('DUP','duplication'))
		if cyto.get((c,s,e))!=None: CYTOB=cyto[(c,s,e)]
		if flags.get((c,s,e,'abparts'))!=None: ABPTS=format(flags[(c,s,e,'abparts')],'.2f')
		if flags.get((c,s,e,'centromere'))!=None: CENTMER=format(flags[(c,s,e,'centromere')],'.2f')
		if flags.get((c,s,e,'segDup'))!=None: SEGD=format(flags[(c,s,e,'segDup')],'.2f')
		if flags.get((c,s,e,'STR'))!=None: STR=format(flags[(c,s,e,'STR')],'.2f')
		if flags.get((c,s,e,'unmapable'))!=None: UNMAP=format(flags[(c,s,e,'unmapable')],'.2f')
		if thouGen.get((c,s,e,cl))!=None: THOUGEN_ID,THOUGEN_OVR = thouGen[(c,s,e,cl)] 
		if genes.get((c,s))!=None: BP1=genes[(c,s)]
		if genes.get((c,e))!=None: BP2=genes[(c,e)]
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
		if QUAL != 'NA' and int(QUAL) < 12 and 'DEL' in cl: fail.append('GQ-FAIL')
		if QUAL != 'NA' and int(QUAL) < 6 and 'DUP' in cl: fail.append('GQ-FAIL')
		if len(fail) > 0: PASS = ','.join(fail)
		INFO = 'SVTYPE={};END={};SVLEN={};STRANDS={};AF={};CYTOBAND={};1000G_ID={};1000G_OVERLAP={};DESCRIPTION={};BP1={};BP2={};GENES={};ABPARTS={};CENTROMERE={};SEGDUP={};STR={};UNMAPABLE={};MEDREFGL={}'.format(TYPE,e,SZ,STRANDS,ALLELE,CYTOB,THOUGEN_ID,THOUGEN_OVR,DESX,BP1,BP2,GENE,ABPTS,CENTMER,SEGD,STR,UNMAP,MEDREF)
		for x in IID:
			if GT.get(((c,s,e,cl),x))!=None: GTS.append(GT[((c,s,e,cl),x)])
		out= '\t'.join(map(str,(c.replace('chr',''),s,CNV_ID,'.','<'+TYPE+'>',QUAL,PASS,INFO,'GT:PE:SR:CN:SQ:GL','\t'.join(GTS)))) 
		VCF.append(out)
		CNV_ID+=1
	date=[]
	date.append(datetime.date.today())
	VCFHEAD=[	'##fileformat=VCFv4.2',
			'##fileDate={}'.format(date[0]),
			'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
			'##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
			'##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency,in the range (0,1)">',
			'##INFO=<ID=CYTOBAND,Number=1,Type=String,Description="cytoBand of the variant, pipe separated if more than one">',
			'##INFO=<ID=1000G_ID,Number=1,Type=String,Description="ID of variant in 1000 Genomes Phase 3 integrated SV callset">',
			'##INFO=<ID=1000G_ID,Number=1,Type=Float,Description="Overlap of variant with 1000 Genomes Phase 3 integrated SV callset, in the range (0,1)">',
			'##INFO=<ID=DESCRIPTION,Number=1,Type=String,Description="Verbose description of SV, particularly important for complex SVs"',
			'##INFO=<ID=BP1,Number=1,Type=String,Description="Location of the first breakpoint, pipe-separated list of all transcripts affected, with comma separated information about the transcript, refGene name, and location with respect to genes"',
			'##INFO=<ID=BP2,Number=1,Type=String,Description="Location of the second breakpoint, pipe-separated list of all transcripts affected, with comma separated information about the transcript, refGene name, and location with respect to genes"',
			'##INFO=<ID=GENES,Number=1,Type=String,Description="Description of genes within this SV locus, pipe-separated list of all transcripts affected, with comma separated information about the transcript, refGene name, and location with respect to genes"',
			'##INFO=<ID=ABPARTS,Number=A,Type=Float,Description="Parts of antibodies overlap, in the range (0,1)">',
			'##INFO=<ID=CENTROMERE,Number=A,Type=Float,Description="Centromere overlap, in the range (0,1)">',
			'##INFO=<ID=SEGDUP,Number=A,Type=Float,Description="Segmental duplication overlap, in the range (0,1)">',
			'##INFO=<ID=STR,Number=A,Type=Float,Description="Short Tandem Repeat overlap, in the range (0,1)">',
			'##INFO=<ID=UNMAPABLE,Number=A,Type=Float,Description="Overlap with regions unmapable using 100bp reads, in the range (0,1)">',
			'##INFO=<ID=MEDREFGL,Number=A,Type=Float,Description="Median Phred-scaled genotype likelihood for individuals genotyped homozygous reference">',
			'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
			'##FORMAT=<ID=PE,Number=1,Type=Float,Description="Ratio of Discordant Paired-Ends to Concordant Paired-Ends">',
			'##FOTMAT=<ID=SR,Number=1,Type=Float,Description="Ratio of Split reads to Concordant Paired-Ends">',
			'##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy Number">',
			'##FORMAT=<ID=SQ,Number=1,Type=Float,Description="Phred-scaled genotype likelihood">',
			'##FORMAT=<ID=GL,Number=2|3,Type=Float,Description="Phred-scaled genotype likelihood; homozygous alt, heterogygous alt, homozygous ref">',
			'##ALT=<ID=DUP,Description="Duplication, if 80% reciprocal overlap with transposable element then TE class, family, and name given separated by colons">',
			'##ALT=<ID=DEL,Description="Deletion, if 80% reciprocal overlap with transposable element then TE class, family, and name given separated by colons">'
		]
	with open('resources/{}.chrom.sizes'.format(gen)) as f:
		for l in f:
			(chrom,size) = l.rstrip('\n').split('\t')
			VCFHEAD.append('##contig=<ID={},length={}>'.format(chrom.replace('chr',''),size))
	outfh = open(OFH,'w')
	outfh.write('\n'.join(VCFHEAD)+'\n')
	outfh.write('\t'.join(('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','\t'.join(IID)))+'\n')
	outfh.write('\n'.join(VCF)+'\n')
	outfh.close()
