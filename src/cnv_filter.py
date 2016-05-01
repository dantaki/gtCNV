#!/usr/env python
import pybedtools as pbed
import bed
def reduce_dict(a):
	"""reduces the masked regions to a dictionary of lists"""
	d={}
	for i in a:
		(c,s,e,cl,tag) = i
		if d.get(tag) == None: d[tag]=[(c,s,e)]
		else: d[tag].append((c,s,e))
	return d
def filter_calls(cnv,gen):
        """add unique annotation to each CNV"""
        tag = 1
        annot_cnv=[]
        flank_pos={}
        start_flank=[]
        end_flank=[]
        for r in cnv:
                (c,s,e,cl) = r
                annot_cnv.append((c,s,e,cl,tag))
                start_flank.append((c,s,s,cl,tag))
                end_flank.append((c,e,e,cl,tag))
                tag+=1
        annot_dict={}
	for i in annot_cnv: annot_dict[i[4]]=(i[0],i[1],i[2],i[3])
	"""add 500bp to each start and end position"""
        start_flank500bp=pbed.BedTool(start_flank).slop(b=500,g='resources/'+gen+'.bedtools.genome')
        end_flank500bp=pbed.BedTool(end_flank).slop(b=500,g='resources/'+gen+'.bedtools.genome')
	start_flank_dict=dict((i[4],(i[1],i[2])) for i in start_flank500bp)
	end_flank_dict=dict((i[4],(i[1],i[2])) for i in end_flank500bp)
	"""mask regions to segmental duplications and unmappable regions"""
        masked_cnv=bed.annot_toList(pbed.BedTool(annot_cnv).subtract(pbed.BedTool('resources/'+gen+'_unmapped.bed')))
        masked_start_flank=bed.annot_toList(start_flank500bp.subtract(pbed.BedTool('resources/'+gen+'_unmapped.bed')))
        masked_end_flank=bed.annot_toList(end_flank500bp.subtract(pbed.BedTool('resources/'+gen+'_unmapped.bed')))
	"""take the union of the two masked call sets"""
        if len(masked_cnv) == 0:
		import sys
		print 'ERROR: CNVs completed masked by excluded regions'
		sys.exit()
	cnv_tag = list(zip(*masked_cnv)[4])
	passed_tag = list(set(cnv_tag))
	masked_cnv_dict=reduce_dict(masked_cnv)
	masked_start_flank_dict=reduce_dict(masked_start_flank)
	masked_end_flank_dict=reduce_dict(masked_end_flank)
	master_cnv={}
	master_cnv_list=[]
	for i in passed_tag:
		windows = (start_flank_dict[i],end_flank_dict[i])
		cnv_span = masked_cnv_dict[i]
		flank_span=None
		if masked_start_flank_dict.get(i) != None and masked_end_flank_dict.get(i) != None:
			flank_span = masked_start_flank_dict[i] + masked_end_flank_dict[i]
		master_cnv[annot_dict[int(i)]]=(cnv_span,flank_span,windows)
		master_cnv_list.append(annot_dict[int(i)])
	final_sorted=[]
	for i in pbed.BedTool(master_cnv_list).sort(): final_sorted.append((i[0],i[1],i[2],i[3])) 
	return final_sorted,master_cnv		
