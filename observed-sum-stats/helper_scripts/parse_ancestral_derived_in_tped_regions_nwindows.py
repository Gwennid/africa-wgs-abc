#Edit 20210531: Update the path.
import gzip
import sys




TPED_f=sys.argv[1]




def parse_info(info):
	d=info[1].split(':')
	return [int(d[0]),int(d[1])]

def parse_anc(info,anc_info,gt):
	gt_set=set(gt).difference(set(['0']))
	if not gt_set.issubset(nt):
		return ''
	if len(gt_set)==1:
		return ''
	if len(gt_set)==2:
		anc=anc_info[1].upper()
		ref=anc_info[2].upper()
		if anc in gt_set and ref in gt_set:
			out_str=' '.join(info)
			if ref==anc:
				for x in gt:
					if x=='0':
						out_str+=' N'
					elif x==anc:
						out_str+=' 0'
					else:
						out_str+=' 1'
				return out_str+'\n'
			else:
				for x in gt:
					if x=='0':
						out_str+=' N'
					elif x==anc:
						out_str+=' 0'
					elif x==ref:
						out_str+=' 1'
					else:
						return ''
				return out_str+'\n'
	return ''


f=open(TPED_f,'r')
parse_pos={}
for l in f:
	d=l.split()
	info=d[:4]
	[chr,pos]=parse_info(info)
	if not chr in parse_pos.keys():
		parse_pos.update({chr:[]})
	parse_pos[chr].append(pos)
f.close()


three_apes_path='/crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/Comparative_datasets/three_ape_state/DIR_per_chr/'

nt=set(['A','C','G','T'])

anc_state={}
for chr in parse_pos.keys():
	anc_state.update({chr:[]})
	#print chr,len(parse_pos[chr])
	three_apes_f=three_apes_path+'chr'+str(chr)+'.txt.gz'
	with gzip.open(three_apes_f,'rb') as apef:
		ape_d=apef.readline().split()
		ape_pos=0
		for pos in parse_pos[chr]:
			out_tuple=(pos,'N','N')
			
			while ape_pos<pos and ape_pos>-1:
				ape_d=apef.readline().split()
				if not ape_d:
					ape_pos=-1
				else:
					[ape_pos,ref,chimp,gor,oran]=ape_d
					ape_pos=int(ape_pos)
			if ape_pos>-1:
				ape_nt=set([chimp.upper(),gor.upper(),oran.upper()])
				if chimp.upper() in nt and len(ape_nt)==1:
					if ape_nt==set([chimp,gor,oran]):
						out_tuple=(pos,chimp,ref)
					else:
						out_tuple=(pos,chimp.lower(),ref)
			anc_state[chr].append(out_tuple)
	apef.close()
	

outf=open(TPED_f[:-len('.tped')]+'_ancestral_state.tped','w')
f=open(TPED_f,'r')
at_chr=0
for l in f:
	d=l.split()
	info=d[:4]
	[chr,pos]=parse_info(info)
	if not chr==at_chr:
		at_chr=chr
		the_anc_l=anc_state[at_chr]
	anc_info=the_anc_l.pop(0)
	if not 'N' in anc_info:
		gt=d[4:]
		out_str=parse_anc(info,anc_info,gt)
		if out_str:
			outf.write(out_str)	
f.close()
outf.close()



