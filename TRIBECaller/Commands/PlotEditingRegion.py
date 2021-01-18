# -*- coding: utf-8 -*-
#!/usr/bin/env python
#
# ------------------------- #
# own Python Modules
# ------------------------- #

from TRIBECaller.TribeCaller import *
from TRIBECaller.GtfReads import *
from TRIBECaller.Plotting.PlotNucleotides import *
from TRIBECaller.Plotting.PlotGene import *
from TRIBECaller.Utilities.Parsers import parse_region
from TRIBECaller.Utilities.utils import *


# ------------------------- #
# Python Modules
# ------------------------- #

import matplotlib.pyplot as plt
from matplotlib import gridspec

def plot_editing_sites(target_path:str,
					   control_path:str,
					   region:str,
					   output_path:str=None,
					   call_editing_sites=True,
					   dpi=300):
	TEC = TribeCaller(target_path,control_path)
	chr_,start,end=parse_region(region)
	if call_editing_sites:
		region,perc,odds,pval=TEC.compute_region_as_percentage(chr_,start,end,call_editing_sites=True)
	else:
		region,perc=TEC.compute_region_as_percentage(chr_,start,end)
	fig,(ax1,ax2)=plt.subplots(2,1,figsize=(16,6))
	plt.subplots_adjust(hspace=1,bottom=.2)
	ax1.set_xbound(region[0],region[-1]+1)
	ax2.set_xbound(region[0],region[-1]+1)
	ax1.get_yaxis().set_ticks([])
	ax2.get_yaxis().set_ticks([])
	ax1.get_xaxis().set_ticks(list(map(lambda x:x+0.5,region)))
	ax2.get_xaxis().set_ticks(list(map(lambda x:x+0.5,region)))
	ax1.set_xlabel("chromosome "+chr_,fontfamily="Arial")
	ax2.set_xlabel("chromosome "+chr_,fontfamily="Arial")
	ax1.set_xticklabels(list(map(lambda x:str(x+1),region)),fontfamily="Arial")
	ax2.set_xticklabels(list(map(lambda x:str(x+1),region)),fontfamily="Arial")   
	ax1.set_title("Experient",fontsize=16,fontweight=600,fontfamily="Arial")  
	ax2.set_title("Control",fontsize=16,fontweight=600,fontfamily="Arial")  

	for ind,(pos,i) in enumerate(zip(region,perc)):
		y1=0
		y2=0
		if call_editing_sites:
			if odds[ind] > 2 and pval[ind] < 0.05 and perc[ind][1][0]>0.8:
				#with np.errstate(invalid='ignore'):
				#	where=(np.array(odds) > 2) & (np.array(pval) < 0.05) & np.array(list(map(lambda x:x[1][0] > 0.8, perc)))
				ax1.axvspan(region[ind],region[ind]+1,color="yellow")
		for n,j in enumerate(i[0]):
			ax1.add_patch(make_nucleotides_elements(NUC[n][0],x=pos,y=y1,height=j,color=NUC[n][1]))
			y1 += j
		for n,j in enumerate(i[1]):
			ax2.add_patch(make_nucleotides_elements(NUC[n][0],x=pos,y=y2,height=j,color=NUC[n][1]))
			y2 += j
	for i in ["left","top","right"]:
		ax1.spines[i].set_visible(False)
		ax2.spines[i].set_visible(False)
	ax1.set_xbound(region[0],region[-1]+1)
	ax2.set_xbound(region[0],region[-1]+1)
	ax1.tick_params(axis='x',labelrotation=45)
	ax2.tick_params(axis='x',labelrotation=45)
	if output_path:
		fig.savefig(output_path,dpi=dpi)
	return


def plot_editing_region(target_path:str,
					   control_path:str,
					   gtf_path:str=None,
					   region:str=None,
					   gene_ensembl_id=None,
					   gene_symbol=None,
					   output_path:str=None,
					   call_editing_sites=True,
					   dpi=300):
	TEC = TribeCaller(target_path,control_path)
	print(GET_CUR_TIME("Fetching GTF file"))
	gr = GtfReads(gtf_path)
	if region:
		result = gr.fetch(region=region)
		chr_,start,end=parse_region(region)
	elif gene_ensembl_id:
		result = gr.fetch(gene_ensembl_id = gene_ensembl_id)
		chr_,start,end=result[0][gr.Flag.Chromosome], min(list(map(lambda i:int(i[gr.Flag.Start]), result))),max(list(map(lambda i:int(i[gr.Flag.End]), result)))
	elif gene_symbol:
		result = gr.fetch(gene_symbol = gene_symbol)
		chr_,start,end=result[0][gr.Flag.Chromosome], min(list(map(lambda i:int(i[gr.Flag.Start]), result))),max(list(map(lambda i:int(i[gr.Flag.End]), result)))
	else:
		raise ValueError("You must provide either a genomic region, gene ensembl id or gene symbol")
	print(GET_CUR_TIME("Computing coverage and editing sites"))
	reg, cov, edi = TEC.call_editing_region_coverage(chr_,start,end)
	plt.subplots_adjust(hspace=1,bottom=.2)
	print(GET_CUR_TIME("Plotting coverage and editing sites"))
	gs = gridspec.GridSpec(4, 1, height_ratios=[3,1,1,1]) 
	fig=plt.Figure(figsize=(10,8))
	ax1 = fig.add_subplot(gs[0])
	ax2 = fig.add_subplot(gs[1])
	ax3 = fig.add_subplot(gs[2])
	ax4 = fig.add_subplot(gs[3])
	y_min = make_gene_elements(result,1,ax1)
	ax2.bar(height=list(map(lambda x:x[0],cov)),x=reg)
	ax4.bar(height=list(map(lambda x:x[1],cov)),x=reg)
	ax3.scatter(x=reg,y=list(map(lambda x:np.random.random() if x else None, edi)), marker='^',color='red')
	for i in ["left","top","right"]:
		ax1.spines[i].set_visible(False)
		ax2.spines[i].set_visible(False)
		ax3.spines[i].set_visible(False)
		ax4.spines[i].set_visible(False)
	ax1.spines['bottom'].set_visible(False)
	ax3.spines['bottom'].set_visible(False)
	ax1.get_yaxis().set_ticks([])
	ax2.get_yaxis().set_ticks([])
	ax3.get_yaxis().set_ticks([])
	ax4.get_yaxis().set_ticks([])
	ax1.get_xaxis().set_ticks([])
	ax2.get_xaxis().set_ticks([])
	ax3.get_xaxis().set_ticks([])
	ax1.set_ybound(y_min-0.2,1.2) if y_min < 0 else ax1.set_ybound(y_min-0.2,1.2) 
	x_min,x_max=ax1.get_xlim()
	ax2.set_xbound(x_min,x_max)
	ax3.set_xbound(x_min,x_max)
	ax4.set_xbound(x_min,x_max)
	ticks_range = PICK(list(range(int(x_min),int(x_max))),(x_max-x_min) // 10)
	ax4.get_xaxis().set_ticks(ticks_range)
	ax4.set_xticklabels(list(map(lambda x:str(x),ticks_range)),fontfamily="Arial")
	ax4.set_xlabel("chromosome "+chr_,fontfamily="Arial",fontweight=600)
	if output_path:
		fig.savefig(output_path,dpi=dpi)
	return 