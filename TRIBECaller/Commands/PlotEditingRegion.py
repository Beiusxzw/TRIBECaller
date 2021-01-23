# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
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
import seaborn as sns

def plot_editing_sites(target_path:str,
					   control_path:str,
					   region:str,
					   output_path:str=None,
					   call_editing_sites=True,
					   a_content_threshold=0.8,
					   dpi=300):
	TEC = TribeCaller(target_path,control_path)
	chr_,start,end=parse_region(region)
	if call_editing_sites:
		region,perc,edi=TEC.compute_region_as_percentage(chr_,start,end,call_editing_sites=True)
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
			if edi[ind]:
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
	reg, cov, edi_f, edi_r, diff = TEC.call_editing_region_coverage(chr_,start,end)
	plt.subplots_adjust(hspace=1,bottom=.2)
	print(GET_CUR_TIME("Plotting coverage and editing sites"))
	gs = gridspec.GridSpec(5, 1, height_ratios=[3,.5,1,.5,.5]) 
	fig=plt.Figure(figsize=(12,8))
	ax1 = fig.add_subplot(gs[0])
	ax2 = fig.add_subplot(gs[1])
	ax3 = fig.add_subplot(gs[2])
	ax4 = fig.add_subplot(gs[3])
	ax5 = fig.add_subplot(gs[4])
	y_min = make_gene_elements(result,1,ax1)
	ax2.bar(height=list(map(lambda x:x[0],cov)),x=reg)
	ax5.bar(height=list(map(lambda x:x[1],cov)),x=reg,color="#FDBE84")
	ax3.scatter(x=reg,y=list(map(lambda x:np.random.random() if x else None, edi_f)), marker='$+$',color='red',s=3)
	ax3.scatter(x=reg,y=list(map(lambda x:np.random.random() if x else None, edi_r)), marker='$-$',color='blue',s=3)
	ax4.bar(x=reg,height=diff)
	#sns.kdeplot(reg,diff,ax=ax4)
	for i in ["left","top","right"]:
		ax1.spines[i].set_visible(False)
		ax2.spines[i].set_visible(False)
		ax3.spines[i].set_visible(False)
		ax5.spines[i].set_visible(False)
	ax1.spines['bottom'].set_visible(False)
	ax3.spines['bottom'].set_visible(False)
	ax1.set_ylabel("Gene",rotation = 90,fontfamily="Arial",fontsize=6)
	ax2.set_ylabel("TRIBE RNA-seq \ncoverage",rotation = 90,fontfamily="Arial",fontsize=6,fontweight=600)
	ax3.set_ylabel("Editing Events",rotation = 90,fontfamily="Arial",fontsize=6,fontweight=600)
	ax4.set_ylabel("Relative S/A \npercentage",rotation = 90,fontfamily="Arial",fontsize=6,fontweight=600)
	ax5.set_ylabel("Control RNA-seq \ncoverage",rotation = 90,fontfamily="Arial",fontsize=6,fontweight=600)
	ax1.get_yaxis().set_ticks([])
	ax2.get_yaxis().set_ticks([])
	ax3.get_yaxis().set_ticks([])
	ax5.get_yaxis().set_ticks([])
	ax1.get_xaxis().set_ticks([])
	ax2.get_xaxis().set_ticks([])
	ax3.get_xaxis().set_ticks([])
	ax4.get_xaxis().set_ticks([])
	ax1.set_ybound(y_min-0.2,1.2) if y_min < -1 else ax1.set_ybound(y_min-0.2,1.2) 
	x_min,x_max=ax1.get_xlim()
	ax2.set_xbound(x_min,x_max)
	ax3.set_xbound(x_min,x_max)
	ax4.set_xbound(x_min,x_max)
	ax5.set_xbound(x_min,x_max)
	ticks_range = PICK(list(range(int(x_min),int(x_max))),(x_max-x_min) // 10)
	ax5.get_xaxis().set_ticks(ticks_range)
	ax5.set_xticklabels(list(map(lambda x:str(x),ticks_range)),fontfamily="Arial")
	ax5.set_xlabel("chromosome "+chr_,fontfamily="Arial",fontweight=600)
	if output_path:
		fig.savefig(output_path,dpi=dpi)
	return 