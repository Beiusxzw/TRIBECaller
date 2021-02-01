# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue

from TRIBECaller.GtfReads import *

from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.transforms import Affine2D
from matplotlib.lines import Line2D

import numpy as np

def make_gene_elements(gene_gtfs, y_pos, ax, enable_arrow=False):
	result = []
	flag=False
	y_pos_min = y_pos
	pos_dict = {}
	ga_esti = (int(gene_gtfs[-1][GtfReads.Flag.End]) - int(gene_gtfs[0][GtfReads.Flag.Start])) // 10
	for i in gene_gtfs:
		mark = '$>$' if i[GtfReads.Flag.Strand] == '+' else '$<$'
		if i[GtfReads.Flag.Feature] == "transcript":
			# The while loop is to calculate the y-position of each transcripts. 
			# If there is no overlap between two transcripts in the same level, 
			# break and the y-position will return.
			if flag:
				y_pos=1
				while y_pos in pos_dict.keys():
					if ((int(i[GtfReads.Flag.Start]) >  pos_dict[y_pos][0]) and (int(i[GtfReads.Flag.End]) < pos_dict[y_pos][1])) or ((int(i[GtfReads.Flag.End]) > pos_dict[y_pos][0]) and (int(i[GtfReads.Flag.Start]) <  pos_dict[y_pos][1])) :
						y_pos -= 0.3
						y_pos_min = min(y_pos_min, y_pos)
					else:
						pos_dict[y_pos] = (min(pos_dict[y_pos][0], int(i[GtfReads.Flag.Start])), max(pos_dict[y_pos][1], int(i[GtfReads.Flag.End])))
						break
				if y_pos not in pos_dict.keys():
					pos_dict[y_pos] = (int(i[GtfReads.Flag.Start]), int(i[GtfReads.Flag.End]))
			else:
				pos_dict[y_pos] = (int(i[GtfReads.Flag.Start]), int(i[GtfReads.Flag.End]))
			gene_length = int(i[GtfReads.Flag.End]) - int(i[GtfReads.Flag.Start])
			cur_end = int(i[GtfReads.Flag.End])
			x = np.arange(int(i[GtfReads.Flag.Start]), int(i[GtfReads.Flag.End]), ga_esti)
			line = ax.plot(x, [y_pos] * len(x), marker=mark, fillstyle='none',color="#B8860B",mew=0.5,lw=.5,ms=5,mfc="#B8860B")
			x = np.arange(int(i[GtfReads.Flag.Start]), int(i[GtfReads.Flag.End]))
			line = ax.plot(x, [y_pos] * len(x), linewidth=1,color="#B8860B")
			if enable_arrow:
				ax.text(int(i[GtfReads.Flag.Start]), y_pos-0.1, i[GtfReads.Flag.Attribute]["gene_name"] + ">" if i[GtfReads.Flag.Strand] == '+' else i[GtfReads.Flag.Attribute]["gene_name"] + "<" , fontfamily="Arial", fontsize=4)
			flag=True
			pos_dict[y_pos] = (int(i[GtfReads.Flag.Start]), int(i[GtfReads.Flag.End]))
		elif i[GtfReads.Flag.Feature] == "exon":
			result.append(Rectangle((int(i[GtfReads.Flag.Start]), y_pos-0.05/2), int(i[GtfReads.Flag.End]) - int(i[GtfReads.Flag.Start]), 0.05, linewidth=1,facecolor='#B8860B'))
		elif i[GtfReads.Flag.Feature] == "CDS":
			result.append(Rectangle((int(i[GtfReads.Flag.Start]), y_pos-0.1/2), int(i[GtfReads.Flag.End]) - int(i[GtfReads.Flag.Start]), 0.1, linewidth=1,facecolor='#B8860B'))
		elif i[GtfReads.Flag.Feature] == "five_prime_utr" or i[GtfReads.Flag.Feature] == "three_prime_utr":
			result.append(Rectangle((int(i[GtfReads.Flag.Start]), y_pos-0.05/2), int(i[GtfReads.Flag.End]) - int(i[GtfReads.Flag.Start]), 0.02, linewidth=1,facecolor='#B8860B'))		
		else: 
			pass
	for i in result:
		ax.add_patch(i)
	return y_pos_min