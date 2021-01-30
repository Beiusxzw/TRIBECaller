
# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue

from TRIBECaller.TribeCaller import *


def compute_coverage(target_path:str,control_path:str,out_prefix:str,g_zip:bool, contig=None,n_threads=None, paired=False, exclude_gap=False):
	TEC = TribeCaller(target_path,control_path)
	if n_threads:
		TEC.run_editing_percentage_par(out_prefix,g_zip=g_zip,contig=contig, n_threads=n_threads,paired=paired,exclude_gap=exclude_gap)
	else:
		TEC.run_editing_percentage(out_prefix,g_zip=g_zip,contig=contig,paired=paired,exclude_gap=exclude_gap)