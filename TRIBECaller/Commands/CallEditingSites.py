# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue

from TRIBECaller.TribeCaller import *


def call_editing_sites(target_path:str,control_path:str,out_prefix:str,g_zip:bool, contig=None, n_threads=None):
	TEC = TribeCaller(target_path,control_path)
	if contig:
		if n_threads:
			TEC.run_par(out_prefix, contig=contig,g_zip=g_zip)
		else:
			TEC.run(out_prefix, contig=contig,g_zip=g_zip)

	else:
		if n_threads:
			TEC.run_par(out_prefix, g_zip=g_zip)
		else:
			TEC.run(out_prefix, g_zip=g_zip)
