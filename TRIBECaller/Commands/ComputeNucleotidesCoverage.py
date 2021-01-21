
# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue

from TRIBECaller.TribeCaller import *


def compute_nucleotides_coverage(target_path:str,control_path:str,out_prefix:str):
	TEC = TribeCaller(target_path,control_path)
	result = TEC.run_editing_percentage()
	res = []
	with open(out_prefix + ".bed", "w+") as f:
		for k,v in result.items():
			for j in v:
				f.write("\t".join([k,str(j[0]),str(j[1]), '\t'.join(list(map(str, j[2]))), '\t'.join(list(map(str, j[3])))])+ "\n")