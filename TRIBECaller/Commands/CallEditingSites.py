# -*- coding: utf-8 -*-

from TRIBECaller.TribeCaller import *


def call_editing_sites(target_path:str,control_path:str,out_prefix:str):
	TEC = TribeCaller(target_path,control_path)
	result = TEC()
	res = []
	with open(out_prefix + ".bed", "w+") as f:
		for k,v in result.items():
			for j in v:
				f.write("\t".join([k,str(j[0]+1),str(j[0]+2),str(j[1]),str(j[2])]) + "\n")