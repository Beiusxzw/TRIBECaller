# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue

from TRIBECaller.TribeCaller import *
from TRIBECaller.TribeCriteria import TribeCriteria

def call_editing_sites(target_path:str,
                       control_path:str,
                       out_prefix:str,
                       g_zip:bool, 
                       criteria,
                       contig=None, 
                       verbose=False,
                       n_threads=None):
    TEC = TribeCaller(target_path,control_path)
    tc = TribeCriteria()
    for k,v in criteria.items():
    	if v:
    		tc.set_args(k,v)
    if n_threads:
        TEC.run_par(tc, out_prefix, contig=contig,g_zip=g_zip,n_threads=n_threads,verbose=verbose)
    else:
        TEC.run(tc, out_prefix, contig=contig,g_zip=g_zip,verbose=verbose)
