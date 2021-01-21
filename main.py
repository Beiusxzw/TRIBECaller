# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
#
# Description: TribeCaller main executable
#
#
# ------------------------- #
# Python Modules
# ------------------------- #

import argparse
import os
import sys

# ------------------------- #
# Own Python Modules
# ------------------------- #

from TRIBECaller._version import __version__
from TRIBECaller.Utilities.utils import *

# ------------------------- #
# Main Function
# ------------------------- #

def main():
	parser = prepare_argparser()
	args = parser.parse_args()
	subcommand = args.subcommand
	PRINT_PRELOGUE()
	PRINT_LOGO()
	PRINT_INFO()
	if subcommand == "callEditingSites":
		from TRIBECaller.Commands.CallEditingSites import call_editing_sites
		print(GET_CUR_TIME("Start running program " + GET_BLUE("callEditingSites")))
		call_editing_sites(args.target,args.control,args.outPrefix)
		print(GET_CUR_TIME("Running completed."))
	elif subcommand == "computeCoverage":
		from TRIBECaller.Commands.ComputeCoverage import compute_coverage
		print(GET_CUR_TIME("Start running program " + GET_BLUE("computeCoverage")))
		compute_coverage(args.target,args.control,args.outPrefix)
		print(GET_CUR_TIME("Running completed."))
	elif subcommand == "plotEditingSite":
		from TRIBECaller.Commands.PlotEditingRegion import plot_editing_sites
		print(GET_CUR_TIME("Start running program " + GET_BLUE("plotEditingSite")))
		dpi = args.dpi if args.dpi else 300
		plot_editing_sites(args.target,args.control,region=args.region,output_path=args.out,call_editing_sites=True,dpi=dpi)
		print(GET_CUR_TIME("Running completed."))
	elif subcommand == "plotEditingRegion":
		from TRIBECaller.Commands.PlotEditingRegion import plot_editing_region
		print(GET_CUR_TIME("Start running program " + GET_BLUE("plotEditingRegion")))
		dpi = args.dpi if args.dpi else 300
		if args.region:
			plot_editing_region(args.target,args.control,gtf_path=args.gtf,region=args.region,output_path=args.out,call_editing_sites=True,dpi=dpi)
		elif args.gene:
			gene = args.gene
			if gene.startswith("ENS"):
				gene = list(map(lambda x:x.strip(), gene.split(','))) if ',' in gene else gene.strip()
				plot_editing_region(args.target,args.control,gtf_path=args.gtf,gene_ensembl_id=gene,output_path=args.out,call_editing_sites=True,dpi=dpi)
			else:
				gene = list(map(lambda x:x.strip(), gene.split(','))) if ',' in gene else gene.strip()
				plot_editing_region(args.target,args.control,gtf_path=args.gtf,gene_symbol=gene,output_path=args.out,call_editing_sites=True,dpi=dpi)
		print(GET_CUR_TIME("Running completed."))
	else:
		raise TypeError("Unknown program name")

def prepare_argparser():
	"""
	Prepare optparser object.
	"""
	desc = "Editing site caller for TRIBE"
	epilog = "For command line options for each command, type COMMAND -h"
	parser = argparse.ArgumentParser(description = desc, epilog = epilog)
	parser.add_argument('-v', '--version', action="version",version="TRIBECaller\t" + __version__)
	subparsers = parser.add_subparsers(dest='subcommand')
	subparsers.required = True

	add_call_editing_sites_parser(subparsers)
	add_plot_editing_site_parser(subparsers)
	add_plot_editing_region_parser(subparsers)
	add_compute_coverage_parser(subparsers)
	return parser

def add_call_editing_sites_parser(subparsers):
	"""
	Add main function 'call editing sites' argument parsers
	"""
	argparser_call_editing_sites = subparsers.add_parser("callEditingSites", help="Main TRIBECaller function",epilog="""Examples:
python main.py callEditingSites -t Experiment_rep1_1_srt.bam -c Control_rep1_1_srt.bam -o testPrefix
	""")

	argparser_call_editing_sites.add_argument('-t', '--target', type=str, help='input target bam file',required=True)
	argparser_call_editing_sites.add_argument('-c', '--control', type=str, help='input control bam file',required=True)
	argparser_call_editing_sites.add_argument('-b', '--binSize', type=int, default=1,help='bin size of calling editing event')
	argparser_call_editing_sites.add_argument('-o', '--outPrefix', type=str, default = "TribeCallerOutput", help='prefix of output bed file',required=True)
	return 

def add_compute_coverage_parser(subparsers):
	"""
	Add main function 'call editing sites' argument parsers
	"""
	argparser_compute_coverage = subparsers.add_parser("computeCoverage", help="Compute coverage and nucleotides content",epilog="""Examples:
python main.py computeCoverage -t Experiment_rep1_1_srt.bam -c Control_rep1_1_srt.bam -o testPrefix
	""")

	argparser_compute_coverage.add_argument('-t', '--target', type=str, help='input target bam file',required=True)
	argparser_compute_coverage.add_argument('-c', '--control', type=str, help='input control bam file',required=True)
	argparser_compute_coverage.add_argument('-o', '--outPrefix', type=str, default = "TribeCallerOutput", help='prefix of output bed file',required=True)
	return 

def add_plot_editing_site_parser(subparsers):
	"""
	Add main function 'plot editing site' argument parsers
	"""
	argparser_plot_editing_sites = subparsers.add_parser("plotEditingSite", help="Plot nucleotides from reads within a genomic region",epilog="""Examples:
python main.py plotEditingSite -t Experiment_rep1_1_srt.bam -c Control_rep1_1_srt.bam --region 1:1214232:1214243
	""")
	argparser_plot_editing_sites.add_argument('-t', '--target', type=str, help='input target bam file',required=True)
	argparser_plot_editing_sites.add_argument('-c', '--control', type=str, help='input control bam file',required=True)
	argparser_plot_editing_sites.add_argument('-r', '--region', type=str, help='genomic region of interests. Format: 1:1072894-1078023',required=True)
	argparser_plot_editing_sites.add_argument('-o', '--out', type=str,  help='prefix of output figure',required=True)
	argparser_plot_editing_sites.add_argument('--dpi', type=int, help='dpi of the output figure')
	return

def add_plot_editing_region_parser(subparsers):
	"""
	Add main function 'plot editing region' argument parsers
	"""
	argparser_plot_editing_sites = subparsers.add_parser("plotEditingRegion", help="Plot editing events within a genomic region",epilog="""Examples:
python main.py plotEditingRegion -t Experiment_rep1_1_srt.bam -c Control_rep1_1_srt.bam --gene POU5F1
	""")
	argparser_plot_editing_sites.add_argument('-t', '--target', type=str, help='input target bam file',required=True)
	argparser_plot_editing_sites.add_argument('-c', '--control', type=str, help='input control bam file',required=True)
	argparser_plot_editing_sites.add_argument('-r', '--region', type=str, help='genomic region of interests. Format: 1:1072894-1078023',required=False)
	argparser_plot_editing_sites.add_argument('-g', '--gene', type=str, help='gene of intersets. For multiple genes, use comma separated input.',required=False)
	argparser_plot_editing_sites.add_argument('-gtf', '--gtf', type=str, help='gtf file to plot the genes.',required=True)
	argparser_plot_editing_sites.add_argument('-o', '--out', type=str,  help='prefix of output figure',required=True)
	argparser_plot_editing_sites.add_argument('--dpi', type=int, help='dpi of the output figure')
	return
	
if __name__ == '__main__':
	__spec__ = None
	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("Exiting due to user interruption")
		sys.exit(0)
	except MemoryError:
		sys.exit("MemoryError occured")
		sys.exit(1)

