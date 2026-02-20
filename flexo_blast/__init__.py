#!/usr/bin/env python3

import argparse
import contextlib
import datetime
import importlib.resources
import logging
import os
import sys
import shutil
import tempfile
import urllib.request
import warnings

from Bio import SeqIO

from .blast import Blast
from .print_final_results import AnnoFasta
from .print_final_results import FinalResults

__author__ = "Laura M. Carroll <lmc297@cornell.edu>"
__version__ = "0.2.0"


@contextlib.contextmanager
def _forward_warnings():
	_showwarning = warnings.showwarning
	try:
		warnings.showwarning = lambda message, category, filename, lineno, file=None, line=None: logging.warning(f"Warning: {message}")
		yield
	finally:
		warnings.showwarning = _showwarning


def run_pipeline(args):
	"""
	run flexo_blast
	"""

	# get path to flexo_blast executable
	flexo_blast_path = os.path.realpath(__file__)
	flexo_blast_path = flexo_blast_path.rpartition("/")[0].strip() + "/"

	# get input and output arguments and prefix to use for output files
	infile = args.input[0]
	prefix = infile.split("/")[-1].strip()
	prefix = ".".join(prefix.split(".")[0:-1])
	output = args.output[0]

	# add a / to output directory name, if one is not already supplied
	if not output.endswith("/"):
		output = output.strip() + "/"

	# if the output directory doesn't yet exist, make it
	if not os.path.isdir(output):
		os.mkdir(output)

	# if a flexo_blast_final_results directory does not exist
	# in the output directory, make it
	if not os.path.isdir(output + "flexo_blast_final_results"):
		os.mkdir(output + "flexo_blast_final_results/")

	# make a logs directory for log files
	if not os.path.isdir(output + "flexo_blast_final_results/logs"):
		os.mkdir(output + "flexo_blast_final_results/logs")

	# save path to final results directory
	final_results_directory = output + "flexo_blast_final_results/"

	# initialize log file
	logging.basicConfig(level = logging.DEBUG, filename = final_results_directory + "logs/" + prefix + ".log", filemode = "a+", format = "%(message)s")
	logging.getLogger().addHandler(logging.StreamHandler())


	# get arguments
	# evalue is stored as a string because it gets passed to blast command as a string
	multi = args.multi
	single = args.single
	pthresh = int(args.pthresh)
	qthresh = int(args.qthresh)
	overlap = float(args.overlap)
	evalue = str(args.evalue)

	# log to file
	now = datetime.datetime.now
	logging.info("You're running FLExo_BLAST!")
	logging.info("You are initializing this run at " + now().strftime("%Y-%m-%d %H:%M"))
	logging.info("You ran the following command: ")
	logging.info(" ".join([str(sa) for sa in sys.argv]))
	logging.info("Report bugs/concerns to Laura M. Carroll <lmc297@cornell.edu>")

	# detect multi-gene exotoxins
	if multi == "True":

		with importlib.resources.path("flexo_blast.seq_multi_db", "multiblast.faa") as multi_path:

			get_multi = Blast(
				task = "tblastn",
				dbseqs = infile,
				fasta = multi_path,
				final_results_directory = final_results_directory,
				prefix = prefix,
				suffix = "multi",
				pthresh = pthresh,
				qthresh = qthresh,
				overlap = overlap,
				evalue = evalue)

			logging.info("Using tblastn to identify multi-gene exotoxins in " + prefix + " at " + now().strftime("%Y-%m-%d %H:%M"))

			multi_results = get_multi.run_blast("tblastn", infile, multi_path, final_results_directory, prefix, "multi", evalue)
			multi_final = get_multi.parse_exo(multi_results, pthresh, qthresh, overlap)

			logging.info("Finished multi-gene exotoxin detection for " + prefix + " at " + now().strftime("%Y-%m-%d %H:%M"))


	# detect single-gene exotoxins
	if single == "True":

		with importlib.resources.path("flexo_blast.seq_single_db", "singleblast.faa") as single_path:

			get_single = Blast(
				task = "tblastn",
				dbseqs = infile,
				fasta = single_path,
				final_results_directory = final_results_directory,
				prefix = prefix,
				suffix = "single",
				pthresh = pthresh,
				qthresh = qthresh,
				overlap = overlap,
				evalue = evalue)

			logging.info("Using tblastn to identify single-gene exotoxins in " + prefix + " at " + now().strftime("%Y-%m-%d %H:%M"))

			single_results = get_single.run_blast("tblastn", infile, single_path, final_results_directory, prefix, "single", evalue)
			single_final = get_single.parse_exo(single_results, pthresh, qthresh, overlap)

			logging.info("Finished single-gene exotoxin detection for " + prefix + " at " + now().strftime("%Y-%m-%d %H:%M"))


	# print results to a final results file
	if multi == "True":

		with importlib.resources.path("flexo_blast.seq_multi_db", "multiblast.faa") as multi_path:
			collect_toxins = AnnoFasta(
				fasta = multi_path)
			multi_dict = collect_toxins.get_annotations(multi_path)

		get_final_results = FinalResults(
			final_results_directory = final_results_directory,
			infile = infile,
			prefix = prefix,
			suffix = "multi",
			results = multi_final,
			genedict = multi_dict)

		get_final_results.print_final_results(final_results_directory, infile, prefix, "multi", multi_final, multi_dict)


	if single == "True":

		with importlib.resources.path("flexo_blast.seq_single_db", "singleblast.faa") as single_path:
			collect_toxins = AnnoFasta(
				fasta = single_path)
			single_dict = collect_toxins.get_annotations(single_path)

		get_final_results = FinalResults(
			final_results_directory = final_results_directory,
			infile = infile,
			prefix = prefix,
			suffix = "single",
			results = single_final,
			genedict = single_dict)

		get_final_results.print_final_results(final_results_directory, infile, prefix, "single", single_final, single_dict)


	for blastdb_ext in ("nsq", "nin", "nhr", "ndb", "not", "ntf", "nto", "njs"):
		blastdb_file = "{}.{}".format(infile, blastdb_ext)
		if os.path.isfile(blastdb_file):
			os.remove(blastdb_file)


	logging.info("")
	logging.info("")
	logging.info("")
	logging.info("FLExo_BLAST finished at " + now().strftime("%Y-%m-%d %H:%M"))
	logging.info("Report bugs/concerns to Laura M. Carroll, lmc297@cornell.edu\n")
	logging.info("Have a nice day!")



def main():

	# FLExo_BLAST arguments

	parser = argparse.ArgumentParser(prog = "flexo_blast", usage = "flexo_blast -i </path/to/genome.fasta> -o </path/to/output/directory/> [other options]")

	parser.add_argument("-i", "--input", help = "Path to input genome in fasta format", nargs = 1, required = True)

	parser.add_argument("-o", "--output", help = "Path to desired output directory", nargs = 1, required = True)

	parser.add_argument("--multi", help = "Optional argument; True or False; perform multi-gene exotoxin detection; default = True", nargs = "?", default = "True")

	parser.add_argument("--single", help = "Optional argument; True or False; perform single-gene exotoxin detection; default = True", nargs = "?", default = "True")

	parser.add_argument("--pthresh", help = "Optional argument; integer from 0 to 100; minimum percent amino acid identity threshold for a gene to be considered present; default = 0", nargs = "?", default = 0)

	parser.add_argument("--qthresh", help = "Optional argument; integer from 0 to 100; minimum percent coverage threshold for a gene to be considered present; default = 80", nargs = "?", default = 80)

	parser.add_argument("--overlap", help = "Optional argument; float from 0 to 1; specify maximum proportion of overlap for overlapping genes to be considered separate genes; genes below this threshold will be considered separate, while those above it will be considered overlapping, and only the top hit will be reported; default=0.7", nargs = "?", default = 0.7)

	parser.add_argument("--evalue", help = "Optional argument; float >= 0; maximum blast e-value for a hit to be saved; default = 1", nargs = "?", default = 1)

	parser.add_argument("--version", action="version", version='%(prog)s {}'.format(__version__), help="Print version")

	args = parser.parse_args()

	run_pipeline(args)

if __name__ == "__main__":

	# run BTyper3

	main()
