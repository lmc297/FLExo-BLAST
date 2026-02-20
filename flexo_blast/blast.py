import os
import subprocess

import pandas as pd
from pandas.errors import EmptyDataError

class Blast:
	"""
	Detect genes with blast

	run_blast
	purpose:
		blast a query against a reference and save results to a text file
	input:
		task = blast algorithm to use (e.g., tblastn)
		dbseqs = blast database fasta file
		fasta = blast query fasta file
		final_results_directory = FLExo_BLAST final results directory
		prefix = genome prefix for output files
		suffix = suffix for output files corresponding to FLExo_BLAST task (e.g., multi, single)
		pthresh = minimum percent identity threshold for blast hits
		qthresh = minimum percent query coverage for blast hits
		overlap = maximum proportion of overlap for overlapping blast hits to be considered separate (hits overlapping below this threshold will be considered separate hits)
		evalue = maximum blast e-value for blast hits

	output:
		path to blast output file of results

	parse_exo
	purpose:
		determine which exotoxins are present, using blast results from run_blast
	input:
		exofile = path to blast results file containing exotoxin hits
		pthresh = minimum percent identity threshold for blast hits
		qthresh = minimum percent query coverage for blast hits
		overlap = maximum proportion of overlap for overlapping genes to be considered separate genes
			Genes below this threshold will be considered separate, while those above it will be considered overlapping, and only the top hit will be reported
	output:
		a sorted list of exotoxins detected in a sequence

	"""

	def __init__(self, task, dbseqs, fasta, final_results_directory, prefix, suffix, pthresh, qthresh, overlap, evalue):

		self.task = task
		self.dbseqs = dbseqs
		self.fasta = fasta
		self.final_results_directory = final_results_directory
		self.prefix = prefix
		self.suffix = suffix
		self.pthresh = pthresh
		self.qthresh = qthresh
		self.overlap = overlap
		self.evalue = evalue

	def run_blast(self, task, dbseqs, fasta, final_results_directory, prefix, suffix, evalue):
		# create output folder if it doesn't exist
		blast_results_dir = os.path.join(final_results_directory, suffix)
		os.makedirs(blast_results_dir, exist_ok=True)

		# build a BLAST database if it doesn't exist yet
		if not os.path.exists(dbseqs + ".nsq"):
			proc = subprocess.run(["makeblastdb", "-in", dbseqs, "-dbtype", "nucl"])
			proc.check_returncode()

		# run the BLAST task
		blast_results = os.path.join(blast_results_dir, "{}_{}.txt".format(prefix, suffix))
		proc = subprocess.run([
			task,
			"-query", fasta,
			"-db", dbseqs,
			"-out", blast_results,
			"-max_target_seqs", "1000000000",
			"-evalue", evalue,
			"-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp sseq",
		])
		proc.check_returncode()

		# return path to results
		return blast_results

	def parse_exo(self, exofile, pthresh, qthresh, overlap):

		try:

			blast_results_file = pd.read_csv(exofile, sep = "\s+", header = None)
			blast_results_file = blast_results_file.sort_values(by = [0], ascending = True) # genome as query and db as db
			genes = blast_results_file.iloc[:,0] # genome as db and db as query

			rangedict = {}
			bitdict = {}
			genes = genes.unique()
			for gene in genes:
				gene_subset = blast_results_file[blast_results_file[0].str.match(gene)] # genome as query and db as db
				gene_subset = gene_subset.sort_values(by = [11], ascending = False) # sort by bitscore
				max_gene = gene_subset.iloc[0,0]
				pid = float(gene_subset.iloc[0,2])
				qid = float(gene_subset.iloc[0,14]) # qcovs
				#qid = 100*(float(gene_subset.iloc[0,3])/float(gene_subset.iloc[0,12])) # original BTyper query coverage value; genome as db and mlst db as query
				glen = float(gene_subset.iloc[0,3])
				gstart = int(gene_subset.iloc[0,8])
				gend = int(gene_subset.iloc[0,9])
				if gstart < gend:
					grange = list(range(gstart, gend + 1))
				elif gend < gstart:
					grange = list(range(gend, gstart + 1))
				if pid >= pthresh and qid >= qthresh:
					rangedict[max_gene] = grange
					bitdict[max_gene] = float(gene_subset.iloc[0,11])


			exo = []
			#comparisons = []
			for key, val in rangedict.items():
				overlap_dict = {}
				test_bits = bitdict[key]
				grange = val
				for key2, val2 in rangedict.items():
					if key != key2: #and key + "VS" + key2 not in comparisons and key2 + "VS" + key not in comparisons:
						#comparisons.append(key + "VS" + key2)
						#comparisons.append(key2 + "VS" + key)
						test_overlap = float(len(list(set(grange) & set(val2))))/float(len(val2))
						og_bits = bitdict[key2]
						if test_overlap > overlap: #and test_bits  > og_bits:
							if key not in overlap_dict.keys():
								overlap_dict[key] = test_bits
							overlap_dict[key2] = og_bits
				if len(overlap_dict.keys()) > 0:
					maxbits = 0
					max_gene = []
					for okey in overlap_dict.keys():
						oval = overlap_dict[okey]
						if oval > maxbits:
							maxbits = oval
							max_gene = [okey]
						elif oval == maxbits:
							max_gene.append(okey)
				else:
					max_gene = [key]
				max_gene = sorted(max_gene)
				for mg in max_gene:
					if mg not in exo:
						exo.append(mg)

		except EmptyDataError:
			exo = []
		return(sorted(exo))