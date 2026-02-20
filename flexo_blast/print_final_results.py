from Bio import SeqIO

class AnnoFasta:
	"""
	Parse FASTA to get genes linked to exotoxins
	
	get_annotations
	purpose: parse FASTA to get annotated genes linked to exotoxins
	input:
	    fasta = path to FASTA file
	output:
	    a dictionary linking exotoxins to genes
	"""
	
	def __init__(self, fasta):
		
		self.fasta = fasta
		
	def get_annotations(self, fasta):
		
		d = {}
		
		with open(fasta, "r") as infile:
			for record in SeqIO.parse(infile, "fasta"):
				seqid = str(record.id).strip()
				toxin = seqid.split("~~~")[1].strip()
				if toxin not in d.keys():
					d[toxin] = []
				d[toxin].append(seqid)
				
		return d
				

class FinalResults:
	"""
	Print final results file

	print_final_results
	purpose: print a final results file containing multi- and single-gene exotoxins detected
	input:
		final_results_directory = FLExo_ final results directory
		infile = path to file used as input
		prefix = prefix used for output files associated with a particular genome
		suffix = suffix added to final results file name 
		results = results to print
		genedict = dictionary linking exotoxins to genes
	output: 
		prints final results file for query genome
		
	"""

	def __init__(self, final_results_directory, infile, prefix, suffix, results, genedict):

		self.final_results_directory = final_results_directory
		self.infile = infile,
		self.prefix = prefix
		self.suffix = suffix
		self.results = results
		self.genedict = genedict

	def print_final_results(self, final_results_directory, infile, prefix, suffix, results, genedict):

		header = ["#filename", "prefix"]
		final_line = [infile, prefix]
		print(results)
		
		for key in sorted(genedict.keys()):
			header.append(key)
			vals = sorted(set(genedict[key]))
			gene_hits = []
			for val in vals:
				print("gene")
				print(val)
				if val in results:
					gene_hits.append(val.split("~~~")[0].strip())
			gene_hits = sorted(set(gene_hits))
			final_cell = str(len(gene_hits)) + "/" + str(len(vals))
			final_cell = final_cell + "(" + ",".join(gene_hits) + ")"
			final_line.append(final_cell)
			
		with open(final_results_directory + prefix + "_" + suffix + "_pa.tsv", "a") as outfile:
			print("\t".join(header), file = outfile)
			print("\t".join(final_line), file = outfile)