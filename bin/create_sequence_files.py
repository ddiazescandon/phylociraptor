#!/usr/bin/env python
# written by Philipp Resl

import os
import sys
import pandas as pd
from Bio import SeqIO
import argparse


if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="create_sequence_files.py", description = """This script will create fasta files for all the buscos from all species with >% of buscos present""", epilog = """written by Philipp Resl""")
pars.add_argument('--busco_table', dest="busco_table", required=True, help="Path to BUSCO table.")
pars.add_argument('--busco_results', dest="busco_results", required=True, help="Results directory containing all BUSCO runs.")
pars.add_argument('--cutoff', dest="cutoff", required=True, help="Percent cutoff % for BUSCO presence. Species below this threshold will be excluded.")
pars.add_argument('--outdir', dest="outdir", required=True, help="Path to output directory.")
pars.add_argument('--minsp', dest="minsp", required=True, help="Minimum number of species which have to be present to keep the sequences.")
pars.add_argument('--type' , dest="type", required=True, help="Type of sequences (aa or nu).")
pars.add_argument('--genome_statistics' , dest="genome_stats", required=True, help="Path to genome statistics output file.")
pars.add_argument('--gene_statistics' , dest="gene_stats", required=True, help="Path to gene statistics output file.")

args=pars.parse_args()

extension=""
if args.type == "nu":
	extension = ".fna"
else:
	extension = ".faa" 


busco_overview = pd.read_csv(args.busco_table, sep="\t")
genomes = os.listdir(args.busco_results)
#print(busco_overview)
print("Settings:")
print("cutoff: ", args.cutoff)
print("minsp: ", args.minsp)
print("type: ", args.type)
print("outdir: ", args.outdir)

species_list = busco_overview.species.tolist()
print("Original number of species:", len(species_list))
#print(species_list)
#first remove species with too low busco coverage
busco_overview = busco_overview.set_index("species")

genome_file = open(args.genome_stats, "w")
for sp in species_list:
	if busco_overview.loc[sp, "percent_complete"] < float(args.cutoff):
		out = sp + "\tFAILED" + "\t" + str(busco_overview.loc[sp, "percent_complete"]) + "\t" + str(float(args.cutoff))
		print(out, file=genome_file)
		busco_overview = busco_overview.drop([sp])
	else:
		out = sp + "\tOK" + "\t" + str(busco_overview.loc[sp, "percent_complete"]) + "\t" + str(float(args.cutoff)) 
		print(out, file=genome_file)
species_list =  list(busco_overview.index)
print("Species remaining after applying cutoff:", len(species_list))
genome_file.close()

#now loop through each busco and extract sequence for each species
buscos = list(busco_overview.columns.values)
buscos.remove("percent_complete")

gene_file = open(args.gene_stats, "w")
for busco in buscos:
	#print("Processing: " + busco)
	numseqs = 0
	outstring = ""
	for species in species_list:
		for genome in genomes: # this loop is to get the correct directory name, it is very unelegant
			#print(args.busco_results+"/"+genome+"/run_busco/single_copy_busco_sequences/"+busco+extension)
			if species == genome:
				try:
					seqfile = open(args.busco_results + "/" + genome + "/run_busco/single_copy_busco_sequences/" + busco + extension, "r")
					for seq_record in SeqIO.parse(seqfile, "fasta"):
						name = ">" +species+"\n"
						#outfile.write(name)
						#outfile.write(str(seq_record.seq)+"\n")
						outstring = outstring + name
						outstring = outstring + str(seq_record.seq) + "\n"
					seqfile.close()
				except: # skip missing buscos
					continue
	if outstring.count(">") >= int(args.minsp):	# only keep sequences if total number is larger than specified cutoff above.		
		print(busco + "\t" + "OK" + "\t" + str(outstring.count(">")) +"\t" + str(int(args.minsp)), file=gene_file)
		outfile = open(args.outdir+"/"+busco+"_all.fas", "w")
		outfile.write(outstring)
		outfile.close()
	else:
		print(busco + "\t" + "FAILED" + "\t" + str(outstring.count(">")) +"\t" + str(int(args.minsp)), file=gene_file)
gene_file.close()
