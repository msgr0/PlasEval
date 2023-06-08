__author__ = 'amane'

from sys import argv
import os
import networkx as nx
import itertools
import copy
import argparse
from bidict import bidict
import compare_sets_v2

#USAGE: python plasmid_comparison_main.py --l left_plasmids.tsv --r right_plasmids.tsv --out out_dir --res out_file

def read_file(filename):
    string = open(filename, "r").read()
    string_list = string.split("\n")
    string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
    return string_list
	
def get_plasmid_details(plasmids, contigs_dict, filename, prefix):
	string_list = read_file(filename)
	plasmids_keys = bidict()
	count = 0
	for string in string_list:
		string = string.split("\t")
		plasmid, contig, length = string[0], string[1], int(string[2])

		if contig not in contigs_dict:
			contigs_dict[contig] = {'length': length, 'L_copies': [], 'R_copies': []}
		
		if plasmid not in plasmids_keys:
			plasmids_keys[plasmid] = count
			plasmids.append([])
			count += 1

		pls_id = plasmids_keys[plasmid]
		plasmids[pls_id].append(contig)
		contigs_dict[contig][prefix+'_copies'].append([contig, pls_id, len(plasmids[pls_id])])
			
	#print(plasmid_dict)	
	return plasmids, contigs_dict, plasmids_keys

if __name__ == "__main__":
	#Parsing arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--l", help="Path to file with 1st set of plasmids")
	parser.add_argument("--r", help="Path to file with 2nd set of plasmids")
	parser.add_argument("--out", help="Path to output dir")
	parser.add_argument("--res", help="Path to output file")

	args = parser.parse_args()

	left_plasmids_file = args.l
	right_plasmids_file = args.r
	output_dir = args.out
	results = args.res

	if not os.path.exists(output_dir):
		os.makedirs(output_dir)	
	results_file = open(os.path.join(output_dir, results), "w")	

	#Reading data and saving it to a dictionary with plasmids as keys and a nested dictionary of contigs as values
	left_plasmids = []
	right_plasmids = []

	contigs_dict = {}
	pls_ids_dict = {'L': {}, 'R': {}}
	#family_dict = {}

	left_plasmids, contigs_dict, pls_ids_dict['L'] = get_plasmid_details(left_plasmids, contigs_dict, left_plasmids_file, 'L')
	right_plasmids, contigs_dict, pls_ids_dict['R'] = get_plasmid_details(right_plasmids, contigs_dict, right_plasmids_file, 'R')

	compare_sets_v2.run_compare_plasmids(left_plasmids, right_plasmids, contigs_dict, pls_ids_dict)

	#print(contigs_dict)

	#print(common_contigs)