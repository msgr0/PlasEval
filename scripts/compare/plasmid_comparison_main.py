__author__ = 'amane'

from sys import argv
import os
import networkx as nx
import itertools
import copy
import argparse

import compare_sets

#USAGE: python plasmid_comparison_main.py --l left_plasmids.tsv --r right_plasmids.tsv --out out_dir --res out_file

def read_file(filename):
    string = open(filename, "r").read()
    string_list = string.split("\n")
    string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
    return string_list
	
def get_plasmid_details(plasmid_dict, filename, prefix):
	string_list = read_file(filename)
	for string in string_list:
		string = string.split("\t")
		plasmid, contig, length = prefix+string[0], string[1], int(string[2])
		if plasmid not in plasmid_dict:
			plasmid_dict[plasmid] = {}
		if contig not in plasmid_dict[plasmid]:	
			plasmid_dict[plasmid][contig] = {}
			plasmid_dict[plasmid][contig]['length'] = length
			plasmid_dict[plasmid][contig]['copies'] = 1
		else:
			plasmid_dict[plasmid][contig]['copies'] += 1	
			
	#print(plasmid_dict)	
	return plasmid_dict

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
	left_plasmids = {}
	right_plasmids = {}
	#family_dict = {}

	left_plasmids = get_plasmid_details(left_plasmids, left_plasmids_file, 'L_')
	right_plasmids = get_plasmid_details(right_plasmids, right_plasmids_file, 'R_')

	compare_sets.run_compare_plasmids(left_plasmids, right_plasmids)

	#print(common_contigs)