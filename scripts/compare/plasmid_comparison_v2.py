__author__ = 'amane'

import os
import argparse
import logging

from bidict import bidict

import compare_sets_v2

from log_errors_utils import (
	check_file,
	create_directory
)

'''
USAGE: 
python plasmid_comparison_main.py --l LEFT_BINS_TSV --r RIGHT_BINS_TSV \
				--out_file OUT_FILE --log_file LOG_FILE
'''

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
	return plasmids, contigs_dict, plasmids_keys

if __name__ == "__main__":
	description = 'PlasEval: Comparing plasmid bins'
	#Parsing arguments
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument("--l", help="Path to file with 1st set of plasmids")
	parser.add_argument("--r", help="Path to file with 2nd set of plasmids")
	parser.add_argument("--out_file", help="Path to output file")
	parser.add_argument("--log_file", help="Path to log file")
	args = parser.parse_args()

	left_plasmids_file = args.l
	right_plasmids_file = args.r
	for in_file in [left_plasmids_file, right_plasmids_file]:
		check_file(in_file)
	output_file = args.out_file
	log_file = args.log_file
	output_dir = os.path.dirname(output_file)
	log_dir = os.path.dirname(log_file)
	create_directory([output_dir, log_dir])
	results_file = open(output_file, "w")

	# Initialize logging
	logging.basicConfig(
		filename=log_file,
		filemode='w',
		level=logging.INFO,
		format='%(name)s - %(levelname)s - %(message)s'
	)

	#Reading data and saving it to a dictionary with plasmids as keys and a nested dictionary of contigs as values
	left_plasmids = []
	right_plasmids = []
	contigs_dict = {}
	pls_ids_dict = {'L': {}, 'R': {}}

	left_plasmids, contigs_dict, pls_ids_dict['L'] = get_plasmid_details(left_plasmids, contigs_dict, left_plasmids_file, 'L')
	right_plasmids, contigs_dict, pls_ids_dict['R'] = get_plasmid_details(right_plasmids, contigs_dict, right_plasmids_file, 'R')

	compare_sets_v2.run_compare_plasmids(left_plasmids, right_plasmids, contigs_dict, pls_ids_dict, results_file)