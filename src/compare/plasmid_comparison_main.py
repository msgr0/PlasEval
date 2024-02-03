__author__ = 'amane'

import os
import argparse
import logging

import pandas as pd
from bidict import bidict

import compare_sets

from log_errors_utils import (
	check_file,
	create_directory
)

'''
USAGE: 
python plasmid_comparison_v2.py --l LEFT_BINS_TSV --r RIGHT_BINS_TSV \
				--out_file OUT_FILE --log_file LOG_FILE
'''
	
def get_plasmid_details(contigs_dict, filename, side):
	'''
	Input: 
		contigs_dict: Key: contig (str), Value: Nested dictionary: 	length (int), 
																	L_copies/R_copies (list of contig copies in plasmid set)
		path to input file
		side ('L' or 'R')
	Returns:
		plasmids: list of list of contig ids
		updated contigs_dict
		plasmids_keys: bidict of plasmid indices <-> plasmid names/ids
	'''
	plasmids = []
	plasmids_keys = bidict()
	count = 0
	pls_ctg_df = pd.read_csv(filename, sep='\t')

	for _, row in pls_ctg_df.iterrows():
		plasmid, contig, length = row['plasmid'], str(row['contig']), row['contig_len']
		if contig not in contigs_dict:
			contigs_dict[contig] = {'length': length, 'L_copies': [], 'R_copies': []}
		if plasmid not in plasmids_keys:
			plasmids_keys[plasmid] = count
			plasmids.append([])
			count += 1
		pls_index = plasmids_keys[plasmid]
		plasmids[pls_index].append(contig)
		contigs_dict[contig][side+'_copies'].append([contig, pls_index, len(plasmids[pls_index])])	
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
	contigs_dict = {}
	pls_ids_dict = {'L': {}, 'R': {}}

	left_plasmids, contigs_dict, pls_ids_dict['L'] = get_plasmid_details(contigs_dict, left_plasmids_file, 'L')
	right_plasmids, contigs_dict, pls_ids_dict['R'] = get_plasmid_details(contigs_dict, right_plasmids_file, 'R')
	compare_sets.run_compare_plasmids(contigs_dict, pls_ids_dict, results_file)