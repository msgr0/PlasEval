__author__ = 'amane'

from sys import argv
import os
import networkx as nx
import itertools
import copy
import argparse

import compare_plasmids

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
			#plasmid_dict[plasmid][contig]['family'] = family
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

	score_db = compare_plasmids.run_compare_plasmids(left_plasmids, right_plasmids)		

	#####################################################
	#	SORTING SCORES AND SAVING RESULTS TO LOG FILE   #
	#####################################################

	#Sorting the combinations in ascending order of the scores
	def sort_by_score(tup):
		tup.sort(key = lambda x: x[1])
		return tup

	score_list = []
	for count in score_db:
		score = score_db[count]['score']
		score_list.append((count, score))

	score_list = sort_by_score(score_list)	

	#Saving results to log file with labeling corresponding to best score noted first
	for tup in score_list:
		count = tup[0]

		results_file.write("Combination "+str(count)+"\n")

		results_file.write("Score:\t\t"+str(score_db[count]['score'])+"\n\n")

		results_file.write("Contribution of lengths to the difference score Left|Right|Common:\t\t")
		for x in score_db[count]['decomposition']:
			results_file.write(str(float(x))+"|")
		results_file.write("\n\n")

		results_file.write("Left plasmids:\n")
		for x in score_db[count]['left_dict']:
			results_file.write(str(x)+":\t")
			for y in list(score_db[count]['left_dict'][x].keys()):
				results_file.write(str(y)+":"+str(score_db[count]['left_dict'][x][y]['length'])+",")
			results_file.write("\n")
		results_file.write("\n")	

		results_file.write("Right plasmids:\n")
		for x in score_db[count]['right_dict']:
			results_file.write(str(x)+":\t")
			for y in list(score_db[count]['right_dict'][x].keys()):
				results_file.write(str(y)+":"+str(score_db[count]['right_dict'][x][y]['length'])+",")
			results_file.write("\n")		
		results_file.write("\n")	

		results_file.write("Left splits:\n")
		for x in score_db[count]['left_splits']:
			for y in x:
				if y != set():
					results_file.write(str(y))
			results_file.write("\n")	
		results_file.write("\n")
		
		results_file.write("Right splits:\n")
		for x in score_db[count]['right_splits']:
			for y in x:
				if y != set():			
					results_file.write(str(y))
			results_file.write("\n")	

		results_file.write("#-------------------------\n\n\n\n")	


