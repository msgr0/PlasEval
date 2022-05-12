from __future__ import division
import os
import argparse
import pandas as pd
import numpy as np

#USAGE: python evaluate_sample.py --pred prediction_file --map mapping_file --out out_dir --res out_file --amb 0/1 --ori 0/1


def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

def parse_line(line):
	line = line.split(",")
	c_id = line[0]
	c_len = line[1]
	t_id = line[2]
	t_len = line[3]
	label = line[6]
	return c_id, c_len, t_id, t_len, label

def create_entry(seq_dict, seq_id, seq_len, map_id, label):
	seq_dict[seq_id] = {}
	seq_dict[seq_id]['label'] = [label]
	seq_dict[seq_id]['length'] = int(seq_len)
	if label == 'ambiguous':
		seq_dict[seq_id]['map'] = [None]
	else:
		seq_dict[seq_id]['map'] = [map_id]	
	return seq_dict				

def update_mapping(seq_dict, seq_id, map_id, label):
	if label == 'ambiguous':
		seq_dict[seq_id]['map'].append(None)
	else:	
		seq_dict[seq_id]['map'].append(map_id)
	return seq_dict

def update_labels(seq_dict, seq_id, label):
	seq_dict[seq_id]['label'].append(label)	
	return seq_dict

if __name__ == "__main__":
	#Parsing arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--pred", help="Path to predictions")
	parser.add_argument("--map", help="Path to contig mapping file")
	parser.add_argument("--out", help="Path to output dir")
	parser.add_argument("--res", help="Path to output file")
	parser.add_argument("--ori", nargs='?', default = 0, help="0/1 to denote if predicted contigs are oriented")
	parser.add_argument("--amb", nargs='?', default = 1, help="0/1 to denote if ambiguous contigs should be considered")
		
	args = parser.parse_args()

	prediction_file = args.pred
	mapping_file = args.map
	output_dir = args.out
	results = args.res
	ori_used = int(args.ori)
	ambiguous = int(args.amb)

	if not os.path.exists(output_dir):
		os.makedirs(output_dir)	
	results_file = open(os.path.join(output_dir, results), "w")	

	###################################
	# Storing input from mapping file #
	###################################
	contig_dict = {}	#Key: contig ID, Values: List of true sequences mapped to, List of true sequences labels, contig length
	truth_dict = {}		#Key: sequence ID, Values: List of contigs mapped to, sequence label, sequence length

	chr_list, pl_list = [], []

	string_list = read_file(mapping_file)
	for line in string_list[1:-1]:
		c_id, c_len, t_id, t_len, label = parse_line(line)
		label = label[:-1]
		if c_id not in contig_dict:
			contig_dict = create_entry(contig_dict, c_id, c_len, t_id, label)
		else:
			contig_dict = update_mapping(contig_dict, c_id, t_id, label)
			contig_dict = update_labels(contig_dict, c_id, label)

		if t_id not in truth_dict:	 
			if label != 'ambiguous':
				truth_dict = create_entry(truth_dict, t_id, t_len, c_id, label)
			if label == 'chromosome':
				chr_list.append((t_id, int(t_len)))
			elif label == 'plasmid':
				pl_list.append((t_id, int(t_len)))
		else:
			truth_dict = update_mapping(truth_dict, t_id, c_id, label)

	######################################
	# Storing input from prediction file #
	######################################
	plasmid_dict = {}	#Key: pred. plasmid ID, Values: 
	n_pred = 0
	string_list = read_file(prediction_file)
	for line in string_list:
		line = line.split(';')
		p_id = line[0]
		plasmid_dict[p_id] = []
		n_pred += 1
		chain = line[1].split(',')
		for c_id in chain:
			if c_id:
				if ori_used == 1:
					c_id = c_id[:-1]	
				plasmid_dict[p_id].append(c_id)	

	#########################################
	# Comparing true and predicted plasmids #
	#########################################
	true_plasmids = []
	for t_id in truth_dict:
		if truth_dict[t_id]['label'][0] == 'plasmid':
			true_plasmids.append(t_id)
	if ambiguous == 1:
		true_plasmids.append('ambiguous')

	pred_plasmids = list(plasmid_dict.keys())

	n_cols = len(true_plasmids)
	n_rows = len(pred_plasmids)

	covered_len_df = pd.DataFrame(np.zeros(shape=(n_rows,n_cols)),columns=true_plasmids,index=pred_plasmids)
	pred_lens = {}

	for p_id in pred_plasmids:
		pred_lens[p_id] = 0
		for c_id in plasmid_dict[p_id]:
			t_id = contig_dict[c_id]['map'][0]
			c_len = contig_dict[c_id]['length']
			pred_lens[p_id] += c_len
			if t_id == None:
				if ambiguous == 1:
					covered_len_df.ix[p_id,'ambiguous'] += c_len
			else:
				label = truth_dict[t_id]['label'][0]
				if label == 'plasmid':
					covered_len_df.ix[p_id,t_id] += c_len			

	#Precision df
	prec_df = covered_len_df.copy()
	for p_id in pred_plasmids:
		prec_df.loc[p_id] = prec_df.loc[p_id].div(pred_lens[p_id])	

	#Recall df
	rec_df = covered_len_df.copy()
	if ambiguous == 1:
		plasmid_cols = covered_len_df[true_plasmids[:-1]]
		rec_df = plasmid_cols.copy()	

	for t_id in rec_df.columns:
		rec_df[t_id] = rec_df[t_id].div(truth_dict[t_id]['length'])

	#####################################
	# Writing evaluation to output file #
	#####################################
	results_file.write("##### general information #####\n")
	results_file.write("No. of reference chromosomes: "+ str(len(chr_list)) +"\n")
	for pair in chr_list:
		t_id = pair[0]
		t_len = pair[1]
		results_file.write(t_id + ": " + str(t_len)+ " nt\n")
	results_file.write("No. of reference plasmids: "+ str(len(pl_list)) +"\n")
	for pair in pl_list:
		t_id = pair[0]
		t_len = pair[1]
		results_file.write(t_id + ": " + str(t_len)+ " nt\n")
	results_file.write("\n")
	
	results_file.write("No. of predicted plasmids: "+ str(len(pred_lens.keys())) +"\n")
	for p_id in pred_lens:
		p_len = pred_lens[p_id]
		results_file.write(p_id + ": " + str(p_len)+ " nt\n")
	results_file.write("\n")	

	results_file.write("> predicted plasmid covers <proportion> of reference plasmid\n")
	rec_df.to_csv(results_file, sep='\t')
	results_file.write("\n")

	results_file.write("> reference plasmid covers <proportion> of predicted plasmid\n")
	prec_df.T.to_csv(results_file, sep='\t')
	results_file.write("\n")

	results_file.write("> in total, how much of predicted plasmid is covered by reference plasmids\n")
	pred_proportions = []
	for p_id in pred_plasmids:
		proportion = prec_df.loc[p_id].sum()
		results_file.write(p_id+"\t"+str(proportion)+"\n")
		pred_proportions.append((p_id, proportion))
	results_file.write("\n")

	results_file.write("> in total, how much of reference plasmid is covered by predicted plasmids\n")	
	ref_proportions = []
	for t_id in rec_df.columns:
		proportion = rec_df[t_id].sum()
		results_file.write(t_id+"\t"+str(proportion)+"\n")
		ref_proportions.append((t_id, proportion))
	results_file.write("\n")		

	threshold = 0.8
	results_file.write("> pairs of predicted and reference plasmids with coverage >= %f in both directions\n" % threshold)
	empty = True
	for p_id in pred_plasmids:
		for t_id in rec_df.columns:
			if prec_df.loc[p_id,t_id] >= threshold and rec_df.loc[p_id, t_id] >= threshold:
				results_file.write(p_id+ " <-> "+t_id+"\n")
				empty = False
	if empty:
		out.write("none\n")
	results_file.write("\n")

	results_file.write("> summary scores\n")
	precision = 0
	recall = 0
	f1 = 0

	covered = 0
	total_len = 0
	for pair in pred_proportions:
		p_id = pair[0]
		prec = pair[1]
		length = pred_lens[p_id]
		covered += prec*length
		total_len += length
	if total_len > 0:
		precision = covered / total_len	
	results_file.write("precision: "+str(precision)+"\n")	

	covered = 0
	total_len = 0
	for pair in ref_proportions:
		t_id = pair[0]
		rec = pair[1]
		length = truth_dict[t_id]['length']
		covered += rec*length
		total_len += length
	if total_len > 0:
		recall = covered / total_len	
	results_file.write("recall: "+str(recall)+"\n")		

	if (precision + recall) != 0:
		f1 = 2*precision*recall / (precision + recall)
	results_file.write("f1_score: "+str(f1)+"\n")	




					
	


			


				


