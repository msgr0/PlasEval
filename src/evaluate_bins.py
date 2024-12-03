#!/usr/bin/python

from __future__ import division
import pandas as pd
import os
from log_errors_utils import (
	check_file,
	create_directory
)

def eval_bins(pred_dict, pls_dict, len_dict, th_len, eval):
	#Following functions are used to compute precision and recall,
	#	for each predicted bin and true plasmid bin respectively
	def create_bin_entry():
		'''
		Input: None
		Returns: Dictionary entry for weighted and unweighted statistics
		'''
		stat_dict = {'wtd': {}, 'unwtd': {}}
		for eval_type in ['wtd', 'unwtd']: 
			stat_dict[eval_type] = {'Val': 0, 'Bin': 'NA', 'Common': 0, 'Total': 0}
		return stat_dict
	
	def get_total_ctgs(ctg_list, len_dict, th_len):
		'''
		Input: 
			List of contigs
			Dictionary of contig lengths
			Length threshold
		Returns:
			Number and total length of contigs in the list
		'''
		n_ctgs = 0
		len_ctgs = 0	
		for ctg in ctg_list:
			ctg_len = len_dict[ctg]
			if ctg_len >= th_len:
				n_ctgs += 1
				len_ctgs += ctg_len
		return n_ctgs, len_ctgs	
	
	def compute_best_bin(stat_dict, ctg_list, opp_bins_dict, len_dict, th_len):
		'''
		Input:
			Dictionary of weighted and unweighted statistics for the bin in question
			List of contigs forming the bin
			Dictionary of bins against which to compute statistics:
				For computing precision: dictionary of true plasmid bins 
				For computing recall: dictionary of predicted bins 
				Format: (Key: Bin id, Value: List of contigs)
			Dictionary of contig lengths
			Length threshold
		Returns:
			Dictionary of weighted and unweighted statistics for the bin in question
			updated with best matched bin details
		'''
		n_ctgs, len_ctgs = stat_dict['unwtd']['Total'], stat_dict['wtd']['Total']
		for bin_id in opp_bins_dict:
			common_ctgs = set(opp_bins_dict[bin_id]).intersection(set(ctg_list))
			ncommon_ctgs, lencommon_ctgs = 0, 0
			for ctg in common_ctgs:
				ctg_len = len_dict[ctg]
				if ctg_len >= th_len:
					ncommon_ctgs += 1
					lencommon_ctgs += ctg_len
			n_stat, len_stat = 0, 0
			if n_ctgs >= 1:
				n_stat = ncommon_ctgs / n_ctgs
				len_stat = lencommon_ctgs / len_ctgs
			if n_stat > stat_dict['unwtd']['Val']:
				stat_dict['unwtd'] = {'Val': n_stat, 'Bin': bin_id, 'Common': ncommon_ctgs, 'Total': n_ctgs}
			if len_stat > stat_dict['wtd']['Val']:
				stat_dict['wtd'] = {'Val': len_stat, 'Bin': bin_id, 'Common': lencommon_ctgs, 'Total': len_ctgs}
		return stat_dict
	
	recall = {}
	precision = {}

	for ref_pls in pls_dict:
		recall[ref_pls] = create_bin_entry()
		nref_ctgs, lenref_ctgs = get_total_ctgs(pls_dict[ref_pls], len_dict, th_len)
		recall[ref_pls]['unwtd']['Total'] = nref_ctgs
		recall[ref_pls]['wtd']['Total'] = lenref_ctgs
		recall[ref_pls] = compute_best_bin(recall[ref_pls], pls_dict[ref_pls], pred_dict, len_dict, th_len)

	for pred_pls in pred_dict:
		precision[pred_pls] = create_bin_entry()
		npred_ctgs, lenpred_ctgs = get_total_ctgs(pred_dict[pred_pls], len_dict, th_len)
		precision[pred_pls]['unwtd']['Total'] = npred_ctgs
		precision[pred_pls]['wtd']['Total'] = lenpred_ctgs
		precision[pred_pls] = compute_best_bin(precision[pred_pls], pred_dict[pred_pls], pls_dict, len_dict, th_len)

	#Following functions are used to compute overall statistics and to write to the output file
	def compute_overall_details(stat_dict, best_match, ovr_dict):
		'''
		Input:
			Dictionary of weighted and unweighted statistics for the bin in question
			Dictionary of best match details
		'''
		best_match['n_stat'] = float("{:.4f}".format(stat_dict['unwtd']['Val']))
		best_match['len_stat'] = float("{:.4f}".format(stat_dict['wtd']['Val']))
		best_match['n_bin'] = stat_dict['unwtd']['Bin']
		best_match['len_bin'] = stat_dict['wtd']['Bin']
		ovr_dict['ovr_n_common'] += stat_dict['unwtd']['Common']
		ovr_dict['ovr_n_total'] += stat_dict['unwtd']['Total']
		ovr_dict['ovr_len_common'] += stat_dict['wtd']['Common']
		ovr_dict['ovr_len_total'] += stat_dict['wtd']['Total']
		return best_match, ovr_dict
	
	def write_best_match_details(eval_file, bin_id, best_match):
		'''
		Input: 
			Output file
			Bin id
			Details of bin matched to bin in question
		Output: None
		'''
		eval_file.write(bin_id + '\t' \
			 	+ str(best_match['n_stat']) + '\t' + best_match['n_bin'] + '\t' \
				+ str(best_match['len_stat']) + '\t' + best_match['len_bin'] +"\n")
		
	def compute_overall_stat(ovr_details):
		'''
		Input: Dictionary with overall details for the statistic (precision or recall) in question
		Output: Weighted and unweighted overall statistic (precision or recall)
		'''
		ovr_n_stat, ovr_len_stat = 0, 0
		if ovr_details['ovr_n_total'] != 0:
			ovr_n_stat = ovr_details['ovr_n_common']/ovr_details['ovr_n_total']
		if ovr_details['ovr_len_total'] != 0:
			ovr_len_stat = ovr_details['ovr_len_common']/ovr_details['ovr_len_total']	
		return ovr_n_stat, ovr_len_stat
	
	eval.write("#Precision: Proportion of correctedly identified contigs for each prediction\n")
	eval.write(">Precision details\n")
	eval.write("#Predicted_bin\tUnwtd_Precision\tUnwtd_Reference_plasmid\tWtd_Precision\tWtd_Reference_plasmid\n")
	ovr_details = {'ovr_n_common': 0, 'ovr_len_common': 0, 'ovr_n_total': 0, 'ovr_len_total': 0}
	for bin_id in precision:
		best_match_details = {'n_stat': None, 'n_bin': None, 'len_stat': None, 'len_bin': None}
		best_match_details, ovr_details = \
			compute_overall_details(precision[bin_id], best_match_details, ovr_details)
		write_best_match_details(eval, bin_id, best_match_details)
	ovr_n_prec, ovr_len_prec = compute_overall_stat(ovr_details)

	eval.write("\n")

	eval.write("#Recall: Proportion of correctedly identified contigs for each reference\n")
	eval.write(">Recall details\n")
	eval.write("#Reference_plasmid\tUnwtd_Recall\tUnwtd_Predicted_bin\tWtd_Recall\tWtd_Predicted_bin\n")
	ovr_details = {'ovr_n_common': 0, 'ovr_len_common': 0, 'ovr_n_total': 0, 'ovr_len_total': 0}
	ovr_n_rec, ovr_len_rec = 0, 0
	for bin_id in recall:
		best_match_details = {'n_stat': None, 'n_bin': None, 'len_stat': None, 'len_bin': None}
		best_match_details, ovr_details = \
			compute_overall_details(recall[bin_id], best_match_details, ovr_details)
		write_best_match_details(eval, bin_id, best_match_details)
	ovr_n_rec, ovr_len_rec = compute_overall_stat(ovr_details)	

	eval.write("\n")

	n_f1, len_f1 = 0, 0
	if (ovr_n_prec + ovr_n_rec) != 0:
		n_f1 = 2*ovr_n_prec*ovr_n_rec / (ovr_n_prec + ovr_n_rec)
	if (ovr_len_prec + ovr_len_rec) != 0:
		len_f1 = 2*ovr_len_prec*ovr_len_rec / (ovr_len_prec + ovr_len_rec)		
	
	eval.write("#Final statistics (Unwtd and Wtd)\n")
	eval.write(">Overall details\n")
	eval.write("#Overall_statistic\tUnwtd_statistic\tWtd_statistic\n")
	eval.write('Precision\t' + str(ovr_n_prec) + '\t' + str(ovr_len_prec) + '\n')
	eval.write('Recall\t' + str(ovr_n_rec) + '\t' + str(ovr_len_rec) + '\n')
	eval.write('F1\t' + str(n_f1)  + '\t' + str(len_f1) + '\n')


def get_bin_details(len_dict, bins_file):
	'''
	Input: 
		path to input file
		len_dict: Key: contig (str), Value: length (int), 
	Returns:
		pls_dict: Key: plasmid id (str), Value: list of contig ids
		updated len_dict
	'''
	pls_ctg_df = pd.read_csv(bins_file, sep='\t')
	pls_dict = {}
	for _, row in pls_ctg_df.iterrows():
		plasmid, contig, length = row['plasmid'], str(row['contig']), row['contig_len']
		len_dict[contig] = length
		if plasmid not in pls_dict and "chromosome" not in plasmid:
			pls_dict[plasmid] = []
		if contig not in set(pls_dict[plasmid]):
			pls_dict[plasmid].append(contig)
	return pls_dict, len_dict

def eval_mode(pred_file, gt_file, min_len, output_file):
	'''
	Reads prediction and ground truth files
	Initializes dictionaries and stores prediction and ground truth bins
	Initializes and populates a dictionary of contig lengths 
	Calls the eval_bins function to compute the precision and recall statistics
	'''
	for in_file in [pred_file, gt_file]:
		check_file(in_file)
	output_dir = os.path.dirname(output_file)
	create_directory([output_dir])
	eval_file = open(output_file, "w")
	#Reading data and saving it to a dictionary with plasmids as keys and a nested dictionary of contigs as values
	len_dict = {}
	pred_dict, len_dict = get_bin_details(len_dict, pred_file)
	gt_dict, len_dict = get_bin_details(len_dict, gt_file)
	eval_bins(pred_dict, gt_dict, len_dict, min_len, eval_file)

			

