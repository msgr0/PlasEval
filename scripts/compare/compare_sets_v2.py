import networkx as nx
import itertools
from collections import defaultdict
import copy
from math import factorial

import logging
import psutil
import time

'''
TODO: 
add many more comments on method: long docstring at the beginning of the file;
revisit variable names and function names; 
maintain consistency in naming variables and functions;
document input and output of every function;
code factorization with appropriate function names.
'''

logger = logging.getLogger(__name__)

def generate_matchings(m, n):
	'''
	Takes the number of copies of a family on both sides and 
	returns a list of matchings where each matching is
	a pair of indices (int) of the contig copies of both sides.
	'''
	matchings = []
	if m >= n:
		pmutns = list(itertools.permutations(list(range(m)), n))
		n_list = list(range(n))
		for pmutn in pmutns:
			matchings.append((list(pmutn), n_list))
	else:
		pmutns = list(itertools.permutations(list(range(n)), m))
		m_list = list(range(m))
		for pmutn in pmutns:
			matchings.append((m_list, list(pmutn)))
	return matchings

def get_matching_positions(ctg_copies, matching):
	'''
	Takes a dictionary of contig copy positions and a matching as input ad
	returns a list of contig positions for each side according to respective indices in the matching 
	'''	
	L, R = ctg_copies['L_copies'], ctg_copies['R_copies']
	lpos, rpos = [], []
	for x in matching[0]:
		lpos.append(L[x])
	for x in matching[1]:
		rpos.append(R[x])	
	return lpos, rpos

#def rename_by_matching(matchings, left_list, right_list, pls_ids):
def rename_by_matching(matching):
	'''
	Takes a matching and renames contigs on both sides according to the indices in the matching
	and returns the list of renamed contigs where 
	each contig is a triple (contig_id, plasmid_id, contig_posn_in_plasmid)
	'''
	reached_contigs = matching.keys()
	left_renamed, right_renamed = [], []
	for contig in reached_contigs:
		M = matching[contig]
		left_ctgs = M[0]
		right_ctgs = M[1]
		for i in range(len(M[0])):
			lpls, lidx = left_ctgs[i][1], left_ctgs[i][2]
			left_renamed.append([contig+'_'+str(i), lpls, lidx])
			rpls, ridx = right_ctgs[i][1], right_ctgs[i][2]
			right_renamed.append([contig+'_'+str(i), rpls, ridx])
	return left_renamed, right_renamed

#Remove contigs unique to a single tool but include them in the score
def remove_unique_contigs(pls_ctg_list, common_contigs, contigs_dict):
	'''
	Takes a dictionary of plasmids and a set of matched contigs and
	returns a new plasmid dictionary consisting of only the matched contigs 
	as well as the cumulative length of removed contigs 
	'''
	unique_ctg_len = 0
	new_list, pls_idx = [], 0
	for pls in pls_ctg_list:
		new_list.append([])
		for c in pls:
			if c in common_contigs:
				new_list[pls_idx].append(c)
			else:
				ctg = c.split('_')[0]
				unique_ctg_len += contigs_dict[ctg]['length']	
		pls_idx += 1		
	return new_list, unique_ctg_len	

#Adds nodes to the bipartite graph G
def add_nodes(G, left_list, right_list, pls_ids):
	G.add_nodes_from(list(set([pls_ids['L'].inv[x[1]] for x in left_list])), bipartite=0)
	G.add_nodes_from(list(set([pls_ids['R'].inv[x[1]] for x in right_list])), bipartite=1)
	return G

def add_edges(G, left_list, right_list, pls_ids):
	'''
	Takes the bipartite graph G, a list of renamed contigs each for both sides and the list of plasmid ids for boths sides,
	adds edges between plasmids if they share a common contig and returns the graph with added edges.
	'''
	def list_ctg_ids_by_pls(pls_ids_dict, ctg_ids_list):
		'''
		For a given dictionary of plasmids and list of contig ids 
		returns the list of contig ids for each plasmid
		'''
		ctg_ids_by_pls = defaultdict(list)
		for x in ctg_ids_list:
			ctg_id, pls_id = x[0], pls_ids_dict.inv[x[1]]
			ctg_ids_by_pls[pls_id].append(ctg_id)
		return ctg_ids_by_pls
	
	left_ctg_ids_by_pls = list_ctg_ids_by_pls(pls_ids['L'], left_list)
	right_ctg_ids_by_pls = list_ctg_ids_by_pls(pls_ids['R'], right_list)
	
	edges_list = []
	for left_plasmid in left_ctg_ids_by_pls:
		for right_plasmid in right_ctg_ids_by_pls:
			left_keys = set(left_ctg_ids_by_pls[left_plasmid])
			right_keys = set(right_ctg_ids_by_pls[right_plasmid])
			if left_keys.intersection(right_keys) != set():
				edges_list.append((left_plasmid, right_plasmid))
	G.add_edges_from(edges_list)
	return G

def modify_partitions(partitions, common):
	'''
	Takes a list of partitions and 
	splits each partition according to intersection with an opposite partition.
	Returns modified list of partitions
	'''
	modified_partitions = []
	for S in partitions:
		if len(S.intersection(common)) != 0 or S.intersection(common) != S:
			modified_partitions.append(S.intersection(common))
			modified_partitions.append(S.difference(common))
		elif S.intersection(common) == S:
			modified_partitions.append(S)				
	return modified_partitions		

def get_partition_cost(partitions, contigs_dict):
	'''
	Takes a list of partitions and dictionary of contig details and
	returns the total length of contigs in the partitions and the cost of partitioning
	'''	
	partition_len = 0
	cost = 0
	for S in partitions:
		S_cost = 0
		for contig in S:
			contig_len = contigs_dict[contig.split('_')[0]]['length']
			partition_len += contig_len
			S_cost += contig_len
			cost = max(cost, S_cost)
	cost = partition_len - cost
	return partition_len, cost

def one_side_cost(side, side_list, opp_list, B, flag, pls_ids, contigs_dict):
	'''
	Computes total cost on one side (splits or joins)
	'''
	def get_ctg_list_by_pls(pls_list, pls_ids, side):
		'''
		Returns list of contigs for each plasmid
		'''
		ctgs_by_pls = {}
		for x in pls_list:
			ctg, pls = x[0], pls_ids[side].inv[x[1]]
			if pls not in ctgs_by_pls:
				ctgs_by_pls[pls] = []
			ctgs_by_pls[pls].append(ctg)
		return ctgs_by_pls
	[s, o] = ['L', 'R'] if flag == 0 else ['R', 'L']
	side_ctgs_by_pls = get_ctg_list_by_pls(side_list, pls_ids, s)
	opp_ctgs_by_pls = get_ctg_list_by_pls(opp_list, pls_ids, o)
	side_len, side_cost = 0, 0
	for node in side:
		partitions = [set(side_ctgs_by_pls[node])]
		for edge in B.edges:
			if edge[flag] == node:
				side_contigs, opp_contigs = set(side_ctgs_by_pls[edge[flag]]), set(opp_ctgs_by_pls[edge[1-flag]])
				common = side_contigs.intersection(opp_contigs)
				partitions = modify_partitions(partitions, common)
		node_len, cost = get_partition_cost(partitions,contigs_dict)
		side_len += node_len
		side_cost += cost
	return side_len, side_cost

#Compute total cost of splits and joins for a particular matching
def compute_match_cost(left_list, right_list, pls_ids, contigs_dict):
	#Create graph with vertices named according to matching and obtain connected components
	B = nx.Graph()
	B = add_nodes(B, left_list, right_list, pls_ids)
	B = add_edges(B, left_list, right_list, pls_ids)
	A = [B.subgraph(c) for c in nx.connected_components(B)]
	left_pls = set([pls_ids['L'].inv[x[1]] for x in left_list])
	n_conn_comp = len(list(A))
	total_cost = 0		
	for i in range(n_conn_comp):
		C = list(A)[i]
		left, right = nx.bipartite.sets(C)	#Split the component according to bipartite sets
		if len(list(right)) != 0 and list(right)[0] in left_pls: #Ensuring proper assignments of bipartite parts
			right,left = left,right
		left_len, left_cost = one_side_cost(left, left_list, right_list, B, 0, pls_ids, contigs_dict)
		right_len, right_cost = one_side_cost(right, right_list, left_list, B, 1, pls_ids, contigs_dict)
	return left_cost, right_cost	

#Compute cost of matching for reached contigs
def compute_current_cost(matching, left_plasmids, right_plasmids, pls_ids, contigs_dict):
	left_renamed, right_renamed = rename_by_matching(matching)
	left_ctg_ids, right_ctg_ids = set(), set()
	for x in left_renamed:
		left_ctg_ids.add(x[0])
	for x in right_renamed:
		right_ctg_ids.add(x[0])
	return compute_match_cost(left_renamed, right_renamed, pls_ids, contigs_dict)

def run_compare_plasmids(left_plasmids, right_plasmids, contigs_dict, pls_ids_dict, results_file):	

	#Common contigs in the solutions obtained from both tools
	left_ctg_ids = set([ctg for ctg in contigs_dict.keys() if len(contigs_dict[ctg]['L_copies']) >= 1])
	right_ctg_ids = set([ctg for ctg in contigs_dict.keys() if len(contigs_dict[ctg]['R_copies']) >= 1])
	common_contigs = left_ctg_ids.intersection(right_ctg_ids)

	max_cost = 0		#Computing upperbound on final_cost
	n_matchings = {}
	max_n_matchings = 1
	for contig in common_contigs:
		m = len(contigs_dict[contig]['L_copies'])
		n = len(contigs_dict[contig]['R_copies'])
		max_cost += m * contigs_dict[contig]['length']
		max_cost += n * contigs_dict[contig]['length']	
		n_matchings[contig] = int(factorial(n)/factorial(n-m)) if n > m else int(factorial(m)/factorial(m-n))
		max_n_matchings *= n_matchings[contig]
	logger.info(f'Maximum possible matchings: {max_n_matchings}')

	start_time = time.time()
	if max_n_matchings <= 10000000:
		### Branch-N-Bound ###
		current_state = {'level': 0, 'total_cost': 0, 'matching': {}, 'cuts_cost': 0, 'joins_cost': 0}
		final_state = {'total_cost': max_cost, 'matching': {}, 'cuts_cost': 0, 'joins_cost': 0}

		contig_list = list(common_contigs)
		sorted_contig_list = sorted(contig_list, key=lambda ctg: n_matchings[ctg])

		count = [0]
		memory_usage = psutil.virtual_memory()

		def recursive_compare(current_state, sorted_contig_list, left_plasmids, right_plasmids, pls_ids_dict, contigs_dict, count):
			nonlocal final_state
			count[0] += 1
			if current_state['level'] < len(sorted_contig_list):				#Compute cost upto this level
				current_contig = sorted_contig_list[current_state['level']]		#Retrieve contig for current level				
				m = len(contigs_dict[current_contig]['L_copies'])
				n = len(contigs_dict[current_contig]['R_copies'])			
				matchings = generate_matchings(m,n)

				for matching in matchings:
					matched_posns = get_matching_positions(contigs_dict[current_contig], matching)
					current_state['matching'][current_contig] = matched_posns
					current_state['cuts_cost'], current_state['joins_cost'] \
						= compute_current_cost(current_state['matching'], left_plasmids, right_plasmids, pls_ids_dict, contigs_dict)	
					current_state['total_cost'] = current_state['cuts_cost'] + current_state['joins_cost']

					if current_state['total_cost'] < final_state['total_cost']:	
						current_state['level'] += 1 
						recursive_compare(current_state, sorted_contig_list, left_plasmids, right_plasmids, pls_ids_dict, contigs_dict, count)
						current_state['level'] -= 1
					del current_state['matching'][current_contig]

			else:
				final_state['total_cost'] = current_state['total_cost']
				final_state['cuts_cost'], final_state['joins_cost'] = current_state['cuts_cost'], current_state['joins_cost']
				final_state['matching'] = copy.deepcopy(current_state['matching'])

		recursive_compare(current_state, sorted_contig_list, left_plasmids, right_plasmids, pls_ids_dict, contigs_dict, count)
		
		end_time = time.time()
		logger.info(f'Total function calls: {count[0]}')
		logger.info(f'Time taken: {end_time - start_time}')
		logger.info(f'Memory usage: {memory_usage.used}')	
		
		total_len, unique_left, unique_right = 0, 0, 0
		for c in contigs_dict:
			l_copies, r_copies = len(contigs_dict[c]['L_copies']), len(contigs_dict[c]['R_copies'])
			ctg_len = contigs_dict[c]['length']
			unique_left += max(l_copies - r_copies, 0) * ctg_len
			unique_right += max(r_copies - l_copies, 0) * ctg_len
			total_len += max(l_copies, r_copies) * ctg_len

		dissimilarity = (unique_left + unique_right + final_state['total_cost']/2)/total_len
		print("Total_ctg_length\t", total_len)
		print("Cuts_cost\t", final_state['cuts_cost'])
		print("Joins_cost\t", final_state['joins_cost'])
		print("Unique_left_ctgs\t", unique_left)
		print("Unique_right_ctgs\t", unique_right)
		print("Dissimilarity_score\t", dissimilarity)
		print("Matching\n", final_state['matching'])	

		logger.info(f'{final_state["matching"]}')

		results_file.write("Total_ctg_length\t" + str(total_len) + "\n")
		results_file.write("Cuts_cost\t" + str(final_state['cuts_cost']) +"\n")
		results_file.write("Joins_cost\t" + str(final_state['joins_cost']) +"\n")
		results_file.write("Unique_left_ctgs\t" + str(unique_left) + "\n")
		results_file.write("Unique_right_ctgs\t" + str(unique_right) + "\n")
		results_file.write("Dissimilarity_score\t" + str(dissimilarity) + "\n")
