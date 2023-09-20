import networkx as nx
import itertools
from collections import defaultdict
import copy
from math import factorial
from collections import defaultdict
#import random
#from bidict import bidict

'''
TODO: 
add many more comments on method: long doctsring at the beginning of the file;
revisit variable names and function names; 
maintain consistency in naming variables and functions;
document input and output of every function;
code factorization with appropriate function names.
'''

#Generate list of permutations by contig family
def generate_permutations(m, n):
	permutations = []
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
	left_ctgs_by_pls = {}
	right_ctgs_by_pls = {}
	
	for x in left_list:
		ctg, pls = x[0], pls_ids['L'].inv[x[1]]
		if pls not in left_ctgs_by_pls:
			left_ctgs_by_pls[pls] = []
		left_ctgs_by_pls[pls].append(ctg)

	for x in right_list:
		ctg, pls = x[0], pls_ids['R'].inv[x[1]]
		if pls not in right_ctgs_by_pls:
			right_ctgs_by_pls[pls] = []
		right_ctgs_by_pls[pls].append(ctg)		

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
	side_ctgs_by_pls = {}
	opp_ctgs_by_pls = {}
	
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
	right_pls = set([pls_ids['R'].inv[x[1]] for x in right_list])

	n_conn_comp = len(list(A))
	total_cost = 0		
	for i in range(n_conn_comp):
		C = list(A)[i]
		left, right = nx.bipartite.sets(C)	#Split the component according to bipartite sets
		#Assign the split sets to the correct tool (as the bipartite split is random)
		temp = []
		if len(list(right)) != 0 and list(right)[0] in left_pls:
			temp = right
			right = left
			left = temp
		if len(list(left)) != 0 and list(left)[0] in right_pls:
			temp = right
			right = left
			left = temp		
		if len(list(left)) == 0:	#If right side plasmid forms a singleton vertex	
			plasmid = list(right)[0]
		elif len(list(right)) == 0:	#If left side plasmid forms a singleton vertex
			plasmid = list(left)[0]
		else:						#If component has vertices from both sides
			left_len, left_cost = one_side_cost(left, left_list, right_list, B, 0, pls_ids, contigs_dict)
			right_len, right_cost = one_side_cost(right, right_list, left_list, B, 1, pls_ids, contigs_dict)
			total_cost += left_cost + right_cost
	return total_cost

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

	print(sorted_contig_list)

		def recursive_compare(current_state, sorted_contig_list, left_plasmids, right_plasmids, pls_ids_dict, contigs_dict, count):
			nonlocal final_state
			count[0] += 1
			if current_state['level'] < len(sorted_contig_list):				#Compute cost upto this level
				current_contig = sorted_contig_list[current_state['level']]		#Retrieve contig for current level				
				m = len(contigs_dict[current_contig]['L_copies'])
				n = len(contigs_dict[current_contig]['R_copies'])			
				matchings = generate_matchings(m,n)

			for perm in permutations:
				matching = get_matching_from_perm(contigs_dict[current_contig], perm)
				#print(matching)
				current_state['matching'][current_contig] = matching
				current_state['cost'] = compute_current_cost(current_state['matching'], left_plasmids, right_plasmids, pls_ids_dict, contigs_dict)	
				print(current_state['level'], sorted_contig_list[current_state['level']], perm, current_state['cost'],final_state['cost'])


				if current_state['cost'] < final_state['cost']:	
					current_state['level'] += 1 
					recursive_compare(current_state, sorted_contig_list, left_plasmids, right_plasmids, pls_ids_dict, contigs_dict)
					current_state['level'] -= 1
				
				del current_state['matching'][current_contig]

		else:
			final_state['cost'] = current_state['cost']
			final_state['matching'] = copy.deepcopy(current_state['matching'])
			print(final_state['cost'])

	recursive_compare(current_state, sorted_contig_list, left_plasmids, right_plasmids, pls_ids_dict, contigs_dict)

	print(final_state['cost'])
	print(final_state['matching'])
