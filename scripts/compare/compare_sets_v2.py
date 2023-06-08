import networkx as nx
#import queue
import itertools
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
			permutations.append((list(pmutn), n_list))
	else:
		pmutns = list(itertools.permutations(list(range(n)), m))
		m_list = list(range(m))
		for pmutn in pmutns:
			permutations.append((m_list, list(pmutn)))
	return permutations

#Obtain matching of contigs using given permutation
def get_matching_from_perm(ctg_copies, perm):
	L, R = ctg_copies['L_copies'], ctg_copies['R_copies']
	lpos, rpos = [], []
	for x in perm[0]:
		lpos.append(L[x])
	for x in perm[1]:
		rpos.append(R[x])	
	return lpos, rpos

#Rename contigs in both dictionaries according to the permutation
def rename_by_pmutn(matchings, left_list, right_list, pls_ids):
	reached_contigs = matchings.keys()
	left_renamed, right_renamed = [], []
	for contig in reached_contigs:
		M = matchings[contig]
		left_ctgs = M[0]
		right_ctgs = M[1]

		for i in range(len(M[0])):
			lpls, lidx = left_ctgs[i][1], left_ctgs[i][2]
			left_renamed.append([contig+'_'+str(i), lpls, lidx])

			rpls, ridx = right_ctgs[i][1], right_ctgs[i][2]
			right_renamed.append([contig+'_'+str(i), rpls, ridx])

	return left_renamed, right_renamed

#Remove contigs unique to a single tool but include them in the score
def remove_unique_contigs(plasmids_dict, common_contigs):
	total_len = 0
	new_dict = copy.deepcopy(plasmids_dict)
	for p in plasmids_dict:
		for c in plasmids_dict[p]:
			if c not in common_contigs:
				new_dict[p].pop(c)
				total_len += plasmids_dict[p][c]['length']
		if len(new_dict[p].keys()) == 0:
			new_dict.pop(p)				
	return new_dict, total_len	

#Adds nodes to the bipartite graph G
def add_nodes(G, left_list, right_list, pls_ids):
	G.add_nodes_from(list(set([pls_ids['L'].inv[x[1]] for x in left_list])), bipartite=0)
	G.add_nodes_from(list(set([pls_ids['R'].inv[x[1]] for x in right_list])), bipartite=1)
	return G

#Adds edges to the bipartite graph G
def add_edges(G, left_list, right_list, pls_ids):
	
	def list_ctg_ids_by_pls(pls_ids_dict ctg_ids_list):
		'''
		For a given dictionary of plasmids and list of contig ids, ...
		'''
		ctg_ids_by_pls = defaultdict(list)
		for x in cg_ids_list:
			ctg_id, pls_id = x[0], pls_ids_dict.inv[x[1]]
			ctg_ids_by_pls[pls_id].append(ctg_id)
		return ctg_ids_by_pls
	
	left_ctg_ids_by_pls = list_ctg_ids_by_pls(pls_ids['L'], left_list)
	right_ctg_ids_by_pls = list_ctg_ids_by_pls(pls_ids['R'], right_list)
	
	edges_list = []
	for left_plasmid in left_ctgs_by_pls:
		for right_plasmid in right_ctgs_by_pls:
			left_keys = set(left_ctgs_by_pls[left_plasmid])
			right_keys = set(right_ctgs_by_pls[right_plasmid])
			if left_keys.intersection(right_keys) != set():
				edges_list.append((left_plasmid, right_plasmid))
	G.add_edges_from(edges_list)
	return G

#Modifies set of splits to set of mutually exclusive partitions  
def modify_partitions(partitions, common):
	modified_partitions = []
	for S in partitions:
		if len(S.intersection(common)) != 0 or S.intersection(common) != S:
			modified_partitions.append(S.intersection(common))
			modified_partitions.append(S.difference(common))
		elif S.intersection(common) == S:
			modified_partitions.append(S)				
	return modified_partitions		

#Obtain cost of a partition
def get_partition_cost(node, partitions, plasmids, contigs_dict):
	total_len = 0
	cost = 0
	for S in partitions:
		S_cost = 0
		for contig in S:
			contig_len = contigs_dict[contig.split('_')[0]]['length']
			total_len += contig_len
			S_cost += contig_len
			cost = max(cost, S_cost)
	cost = total_len - cost
	return total_len, cost

#Compute total cost on one side (splits or joins)
def one_side_cost(side, side_list, opp_list, B, flag, pls_ids, contigs_dict):
	side_ctgs_by_pls = {}
	opp_ctgs_by_pls = {}
	
	# TO DO: Code factorization
	[s, o] = ['L', 'R'] if flag == 0 else ['R', 'L']
	for x in side_list:
		ctg, pls = x[0], pls_ids[s].inv[x[1]]
		if pls not in side_ctgs_by_pls:
			side_ctgs_by_pls[pls] = []
		side_ctgs_by_pls[pls].append(ctg)

	for x in opp_list:
		ctg, pls = x[0], pls_ids[o].inv[x[1]]
		if pls not in opp_ctgs_by_pls:
			opp_ctgs_by_pls[pls] = []
		opp_ctgs_by_pls[pls].append(ctg)	

	side_len, side_cost = 0, 0

	for node in side:
		partitions = [set(side_ctgs_by_pls[node])]
		for edge in B.edges:
			if edge[flag] == node:
				side_contigs = set(side_ctgs_by_pls[edge[flag]])
				opp_contigs = set(opp_ctgs_by_pls[edge[1-flag]])
				common = side_contigs.intersection(opp_contigs)
				partitions = modify_partitions(partitions, common)
		node_len, cost = get_partition_cost(node, partitions, side_list, contigs_dict)
		side_len += node_len
		side_cost += cost
	return side_len, side_cost

#Compute total cost of splits and joins for a particular matching
def compute_match_cost(left_list, right_list, pls_ids, contigs_dict):
	#Create graph with vertices named according to permutation and obtain connected components
	B = nx.Graph()
	B = add_nodes(B, left_list, right_list, pls_ids)
	B = add_edges(B, left_list, right_list, pls_ids)
	A = [B.subgraph(c) for c in nx.connected_components(B)]

	n_conn_comp = len(list(A))
	total_cost = 0		

	for i in range(n_conn_comp):
		C = list(A)[i]
		left, right = nx.bipartite.sets(C)	#Split the component according to bipartite sets
		if len(list(right)) != 0 and list(right)[0] in left_pls: #Ensuring proper assignments of bipartite parts
			right,left = left,right
		left_len, left_cost = one_side_cost(left, left_list, right_list, B, 0, pls_ids, contigs_dict)
		right_len, right_cost = one_side_cost(right, right_list, left_list, B, 1, pls_ids, contigs_dict)
		total_cost += left_cost + right_cost
	return total_cost

#Compute cost of matching for reached contigs
def compute_current_cost(matching, left_plasmids, right_plasmids, pls_ids, contigs_dict):
	left_renamed, right_renamed = rename_by_pmutn(matching, left_plasmids, right_plasmids, pls_ids)
	left_ctg_ids, right_ctg_ids = set(), set()
	for x in left_renamed:
		left_ctg_ids.update(x[0])
	for x in right_renamed:
		right_ctg_ids.update(x[0])

	#common_contigs = left_ctg_ids.intersection(right_ctg_ids)	
	#Delete contigs unique to a single tool
	#left_renamed, left_only_len = remove_unique_contigs(left_renamed, common_contigs)
	#right_renamed, right_only_len = remove_unique_contigs(right_renamed, common_contigs)	
	return compute_match_cost(left_renamed, right_renamed, pls_ids, contigs_dict)

def run_compare_plasmids(left_plasmids, right_plasmids, contigs_dict, pls_ids_dict):
	#left_keys = pls_ids_dict['L'].keys()
	#right_keys = pls_ids_dict['R'].keys()

	#Common contigs in the solutions obtained from both tools
	left_ctg_ids = set([ctg for ctg in contigs_dict.keys() if len(contigs_dict[ctg]['L_copies']) >= 1])
	right_ctg_ids = set([ctg for ctg in contigs_dict.keys() if len(contigs_dict[ctg]['R_copies']) >= 1])
	common_contigs = left_ctg_ids.intersection(right_ctg_ids)

	max_cost = 0		#Computing upperbound on final_cost
	n_permutations = {}
	for contig in common_contigs:
		m = len(contigs_dict[contig]['L_copies'])
		n = len(contigs_dict[contig]['R_copies'])
		max_cost += m * contigs_dict[contig]['length']
		max_cost += n * contigs_dict[contig]['length']	
		n_permutations[contig] = int(factorial(n)/factorial(n-m)) if n > m else int(factorial(m)/factorial(m-n))

	### BNB ###
	current_state = {'level': 0, 'cost': 0, 'matching': {}}
	final_state = {'cost': max_cost, 'matching': {}}

	contig_list = list(common_contigs)
	sorted_contig_list = sorted(contig_list, key=lambda ctg: n_permutations[ctg])

	#print(sorted_contig_list)

	def recursive_compare(current_state, sorted_contig_list, left_plasmids, right_plasmids, pls_ids_dict, contigs_dict):
		nonlocal final_state
		if current_state['level'] < len(sorted_contig_list):				#Compute cost upto this level
			current_contig = sorted_contig_list[current_state['level']]		#Retrieve contig for current level

			print(current_state['level'], sorted_contig_list[current_state['level']])
			
			m = len(contigs_dict[current_contig]['L_copies'])
			n = len(contigs_dict[current_contig]['R_copies'])			
			permutations = generate_permutations(m,n)

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
	
# Proposed solution avoiding nonlocal variables
	
def run_compare_plasmids(left_plasmids, right_plasmids, contigs_dict, pls_ids_dict):
	
	#Common contigs in the solutions obtained from both tools
	left_ctg_ids = set([ctg for ctg in contigs_dict.keys() if len(contigs_dict[ctg]['L_copies']) >= 1])
	right_ctg_ids = set([ctg for ctg in contigs_dict.keys() if len(contigs_dict[ctg]['R_copies']) >= 1])
	common_contigs = left_ctg_ids.intersection(right_ctg_ids)

	max_cost = 0		#Computing upperbound on final_cost
	n_permutations = {}
	for contig in common_contigs:
		m = len(contigs_dict[contig]['L_copies'])
		n = len(contigs_dict[contig]['R_copies'])
		max_cost += m * contigs_dict[contig]['length']
		max_cost += n * contigs_dict[contig]['length']	
		n_permutations[contig] = int(factorial(n)/factorial(n-m)) if n > m else int(factorial(m)/factorial(m-n))

	### BNB ###
	current_state = {'level': 0, 'cost': 0, 'matching': {}}
	final_state = {'cost': max_cost, 'matching': {}}

	contig_list = list(common_contigs)
	sorted_contig_list = sorted(contig_list, key=lambda ctg: n_permutations[ctg])

	def recursive_compare(current_state, sorted_contig_list, left_plasmids, right_plasmids, pls_ids_dict, contigs_dict, best_cost):
		best_sol_in_subtree = None
		if current_state['level'] < len(sorted_contig_list):				#Compute cost upto this level
			current_contig = sorted_contig_list[current_state['level']]		#Retrieve contig for current level

			print(current_state['level'], sorted_contig_list[current_state['level']])
			
			m = len(contigs_dict[current_contig]['L_copies'])
			n = len(contigs_dict[current_contig]['R_copies'])			
			permutations = generate_permutations(m,n)

			for perm in permutations:
				matching = get_matching_from_perm(contigs_dict[current_contig], perm)
				current_state['matching'][current_contig] = matching
				current_state['cost'] = compute_current_cost(current_state['matching'], left_plasmids, right_plasmids, pls_ids_dict, contigs_dict)	
				print(current_state['level'], sorted_contig_list[current_state['level']], perm, current_state['cost'],final_state['cost'])

				if current_state['cost'] < final_state['cost']:	
					current_state['level'] += 1 
					best_sol_in_subtree = recursive_compare(current_state, sorted_contig_list, left_plasmids, right_plasmids, pls_ids_dict, contigs_dict, best_cost)
					current_state['level'] -= 1
				
				del current_state['matching'][current_contig]

		else:
			return copy.deepcopy(current_state)
			print(current_state['cost'])

	optimal_solution = recursive_compare(current_state, sorted_contig_list, left_plasmids, right_plasmids, pls_ids_dict, contigs_dict, max_cost)

