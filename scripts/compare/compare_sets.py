import networkx as nx
import itertools
import copy
import random

#Function to create a dictionary of contigs in order to record lengths, sequences and number of copies for each.
def create_contigs_dict(plasmids_dict):
	contigs_dict = {}
	for plasmid in plasmids_dict:
		for contig in plasmids_dict[plasmid]:
			if contig not in contigs_dict:
				contigs_dict[contig] = {}
				contigs_dict[contig]['length'] = plasmids_dict[plasmid][contig]['length']
				contigs_dict[contig]['copies'] = 0
			contigs_dict[contig]['copies'] += plasmids_dict[plasmid][contig]['copies']
	return contigs_dict			

#Generate list of permutations by contig family
def generate_permutations(contig, m, n, permutations):
	permutations[contig] = []
	if m >= n:
		pmutns = list(itertools.permutations(list(range(m))))
		n_list = list(range(n))
		for pmutn in pmutns:
			permutations[contig].append((list(pmutn), n_list))
	else:
		pmutns = list(itertools.permutations(list(range(n))))
		m_list = list(range(m))
		for pmutn in pmutns:
			permutations[contig].append((m_list, list(pmutn)))
	return permutations

#Update the name of contig
def update_contig_name(plasmids_dict, contig, idxs, new_dict, pls_ids):
	count = 0

	for plasmid in pls_ids:
		if contig in plasmids_dict[plasmid].keys():
			copies = plasmids_dict[plasmid][contig]['copies']
			length = plasmids_dict[plasmid][contig]['length']
			if plasmid not in new_dict:
				new_dict[plasmid] = {}
			for i in range(copies):
				idx = idxs[count]
				contig_new = contig+'_'+str(idx)
				new_dict[plasmid][contig_new] = {}
				new_dict[plasmid][contig_new]['length'] = length
				new_dict[plasmid][contig_new]['copies'] = 1
				count += 1			
	return new_dict

#Rename contigs in both dictionaries according to the permutation
def rename_by_pmutn(matchings, left_dict, right_dict, lpls, rpls):
	reached_contigs = matchings.keys()
	ld, rd = {}, {}
	for contig in reached_contigs:
		M = matchings[contig]
		ll = M[0]
		rl = M[1]
		ld = update_contig_name(left_dict, contig, ll, ld, lpls)
		rd = update_contig_name(right_dict, contig, rl, rd, rpls)
	return ld, rd

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
def add_nodes(G, left_dict, right_dict):
	G.add_nodes_from(list(left_dict.keys()), bipartite=0)
	G.add_nodes_from(list(right_dict.keys()), bipartite=1)
	return G

#Adds edges to the bipartite graph G
def add_edges(G, left_dict, right_dict):
	edges_list = []
	for left_plasmid in left_dict:
		for right_plasmid in right_dict:
			left_keys = set(left_dict[left_plasmid].keys())
			right_keys = set(right_dict[right_plasmid].keys())
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
			#if S in modified_partitions:
			#	modified_partitions.remove(S)
		elif S.intersection(common) == S:
			modified_partitions.append(S)				
	return modified_partitions		

#Obtain cost of a partition
def get_partition_cost(node, partitions, plasmids):
	total_len = 0
	cost = 0
	for S in partitions:
		S_cost = 0
		for contig in S:
			contig_len = plasmids[node][contig]['length']
			total_len += contig_len
			S_cost += contig_len
			cost = max(cost, S_cost)
	cost = total_len - cost
	return total_len, cost

#Compute total cost on one side (splits or joins)
def one_side_cost(side, side_dict, opp_dict, B, flag):
	side_len, side_cost = 0, 0
	for node in side:
		partitions = [set(side_dict[node].keys())]
		for edge in B.edges:
			if edge[flag] == node:
				side_contigs = set(side_dict[edge[flag]].keys())
				opp_contigs = set(opp_dict[edge[1-flag]].keys())
				common = side_contigs.intersection(opp_contigs)
				partitions = modify_partitions(partitions, common)
		node_len, cost = get_partition_cost(node, partitions, side_dict)
		side_len += node_len
		side_cost += cost
	return side_len, side_cost

#Compute total cost of splits and joins for a particular matching
def compute_current_cost(left_dict, right_dict):
	#Create graph with vertices named according to permutation and obtain connected components
	B = nx.Graph()
	B = add_nodes(B, left_dict, right_dict)
	B = add_edges(B, left_dict, right_dict)
	A = [B.subgraph(c) for c in nx.connected_components(B)]
	n_conn_comp = len(list(A))
	total_cost = 0						
	for i in range(n_conn_comp):
		C = list(A)[i]
		left, right = nx.bipartite.sets(C)	#Split the component according to bipartite sets
		#Assign the split sets to the correct tool (as the bipartite split is random)
		temp = []
		if len(list(right)) != 0 and list(right)[0] in left_dict.keys():
			temp = right
			right = left
			left = temp
		if len(list(left)) != 0 and list(left)[0] in right_dict.keys():
			temp = right
			right = left
			left = temp					
		if len(list(left)) == 0:	#If right side plasmid forms a singleton vertex	
			plasmid = list(right)[0]
		elif len(list(right)) == 0:	#If left side plasmid forms a singleton vertex
			plasmid = list(left)[0]
		else:						#If component has vertices from both sides
			left_len, left_cost = one_side_cost(left, left_dict, right_dict, B, 0)
			right_len, right_cost = one_side_cost(right, right_dict, left_dict, B, 1)
			total_cost += left_cost + right_cost
	return total_cost


def recursive_compare(current_level, final_matching, final_cost, current_matching, current_cost, permutations, contig_list, left_plasmids, right_plasmids, lpls, rpls, test_flag):
	if current_level >= len(contig_list):
		final_cost = current_cost
		final_matching = copy.deepcopy(current_matching)
		print(final_cost)
		test_flag += 1
		print(test_flag)
		#return final_cost, test_flag
		#print(final_matching)
		#return current_matching, final_cost
	else: 		#Compute cost upto this level
		test_flag += 1
		print("\n") 
		print(current_level, contig_list[current_level], test_flag)
		current_contig = contig_list[current_level]
		for i in range(len(permutations[current_contig])):
			current_matching[current_contig] = permutations[current_contig][i]
			ld, rd = rename_by_pmutn(current_matching, left_plasmids, right_plasmids, lpls, rpls)
			left_keys, right_keys = set(), set()
			for x in ld:
				left_keys.update(set(ld[x].keys()))
			for x in rd:
				right_keys.update(set(rd[x].keys()))
			common_contigs = left_keys.intersection(right_keys)	
			#Delete contigs unique to a single tool
			ld, left_only_len = remove_unique_contigs(ld, common_contigs)
			rd, right_only_len = remove_unique_contigs(rd, common_contigs)	
			current_cost = compute_current_cost(ld, rd)
			if current_cost < final_cost:
				print(current_cost, final_cost)
				current_level += 1 
				recursive_compare(current_level, final_matching, final_cost, current_matching, current_cost, permutations, contig_list, \
					      left_plasmids, right_plasmids, lpls, rpls, test_flag)

def run_compare_plasmids(left_plasmids, right_plasmids):
	#Creating a dictionary with contigs as keys from dictionary of plasmids and contig lengths and sequences as values.
	left_contigs = create_contigs_dict(left_plasmids)
	right_contigs = create_contigs_dict(right_plasmids)

	lpls = list(left_plasmids.keys())
	rpls = list(right_plasmids.keys())

	#Common contigs in the solutions obtained from both tools
	common_contigs = []
	left_keys = set(left_contigs.keys())
	right_keys = set(right_contigs.keys())
	common_contigs = left_keys.intersection(right_keys)

	permutations = {}
	final_cost = 0
	for contig in common_contigs:
		m = left_contigs[contig]['copies']
		n = right_contigs[contig]['copies']
		final_cost += m * left_contigs[contig]['length']
		final_cost += n * right_contigs[contig]['length']

		permutations = generate_permutations(contig, m, n, permutations)

	### BNB ###
	current_matching = {}
	final_matching = {}
	contig_list = list(common_contigs)

	current_level = 0
	current_cost = 0

	test_flag = 0
	recursive_compare(current_level, final_matching, final_cost, current_matching, current_cost, permutations, contig_list, \
					      left_plasmids, right_plasmids, lpls, rpls, test_flag)



