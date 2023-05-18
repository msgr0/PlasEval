import networkx as nx
import queue
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
def rename_by_pmutn(matchings, left_dict, right_dict, pls_ids):
	reached_contigs = matchings.keys()
	left_renamed, right_renamed = {}, {}
	for contig in reached_contigs:
		M = matchings[contig]
		left_ctgs = M[0]
		right_ctgs = M[1]
		left_renamed = update_contig_name(left_dict, contig, left_ctgs, left_renamed, pls_ids['left'])
		right_renamed = update_contig_name(right_dict, contig, right_ctgs, right_renamed, pls_ids['right'])
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
def compute_match_cost(left_dict, right_dict):
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

#Compute cost of matching for reached contigs
def compute_current_cost(matching, left_plasmids, right_plasmids, pls_ids):
	left_renamed, right_renamed = rename_by_pmutn(matching, left_plasmids, right_plasmids, pls_ids)
	left_ctg_ids, right_ctg_ids = set(), set()
	for x in left_renamed:
		left_ctg_ids.update(set(left_renamed[x].keys()))
	for x in right_renamed:
		right_ctg_ids.update(set(right_renamed[x].keys()))
	common_contigs = left_ctg_ids.intersection(right_ctg_ids)	
	#Delete contigs unique to a single tool
	left_renamed, left_only_len = remove_unique_contigs(left_renamed, common_contigs)
	right_renamed, right_only_len = remove_unique_contigs(right_renamed, common_contigs)	
	return compute_match_cost(left_renamed, right_renamed)


def run_compare_plasmids(left_plasmids, right_plasmids):
	#Creating a dictionary with contigs as keys from dictionary of plasmids and contig lengths and sequences as values.
	left_contigs = create_contigs_dict(left_plasmids)
	right_contigs = create_contigs_dict(right_plasmids)

	left_pls_ids = list(left_plasmids.keys())
	right_pls_ids = list(right_plasmids.keys())
	pls_ids = {'left': left_pls_ids, 'right': right_pls_ids}

	#Common contigs in the solutions obtained from both tools
	common_contigs = []
	left_ctg_ids = set(left_contigs.keys())
	right_ctg_ids = set(right_contigs.keys())
	common_contigs = left_ctg_ids.intersection(right_ctg_ids)

	permutations = {}
	max_cost = 0		#Computing upperbound on final_cost
	for contig in common_contigs:
		m = left_contigs[contig]['copies']
		n = right_contigs[contig]['copies']
		max_cost += m * left_contigs[contig]['length']
		max_cost += n * right_contigs[contig]['length']

		permutations = generate_permutations(contig, m, n, permutations)

	### BNB ###
	current_state = {'level': 0, 'cost': 0, 'matching': {}}
	final_state = {'cost': max_cost, 'matching': {}}

	contig_list = list(common_contigs)
	sorted_contig_list = sorted(contig_list, key=lambda ctg:len(permutations[ctg]))

	def recursive_compare(current_state, permutations, sorted_contig_list, left_plasmids, right_plasmids, pls_ids):
		nonlocal final_state
		if current_state['level'] < len(sorted_contig_list):	#Compute cost upto this level
			current_contig = sorted_contig_list[current_state['level']]		#Retrieve contig for current level
			print(current_state['level'], sorted_contig_list[current_state['level']])
			for perm in permutations[current_contig]:
				current_state['matching'][current_contig] = perm
				current_state['cost'] = compute_current_cost(current_state['matching'], left_plasmids, right_plasmids, pls_ids)	
				print(current_state['level'], sorted_contig_list[current_state['level']], perm, current_state['cost'],final_state['cost'])


				if current_state['cost'] < final_state['cost']:	
					current_state['level'] += 1 
					recursive_compare(current_state, permutations, sorted_contig_list, left_plasmids, right_plasmids, pls_ids)
					current_state['level'] -= 1
				
				del current_state['matching'][current_contig]

		else:
			final_state['cost'] = current_state['cost']
			final_state['matching'] = copy.deepcopy(current_state['matching'])
			print(final_state['cost'])
	
		
			


	recursive_compare(current_state, permutations, sorted_contig_list, left_plasmids, right_plasmids, pls_ids)

	print(final_state['cost'])
	print(final_state['matching'])
