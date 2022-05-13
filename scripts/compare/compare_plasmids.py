import networkx as nx
import itertools
import copy

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

#Function to rename copies of the same contig in order to differentiate them from each other
def rename_contigs(contig, plasmids_dict):
	count = 0
	for plasmid in plasmids_dict:
		if contig in plasmids_dict[plasmid].keys():
			length = plasmids_dict[plasmid][contig]['length']
			copies = plasmids_dict[plasmid][contig]['copies']
			for i in range(copies):
				plasmids_dict[plasmid][contig+'_'+str(count)] = {}
				plasmids_dict[plasmid][contig+'_'+str(count)]['length'] = length
				plasmids_dict[plasmid][contig+'_'+str(count)]['copies'] = 1
				count += 1
			plasmids_dict[plasmid].pop(contig)
	return plasmids_dict

#Obtain list of permutations 
def get_permutations(copies, count):
	return list(itertools.permutations(copies, count))

#Obtain list of combinations of permutations
def get_combinations(pmutn_list):
	return list(itertools.product(*pmutn_list, repeat=1))

#Update the name of contig
def update_contig_name(plasmids_dict, contig, idxs):
	count = 0
	for plasmid in plasmids_dict:
		if contig in plasmids_dict[plasmid].keys():
			copies = plasmids_dict[plasmid][contig]['copies']
			length = plasmids_dict[plasmid][contig]['length']
			for i in range(copies):
				idx = idxs[count]
				contig_new = contig+'_'+str(idx)
				plasmids_dict[plasmid][contig_new] = {}
				plasmids_dict[plasmid][contig_new]['length'] = length
				plasmids_dict[plasmid][contig_new]['copies'] = 1
				count += 1			
			plasmids_dict[plasmid].pop(contig)
	return plasmids_dict

#Rename contigs in both dictionaries according to the permutation
def rename_by_pmutn(pmutn, left_dict, right_dict, ref_array):
	renamed_contigs = []
	for item in ref_array:
		contig, m, n = item[0], item[1], item[2]
		k = ref_array.index(item)
		idxs = pmutn[k]

		if m < n:
			left_dict = update_contig_name(left_dict, contig, idxs)
		else:
			right_dict = update_contig_name(right_dict, contig, idxs)
	return left_dict, right_dict

#Create database for scores with iteration number for a labeling as the key. 
#The values are nested dictionaries consisting of scores, plasmid splits and dictionaries corresponding to the labeling.  
def create_score_db(iter_no, score_db):
	score_db[iter_no] = {}
	score_db[iter_no]['score'] = 0 #default
	score_db[iter_no]['left_splits'] = []
	score_db[iter_no]['right_splits'] = []
	score_db[iter_no]['left_dict'] = []
	#score_db[iter_no]['left_lengths'] = []
	score_db[iter_no]['right_dict'] = []
	#score_db[iter_no]['right_lengths'] = []
	score_db[iter_no]['decomposition'] = []
	return score_db

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

#Obtain length of plasmid as sum of contig lengths
def get_plasmid_length(plasmid, plasmids_dict):
	length = 0
	for contig in plasmids_dict[plasmid]:
		contig_len = plasmids_dict[plasmid][contig]['length']
		length += contig_len
	return length		

#Weighted score of all component scores
def get_overall_score(component_scores):
	cost = 0
	length = 0
	for comp in component_scores:
		cost += (comp[0]*comp[1])
		length += comp[1]
	score = cost/length
	return score

def run_compare_plasmids(HyAsP_plasmids, MOBsuite_plasmids):
	#Creating a dictionary with contigs as keys from dictionary of plasmids and contig lengths and sequences as values.
	HyAsP_contigs = create_contigs_dict(HyAsP_plasmids)
	MOBsuite_contigs = create_contigs_dict(MOBsuite_plasmids)

	#Common contigs in the solutions obtained from both tools
	common_contigs = []
	HyAsP_keys = set(HyAsP_contigs.keys())
	MOBsuite_keys = set(MOBsuite_contigs.keys())
	common_contigs = HyAsP_keys.intersection(MOBsuite_keys)

	#Form a list of combinations of renamed contigs
	#1. Iterate through list of contigs common to both solutions
	#2. For each common contig, compare number of copies of contig in both solutions
	#3. Fix the order of indexing in tool with more copies
	#4.1. Obtain list of permuations of contig indices in other tool
	#4.2. Obtain all the combinations of indices for all contigs 
	ref_array = []	#reference array of tuples (contig_id, no. of copies in HyAsP, no. of copies in MOBsuite)
	pmutn_list = []	#list of permutations of indices for each contig, ordered according to the contigs in reference array
	for contig in common_contigs:
		m = HyAsP_contigs[contig]['copies']
		n = MOBsuite_contigs[contig]['copies']
		if m >= n:
			HyAsP_plasmids = rename_contigs(contig, HyAsP_plasmids)
			ref_array.append((contig, m, n))
			copies = list(range(m))
			pmutns = get_permutations(copies, n)
			pmutn_list.append(pmutns)
		else:
			MOBsuite_plasmids = rename_contigs(contig, MOBsuite_plasmids)
			ref_array.append((contig, m, n))
			copies = list(range(n))
			pmutns = get_permutations(copies, m)
			pmutn_list.append(pmutns)	
	combinations = get_combinations(pmutn_list)

	#####################################################
	#	CONSTRUCTING THE GRAPH, COMPUTING THE SPLITS    #
	#	AND JOINS FOR EACH COMPONENT, COMPUTING THE     #
	#	SCORE FOR EACH COMPONENT AND OVERALL SCORE   	# 
	#####################################################
	score_db = {}	#Keep track of scores, splits and contig names for each combination/labeling
	count = 0
	for x in combinations:
		
		print("Combination no: ", count)
		score_db = create_score_db(count, score_db)

		left_dict = copy.deepcopy(HyAsP_plasmids)
		right_dict = copy.deepcopy(MOBsuite_plasmids)	
		left_dict, right_dict = rename_by_pmutn(x, left_dict, right_dict, ref_array)
		score_db[count]['left_dict'] = left_dict
		score_db[count]['right_dict'] = right_dict

		left_keys, right_keys = set(), set()
		for x in left_dict:
			left_keys.update(set(left_dict[x].keys()))
		for x in right_dict:
			right_keys.update(set(right_dict[x].keys()))
		common_contigs = left_keys.intersection(right_keys)	
	
		#Delete contigs unique to a single tool
		left_dict, left_only_len = remove_unique_contigs(left_dict, common_contigs)
		right_dict, right_only_len = remove_unique_contigs(right_dict, common_contigs)

		component_scores= []	#Keep track of scores for each component. 
								#If component involves perfect match between plasmids, score should be 0.
		component_scores.append((1, left_only_len))
		component_scores.append((1, right_only_len))
		score_db[count]['decomposition'].append(left_only_len)
		score_db[count]['decomposition'].append(right_only_len)								


		#Create graph with vertices named according to permutation and obtain connected components
		B = nx.Graph()
		B = add_nodes(B, left_dict, right_dict)
		B = add_edges(B, left_dict, right_dict)
		A = [B.subgraph(c) for c in nx.connected_components(B)]
		n_conn_comp = len(list(A))
							
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
				score_db[count]['right_splits'].append([(set(right_dict[plasmid].keys()))])
				length = get_plasmid_length(plasmid, right_dict)
				component_scores.append((1, length))
			elif len(list(right)) == 0:	#If left side plasmid forms a singleton vertex
				plasmid = list(left)[0]
				score_db[count]['left_splits'].append([(set(left_dict[plasmid].keys()))])
				length = get_plasmid_length(plasmid, left_dict)
				component_scores.append((1, length))
			else:						#If component has vertices from both sides
				left_len, right_len = 0, 0
				left_cost, right_cost = 0, 0 

				for node in left:
					partitions = [set(left_dict[node].keys())]
					for edge in B.edges:
						if edge[0] == node:
							left_contigs = set(left_dict[edge[0]].keys())
							right_contigs = set(right_dict[edge[1]].keys())
							common = left_contigs.intersection(right_contigs)
							partitions = modify_partitions(partitions, common)
							score_db[count]['left_splits'].append(partitions)
					node_len, cost = get_partition_cost(node, partitions, left_dict)
					left_len += node_len
					left_cost += cost

				for node in right:	
					partitions = [set(right_dict[node].keys())]
					for edge in B.edges:
						if edge[1] == node:
							left_contigs = set(left_dict[edge[0]].keys())
							right_contigs = set(right_dict[edge[1]].keys())
							common = right_contigs.intersection(left_contigs)
							partitions = modify_partitions(partitions, common)			
							score_db[count]['right_splits'].append(partitions)
					node_len, cost = get_partition_cost(node, partitions, right_dict)			
					right_len += node_len
					right_cost += cost	

				normalized_cost = (left_cost/left_len + right_cost/right_len)/2

				component_scores.append((normalized_cost, left_len))	

		print(component_scores)
		weighted_len = 0
		for pair in component_scores[2:]:
			weighted_len += pair[0]*pair[1]
		score_db[count]['decomposition'].append(weighted_len)		
		score = get_overall_score(component_scores)	
		score_db[count]['score'] = score	
		count+=1

		print("\n")	

	return score_db		