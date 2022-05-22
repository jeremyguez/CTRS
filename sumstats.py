import numpy as np
import tskit, os, math, shutil
import pandas as pd
import psutil
import sys
import dadi


### DISCLAIMER: this is the code used for the computations in the submitted manuscript
### More efficient versions of many of these are under developpement in tskit package :
### https://github.com/tskit-dev/tskit/discussions/2245


def path_length (tree, parent, child):
    
    counter = 0
    
    c = child
    p = c
    
    while p != parent:
        p = tree.parent(c)
        c = p
        counter += 1
        
    return(counter)


# Computes indices that need looping on nodes

def nodes_sumstats (tree_seq):
	
	nodes_num_list = [] # Nodes count
	
	counter_poly_list = [] # Polytomies count
	
	# Colless index
	sum_diff_list = []
	mean_diff_list = []
	var_diff_list = []
	
	# B1
	b1_list = []
	b1_norm_list = []
	
	for tree in tree_seq.trees():	
		
		counter = 0 # Polytomies count
		diff_list = [] # Colless index
		max_list = [] # B1
	
		for n in tree.nodes():
			
			# Polytomies count
			if len(tree.children(n)) > 2:
				counter+=1
			
			
			# Colless index
			if n != tree.root and tree.is_leaf(n) == False :
				leaves_list = []
				for c in tree.children(n):
					leaves_list.append(tree.num_samples(c))	
				if len(leaves_list) > 0:
					diff_list.append(max(leaves_list)-min(leaves_list))
			
			
			# B1
			if n != tree.root and tree.is_leaf(n) == False :
				depth_list =  []
				for leaf in tree.leaves(n):
					depth_list.append(path_length(tree, n, leaf))
				max_list.append(max(depth_list))
		
		
		nodes_num_list.append(len(list(tree.nodes()))-(len(list(tree.leaves()))+1))
		
		counter_poly_list.append(counter) # Polytomies count
		
		sum_diff_list.append(np.sum(diff_list)) # Colless index
		mean_diff_list.append(np.mean(diff_list)) # Colless index
		var_diff_list.append(np.var(diff_list)) # Colless index
				
		b1_list.append(np.sum([1/m for m in max_list])) # B1
		b1_norm_list = [b1 / n_nodes for b1,n_nodes in zip(b1_list,nodes_num_list)]
		
		
		
	return((nodes_num_list, counter_poly_list, sum_diff_list, mean_diff_list, var_diff_list, b1_list, b1_norm_list))



# Computes indices that need looping on leaves

def leaves_sumstats (tree_seq, sample_size, vec_bins_hist2):
	

	b2_list = []
	
	sackin_sum =  []
	sackin_mean = []
	sackin_var = []
	diff_max_min = []
	depth_hist_densities = []
	
	counter = 0
	
	for tree in tree_seq.trees():
			
		depth_list = []
		prob_list = []
		
		counter += 1
		
		for l in tree.leaves():
			
			depth_list.append(tree.depth(l))
			
			random_walk = []
			p = l
			
			while p != tree.root:
				p = tree.parent(l)
				l = p
				random_walk.append(1/tree.num_children(p))
			 
			prob_list.append(np.prod(random_walk))
			del random_walk
			
		depth_list_normalized = [val / math.log2(sample_size*2) for val in depth_list]
		
		sackin_sum.append(np.sum(depth_list)) # Sackin
			
		sackin_mean.append(np.mean(depth_list_normalized)) # Sackin mean
		
		sackin_var.append(np.var(depth_list_normalized)) # Sackin var
		
		depth_hist = np.histogram(depth_list_normalized, bins = vec_bins_hist2, density=False)
		depth_hist_densities.append(depth_hist[0]/(sample_size*2)) # Histogram of depths
		
		diff_max_min.append(np.quantile(depth_list_normalized,0.5)) # New index
			
		b2_list.append(-np.sum([prob*np.log2(prob) for prob in prob_list])) # B2
		
		
	return((sackin_sum, sackin_mean, sackin_var, diff_max_min, depth_hist_densities, b2_list))
		
		
# Computes Ib (based on a script from Brandenburg et al. (2012))
		
def Ib (tree_seq, current_scenario, current_replicate, generation):
	
	# 1) Newick conversion
		
	newick_string = ""
	for tree in tree_seq.aslist():
		newick_string = newick_string + str(tree.newick()) + "\n"
	
	newick_name = "simu_scenario_" + str(current_scenario) + "/simupop" + str(current_replicate) + "-" + str(generation) + ".newick"
	with open(newick_name, 'w') as f:
		f.write(newick_string)


	# 2) JT conversion
	
	newick_name_out =  str(newick_name) + ".out"
	os.system("../format_newick.py " + str(newick_name) + " " + str(newick_name_out))


	# 3) Imbalance computation
	
	newick_name_results =  str(newick_name) + ".results"
	
		
	os.system("../CalcIPrimeV2_linux64.exe -f " + str(newick_name_out) + " > " + str(newick_name_results))
	
	data = pd.read_table(newick_name_results)
	
	return(data['IprimeNonBin2N'])



def Ib_simple (tree_seq, current_replicate):
	
	# 1) Newick conversion
		
	newick_string = ""
	for tree in tree_seq.aslist():
		newick_string = newick_string + str(tree.newick()) + "\n"
	
	dir_name = "compute_JT_expectation"
	
	try:
		shutil.rmtree(dir_name)
	except:
		pass
	
	os.mkdir(dir_name)
	
	newick_name = dir_name + "/simupop" + str(current_replicate) + ".newick"
	with open(newick_name, 'w') as f:
		f.write(newick_string)


	# 2) JT conversion
	
	newick_name_out =  str(newick_name) + ".out"
	os.system("./format_newick.py " + str(newick_name) + " " + str(newick_name_out))


	# 3) Imbalance computation
	
	newick_name_results =  str(newick_name) + ".results"
	
		
	os.system("./CalcIPrimeV2_linux64.exe -f " + str(newick_name_out) + " > " + str(newick_name_results))
	
	data = pd.read_table(newick_name_results)
	
	return(data['IprimeNonBin2N'])


### Function for inferring demographic parameters with dadi

def dadi_inference(ts):

	# Compute SFS on the previously loaded tree sequence ts using tskit
	sfs = ts.allele_frequency_spectrum(mode="site", windows=None, polarised=True, span_normalise=False)
				
	pts = [20]

	data = dadi.Spectrum(sfs)
	func = dadi.Demographics1D.two_epoch
	func_ex = dadi.Numerics.make_extrap_log_func(func)


	upper_bound = [100, 10]
	lower_bound = [0.01, 0]

	test_retry = []

	for test in range(0,2):
		p0 = [random.uniform(1,5),random.uniform(1,2)]
		p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)
		popt = dadi.Inference.optimize_log_fmin(p0, data, func_ex, pts, lower_bound=lower_bound, upper_bound=upper_bound)
		if popt[0] > 1/99 and popt[0] < 99 and popt[1] < 9.9 and popt[1] > 0.01:
			test_retry.append(popt)

	if not test_retry:
		test_retry.append(popt)
		
	popt_kept = np.median(test_retry, axis=0)

	model = func_ex(popt_kept, [sample_size*2], pts)
	theta = dadi.Inference.optimal_sfs_scaling(model, data)

	ancient_pop_size_dadi = theta/(4*mu*L)
		
	return((popt_kept[0], ancient_pop_size_dadi))

