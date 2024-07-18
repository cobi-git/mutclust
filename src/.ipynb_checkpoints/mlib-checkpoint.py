from math import ceil
from os.path import join
import sys
import numpy as np

## hyperparams
# min_ps=11	# weakest clade mut's per_sum=11.52 (next 26.87)
# maxeps=1000
# eps_scaler_const=1.5
# exp_dim=30.	# the eps_scaler diminishing effect (higher the exp_dim the lesser the es penalized)
# min_clength=10	# minimum length of cluster

from src.utils import mutation_filtering

POS = 'Position'
FREQ = 'Frequency'
PER = 'Percentage'
ENT = 'Entropy'

HSCORE = 'H-score'
HSCORE_SUM = 'H-score_sum'
HSCORE_AVR = 'H-score_avr'

PER_SUM = 'per_sum'
ENT_SUM = 'ent_sum'
PER_AVR = 'per_avr'
ENT_AVR = 'ent_avr'

MAX_DEPS = 50
MAX_EPS = 1000
MIN_EPS= 10

EPSILON = 5
EPSILON_SCALING_FACTOR = 10

DIMINISHING_FACTOR = 3

MIN_CLUSTER_LENGTH = 10

CCM_MIN_FREQUENCY = 1000
CCM_MIN_PERCENTAGE = 0.01
CCM_MIN_ENTROPY = 0.45

CCM_MIN_PERCENTAGE_SUM = 0.01
CCM_MIN_ENTROPY_SUM = 0.05

CCM_MIN_PERCENTAGE_AVR = 0.001
CCM_MIN_ENTROPY_AVR = 0.01

CCM_MIN_HSCORE_SUM = 0.05
CCM_MIN_HSCORE_AVR = 0.01
CCM_MIN_HSCORE = 0.03

MIN_MUTATIONS = 5
MIN_CLUSTERS = 5

class Cluster(object):
	def __init__(self, cid, idx, f, af):
		self.cid=cid
		self.cpts=[idx]
		self.totf=[f] 	# total num of mut freq in cluster
		self.af=af 		# average mut freq in cluster
		self.length=0 	# length of cluster
		self.size=1		# num of muts in cluster
		self.seq=''		# sequence of cluster

	def add_cpt(self, idx, f):
		self.cpts.append(idx)
		self.totf.append(f)
		self.af=np.mean(self.totf)	# update average mut freq of cluster
		self.length=self.cpts[-1]-self.cpts[0]	# updated cluster length(bp)
		self.size+=1				# update cluster size

	def __str__(self):
		return '%d\t%d:%d\t%d\t%d\t%.3f\t%s'%(self.cid, self.cpts[0], self.cpts[-1], self.length, self.size, self.af, self.seq)

def init(d, info):
	print('\n--- Configurations ---')
	print('Input data: \'%s\' %s'%(info.fin, d.shape))
	print('Output dir: \'%s\''%(info.outdir))

	print('Parameters:')
	print('  Min Eps=%d'%(info.eps))
	print('  Max Eps=%d'%(info.maxeps))
	print('  Min per_sum=%.1f' % (info.min_persum))
	print('  Eps scaling factor=%.1f'%(info.eps_scaler_const))
	print('  Expansion diminishing factor=%d' % (info.es_control_const))
	print('  Min cluster length=%d' % (info.min_cluster_length))
	print('----------------------\n')

class get_eps_stats(object):
	def __init__(self, idx, pos, df, lr_index, lr_distance, es):		
		self.idx=idx
		self.i=pos
		self.i_per=df.loc[idx, PER]
		self.i_freq=df.loc[idx, FREQ]
		self.i_ent=df.loc[idx, ENT]
		self.i_hscore = df.loc[idx, HSCORE]
		
		self.l_dist= lr_distance[0]
		self.r_dist = lr_distance[1]
		ccm_df = df.loc[lr_index[0]:lr_index[1] + 1, :]

		self.length = len(ccm_df)
		self.l_pos = df.loc[lr_index[0], POS]
		self.r_pos = df.loc[lr_index[1], POS]
		
		self.mut_n = len(ccm_df[ccm_df[HSCORE]>0])
		self.eps_scaler = es 
		
		self.freq_sum= ccm_df[FREQ].sum()
		self.freq_avr=self.freq_sum/self.length
		
		self.per_sum= ccm_df[PER].sum()
		self.per_avr=self.per_sum/self.length
		
		self.ent_sum= ccm_df[ENT].sum()
		self.ent_avr=self.ent_sum/self.length

		self.hscore_sum=ccm_df[HSCORE].sum()
		self.hscore_avr=self.hscore_sum/self.length

	def to_list(self):
		param_list=[self.idx, self.i, self.i_freq, self.i_per, self.length, self.freq_sum, self.freq_avr, self.per_sum, self.per_avr, self.ent_sum, self.ent_avr, self.eps_scaler, self.l_dist, self.r_dist]
		return param_list

	def to_dict(self):
		param_dict={'index':self.idx, POS:self.i, FREQ:self.i_freq, PER:self.i_per,ENT:self.i_ent,HSCORE:self.i_hscore, 'length':self.length, 'freq_sum':self.freq_sum, 'freq_avr':self.freq_avr,
					PER_SUM:self.per_sum, PER_AVR:self.per_avr, ENT_SUM:self.ent_sum, ENT_AVR:self.ent_avr, HSCORE_SUM:self.hscore_sum, HSCORE_AVR:self.hscore_avr, 
					'eps_scaler':self.eps_scaler, 'left_distance':self.l_dist, 'right_distance':self.r_dist, 'l_pos':self.l_pos, 'r_pos':self.r_pos, 'mut_n':self.mut_n}
		return param_dict

def get_eps_region(mutclust_input_df, cur_index, info):
	dist_l = dist_r = 0	# left and right distance from i
	cur_l_index =  cur_index - 1
	cur_r_index = cur_index + 1
	cur_hscore = mutclust_input_df.loc[cur_index, HSCORE]

	eps_scaler = ceil(EPSILON_SCALING_FACTOR * cur_hscore)
	ldeps = rdeps = eps_scaler*EPSILON

	if cur_l_index < 0:
		cur_l_index = cur_index
	if cur_r_index >= mutclust_input_df.shape[0]-1:
		cur_r_index = cur_index
	
	# search left boundary
	if ldeps > info.maxeps: 
		ldeps = info.maxeps # max at maxeps
		
	while dist_l < ldeps and cur_l_index >= 0: 
		ld = mutclust_input_df.loc[cur_index, POS] - mutclust_input_df.loc[cur_l_index, POS]
		if ld > ldeps:
			break
		cur_l_index -= 1
		dist_l = ld
	cur_l_index += 1

	# search right boundary
	if rdeps > info.maxeps:
		rdeps = info.maxeps	# max at maxeps

	while dist_r < rdeps and cur_r_index < mutclust_input_df.shape[0]:
		# print(rdeps)
		rd = mutclust_input_df.loc[cur_r_index, POS] - mutclust_input_df.loc[cur_index, POS]
		if rd > rdeps:
			break
		cur_r_index += 1
		dist_r = rd
	cur_r_index -= 1 
	
	return [cur_l_index, cur_r_index], [cur_index-cur_l_index, cur_r_index-cur_index], eps_scaler

def expand_cluster(ccm_idx, total_mutation_info_list, info):
	left_cur_dist = right_cur_dist = 0			# left and right distance from pos
	left_cur_index = ccm_idx - 1				# left moving index
	right_cur_index = ccm_idx + 1				# right moving index
	mut_n = len(total_mutation_info_list) 
	if right_cur_index >= mut_n:
		right_cur_index = ccm_idx

	es_l = es_r = total_mutation_info_list[ccm_idx]['eps_scaler']		# ccm's pre-computed eps scaler (es)(per/1.5)
	left_max_dist = total_mutation_info_list[ccm_idx]['left_distance']
	right_max_dist = total_mutation_info_list[ccm_idx]['right_distance']		# ccm's pre-computed deps length (dist)

	# expand left
	while left_cur_dist<left_max_dist and left_cur_index>=0:
		ld = total_mutation_info_list[ccm_idx][POS] - total_mutation_info_list[left_cur_index][POS]
		if ld > left_max_dist:
			break
		left_cur_dist = ld
		
		# decrease deps in respect to es of cur_l
		delta_es = es_l - total_mutation_info_list[left_cur_index]['eps_scaler']	# delta(eps)=es_cur-es_ccm
		es_l = es_l - (delta_es) / info.es_control_const# diminish es by delta(eps)/exp_dim delta/30
		mut_deps = info.eps * es_l		

		if mut_deps > 0:
			left_max_dist = mut_deps
		else:
			break
		left_cur_index -= 1

	# expand right
	while right_cur_dist<right_max_dist and right_cur_index < mut_n:
		rd = total_mutation_info_list[right_cur_index][POS] - total_mutation_info_list[ccm_idx][POS]
		if rd>right_max_dist:
			break
		right_cur_dist=rd
		# decrease deps in respect to es of cur_r
		delta_es = es_r - total_mutation_info_list[right_cur_index]['eps_scaler']	# delta(eps)=eps_i-eps_curl
		es_r = es_r - (delta_es)/info.es_control_const	# diminish es by delta(eps)/exp_dim (default: 30)
		mut_deps = info.eps * es_r

		if mut_deps > 0:
			right_max_dist = mut_deps
		else:
			break
		right_cur_index+=1

	if right_cur_index == mut_n:
		right_cur_index-=1
	if left_cur_index < 0:
		left_cur_index = 0 

	ret_dict = { 'length': total_mutation_info_list[right_cur_index][POS] - total_mutation_info_list[left_cur_index][POS] + 1,
				'ccm_position':ccm_idx,
				'mut_positions': sorted([a[POS] for a in total_mutation_info_list[left_cur_index:right_cur_index+1] if a[HSCORE] > 0])}
	ret_dict['left_position'] = ret_dict['mut_positions'][0]
	ret_dict['right_position'] = ret_dict['mut_positions'][-1]

	return ret_dict

def dynaclust(total_mutation_info_list, ccm_index_list, info, tag, i):
	print('Perfoming dynamic clustering...'),
	cluster_list=[]
	for ccm_idx in ccm_index_list:
		ret_dict = expand_cluster( ccm_idx, total_mutation_info_list, info)
		cluster_list.append(ret_dict)
		sys.stdout.write('\r{0} clusters found'.format(len(cluster_list)))
		#print(ret_dict)

	print()

	#merging clusters
	print('merging clusters...')
	merged_clusters=[]

	cluster_list = sorted(cluster_list, key=lambda x: x['left_position'])
	clst_n = len(cluster_list)
	i = 0
	while i < clst_n:
		lpos=cluster_list[i]['left_position']
		rpos=cluster_list[i]['right_position']
		mut_list=cluster_list[i]['mut_positions']
		j = i+1
		while j < clst_n:
			if rpos < cluster_list[j]['left_position']: # cluster not in cluster next 
				i = j
				break
			else: # cluster in cluster 
				lpos = min(lpos, cluster_list[j]['left_position'])
				rpos = max(rpos, cluster_list[j]['right_position'])
				mut_list.extend(cluster_list[j]['mut_positions'])
				mut_list = [a for a in set(mut_list)]
				mut_list.sort()
				j+=1
				if j >= clst_n:
					i=j
					break
	
		if j>= clst_n:
			i+=1

		if len(mut_list) >= MIN_MUTATIONS: # saving merged cluster
			ret_dict = { 'left_position':lpos,
						'right_position':rpos,
							'length': rpos-lpos+1,
						'mut_positions':','.join([str(a) for a in mut_list])}
			merged_clusters.append(ret_dict)
	print('merged_clusters : %d'%len(merged_clusters))

	## write clusters to file
	cl_outf=open('%s/clusters_%s.txt'%(info.outdir,str(tag)), 'w')
	header=merged_clusters[0].keys()
	cl_outf.write('%s\n'%('\t'.join(header)))
	for cluster_info_dict in merged_clusters:		
		cl_outf.write('%s\n'%('\t'.join([str(x) for x in cluster_info_dict.values()])))
	cl_outf.close()

	return merged_clusters

def get_candidate_core_mutations(mutclust_input_df, info, tag, i):
	print('Searching candidate core mutations...'),
	total_mutInfo_list=[]
	ccm_index_list=[]
	total_index_list=[] #추가
	#mutation filtering
	filtered_mutclust_input_df = mutation_filtering(mutclust_input_df)
	for index, pos in enumerate(filtered_mutclust_input_df[POS]):#position
		lr_index, lr_distance, eps_scaler = get_eps_region(filtered_mutclust_input_df, index, info)
		mut_info = get_eps_stats(index, pos, filtered_mutclust_input_df, lr_index, lr_distance, eps_scaler)
		total_mutInfo_list.append(mut_info.to_dict())	# archive eps stats of all mutations (reference data)
        
		total_index_list.append(index) #추가
        
		if mut_info.mut_n < MIN_MUTATIONS:
			continue

		if mut_info.hscore_sum< CCM_MIN_HSCORE_SUM:
			continue
		
		if mut_info.hscore_avr< CCM_MIN_HSCORE_AVR:
			continue

		if mut_info.i_hscore < CCM_MIN_HSCORE:
			continue

		ccm_index_list.append(index)
		#print(str(len(ccm_index_list))+ ' ccms- pos:' + str(mut_info.i) + ' hscore: ' + str(mut_info.i_hscore) + ' hscore_sum: ' + str(mut_info.hscore_sum) + ' length:' + str(mut_info.length) + ' mut_n: ' + str(mut_info.mut_n))
   
	#print_mutatoin_by_gene(filtered_mutclust_input_df.loc[ccm_index_list], info.mutation_info_outdir)
	
# temp
	with open('%s/total_results_%s.tsv'%(info.outdir, tag), 'w') as outf:
		header='\t'.join(total_mutInfo_list[0].keys()) + '\n'
		outf.write(header)
		for index in total_index_list:
			outf.write('%s\n'%('\t'.join([str(x) for x in total_mutInfo_list[index].values()])))
		outf.close()
# temp end

	# write ccm data to file
	with open('%s/ccm_results_%s.tsv'%(info.outdir, tag), 'w') as outf:
		header='\t'.join(total_mutInfo_list[0].keys()) + '\n'
		outf.write(header)
		for index in ccm_index_list:
			outf.write('%s\n'%('\t'.join([str(x) for x in total_mutInfo_list[index].values()])))
		outf.close()

	total_mutInfo_list = np.asarray(total_mutInfo_list)
	return total_mutInfo_list, ccm_index_list
