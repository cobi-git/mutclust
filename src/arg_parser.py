import argparse
import shutil
from os.path import join, exists, isdir, dirname
from os import mkdir

from Bin.mlib import DIMINISHING_FACTOR, EPSILON, EPSILON_SCALING_FACTOR, MAX_EPS, MIN_CLUSTER_LENGTH, CCM_MIN_PERCENTAGE_SUM


class args_info:
	args = {}
	fin=str()
	ref=str()
	outdir=str()

	# hyper parameters
	eps=int()
	maxeps=int()
	min_persum=int()	# weakest clade mut's per_sum=11.52 (next 26.87)
	eps_scaler_const=float()
	es_control_const=float()	# the eps_scaler diminishing effect (higher the exp_dim the lesser the es penalized)
	min_cluster_length=int()	# minimum length of cluster

def set_env(input=None, reference=None, output=None):
	info = args_info()
	parser = argparse.ArgumentParser(prog="cluster.py")

	# required arguements
	parser.add_argument('-f', '--input_file', type=str, default='/data3/projects/2020_MUTCLUST/Data/Rawdata/COVID19/nucleotide_data/mutclust_input_data.txt', help='mutation frequency data file', required=False)
	parser.add_argument('-r', '--ref', type=str, default='/data3/projects/2020_MUTCLUST/Data/Rawdata/COVID19/nucleotide_data/new_reference.fasta', help='the reference genome', required=False)
	# parameters
	parser.add_argument('-e', '--eps', type=int, default=EPSILON, help='width of window (epsilon)', required=False)
	parser.add_argument('-maxeps', '--maxeps', type=int, default=MAX_EPS, help='maximum eps', required=False)
	parser.add_argument('-minps', '--min_ps', type=int, default=CCM_MIN_PERCENTAGE_SUM, help='minimum per_sum', required=False)
	parser.add_argument('-es', '--eps_scaler_const', type=float, default=EPSILON_SCALING_FACTOR, help='eps scaling factor', required=False)
	parser.add_argument('-exd', '--exp_dim', type=float, default=DIMINISHING_FACTOR, help='cluster expansion es diminishing factor', required=False)
	parser.add_argument('-minl', '--min_clength', type=int, default=MIN_CLUSTER_LENGTH, help='minimum cluster length', required=False)
		
	args = vars(parser.parse_args())
	info.args = args
	if input == None:
		info.fin = args['input_file']
	else:
		info.fin = input

	if output == None:
		info.outdir = args['outdir']
	else:
		info.outdir = output
		info.plot_outdir = output
		info.mutation_info_outdir = output

	if not exists(input):
		print('there is no input file!')
		exit()

	#if isdir(info.outdir):
	#	shutil.rmtree(info.outdir)
	#mkdir(info.outdir)
	
	info.eps = args['eps']
	info.maxeps = args['maxeps']
	info.min_persum = args['min_ps']
	info.maxeps = args['maxeps']
	info.eps_scaler_const = args['eps_scaler_const']
	info.es_control_const = args['exp_dim']
	info.min_cluster_length = args['min_clength']

	return info

