import numpy as np
import sympy as sp
import sys
import os
import random
from math import log, sqrt
from collections import Counter, defaultdict
from itertools import combinations
from math import e,log
from boxplot import make_boxplot
from GaussJordan import GaussJordanElimination, augment, PRODUCT_OPERATIONS
from fractions import Fraction
import matplotlib.pyplot as plt
from utilities import max_based, LEAKDEF,CPT

import subprocess as sub


def max_based2(or_cpt):
	leak = or_cpt[-1]
	generator = (list(int(j) for j in bin(i)[2:].zfill(4)) for i in range(15,-1,-1))
	def fun(g, or_cpt, leak):
		name = " ".join([("True" if gi==1 else "False") for gi in g]) + " True"
		leak_exp = sum(g)-1
		value =  1.0 - reduce(lambda x,y: x*y, [1.0 - a*b for a,b in zip(g, or_cpt)][:4])/((1.0 - leak)**leak_exp)
		#TODO: FIX - THIS LEADS TO NEGATIVE VALUES IN HELLINGER DISTANCE
		#value = 1.0 - reduce(lambda x,y: x*y, [1.0 - a*b for a,b in zip(g, or_cpt)][:4])/(1.0 - leak)
		#TODO: ugly fix, doesn't count LEAK at all! yet hellinger works
		#value =  1.0 - reduce(lambda x,y: x*y, [1.0 - a*b for a,b in zip(g, or_cpt)][:4])
		return "%s %f" % (name, value)
		#return (name, value)
	return "\n".join(list(fun(g, or_cpt, leak) for g in generator))

def main():
	def eucl_dist(A, B):
		assert not any(((a < 0.0 or a > 1.0) for a in A)), "Error in set A: %s"%(A)
		assert not any(((b < 0.0 or b > 1.0) for b in B)), "Error in set B: %s"%(B)
		assert len(A) == len(B), "sets differ in size %d vs %d" % (len(A), len(B))

		return sqrt(sum((A[i] - B[i])**2 for i in range(len(A))))

	def kl_dist(P, Q):
		assert len(P) == len(Q)
		suma = 0.0;
		for i in range(len(P)):
	
			q = Q[i]
			if Q[i] < 10e-7:
				q = 10e-7
	
			suma += P[i]*log((P[i]/q),e);
		return suma

	def hellinger_dist(P, Q):
		assert len(P) == len(Q)
		assert all((Q[i]>=0 for i in xrange(len(Q)))), "%s"%(Q)
		suma = sum(( (sqrt(P[i]) - sqrt(Q[i]))**2 for i in range(len(P)) ))
		return sqrt(suma) * 0.7071067811865475  # * 1./sqrt(2)

# OLD PREGENERATED DATA
#	filenames = [ "CancerOR_1k.txt", "CancerOR_2k.txt", "CancerOR_3k.txt", "CancerOR_4k.txt", "CancerOR_1100.txt", "CancerOR_1200.txt", "CancerOR_1300.txt", "CancerOR_1400.txt", "CancerOR_1500.txt", "CancerOR_1600.txt",
#		"CancerOR_1800.txt", "CancerOR_1900.txt", "CancerOR_800.txt", "CancerOR_2500.txt", "CancerOR_2400.txt", "CancerOR_2300.txt", "CancerOR_2200.txt", "CancerOR_2100.txt", "CancerOR_3100.txt",
#		"CancerOR_3300.txt", "CancerOR_3500.txt", "CancerOR_3700.txt", "CancerOR_3900.txt"] 
#	nums = [1000, 2000, 3000, 4000, 1100, 1200, 1300, 1400, 1500, 1600, 1800, 1900, 800, 2500, 2400, 2300, 2200, 2100, 3100, 3300, 3500, 3700, 3900]


	#dist_func = hellinger_dist
	dist_func = eucl_dist

	cpp_generator = "./NoisyMAXSmile/generator"
	cpp_learner = "./NoisyMAXSmile/smile_learner"
	py_learner = "./noisyMAX.py"
	naive_learner = "./naiveMAX.py"
	def exec_command(command):
		p = os.popen(command, "r")
		output = ""
		while 1:
			line = p.readline()
			if not line:
				break
			output += line
		return output[:-1]
	
	def makedict(text):
		r_dict = {}
		for line in text.split('\n')[1:]:
			cut = line.split(' ')
			val = float(cut[-1])
			case = " ".join(cut[:-1])
			r_dict[case] = val
		return r_dict
	
	def distance(A,B, dist_f):
		valA = []
		valB = []
		for key, value in A.iteritems():
			if key.endswith("True"):
				valA.append(value)
				valB.append(B[key])
		return dist_f(valA, valB)
				

	#print makedict(max_based2(real_or))

	#real_or = [ [0.61, 0, 0.25, 0, 0.15, 0, 0.04, 0, 0.01,],
	#		[0.39, 1, 0.75, 0, 0.85, 0, 0.96, 0, 0.99,], ]
	#real_or_dim = 5
	#network = "./src/CancerOR.xdsl"
	#network_cpt = "./src/CancerOR_CPT.xdsl"

	real_or = [ [0.61, 0, 0.25, 0, 0.15, 0, 0.04, 0, 0.35, 0, 0.89, 0, 0.01,],
			[0.39, 1, 0.75, 1, 0.85, 1, 0.96, 1, 0.65, 1, 0.11, 1, 0.99,], ]
	real_or_dim = 7
	network = "./src/OR_6.xdsl"
	network_cpt = "./src/OR_6_CPT.xdsl"


	print real_or
	parent_dims = [2]*(real_or_dim - 1)
	#labels = ["Smoker","Genetic","CoalWorker","BadDiet","LungCancer"]
	labels = ["P%d"%(i) for i in range(1,7)] + ['C1',]
	real_or_cpt = CPT(real_or, parent_dims=parent_dims, network_type=CPT.TYPE_NOISY_MAX, labels=labels, states = [["True","False"]]*real_or_dim).max_based()
	#max_based_real = max_based(real_or, leakdef = LEAKDEF.HENRION)
	#print max_based_real
	#real_out = dict(real_or_cpt.print_raw())
	real_out = makedict(real_or_cpt.print_raw())
	#print real_out

	#sety = zip(filenames, nums)
	#dataset_file = "data/"+filename

	#nums = range(800,5001,500)
	nums = [1000, 3000, 10000, 50000, 100000]
	#nums = [100, 1000, 10000, 100000]
	#nums = [2000, 2500]
	#nums = [500, 1000 , 1500, 2000, 2500]
	labels = [ "%d records"%(num, )  for num in nums]
	#labels = ["1k records", "3k records", "10k records", "50k records", "100k records"]
	#boxNames = ["Gauss-Jordan elimination", "SMILE Noisy-MAX fitting"]
	#boxColors = ['darkkhaki', 'forestgreen']

	reps = 100
	#sety = xrange(1000,5000,50)
	sety = zip("A"*len(nums), nums)
	dataset_file = "tmp/tmp_data.txt"

	#network_generator_seed = 3333  # used for 1k-100k GJ vs naive
	#network_generator_seed = 123456  # used for 500 - 2500 GJ vs naive
	#network_generator_seed = 44444  # used for 500 - 2500 GJ vs smile
	#network_generator_seed = 444445  # used for 1k -100k GJ vs GJFit
	#network_generator_seed = 5555  # used for 500-2500 all
	#network_generator_seed = 3355  # used for 1k-100k all
	#network_generator_seed = 3357  # used for 1k-100k all
	network_generator_seed = 123412333

	rnd = random.Random()
	rnd.seed(network_generator_seed)
	data = []
	for filename, i in sety:
		subdata_gj = []
		subdata_gjf = []
		subdata_genie = []
		subdata_naive = []
		for t in xrange(reps):
			seed = rnd.randint(1,10**6)
			gen_command = "%s %d %s %d > %s" % (cpp_generator, i, network, seed, dataset_file)
			#learn_command_em = "%s %s %s EM" %(cpp_learner, dataset_file, network_cpt)
			learn_command_genie = "%s %s %s Smile" %(cpp_learner, dataset_file, network_cpt)
			learn_command_GJ = "python %s %s GJ" %(py_learner, dataset_file)
			learn_command_GJFit = "python %s %s GJFit" %(py_learner, dataset_file)
			learn_command_naive = "python %s %s" %(naive_learner, dataset_file)


			print "s = %d, t = %d"%(i,t,)
			exec_command(gen_command)
			#em_out = makedict(exec_command(learn_command_em))

			GJ_out = makedict(exec_command(learn_command_GJ))
			#print real_out, GJ_out
			d_r_gj =  distance(real_out, GJ_out, dist_func)
			print "GJ", d_r_gj

			GJFit_out = makedict(exec_command(learn_command_GJFit))
			d_r_gjf =  distance(real_out, GJFit_out, dist_func)
			print "GJFit", d_r_gjf

			print learn_command_genie
			command_out = exec_command(learn_command_genie)
			genie_out = makedict(command_out)
			d_r_g =  distance(real_out, genie_out, dist_func)
			print "Genie", d_r_g

			naive_out = makedict(exec_command(learn_command_naive))
			d_r_naive =  distance(real_out, naive_out, dist_func)
			print "Naive", d_r_naive

			subdata_gj.append(d_r_gj)
			subdata_gjf.append(d_r_gjf)
			subdata_genie.append(d_r_g)
			subdata_naive.append(d_r_naive)
			
		data.append(subdata_gj)
		data.append(subdata_gjf)
		data.append(subdata_genie)
		data.append(subdata_naive)

	boxNames = ["Gauss-Jordan", "Gauss-Jordan (confidence variant)", "Genie", "Naive"]
	boxColors = ['darkkhaki', 'indianred', 'forestgreen', 'royalblue' ]

	make_boxplot(data, labels, boxColors, boxNames, title="Euclidian distance for network %s"%(network,))

if __name__ == "__main__":
	main()

