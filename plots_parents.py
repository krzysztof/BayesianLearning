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

import subprocess as sub


def max_based(max_cpt):
	leak = max_cpt[-1]
	generator = (list(int(j) for j in bin(i)[2:].zfill(4)) for i in range(15,-1,-1))
	def fun(g, max_cpt, leak):
		name = " ".join([("True" if gi==1 else "False") for gi in g]) + " True"
		leak_exp = sum(g)-1
		#value =  1.0 - reduce(lambda x,y: x*y, [1.0 - a*b for a,b in zip(g, max_cpt)][:4])/((1.0 - leak)**leak_exp)
		#TODO: FIX - THIS LEADS TO NEGATIVE VALUES IN HELLINGER DISTANCE
		#value = 1.0 - reduce(lambda x,y: x*y, [1.0 - a*b for a,b in zip(g, max_cpt)][:4])/(1.0 - leak)
		#TODO: ugly fix, doesn't count LEAK at all! yet hellinger works
		value =  1.0 - reduce(lambda x,y: x*y, [1.0 - a*b for a,b in zip(g, max_cpt)][:4])
		return value
		#return (name, value)
	return list(fun(g, max_cpt, leak) for g in generator)

def max_based2(max_cpt):
	leak = max_cpt[-1]
	generator = (list(int(j) for j in bin(i)[2:].zfill(4)) for i in range(15,-1,-1))
	def fun(g, max_cpt, leak):
		name = " ".join([("True" if gi==1 else "False") for gi in g]) + " True"
		leak_exp = sum(g)-1
		#value =  1.0 - reduce(lambda x,y: x*y, [1.0 - a*b for a,b in zip(g, max_cpt)][:4])/((1.0 - leak)**leak_exp)
		#TODO: FIX - THIS LEADS TO NEGATIVE VALUES IN HELLINGER DISTANCE
		value = 1.0 - reduce(lambda x,y: x*y, [1.0 - a*b for a,b in zip(g, max_cpt)][:4])/(1.0 - leak)
		#TODO: ugly fix, doesn't count LEAK at all! yet hellinger works
		#value =  1.0 - reduce(lambda x,y: x*y, [1.0 - a*b for a,b in zip(g, max_cpt)][:4])
		return "%s %f" % (name, value)
		#return (name, value)
	return "\n".join(list(fun(g, max_cpt, leak) for g in generator))

def main():
	def eucl_dist(A, B):
		return sqrt(sum((A[i] - B[i])**2 for i in range(len(A)) ))

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

	filenames = [ "CancerOR_1k.txt", "CancerOR_2k.txt", "CancerOR_3k.txt", "CancerOR_4k.txt", "CancerOR_1100.txt", "CancerOR_1200.txt", "CancerOR_1300.txt", "CancerOR_1400.txt", "CancerOR_1500.txt", "CancerOR_1600.txt",
"CancerOR_1800.txt", "CancerOR_1900.txt", "CancerOR_800.txt", "CancerOR_2500.txt", "CancerOR_2400.txt", "CancerOR_2300.txt", "CancerOR_2200.txt", "CancerOR_2100.txt", "CancerOR_3100.txt",
"CancerOR_3300.txt", "CancerOR_3500.txt", "CancerOR_3700.txt", "CancerOR_3900.txt"] 
	nums = [1000, 2000, 3000, 4000, 1100, 1200, 1300, 1400, 1500, 1600, 1800, 1900, 800, 2500, 2400, 2300, 2200, 2100, 3100, 3300, 3500, 3700, 3900]

	nums = range(800, 5001, 100) #+ [5000, 10000, 100000]
	filenames = [ "CancerOR_%d.txt"%(n) for n in nums]

	real_or = (0.61, 0.25, 0.15, 0.04, 0.01,)
	X_plot = []
	Y_plot = []

	#dist_func = hellinger_dist
	dist_func = eucl_dist

	cpp_generator = "./NoisyMAXSmile/generator"
	cpp_learner = "./NoisyMAXSmile/smile_learner"
	py_learner = "./noisyMAX.py"
	naive_learner = "./naiveMAX.py"
	#network = "./src/CancerOR.xdsl"
	#network_cpt = "./src/CancerOR_CPT.xdsl"

	network = "./src/CancerOR2.xdsl"
	network_cpt = "./src/CancerOR2_CPT.xdsl"
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
		for line in text.split('\n'):
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
	real_out = makedict(max_based2(real_or))


	#sety = zip(filenames, nums)
	#dataset_file = "data/"+filename

	#nums = range(800,5001,500)
	nums = [1000, 3000, 10000, 50000, 100000]
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
	network_generator_seed = 3355  # used for 1k-100k all

	rnd = random.Random()
	rnd.seed(network_generator_seed)
	data = []
	for filename, i in sety:
		subdata_gj = []
		subdata_gjf = []
		subdata_smile = []
		subdata_naive = []
		for t in xrange(reps):
			seed = rnd.randint(1,10**6)
			gen_command = "%s %d %s %d > %s" % (cpp_generator, i, network, seed, dataset_file)
			#learn_command_em = "%s %s %s EM" %(cpp_learner, dataset_file, network_cpt)
			learn_command_smile = "%s %s %s Smile" %(cpp_learner, dataset_file, network_cpt)
			learn_command_GJ = "python %s %s GJ" %(py_learner, dataset_file)
			learn_command_GJFit = "python %s %s GJFit" %(py_learner, dataset_file)
			learn_command_naive = "python %s %s" %(naive_learner, dataset_file)


			print "s = %d, t = %d"%(i,t,)
			exec_command(gen_command)
			#em_out = makedict(exec_command(learn_command_em))

			GJ_out = makedict(exec_command(learn_command_GJ))
			d_r_gj =  distance(real_out, GJ_out, dist_func)
			print "GJ", d_r_gj

			GJFit_out = makedict(exec_command(learn_command_GJFit))
			d_r_gjf =  distance(real_out, GJFit_out, dist_func)
			print "GJFit", d_r_gjf

			smile_out = makedict(exec_command(learn_command_smile))
			d_r_g =  distance(real_out, smile_out, dist_func)
			print "SMILE", d_r_g

			naive_out = makedict(exec_command(learn_command_naive))
			d_r_naive =  distance(real_out, naive_out, dist_func)
			print "Naive", d_r_naive

			subdata_gj.append(d_r_gj)
			subdata_gjf.append(d_r_gjf)
			subdata_smile.append(d_r_g)
			subdata_naive.append(d_r_naive)
			
		data.append(subdata_gj)
		data.append(subdata_gjf)
		data.append(subdata_smile)
		data.append(subdata_naive)

	boxNames = ["Gauss-Jordan", "Gauss-Jordan (confidence variant)", "SMILE", "Naive"]
	boxColors = ['darkkhaki', 'indianred', 'forestgreen', 'royalblue' ]

	make_boxplot(data, labels, boxColors, boxNames )

		#X_plot.append(i)
		#Y_plot.append((
		#	("De_orig_hill", d_r_g), 
		#	("De_orig_gj", d_r_gj),
		#	("De_orig_gj_fit", d_r_gjf),
		#	("De_orig_naive", d_r_naive),
		#	#("De_orig_hill", dist_func(real_cpt, CancerOR_hill_cpt)), 
		#	#("De_orig_gj", dist_func(real_cpt, GJElim_cpt)),
		#	#("De_orig_gj_fit", dist_func(real_cpt, GJElim_fit_cpt)),
		#))
	#plot_data = sorted(zip(X_plot, Y_plot), key=lambda x:x[0])

	#colors = (
	#		(0,'r'),
	#		#(1,'g'),
	#		#(2,'b'),
	#		(3,'y'),
	#		)
	#for n, c in colors:
	#	Y = [y[n][1] for x, y in plot_data] 
	#	X = [x for x, y in plot_data] 
	#	plt.plot(X, Y, color=c, marker="o")

	#	#plt.title("Hellinger distance (R)")

	#plt.title("Euclidian distance")
	#plt.xlabel("Number or records")
	#plt.ylabel("Distance")
	#plt.show()

		
	#for k, v in sorted(GenieTests.items(), key=lambda i:i[1]['num']):
	#	GJElim_distribution, GJElim_fit_distribution = mainGJ(k)

	#	CancerOR_em = v["em"]
	#	CancerOR_hill = v["hill"]

	#	CancerOR_hill_cpt = max_based(CancerOR_hill)
	#	GJElim_cpt = max_based(GJElim_distribution)
	#	GJElim_fit_cpt = max_based(GJElim_fit_distribution)
	#	real_cpt = max_based(real_or)
	#	

	#	print (":"*10)+"Running test for %s"%(k) + (":"*10)
	#	print "real_OR", real_or
	#	print "GJ_OR", GJElim_distribution
	#	print "GJFit_OR", GJElim_fit_distribution
	#	print "hill_OR:", CancerOR_hill
	#	print ""

	#	print "CPT(real_OR -> cpt)%s\n"%(real_cpt,)
	#	print "CPT(EM from data)%s\n"%(CancerOR_em,)
	#	print "CPT(hill -> cpt)%s\n"%(CancerOR_hill_cpt,)
	#	print "CPT(GJ -> cpt) %s\n"%(GJElim_cpt,)
	#	print "CPT(GJFit -> cpt) %s\n"%(GJElim_fit_cpt,)


	#	print "De(EM vs max_based(real_OR))", dist_func(CancerOR_em, real_cpt,)
	#	print "De(EM vs max_based(hill_OR))", dist_func(CancerOR_em, CancerOR_hill_cpt)
	#	print "De(EM vs max_based(GJ_OR))", dist_func(CancerOR_em, GJElim_cpt)
	#	print "De(EM vs max_based(GJFit_OR))", dist_func(CancerOR_em, GJElim_fit_cpt)

	#	print "De(max_based(orig_OR) vs max_based(hill_OR))", dist_func(real_cpt, CancerOR_hill_cpt)
	#	print "De(max_based(orig_OR) vs max_based(GJ_OR))", dist_func(real_cpt, GJElim_cpt)
	#	print "De(max_based(orig_OR) vs max_based(GJFit_OR))", dist_func(real_cpt, GJElim_fit_cpt)
	#	print "De(orig_OR vs hill_OR):", dist_func(real_or, CancerOR_hill)
	#	print "De(orig_OR vs GJ_OR):", dist_func(real_or, GJElim_distribution)
	#	print "De(orig_OR vs GJFit_OR):", dist_func(real_or, GJElim_fit_distribution)

	#	X_plot.append(v["num"])
	#	Y_plot.append((
	#		("De_orig_hill", dist_func(real_cpt, CancerOR_hill_cpt)), 
	#		("De_orig_gj", dist_func(real_cpt, GJElim_cpt)),
	#		("De_orig_gj_fit", dist_func(real_cpt, GJElim_fit_cpt)),
	#	))
	#plot_data = sorted(zip(X_plot, Y_plot), key=lambda x:x[0])

	#colors = ((0,'r'),(1,'g'),(2,'b'))
	#for n, c in colors:
	#	Y = [y[n][1] for x, y in plot_data] 
	#	X = [x for x, y in plot_data] 
	#	plt.plot(X, Y, color=c, marker="o")

	#	#plt.title("Hellinger distance (R)")
	#	plt.title("Hellinger distance (S)")
	#	plt.xlabel("Number or records")
	#	plt.ylabel("Distance")
	#plt.show()

if __name__ == "__main__":
	main()

