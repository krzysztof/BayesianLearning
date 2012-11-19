# Learning of Canonical networks
# author: Krzysztof Nowak
from collections import Counter

from math import e,log, sqrt


import sys

import numpy as np

def read_data(filename):

	file = open(filename)

	data = [list(d.replace('\r\n','').replace('True','1').replace('False','0').split(' ')) for d in file.readlines()]
	data = data[1:]
	for i in xrange(len(data)):
		for j in xrange(len(data[i])):
			data[i][j] = int(data[i][j])
	data = [tuple(d) for d in data]


	full_pairs =  Counter(data)
	#full_pairs = sorted(full_pairs.items(), key=lambda p:p[1], reverse=True)

	data_params = [d[:-1] for d in data]

	param_pairs = Counter(data_params)
	#param_pairs = sorted(param_pairs.items(), key=lambda p:p[1], reverse=True)

	# amount of parameters (without children)
	param_size = len(data_params[0])

	leak_true = tuple([0]*param_size + [1]) # e.g. 0 0 0 1
	leak_false= tuple([0]*(param_size+1)) # e.g 0 0 0 0
	leak_denom = tuple([0]*param_size) #e.g. 0 0 0

	leak = float(full_pairs[leak_true]) / param_pairs[leak_denom]

	#remove leaks
	del full_pairs[leak_true]
	del full_pairs[leak_false]
	del param_pairs[leak_denom]

	return full_pairs, param_pairs, leak

def choose_greedy(full, params):
	param_pairs = sorted(params.items(), key=lambda p:p[1], reverse=True)

	checked = [False]*len(param_pairs[0][0])

	chosen = []
	N = range(len(param_pairs[0][0]))
	n = 0
	for param, count in param_pairs:
		if any((checked[i] == False and param[i] == True for i in range(len(param)))):
			chosen.append((param, count))
			for i in range(len(param)):
				checked[i] = checked[i] or param[i]
			n+=1
		if n==N:
			return chosen

	return chosen

def choose_exact(full, params):
	param_pairs = params.items()

	from itertools import combinations

	all_combinations = list(combinations(range(len(param_pairs)), 4))
	final_combinations = []
	for comb in all_combinations:
		vectors = [param_pairs[i][0] for i in comb]
		A = np.array(vectors)
		A = np.transpose(A)
		detA = np.linalg.det(A)
		if detA > 10e-7:
			final_combinations.append([comb,sum((param_pairs[i][1] for i in comb))])
	final_combinations = sorted(final_combinations, key=lambda p:p[1], reverse=True)

	chosen = [param_pairs[i] for i in final_combinations[0][0]]
	#print chosen

	return chosen


def choose_for_each(full, params, p_num):
	"""
	p_num - parameter index that's to be optimized
	"""
	param_pairs = sorted(params.items(), key=lambda p:p[1], reverse=True)
	for k in param_pairs:
		print k
	print
	print
	print
	elite = filter(lambda item: item[0][p_num]==1, param_pairs)
	rest = filter(lambda item: item[0][p_num]==0, param_pairs)
	param_pairs = elite+rest
	for k in param_pairs:
		print k

	checked = [False]*len(param_pairs[0][0])

	chosen = []
	N = range(len(param_pairs[0][0]))
	n = 0
	for param, count in param_pairs:
		if any((checked[i] == False and param[i] == True for i in range(len(param)))):
			chosen.append([param, count])
			for i in range(len(param)):
				checked[i] = checked[i] or param[i]
			n+=1
		if n==N:
			return chosen

	return chosen

# greedy choice -> maximize amount of records in a single equation set
#choose_equation = choose_for_each

def choose_lazy(full, params):
	N = len(params.items()[0][0])
	chosen = []
	v = [1]+[0]*(N-1)
	for i in range(N):
		chosen.append((tuple(v), params[tuple(v)]))
		v = [v[-1]]+v[:-1]
	return chosen

def apply_leak(equation_set, leak):
	eq_set = list(equation_set)
	for i in range(len(eq_set)):
		eq = eq_set[i][0]
		v = eq_set[i][1]
		N = sum(eq)
		eq_set[i][1] = 1 - (1-v)/((1-leak)**(N-1))
	return eq_set

def prepare_constants(full, params, eq_set_orig):
	eq_set = list(eq_set_orig)
	for i in xrange(len(eq_set)):
		#v = log(1.0-float(full[tuple(list(eq_set[i][0]) + [1])]) / params[eq_set[i][0]])
		v = float(full[tuple(list(eq_set[i][0]) + [1])]) / params[eq_set[i][0]]
		eq_set[i] = list(eq_set[i])[:-1] + [1.0 - v]
	return eq_set

def solve(eq_set_orig, leak):
	eq_set = list(eq_set_orig)

	for i in xrange(len(eq_set)):
		eq_set[i][1] = log(eq_set[i][1])

#	for eq in eq_set:
#		print eq

	#horizontal vectors
	vectors = []
	for i in range(N):
		vect = []
		for eq in eq_set:
			vect.append(eq[0][i])
		vectors.append(vect)
		#for j in range(len(eq_set)):
			#vectors.append([x for x in eq_set[j][0][i]])

	b = [eq_set[i][1] for i in range(len(eq_set))]

#	print "A:"
#	for v in vectors:
#		print v
#	print "b:",
#	print b

	A = np.array(vectors)
	detA = np.linalg.det(A)
	solution = []
	for i in range(N):
		p = np.array(vectors[:i] + [b] + vectors[i+1:])
		detP = np.linalg.det(p)
		solution.append(1.0-e**(detP/detA))
	#print "p1:",(detp1/detA)
	#print 'detA', np.linalg.det(A)
	#print p1
	#print 'detAx', np.linalg.det(p1)
	return tuple(solution)

	#B = np.array(b)
	#X = np.linalg.solve(A, B)
	#print X
	##X = [1.0 - e**x for x in np.linalg.solve(A,B)]
	#print X


def vect_dist(a,b):
	return tuple((a[i]-b[i] for i in range(len(a))))

def euclidian_dist(a,b):
	assert len(a) == len(b)
	sum = 0.0;
	for i in range(len(a)):
		sum+= (a[i] - b[i])**2
	return sqrt(sum)

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
	assert all((Q[i]>=0 for i in xrange(len(Q))))
	suma = sum(( (sqrt(P[i]) - sqrt(Q[i]))**2 for i in range(len(P)) ))
	return sqrt(suma) * 0.7071067811865475  # * 1./sqrt(2)

def main():
	#Network1 - cancer (5 nodes)
	ORIG = (0.61, 0.25, 0.15, 0.04)

	global N
	N = len(ORIG)
	full, params, leak  = read_data(sys.argv[1])
	eq_set = choose_greedy(full, params)
	eq_set = prepare_constants(full, params, eq_set)
	eq_set = apply_leak(eq_set, leak)
	X = solve(eq_set, leak)

	eq_set_lazy = choose_lazy(full, params)
	eq_set_lazy = prepare_constants(full, params, eq_set_lazy)
	X_lazy = solve(eq_set_lazy, leak)

	print ORIG
	print 'new_v:', X
	print 'old_v:', X_lazy
	#print 'new_dist', tuple((abs(x) for x in vect_dist(ORIG,X)))
	#print 'old_dist', tuple((abs(x) for x in vect_dist(ORIG,X_lazy)))
	print 'new_MAX', max(tuple((abs(x) for x in vect_dist(ORIG,X))))
	print 'new_SUM', sum(tuple((abs(x) for x in vect_dist(ORIG,X))))
	print 'new_AVG', sum(tuple((abs(x) for x in vect_dist(ORIG,X))))/len(X)
	print 'old_MAX', max(tuple((abs(x) for x in vect_dist(ORIG,X_lazy))))
	print 'old_SUM', sum(tuple((abs(x) for x in vect_dist(ORIG,X_lazy))))
	print 'old_AVG', sum(tuple((abs(x) for x in vect_dist(ORIG,X_lazy))))/len(X)
	print 'new_euclidian_dist', euclidian_dist(ORIG,X)
	print 'old_euclidian_dist', euclidian_dist(ORIG,X_lazy)

	print 'new_kl_dist', kl_dist(ORIG,X)
	print 'old_kl_dist', kl_dist(ORIG,X_lazy)

	print 'new_hellinger_dist', hellinger_dist(ORIG,X)
	print 'old_hellinger_dist', hellinger_dist(ORIG,X_lazy)

def main2():
	ORIG = (0.61, 0.25, 0.15, 0.04)

	global N
	N = len(ORIG)
	full, params, leak  = read_data(sys.argv[1])

	eq_set_lazy = choose_lazy(full, params)
	eq_set_lazy = prepare_constants(full, params, eq_set_lazy)
	X_lazy = solve(eq_set_lazy, leak)

	X = []
	for k in range(N):
		eq_set = choose_for_each(full, params, k)
		print "K: %s, set:%s"%(k,eq_set)

		#X.append()
		#eq_set = prepare_constants(full, params, eq_set)
		#eq_set = apply_leak(eq_set, leak)
		if len(eq_set)==4:
			Xk = solve(eq_set, leak)
			print Xk
			#X.append[Xk[k]]

	print ORIG
	print 'new_v:', X
	print 'old_v:', X_lazy
	#print 'new_dist', tuple((abs(x) for x in vect_dist(ORIG,X)))
	#print 'old_dist', tuple((abs(x) for x in vect_dist(ORIG,X_lazy)))
	print 'new_MAX', max(tuple((abs(x) for x in vect_dist(ORIG,X))))
	print 'new_SUM', sum(tuple((abs(x) for x in vect_dist(ORIG,X))))
	print 'new_AVG', sum(tuple((abs(x) for x in vect_dist(ORIG,X))))/len(X)
	print 'old_MAX', max(tuple((abs(x) for x in vect_dist(ORIG,X_lazy))))
	print 'old_SUM', sum(tuple((abs(x) for x in vect_dist(ORIG,X_lazy))))
	print 'old_AVG', sum(tuple((abs(x) for x in vect_dist(ORIG,X_lazy))))/len(X)
	print 'new_euclidian_dist', euclidian_dist(ORIG,X)
	print 'old_euclidian_dist', euclidian_dist(ORIG,X_lazy)

	print 'new_kl_dist', kl_dist(ORIG,X)
	print 'old_kl_dist', kl_dist(ORIG,X_lazy)

	print 'new_hellinger_dist', hellinger_dist(ORIG,X)
	print 'old_hellinger_dist', hellinger_dist(ORIG,X_lazy)

def method3():
	ORIG = (0.61, 0.25, 0.15, 0.04)

	global N
	N = len(ORIG)
	full, params, leak  = read_data(sys.argv[1])
	eq_set = choose_exact(full, params)
	eq_set = prepare_constants(full, params, eq_set)
	eq_set = apply_leak(eq_set, leak)
	X = solve(eq_set, leak)

	eq_set_lazy = choose_lazy(full, params)
	eq_set_lazy = prepare_constants(full, params, eq_set_lazy)
	X_lazy = solve(eq_set_lazy, leak)

	print ORIG
	print 'new_v:', X
	print 'old_v:', X_lazy
	#print 'new_dist', tuple((abs(x) for x in vect_dist(ORIG,X)))
	#print 'old_dist', tuple((abs(x) for x in vect_dist(ORIG,X_lazy)))
	print 'new_MAX', max(tuple((abs(x) for x in vect_dist(ORIG,X))))
	print 'new_SUM', sum(tuple((abs(x) for x in vect_dist(ORIG,X))))
	print 'new_AVG', sum(tuple((abs(x) for x in vect_dist(ORIG,X))))/len(X)
	print 'old_MAX', max(tuple((abs(x) for x in vect_dist(ORIG,X_lazy))))
	print 'old_SUM', sum(tuple((abs(x) for x in vect_dist(ORIG,X_lazy))))
	print 'old_AVG', sum(tuple((abs(x) for x in vect_dist(ORIG,X_lazy))))/len(X)
	print 'new_euclidian_dist', euclidian_dist(ORIG,X)
	print 'old_euclidian_dist', euclidian_dist(ORIG,X_lazy)

	print 'new_kl_dist', kl_dist(ORIG,X)
	print 'old_kl_dist', kl_dist(ORIG,X_lazy)

	print 'new_hellinger_dist', hellinger_dist(ORIG,X)
	print 'old_hellinger_dist', hellinger_dist(ORIG,X_lazy)


def test():
	full, params, leak  = read_data(sys.argv[1])
	#print leak

	param_pairs = sorted(params.items(), key=lambda p:p[1], reverse=True)
	for k,v in param_pairs:
		print k,v

if __name__ == "__main__":
	#main2()
	main()
	#method3()
	#test()