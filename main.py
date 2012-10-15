# Learning of Canonical networks
# author: Krzysztof Nowak
from collections import Counter

from math import e,log

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

	leak = float(full_pairs[(0,0,0,0,1)]) / param_pairs[(0,0,0,0)]

	#remove leaks
	del full_pairs[(0,0,0,0,1)]
	del full_pairs[(0,0,0,0,0)]
	del param_pairs[(0,0,0,0)]

	return full_pairs, param_pairs, leak

def choose_equation(full, params):
	param_pairs = sorted(params.items(), key=lambda p:p[1], reverse=True)

	checked = [False]*len(param_pairs[0][0])

	chosen = []
	for param, count in param_pairs:
		if any((checked[i] == False and param[i] == True for i in range(len(param)))):
			chosen.append((param, count))
			for i in range(len(param)):
				checked[i] = checked[i] or param[i]
		if all(checked):
			return chosen
	return chosen

def main():
	full, params, leak  = read_data("Network1.txt")
	eq_set = choose_equation(full, params)

	# calculate values on right
	for i in xrange(len(eq_set)):
		#v = log(1.0-float(full[tuple(list(eq_set[i][0]) + [1])]) / params[eq_set[i][0]])
		v = float(full[tuple(list(eq_set[i][0]) + [1])]) / params[eq_set[i][0]]
		eq_set[i] = list(eq_set[i])[:-1] + [1.0 - v]

	N = len(eq_set[0][0])

	#horizontal vectors
	vectors = []
	for i in range(N):
		vect = []
		for eq in eq_set:
			vect.append(eq[0][i])
		vectors.append(vect)
		#for j in range(len(eq_set)):
			#vectors.append([x for x in eq_set[j][0][i]])
	for eq in eq_set:
		print eq

	b = [eq_set[i][1] for i in range(len(eq_set))]

	print "A:"
	for v in vectors:
		print v
	print "b:",
	print b

	A = np.array(vectors)
	B = np.array(b)
	X = [1.0 - e**x for x in np.linalg.solve(A,B)]
	print X



def test():
	a = np.array([[3,1], [1,2]])
	b = np.array([9,8])
	x = np.linalg.solve(a,b)
	print x


if __name__ == "__main__":
	main()
	#test()