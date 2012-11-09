import numpy as np
import sys
from collections import Counter
from itertools import combinations
from math import e,log

def pt(arr):
	for a in arr:
		print a

def read_data(filename, delimiter):
	"""
	rtype: dict
	return: {'header': ['Smoker','BadDiet','Cancer'], 'data':[ ['True','False','Yes'],...], 'range':{'Smoker':['True','False'],...}}
	"""

	data = tuple([ tuple(line.strip().split(delimiter)) for line in open(filename).readlines() ])
	header, data = data[0], data[1:]
	domain = {}
	for i in range(len(header)):
		domain[header[i]] = tuple( set([ data[j][i] for j in range(len(data)) ]) )

	return {'header': header, 'data':data, 'domain': domain }

def meets_rules(record, rules):
	for i in range(len(record)):
		if not record[i] in rules[i]:
			return False
	return True

def filter_data(data, rules):
	return [data[i] for i in range(len(data)) if meets_rules(data[i], rules)]

def get_positives(data):
	domain = data['domain']
	negatives = data['negatives']
	header = data['header']
	return tuple([tuple(set(domain[header[i]]) - set([negatives[i]])) for i in range(len(negatives))])

def get_leak(data):
	negatives = data['negatives']
	header = data['header']
	domain = data['domain']
	child_node = data['child_node']

	rules = tuple([set([i]) for i in negatives] + [set(domain[header[child_node]])])
	data_total = filter_data(data['data'], rules)

	LEAK = []
	for i in range(len(domain[header[child_node]])):
		rules1 = tuple( [set( [j] ) for j in negatives] + [set( [domain[header[child_node]][i]]) ] )
		data_i = filter_data(data['data'], rules1)
		LEAK.append( (domain[header[child_node]][i], len(data_i)/float(len(data_total)) ) )
	return LEAK

def get_greedy_params(data):
	negatives = data['negatives']
	header = data['header']
	domain = data['domain']
	positives = data['positives']
	child_node = data['child_node']
	counter = data['counter']

	params = []
	for i in range(len(positives)): # for each parameter
		for j in range(len(positives[i])): # for each value
			total = 0
			for k in range(len(domain[header[child_node]])):
				rules_i = tuple(negatives[:i]) + tuple([positives[i][j]]) + tuple(negatives[i+1:]) + tuple([domain[header[child_node]][k]])
				total += counter[rules_i]

			for k in range(len(domain[header[child_node]])):
				rules_i = tuple(negatives[:i]) + tuple([positives[i][j]]) + tuple(negatives[i+1:]) + tuple([domain[header[child_node]][k]])
				value = counter[rules_i]/float(total)
				params.append((header[i], positives[i][j], domain[header[child_node]][k], value ))
				#print rules_i, counter[rules_i], total, counter[rules_i]/float(total)
	return params

def encode_eq(data, eq):
	"""
	@param eq: (('A','B','C'), 333)
	@return: ((0,1,0,1,1,0), 333)
	"""
	positives = data['positives']
	equation_params = eq[0]
	final = []
	for i in range(len(positives)):
		tmp = [0]*(len(positives[i]))
		if equation_params[i] in positives[i]:
			tmp[positives[i].index(equation_params[i])] = 1
		final.extend(tmp)
	return tuple([tuple(final), eq[1]])

def decode_eq_single(data, eq):
	"""
	assumes eq = [0,0, .. ,1, .., 0, 0]
	"""
	new_eq = decode_eq(data, eq)
	positives = data['positives']
	for i in range(len(new_eq)):
		if new_eq[i] in positives[i]:
			return new_eq[i]

def decode_eq_header(data, eq):
	positives = data['positives']
	p_counts = sum([[i]*len(positives[i]) for i in range(len(positives))],[])
	index = eq.index(1)
	return data['header'][p_counts[index]]

def decode_eq(data, eq):
	"""
	@param eq: (0,1,1,1,0,0,0)
	@return: ('A','B','C')
	"""
	positives = data['positives']
	negatives = data['negatives']
	final = []
	p_counts = [len(p) for p in positives]
	suma = 0
	for i in range(len(p_counts)):
		eq_part = eq[suma:suma+p_counts[i]]
		if 1 in eq_part:
			final.append(positives[i][eq_part.index(1)])
		else:
			final.append(negatives[i])
		suma+=p_counts[i]
	return tuple(final)

def equation_valid(encoded_eq):
	vectors = [eq[0] for eq in encoded_eq]
	A = np.array(vectors)
	detA = np.linalg.det(A)
	return detA > 10e-7


def choose_equations(sorted_counts, n, data):
	# remove child params and get equation only
	# TODO: HARDCODED -1 as child node!
	# ( (('A','B'), 33), ...)
	sorted_equations = tuple([(eq[:-1],cnt) for eq,cnt in sorted_counts])

	# ( ((1,0,0,1,1),33), ...)
	encoded_eq = [encode_eq(data, eq) for eq in sorted_equations]

	#for eq in encoded_eq:
		#print eq

	t = n
	chosen = encoded_eq[:n]
	found = False
	while(t <= len(encoded_eq) and not found):
		for comb in combinations(encoded_eq[:t], n):
			if equation_valid(comb):
				chosen = comb
				found = True
				break;
		t+=1
	return chosen

	##weights = [cnt for eq,cnt in sorted_counts]
	##print weights
	##print equations

	#encoded_eq = [encode_eq(data, eq) for eq in equations]
	##for eq in encoded_eq:
	#	#print eq, decode_eq(data,eq[0])

	#vectors = [eq[0] for eq in encoded_eq]
	##print vectors
	#A = np.array(vectors)
	#detA = np.linalg.det(A)
	##print A, detA


def get_smart_params(data):
	negatives = data['negatives']
	header = data['header']
	domain = data['domain']
	positives = data['positives']
	child_node = data['child_node']

	n = sum(( len(d) for d in positives ))

	flat_list = sum(positives,tuple())
	#data['p_codes'] = dict((k,i) for i in range(n)))

	#sorted_counts = sorted(counter.items(), key=lambda i:i[1])
	#for s in sorted_counts:
		#print s
	#print header
	#print positives
	#print sorted_counts

	# make new counter
	counter = Counter(data['data'])

	#delete leak
	for child_param in domain[header[child_node]]:
		assert (tuple(negatives) + tuple([child_param])) in counter, "LEAK combination not found in counter!"
		del counter[tuple(negatives) + tuple([child_param])]

	# CHOOSE FOR MALIGNANT ONLY
#	for child_param in domain[header[child_node]]:
	data_malignant = [d for d in counter.iteritems() if d[0][child_node] == 'Malignant']
	sorted_counts_mali = sorted(data_malignant, key=lambda i:i[1], reverse=True)
	sorted_counts_mali_param = []
	for eq, cnt in sorted_counts_mali:
		params = eq[:-1]
		suma = 0
		for child_param in domain[header[child_node]]:
			suma+= counter[params+tuple([child_param])]
		sorted_counts_mali_param.append(tuple([eq, cnt/float(suma)]))

	#pt(sorted_counts_mali_param)

	#sorted_counts_mali = [(d[0], su) for d in sorted_counts_mali]

	# how many parameters
	n = sum(( len(d) for d in positives ))

	#for s in sorted_counts_mali:
		#print s
	choice = choose_equations(sorted_counts_mali_param, n, data)
	pt(choice)
	#print positives
	vectors = [eq[0] for eq in choice]
	b = [log(1.0-eq[1]) for eq in choice]

	vectors = zip(*vectors)
	A = np.array(vectors)
	detA = np.linalg.det(A)

	Y = []
	for i in range(n):
		p = np.array(vectors[:i] + [b] + vectors[i+1:])
		detP = np.linalg.det(p)
		Y.append(1.0-e**(detP/detA))


	solution = []
	for i in range(len(Y)):
		new_eq = [0]*i+[1]+[0]*(len(Y)-1-i)
		solution.append([decode_eq_header(data,new_eq), decode_eq_single(data,new_eq), Y[i]])
	pt(solution)

	#return tuple(solution)



def main():
	data = read_data(sys.argv[1], ' ')

	# TODO: HARDCODED child_node and negatives!
	data['negatives'] = ['False','False','False','Medium']
	data['positives'] = get_positives(data)
	data['counter'] = Counter(data['data'])
	data['child_node'] = 4

	LEAK  = get_leak(data)
	GREEDY_PARAMS = get_greedy_params(data)
	#for P in PARAMS:
	#	print P
	#for L in LEAK:
	#	print L
	SMART_PARAMS = get_smart_params(data)

def test_meets_rules():
	data = [['a','b','c'], ['a','b','z']]
	rules = (
		set(['a']),
		set(['b']),
		set(['c']),
	)
	print filter_data(data,rules)



if __name__ == "__main__":
	#test_meets_rules()
	main()
