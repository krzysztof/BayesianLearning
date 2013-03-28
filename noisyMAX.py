import numpy as np
import sympy as sp
import sys
import random
from math import log
from collections import Counter, defaultdict
from itertools import combinations
from math import e,log
from GaussJordan import GaussJordanElimination, augment, PRODUCT_OPERATIONS
from fractions import Fraction

def pt(arr):
	for a in arr:
		print a

def read_data(filename, delimiter):
	"""
	rtype: dict
	return: {'header': ['Smoker','BadDiet','Cancer'], 'data':[ ['True','False','Yes'],...], 'domain':{'Smoker':['True','False'],...}}
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
	#print negatives
	return tuple([tuple(set(domain[header[i]]) - set([negatives[i]])) for i in range(len(negatives))])

def get_leak(data):
	negatives = data['negatives']
	header = data['header']
	domain = data['domain']
	child_node = data['child_node']

	rules = tuple([set([i]) for i in negatives] + [set(domain[header[child_node]])])
	data_total = filter_data(data['data'], rules)

	LEAK = {}
	for i in range(len(domain[header[child_node]])):
		rules1 = tuple( [set( [j] ) for j in negatives] + [set( [domain[header[child_node]][i]]) ] )
		data_i = filter_data(data['data'], rules1)
		LEAK[domain[header[child_node]][i]] = len(data_i)/float(len(data_total))
	return LEAK

def get_greedy_params(data, LEAK):
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
				params.append(( domain[header[child_node]][k], header[i], positives[i][j], value ))
				#print rules_i, counter[rules_i], total, counter[rules_i]/float(total)
	return [p for p in params if p[0]!=data['negative_child']]

def encode_eq(data, eq):
	"""
	@param eq: (('A','B','C'), 333)
	@return: ((0, 1, 0, 1, 1, 0), 333)
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
	@param eq: list of type [0,0, .. ,1, .., 0, 0]
	@type eq: list
	@return: named parameter 'True'
	@rtype: string
	"""
	new_eq = decode_eq(data, eq)
	positives = data['positives']
	for i in range(len(new_eq)):
		if new_eq[i] in positives[i]:
			return new_eq[i]

def decode_eq_header(data, eq):
	"""
	@param eq: list of type [0,0, .. ,1, .., 0, 0]
	@type eq: list
	@return: header name e.g. 'Smoker'
	@rtype: string
	"""
	positives = data['positives']
	p_counts = sum([[i]*len(positives[i]) for i in range(len(positives))],[])
	index = eq.index(1)
	return data['header'][p_counts[index]]

def decode_eq(eq, data):
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

def choose_equations(equations, n, data):
	# Assuming that equations are ordered by preference
	#
	# ( (('A','B'), 33), ...)
	#sorted_equations = tuple([(eq[:-1],cnt) for eq,cnt in sorted_counts])

	# ( ((1,0,0,1,1),33), ...)
	encoded_eq = [encode_eq(data, eq) for eq in equations]

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


def get_smart_params(data, LEAK):
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
	# TODO: HARDCODED - assuming that -1 is the only child, and is last in the column)
	counter = Counter([d[:-1] for d in data['data']])

	sorted_counts = sorted(counter.items(), key=lambda i:i[1], reverse=True)
	eq_choice = choose_equations(sorted_counts, n, data)


	#delete leak
	#for child_param in domain[header[child_node]]:
		#assert (tuple(negatives) + tuple([child_param])) in counter, "LEAK combination not found in counter!"
		#del counter[tuple(negatives) + tuple([child_param])]

	final_solution = []
	for child_param_single in domain[header[child_node]]:
		if child_param_single == data['negative_child']:
			continue
		leak = LEAK[child_param_single]

		data_single_child_param = [d for d in counter.iteritems() if d[0][child_node] == child_param_single]
		sorted_counts_single = sorted(data_single_child_param, key=lambda i:i[1], reverse=True)
		sorted_counts_single_param = []
		for eq, cnt in sorted_counts_single:
			params = eq[:-1]
			suma = 0
			for child_param in domain[header[child_node]]:
				suma+= counter[params+tuple([child_param])]
			sorted_counts_single_param.append(tuple([eq, cnt/float(suma)]))

		#pt(sorted_counts_mali_param)

		#sorted_counts_mali = [(d[0], su) for d in sorted_counts_mali]

		# how many parameters
		n = sum(( len(d) for d in positives ))

		#for s in sorted_counts_mali:
			#print s
		choice = choose_equations(sorted_counts_single_param, n, data)

		print child_param_single
		pt(choice)

		vectors = [eq[0] for eq in choice]
		b = [1.0 - (1.0 - eq[1]) * (1.0 - leak) ** (sum(eq[0]) - 1) for eq in choice]
		b = [log(1.0 - bb) for bb in b]

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
			solution.append([child_param_single, decode_eq_header(data,new_eq), decode_eq_single(data,new_eq), Y[i]])
		#pt(solution)
		final_solution.extend(solution)
	return final_solution

	#return tuple(solution)



def main():
	data = read_data(sys.argv[1], ' ')

	# TODO: HARDCODED child_node and negatives!
	data['negatives'] = ['False','False','False','Medium']

	#simple cancer max
	#data['negatives'] = ['False']*2
	# noisy or network
	#data['negatives'] = ['False']*4
	data['negative_child'] = 'No'
	data['positives'] = get_positives(data)
	data['counter'] = Counter(data['data'])
	data['child_node'] = -1

	LEAK  = get_leak(data)
	REAL = (
		('Malignant', 'Smoker', 'TwoPacks', 0.23),
		('Malignant', 'Smoker', 'OnePack', 0.11),
		('Malignant', 'Genetic', 'True', 0.25),
		('Benign', 'Smoker', 'TwoPacks', 0.25),
		('Benign', 'Smoker', 'OnePack', 0.14),
		('Benign', 'Genetic', 'True', 0.55),
		#('No', 'Smoker', 'TwoPacks', 0.52),
		#('No', 'Smoker', 'OnePack', 0.75),
		#('No', 'Genetic', 'True', 0.2),
	)
	REAL_LEAK = (
		('Malignant', 0.005),
		('Benign', 0.015),
		('No', 0.98),
	)
	GREEDY_PARAMS = get_greedy_params(data, LEAK)
	#for P in PARAMS:
	#	print P
	for L in LEAK.items():
		print L
	SMART_PARAMS = get_smart_params(data, LEAK)
	#print "GREEDY"
	#pt(GREEDY_PARAMS)
	#print "SMART"
	#pt(SMART_PARAMS)
	REALmap = sorted([("|".join(v[:-1]), v[-1] ) for v in REAL], key=lambda i:i[0])
	GREEDYmap = sorted([("|".join(v[:-1]), v[-1] ) for v in GREEDY_PARAMS], key=lambda i:i[0])
	SMARTmap = sorted([("|".join(v[:-1]), v[-1] ) for v in SMART_PARAMS], key=lambda i:i[0])
	#for i in SMARTmap:
		#print i
	print "NAME REAL GREEDY SMART"
	for i in range(len(REALmap)):
		print REALmap[i][0], REALmap[i][1], GREEDYmap[i][1], SMARTmap[i][1]
		#print REALmap[i], REALmap[i], GREEDYmap[i], SMARTmap[i]

def match_by_column(data_counter, col_idx):
    """
        Merge data by column
    """
    new_counter = {}
    for k, v in data_counter.iteritems():
        merged_key = k[:col_idx] + k[col_idx+1:]
        if merged_key not in new_counter:
            new_counter[merged_key] = defaultdict(lambda: 0)
        new_counter[merged_key][k[col_idx]] = v
    return new_counter
    #new_counter = {}
    #for k, v in data_counter.iteritems():
def binarize(data_counter):
    binary_data = {}

    def binname(name):
        return {'False':0, 'True':1}[name]

    for k, v in data_counter.iteritems():
        new_key = tuple([binname(name) for name in k])
        denom = v["True"] + v["False"]
        binary_data[new_key] = (v["True"] / float(denom), denom)
    return binary_data
def evalf(eq, s, orig_subs):
    for k, v in s.items():
        locals()[k] = v
    def c(val):
        return float(1.0 - 2**val)
    if type(eval(eq)) != float:
        print eval(eq), c(eval(eq).evalf(subs=orig_subs))
    else:
        print eq, c(eval(eq))

def mainGJ():
    """ Main execution using GaussJordan elimination"""
    data = read_data(sys.argv[1], ' ')
    c = Counter(data['data'])
    new_counter = match_by_column(c, 4)

    binary_data = binarize(new_counter)

    items = sorted(binary_data.items(), key = lambda x: x[1][1], reverse=True)

    def leak_exponent(k):
        return (-sum(k)+1,)
        #return (1,)
        #return ()

    log_base = 2

    A_vect = [leak_exponent(k) + k for k, v in items if v[0] not in(1.0, 0.0)]
    A = np.array(A_vect) * Fraction(1,1)

    b_vect = [ v[0] for k, v in items if v[0] not in (1.0, 0.0) ]
    b_vect = [ log(1.0 - b, log_base) for b in b_vect]
    b_cnt = [ v[1] for k, v in items if v[0] not in (1.0, 0.0) ]

    for i in xrange(A.shape[0]):
        print A_vect[i], b_vect[i], b_cnt[i]

    b = np.array([sp.Symbol('b%d'%(i), real=True) for i in range(0, A.shape[0])])
    subs = dict(zip(b,b_vect))
    #for i in sorted([ (int(str(k)[1:]),v) for k, v in subs.items()], key=lambda x:x[0]):
        #print i
    
    A2, b2 = GaussJordanElimination(A, b)
    b3 = [1.0 - float(log_base**b.evalf(subs=subs)) for b in b2]
    subs_str = tuple([(str(k), v) for k, v in subs.iteritems()]) + tuple([("r%d"%i, b2[i]) for i in range(len(b2)) ])
    subs_str = dict(subs_str)
    print augment([A2, b2, b3])
    inp = 0
    while True:
        inp = raw_input()
        if inp =='0':
            break
        evalf(inp, subs_str, subs)


    #for k,v in items:
    #    print k,v

    #for k,v in c.iteritems():
    #    for i in k:
    #        if i == 'True':
    #            print 1,
    #        else:
    #            print 0,
    #    print v

if __name__ == "__main__":
	#test_meets_rules()
	#main()
    mainGJ()
