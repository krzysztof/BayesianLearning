import numpy as np
import sympy as sp
import sys
import random
from math import log, sqrt
from collections import Counter, defaultdict
from itertools import combinations
from math import e,log
from GaussJordan import GaussJordanElimination, augment, PRODUCT_OPERATIONS
from fractions import Fraction
import matplotlib.pyplot as plt

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
		rules1 = tuple([set( [j] ) for j in negatives] + [set( [domain[header[child_node]][i]]) ] )
		data_i = filter_data(data['data'], rules1)
		LEAK[domain[header[child_node]][i]] = len(data_i) / float(len(data_total))
	return LEAK


def get_greedy_params(data, LEAK):
	negatives = data['negatives']
	header = data['header']
	domain = data['domain']
	positives = data['positives']
	child_node = data['child_node']
	counter = data['counter']

	params = []
	for i in range(len(positives)):  # for each parameter
		for j in range(len(positives[i])):  # for each value
			total = 0
			for k in range(len(domain[header[child_node]])):
				rules_i = tuple(negatives[:i]) + tuple([positives[i][j]]) + tuple(negatives[i+1:]) + tuple([domain[header[child_node]][k]])
				total += counter[rules_i]

			for k in range(len(domain[header[child_node]])):
				rules_i = tuple(negatives[:i]) + tuple([positives[i][j]]) + tuple(negatives[i+1:]) + tuple([domain[header[child_node]][k]])
				value = counter[rules_i] / float(total)
				params.append((domain[header[child_node]][k], header[i], positives[i][j], value ))
				#print rules_i, counter[rules_i], total, counter[rules_i]/float(total)
	return [p for p in params if p[0] != data['negative_child']]


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


def binarize(data_counter):
	binary_data = {}

	def binname(name):
		return {'False':0, 'True':1}[name]

	for k, v in data_counter.iteritems():
		new_key = tuple([binname(name) for name in k])

		#FIXING 1.0 and 0.0
		if v["True"] == 0:
			v["True"] += 1
		if v["False"] == 0:
			v["False"] += 1
		denom = v["True"] + v["False"]
		nominator = v["True"]
		binary_data[new_key] = (v["True"] / float(denom), nominator, denom)
	return binary_data


def max_based(max_cpt):
	leak = max_cpt[-1]
	generator = (list(int(j) for j in bin(i)[2:].zfill(4)) for i in range(15,-1,-1))
	def fun(g, max_cpt, leak):
		name = " ".join([("True" if gi==1 else "False") for gi in g]) + " True"
		leak_exp = sum(g)-1
		value =  1.0 - reduce(lambda x,y: x*y, [1.0 - a*b for a,b in zip(g, max_cpt)][:4])/((1.0 - leak)**leak_exp)
		#TODO: FIX - THIS LEADS TO NEGATIVE VALUES IN HELLINGER DISTANCE
		#value = 1.0 - reduce(lambda x,y: x*y, [1.0 - a*b for a,b in zip(g, max_cpt)][:4])/(1.0 - leak)
		#TODO: ugly fix, doesn't count LEAK at all! yet hellinger works
		#value =  1.0 - reduce(lambda x,y: x*y, [1.0 - a*b for a,b in zip(g, max_cpt)][:4])
		return (name, value)
	return list(fun(g, max_cpt, leak) for g in generator)


def evalf(eq, s, orig_subs):
	for k, v in s.items():
		locals()[k] = v
	def c(val):
		return float(1.0 - 2**val)
	if type(eval(eq)) != float:
		print eval(eq), c(eval(eq).evalf(subs=orig_subs))
	else:
		print eq, c(eval(eq))


def get_or_default(kwargs, key, default):
	return kwargs[key] if key in kwargs else default


def mainGJ(filename, **kwargs):  #TODO: FIX THE KWARGS !!!
	""" Main execution using GaussJordan elimination"""

	DEBUG = get_or_default(kwargs, 'DEBUG', False)
	data = read_data(filename, ' ')
	c = Counter(data['data'])
	new_counter = match_by_column(c, 4)
#	for k,v in sorted(new_counter.iteritems(), key=lambda x: "".join(x[0]), reverse=True):
#		print k,v.items(), v['True']/float(v['False']+ v['True'])
#
#	for i in range(5):
#		priora = [d for d in data['data'] if d[i] =='True']
#		print i,len(priora), len(data['data']), float(len(priora))/len(data['data'])
   #return

	binary_data = binarize(new_counter)

	items = sorted(binary_data.items(), key = lambda x: x[1][2], reverse=True)

	def leak_exponent(k):
		#return (-sum(k)+1,)
		return (1,)
		#return ()

	log_base = 2

	A_vect = [k + leak_exponent(k) for k, v in items if v[0] not in(1.0, 0.0)]
	A = np.array(A_vect) * Fraction(1,1)

	b_vect = [v[0] for k, v in items if v[0] not in (1.0, 0.0) ]
	b_vect = [log(1.0 - b, log_base) for b in b_vect]

	b_cnt = [(v[1], v[2]) for k, v in items if v[0] not in (1.0, 0.0) ]

	#A_vect = [leak_exponent(k) + k for k, v in items]
	#A = np.array(A_vect) * Fraction(1,1)
	#b_vect = [ v[0] for k, v in items]
	#b_vect = [ log(1.0 - b, log_base) for b in b_vect]
	#b_cnt = [ v[1] for k, v in items ]

	if DEBUG:
		for i in xrange(A.shape[0]):
			print "b%d"%i, A_vect[i], b_vect[i], b_cnt[i]
	
	#b = np.array([sp.Symbol('b%d'%(i), real=True) for i in range(0, A.shape[0])])
	b = np.array(sp.symbols('b0:%d' % A.shape[0]))
	subs = dict(zip(b,b_vect))
	subs_cnt = dict(zip(b,b_cnt))
	#for i in sorted([ (int(str(k)[1:]),v) for k, v in subs.items()], key=lambda x:x[0]):
		#print i
	
	A2, b2 = GaussJordanElimination(A, b)
	b3 = [1.0 - float(log_base**b.evalf(subs=subs)) for b in b2]

	subs_str = tuple([(str(k), v) for k, v in subs.iteritems()]) + tuple([("r%d"%i, b2[i]) for i in range(len(b2)) ])
	subs_str = dict(subs_str)

	if DEBUG:
		print augment([A2, b2, b3])

	nonzero_i = (i for i in range(A2.shape[0]) if any(j!=0 for j in A2[i]))
	zero_i = (i for i in range(A2.shape[0]) if all(j==0 for j in A2[i]))
	nonzero_v = list((A2[i], b2[i]) for i in nonzero_i)
	zero_v = list((A2[i], b2[i]) for i in zero_i)

	def product(l):
		return reduce(lambda x,y:x*y, l)

	def _min_fitness(b_val, b_subs_cnt_orig):
		b_subs_cnt = dict( (k,v[1]) for k, v in b_subs_cnt_orig.iteritems())
		total = sum(b_subs_cnt.values())
		coeff = [(b.args if b.args else (1, b)) for b in (b_val.args if not type(b_val)==sp.Symbol else [b_val])]
		min_c = min(b_subs_cnt[c[1]] for c in coeff)
		return min_c/float(total)

	def _avg_fitness(b_val, b_subs_cnt_orig):
		b_subs_cnt = dict( (k,v[1]) for k, v in b_subs_cnt_orig.iteritems())
		total = sum(b_subs_cnt.values())
		coeff = [(b.args if b.args else (1, b)) for b in (b_val.args if not type(b_val)==sp.Symbol else [b_val])]
		#print coeff
		return sum(b_subs_cnt[s[1]]/float(total) for s in coeff)/ float(sum(abs(s) for s,_ in coeff))
		#return sum(abs(s[0])*(b_subs_cnt[s[1]]/float(total)) for s in coeff) / sum(b_subs_cnt[s[1]]/float(total) for s in coeff)
		#return 1

	def _max_count_fitness(b_val, b_subs_cnt_orig):
		b_subs_cnt = dict( (k,v[1]) for k, v in b_subs_cnt_orig.iteritems())
		total = sum(b_subs_cnt.values())
		coeff = [(b.args if b.args else (1, b)) for b in (b_val.args if not type(b_val)==sp.Symbol else [b_val])]
		return sum(b_subs_cnt[s[1]]/abs(s[0]) for s in coeff) / float(total)
	
	def _pu(x,n,c):
		n = float(n)
		x = float(x)
		c = float(c)
		sqr = sqrt(((x/n)*(1.0-x/n))/n)
		return c*sqr
		#return x/n-Ualph*sqr,x/n+Ualph*sqr

	def _pu_fitness(b_val, b_subs_cnt):
		#total = sum(b_subs_cnt.values())
		coeff = [(b.args if b.args else (1, b)) for b in (b_val.args if not type(b_val)==sp.Symbol else [b_val])]
		#return 1.0 - max(b_subs_cnt[b][0]/float(b_subs_cnt[b][1]) - _pu(b_subs_cnt[b][0], b_subs_cnt[b][1], 1.65)[0] for c, b in coeff)
		#return 1.0 - max(b_subs_cnt[b][0]/float(b_subs_cnt[b][1]) - abs(c)*_pu(b_subs_cnt[b][0], b_subs_cnt[b][1], 1.65) for c, b in coeff)
		return 1.0 - max(abs(c)*_pu(b_subs_cnt[b][0], b_subs_cnt[b][1], 1.65) for c, b in coeff)
		
	#fitness = _min_fitness
	#fitness = _avg_fitness
	fitness = _pu_fitness

	#BELOW: poor fitness!
	#fitness = _max_count_fitness

	solutions = []
	for i in nonzero_v:
		for zv in ([(0,0)] + zero_v):
			for coeff in [2, 1,-1, -2]:
				expr = (i[1] + coeff*zv[1])
				fit = fitness(expr, subs_cnt)
				#print i[0], " [",coeff,"]", zv[0], "expr:",expr, "value:",float(1.0 - log_base**expr.evalf(subs=subs)), "fitness:", fit
				solutions.append((i[0],'V' if type(zv[0])!=int else '0', coeff, zv[1],"EXPR:", expr, float(1.0 - log_base ** expr.evalf(subs=subs)), fit))
				if type(zv[0]) == int:
					break

	GJElim_fit_distribution = []
	for i in range(5):
		solutions_filtered = [s for s in sorted(solutions, key= lambda x: x[-1], reverse=True) if s[0][i] == 1][:5]
		GJElim_fit_distribution.append(solutions_filtered[0][-2])
		suma = sum(s[-1]*s[-2] for s in solutions_filtered)
		if DEBUG:
			for s in solutions_filtered:
				print s
			print suma / sum(s[-1] for s in solutions_filtered)
			print ""

	if DEBUG:
		print augment([A2, b2, b3])

	GJElim_distribution = []
	for i in range(5):
		for j in range(A2.shape[0]):
			if A2[j][i] == 1:
				GJElim_distribution.append(b3[j])
				break
	GJElim_distribution = [(d if d>0 else 10e-5) for d in GJElim_distribution]
	GJElim_fit_distribution = [(d if d>0 else 10e-5) for d in GJElim_fit_distribution]

	return GJElim_distribution, GJElim_fit_distribution
#	for i in range(5):
#		genie_dist.append(float(raw_input()))
#	print "eucl_Genie:", eucl_dist(real_or, genie_dist)
	#print "eucl_PU:", eucl_dist(real_or, GJElim_distribution)
	
#	inp = 0
#	while True:
#		break #FIXME TO ENTER TYPE MODE
#		inp = raw_input()
#		if inp =='0':
#			break
#		evalf(inp, subs_str, subs)

	#for k,v in items:
	#	print k,v

	#for k,v in c.iteritems():
	#	for i in k:
	#		if i == 'True':
	#			print 1,
	#		else:
	#			print 0,
	#	print v

def RunTests():
	

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

	CancerOR_2k_em = (0.7637068, 0.7538612499999999,0.8055532765957447,0.6244346428571429,0.6712355999999999,0.41795375,0.6472175444015444,0.6226004694835681,0.84853,0.34221875,0.3025307189542484,0.3225763358778626,0.06405333333333332,0.1655,0.05368442211055276,0.006390658174097664, 0.2362932,0.24613875,0.1944467234042553,0.3755653571428571,0.3287644,0.58204625,0.3527824555984556,0.3773995305164319,0.15147,0.65778125,0.6974692810457517,0.6774236641221374,0.9359466666666667,0.8344999999999999,0.9463155778894472,0.9936093418259023)
	CancerOR_2k_em = CancerOR_2k_em[:16]
	CancerOR_2k_2_em = (0.9405523906542476,0.818769967985576,0.6854689515035244,0.7366789812195553,0.8553335209956948,0.8529915056865248,0.5980989195000936,0.5731359203156429,0.2529192597950442,0.4992795021053976,0.2695247596004046,0.2854046027911718,0.4664779287831682,0.04633065152582788,0.06117990592057187,0.01043861815928428)

	CancerOR_4k_em = (0.7567733324386966,0.9215699715581367,0.751491194077176,0.7443145406874079,0.6699447273165121,0.6121346838796311,0.6581018686820314,0.6059955010892189,0.7032269309025409,0.338461383470633,0.2962101658385292,0.2806334738238965,0.3211246951152935,0.1223849393059023,0.05147109606118423,0.006550251569634347)
	CancerOR_3k_em_error = (0.8019937785359899,0.7770205376906139,0.754325814963958,0.7233405229512437,0.7149023113299757,0.6789447884250167,0.6462679717748254,0.6016540446452869,0.5029690289303663,0.4402817352869736,0.3833139289646046,0.3055353134035032,0.2843538955409726,0.194094092618776,0.1120694481857673,8.052500713440214e-005)
	CancerOR_3k_em = (0.8823862344850944,0.7647521189870088,0.7723250621112042,0.7001730537939015,0.8084053891488456,0.6658808209441537,0.6485117400204923,0.6398550757663786,0.7112086852239463,0.6055918506677425,0.3241390883654872,0.2288411100961994,0.02318613836861917,0.1452788122261852,0.05922565811074033,0.009985793609668024)
	CancerOR_2500_em = (0.7637068,0.9179537499999999,0.6408119603960396,0.7612991176470588,0.5616474666666667,0.534363,0.6201787216494845,0.6177457692307693,0.4647066666666667,0.4737750000000001,0.261669387755102,0.2608487654320987,0.03202666666666666,0.2698125,0.05247394468704512,0.01660033167495854)


	CancerOR_1500_em = (0.8019663872396284,0.9442440053389744,0.7500806648057734,0.7086380938946789,0.7429781278371519,0.4197279907295194,0.6022090574165052,0.574126978883834,0.2514923295821444,0.7201325138312157,0.2337176282554325,0.1990732125393137,0.1422230444484985,0.09707425935402761,0.04754259670217673,0.01554425006478532)

	CancerOR_1k_em = (0.8850064394263117,0.7441086909716402,0.6416538144375366,0.6681303438925033,0.6580000105816407,0.897637756386483,0.61273001713991,0.5978799943882992,0.5008748877009172,0.2223284118737788,0.3645939109286032,0.3392866952333448,0.2013551920495993,0.02228025261973483,0.01808960944272913,0.004424969035062378)
	CancerOR_800_em = (0.8150984406330213,0.900465592958686,0.7869781540687645,0.5925330902938479,0.881254718913275,0.7443127860601015,0.595216983732205,0.5823486011700594,0.7245970127555574,0.2034954910425953,0.234047377333469,0.2698262234734214,0.3231339679101403,0.1191647990404533,0.03391271348139648,0.0158733473100359)
	CancerOR_2100_em = (0.7570173363776163,0.7570088499116994,0.7289893500302636,0.7604699543092656,0.6074009849394084,0.822196745039546,0.6534453405632025,0.5935653253878364,0.1426005854636031,0.1069454429172962,0.2299515374901263,0.3116668289266736,0.430217821952639,0.1660835166713437,0.04568900107205483,0.01098097146474772)
	CancerOR_2200_em = (0.7408792847396678,0.8637452282160795,0.6940910577003088,0.7319484646594838,0.8593205006360379,0.9013677052089797,0.6205762049898314,0.6317292575793824,0.2709233362070374,0.1166258867482374,0.3084225862411021,0.2247152467274453,0.1726217781686412,0.04186195797902215,0.05053463118164066,0.01200068092559615)
	CancerOR_2300_em = (0.5982492621147051,0.775171843394545,0.7135442297647495,0.7441044948831113,0.4386315464215321,0.5462228210024325,0.6023450667169213,0.5450352374766989,0.2838037095998791,0.2430992885579457,0.316229930851719,0.2376007108321504,0.5358641423756911,0.01544092131684935,0.04232534313519233,0.007007457956962515)
	CancerOR_2400_em = (0.9192291441837235,0.8733779402550531,0.6849574073099998,0.7212057403831821,0.5424579160066027,0.8997699916548303,0.6792241628860898,0.6437831230831829,0.07710693408047706,0.09832152195804919,0.2803334305266786,0.2524084909385175,0.1812178122280657,0.005766519019329675,0.05388727520122803,0.0197834042895137)
	CancerOR_3100_em = (0.5534296304766789,0.5862992614654351,0.7191251078285286,0.6922119279112432,0.4322773922369548,0.929877856047702,0.6461676578362225,0.56074732698091,0.05052469964710422,0.3912758985849051,0.2793336165729448,0.2740185117209741,0.1071753981770754,0.009716130759596403,0.05830104696032093,0.01546054113871081)
	CancerOR_3300_em = (0.9026360430457321,0.8924484138015509,0.7213883572687733,0.6893225958243857,0.9300349413607905,0.8454283644480688,0.660267206413334,0.6162799685054158,0.6103156190567145,0.7942664177099915,0.2964373350009358,0.2502645844452496,0.6393980489799679,0.01412096478574662,0.06123811439481969,0.01432299894156858)

	CancerOR_3500_em = (0.8396287833410943,0.8396151660409161,0.6891223512220257,0.66964830254456,0.7010543443556048,0.2680669177666602,0.648336977955042,0.6231526945916434,0.8078906076517771,0.2314971817097771,0.2849870744712538,0.2592976330303004,0.06636778182627838,0.1212826383988348,0.05726435633583268,0.0111323405277349)
	CancerOR_3700_em = (0.9377407672873022,0.9590366616358816,0.7348345595540176,0.7354445763583706,0.7030340622969287,0.6524834236408239,0.6396099342696007,0.6026179103495799,0.6375850054839433,0.2763605568400151,0.2729283453331389,0.2602508914307488,0.1983266102024267,0.009616035311797808,0.04368074114389275,0.005800491812169645)
	CancerOR_3900_em = (0.9329594439426322,0.8484684447119738,0.7079585165241424,0.661628873589612,0.7342635528034799,0.7573171502392182,0.6248919978469948,0.6158116047893651,0.3104000601683837,0.2857239299998058,0.3101011423536778,0.2562347939045284,0.0388284900401115,0.133135281493882,0.04303684842050223,0.008849602624313131)

	CancerOR_1100_em = (0.7824243699388535,0.3901013423673716,0.6904253620843223,0.5488086546780142,0.426870090316565,0.7044934157471144,0.7104796831180581,0.5735038013186814,0.2669466955776569,0.1639701894579935,0.3211913973525691,0.2396852251316635,0.02111643417502209,0.3654830621245268,0.05681898374421952,0.01953720853075903)
	CancerOR_1200_em = (0.6368549928372216,0.6350901137208117,0.8420292821784607,0.7466416719925483,0.2929040178661636,0.6459487672985403,0.5809733182167276,0.6054065226878272,0.2461739038664071,0.2425103207020249,0.3417628529304773,0.3027888121738482,0.02336821135941388,0.03400767462549593,0.04929108004958355,0.01894036133367874)
	CancerOR_1300_em = (0.6815653906364617,0.6815608775501265,0.7281951199500025,0.7185213986307627,0.518455318305393,0.2776798266199361,0.6678111322331951,0.6639219342683596,0.2985485992628626,0.2985386577846487,0.3047031948847608,0.3329630714485696,0.2041110076804497,0.01027057849879115,0.0701362295336862,0.006780292870456975)
	CancerOR_1400_em = (0.8483746477811397,0.8272991634458835,0.8126676926044322,0.7774897257744511,0.5248069554853283,0.3788166502048111,0.5528719204718661,0.6177465667811998,0.3788326468695611,0.2924921857995902,0.3368405460650412,0.269068415465145,0.02565193744809742,0.00177257782358431,0.05209050532587214,0.01009061919169317)
	CancerOR_1600_em = (0.8648236904740543,0.7189471616228647,0.7734860026541823,0.6199806569936621,0.8622453705731462,0.6425949106694898,0.6631447710155379,0.6216383831442848,0.456945419170511,0.1729673708165105,0.2997881758797291,0.2562722533128288,0.339708935501979,3.777820988215552e-005,0.05060663092362704,0.007772230190587289)
	CancerOR_1700_em = (0.3875504814270279,0.7170907844297954,0.7317928412851302,0.700009831288094,0.8712930732911424,0.8380945684083765,0.6078522353142901,0.615584695379913,0.1719702655073553,0.05823912203236359,0.2615220011045646,0.2757482743299367,0.04151510375068712,0.01846732998118546,0.06047491933846234,0.0203045992989147)
	CancerOR_1800_em = (0.7099369405043438,0.9015601326028488,0.7124680016248253,0.6426787234599686,0.6680303764169048,0.8887265563544107,0.6064311113114125,0.5857579927874907,0.08691889075763964,0.05295533767411701,0.2365717332909318,0.2415409493960194,0.2295701398546997,0.02460381259587108,0.04224429159451549,0.002341770520198471)
	CancerOR_1900_em = (0.1845606481393627,0.9828754886000373,0.7027078097553922,0.6164570629413035,0.9818700989484301,0.6963630589505946,0.5846793176200106,0.6559639803246075,0.9099963294828243,0.03802632449920801,0.3355191970765169,0.3178183454010196,0.1653535699049097,0.7650916976958262,0.03597136158788174,0.006928560338436034)


	CancerOR_4k_hill = (0.6016219652511965,0.3054793871261965,0.1940291918136966,0.1119979418136965,8.052500713440214e-005)
	CancerOR_3k_hill_error = (0.6015543254109895,0.3054117472859899,0.19408362228599,0.1120523722859902, 8.052500713440214e-005)
	CancerOR_3k_hill = (0.6391001016725941,0.3836069376100944,0.1852426797975946,0.01593115636009479,9.809829716789231e-005)
	CancerOR_2k_hill = (0.5392194953124999,0.3750349249999999,0.1113630499999999,0.1012312140624999, 4.300192409767423e-005)
	CancerOR_2k_2_hill = (0.664307273466747,0.221436179716748,0.2382818828417481,0.07116762502924823,6.264159678437409e-005)
	CancerOR_1500_hill = (0.5927378716146284,0.3105113091146283,0.1621958794271283,8.650442712831996e-005,4.132037728532545e-005)
	CancerOR_1k_hill = (0.6358609316138117,0.2808804628638116,0.1014371034888116,0.08141757223881163,3.043778506239825e-005)
	CancerOR_800_hill = (0.6310164093830212,0.2189070343830212,0.2254988312580212,0.1317488312580212,4.206685035978985e-006)

	CancerOR_2500_hill = (0.5753523078125,0.3166853156249999,0.1564069953124999,3.492499999990795e-005,0.007322987924958446)
	CancerOR_2100_hill = (0.6445905785651162,0.07903882075261626,0.09075757075261626,0.04913159419011626,0.08434522927724775)
	CancerOR_2200_hill = (0.7134134644271678,0.1638529175521678,0.0605814331771678,0.0870706909896678,0.001258493425596252)
	CancerOR_2300_hill = (0.5217111761772051,0.2010324652397051,0.09300023867720508,0.1397531683647051,0.009570934519462448)
	CancerOR_2400_hill = (0.7081695738712235,0.1404205504337235,3.969105872436707e-005,0.03421937855872326,0.03882637303951364)
	CancerOR_3100_hill = (0.5826044351641789,0.180871036726679,8.490391417903531e-005,8.490391417914633e-005,0.06172518957621043)
	CancerOR_3300_hill = (0.6464104571082321,0.2899651446082322,0.3964104571082321,0.1151604571082322,4.077237906852105e-005)
	CancerOR_3500_hill = (0.5361619864660943,0.3371873770910941,0.06728991615359414,0.1314989005285941,2.394209023492966e-005)
	CancerOR_3700_hill = (0.6514858844748021,0.3041958454123021,0.1129116657248021,0.09472318916230205,6.318712466968357e-005)
	CancerOR_3900_hill = (0.6729496783176321,0.2562016314426323,0.09506881894263242,0.01010788144263253,0.001525383874313113)
	CancerOR_1100_hill = (0.5182642136888532,0.1232446824388535,0.02497808087635345,0.004836479313853559,0.113897560093259)
	CancerOR_1200_hill = (0.5460346803372216,0.2838276490872217,1.417252472157848e-005,1.417252472157848e-005,0.02052727539617871)
	CancerOR_1300_hill = (0.5118051473416034,0.2874399129666034,8.639734160342805e-005,0.1220346395291034,0.007004519124565345)
	CancerOR_1400_hill = (0.5702984759061397,0.3458111712186397,0.0001080462186396902,0.03807191340613969,8.085356669318422e-005)
	CancerOR_1600_hill = (0.6571971783384363,0.1263133892759363,0.05539053771343605,0.2050487408384363,1.212377239223628e-005)
	CancerOR_1700_hill = (0.6489030204895278,0.1153336845520278,9.930955202774872e-005,0.01779950486452775,0.04911319304891471)
	CancerOR_1800_hill = (0.6768558858168438,0.09763225300434375,0.02536662800434375,0.04465373737934375,0.03591110645769724)
	CancerOR_1900_hill = (0.5424708043893627,0.1623438512643627,0.3456934606393627,0.005727640326862682,0.09994613846343603)


	GenieTests = {"CancerOR_1k.txt": {"em": CancerOR_1k_em, "hill": CancerOR_1k_hill, "num":1000},
			"CancerOR_2k.txt": {"em": CancerOR_2k_em, "hill": CancerOR_2k_hill, "num":2000},
			#"CancerOR_2k_2.txt":{"em": CancerOR_2k_2_em, "hill": CancerOR_2k_2_hill, "num":2000 },
			"CancerOR_3k.txt":{"em": CancerOR_3k_em, "hill": CancerOR_3k_hill, "num":3000}, 
			"CancerOR_4k.txt":{"em": CancerOR_4k_em, "hill": CancerOR_4k_hill, "num":4000},
			"CancerOR_1100.txt":{"em": CancerOR_1100_em, "hill": CancerOR_1100_hill, "num":1100},
			"CancerOR_1200.txt":{"em": CancerOR_1200_em, "hill": CancerOR_1200_hill, "num":1200},
			"CancerOR_1300.txt":{"em": CancerOR_1300_em, "hill": CancerOR_1300_hill, "num":1300},
			"CancerOR_1400.txt":{"em": CancerOR_1400_em, "hill": CancerOR_1400_hill, "num":1400},
			"CancerOR_1500.txt":{"em": CancerOR_1500_em, "hill": CancerOR_1500_hill, "num":1500},
			"CancerOR_1600.txt":{"em": CancerOR_1600_em, "hill": CancerOR_1600_hill, "num":1600},
			#"CancerOR_1700.txt":{"em": CancerOR_1700_em, "hill": CancerOR_1700_hill, "num":1700}, #WYSYPUJE SIE BO WSZYSTKIE WEKTORY Z x2 sa 1.0 lub 0.0
			"CancerOR_1800.txt":{"em": CancerOR_1800_em, "hill": CancerOR_1800_hill, "num":1800},
			"CancerOR_1900.txt":{"em": CancerOR_1900_em, "hill": CancerOR_1900_hill, "num":1900},
			"CancerOR_800.txt":{"em": CancerOR_800_em, "hill": CancerOR_800_hill, "num":800},
			"CancerOR_2500.txt":{"em": CancerOR_2500_em, "hill": CancerOR_2500_hill, "num":2500},
			"CancerOR_2400.txt":{"em": CancerOR_2400_em, "hill": CancerOR_2400_hill, "num":2400},
			"CancerOR_2300.txt":{"em": CancerOR_2300_em, "hill": CancerOR_2300_hill, "num":2300},
			"CancerOR_2200.txt":{"em": CancerOR_2200_em, "hill": CancerOR_2200_hill, "num":2200},
			"CancerOR_2100.txt":{"em": CancerOR_2100_em, "hill": CancerOR_2100_hill, "num":2100},
			"CancerOR_3100.txt":{"em": CancerOR_3100_em, "hill": CancerOR_3100_hill, "num":3100},
			"CancerOR_3300.txt":{"em": CancerOR_3300_em, "hill": CancerOR_3300_hill, "num":3300},
			"CancerOR_3500.txt":{"em": CancerOR_3500_em, "hill": CancerOR_3500_hill, "num":3500},
			"CancerOR_3700.txt":{"em": CancerOR_3700_em, "hill": CancerOR_3700_hill, "num":3700},
			"CancerOR_3900.txt":{"em": CancerOR_3900_em, "hill": CancerOR_3900_hill, "num":3900},
			}

	real_or = (0.61, 0.25, 0.15, 0.04, 0.01,)
	X_plot = []
	Y_plot = []
	dist_func = hellinger_dist
	for k, v in sorted(GenieTests.items(), key=lambda i:i[1]['num']):
		GJElim_distribution, GJElim_fit_distribution = mainGJ(k)

		CancerOR_em = v["em"]
		CancerOR_hill = v["hill"]

		CancerOR_hill_cpt = max_based(CancerOR_hill)
		GJElim_cpt = max_based(GJElim_distribution)
		GJElim_fit_cpt = max_based(GJElim_fit_distribution)
		real_cpt = max_based(real_or)
		

		print (":"*10)+"Running test for %s"%(k) + (":"*10)
		print "real_OR", real_or
		print "GJ_OR", GJElim_distribution
		print "GJFit_OR", GJElim_fit_distribution
		print "hill_OR:", CancerOR_hill
		print ""

		print "CPT(real_OR -> cpt)%s\n"%(real_cpt,)
		print "CPT(EM from data)%s\n"%(CancerOR_em,)
		print "CPT(hill -> cpt)%s\n"%(CancerOR_hill_cpt,)
		print "CPT(GJ -> cpt) %s\n"%(GJElim_cpt,)
		print "CPT(GJFit -> cpt) %s\n"%(GJElim_fit_cpt,)


		print "De(EM vs max_based(real_OR))", dist_func(CancerOR_em, real_cpt,)
		print "De(EM vs max_based(hill_OR))", dist_func(CancerOR_em, CancerOR_hill_cpt)
		print "De(EM vs max_based(GJ_OR))", dist_func(CancerOR_em, GJElim_cpt)
		print "De(EM vs max_based(GJFit_OR))", dist_func(CancerOR_em, GJElim_fit_cpt)

		print "De(max_based(orig_OR) vs max_based(hill_OR))", dist_func(real_cpt, CancerOR_hill_cpt)
		print "De(max_based(orig_OR) vs max_based(GJ_OR))", dist_func(real_cpt, GJElim_cpt)
		print "De(max_based(orig_OR) vs max_based(GJFit_OR))", dist_func(real_cpt, GJElim_fit_cpt)
		print "De(orig_OR vs hill_OR):", dist_func(real_or, CancerOR_hill)
		print "De(orig_OR vs GJ_OR):", dist_func(real_or, GJElim_distribution)
		print "De(orig_OR vs GJFit_OR):", dist_func(real_or, GJElim_fit_distribution)

		X_plot.append(v["num"])
		Y_plot.append((
			("De_orig_hill", dist_func(real_cpt, CancerOR_hill_cpt)), 
			("De_orig_gj", dist_func(real_cpt, GJElim_cpt)),
			("De_orig_gj_fit", dist_func(real_cpt, GJElim_fit_cpt)),
		))
	plot_data = sorted(zip(X_plot, Y_plot), key=lambda x:x[0])

	colors = ((0,'r'),(1,'g'),(2,'b'))
	for n, c in colors:
		Y = [y[n][1] for x, y in plot_data] 
		X = [x for x, y in plot_data] 
		plt.plot(X, Y, color=c, marker="o")

		#plt.title("Hellinger distance (R)")
		plt.title("Hellinger distance (S)")
		plt.xlabel("Number or records")
		plt.ylabel("Distance")
	plt.show()

def mainOutData():
	data = read_data(sys.argv[1], ' ')
	c = Counter(data['data'])
	new_counter = match_by_column(c, 4)

	for i in range(4):
		priora = [d for d in data['data'] if d[i] =='True']
		#print i,len(priora), len(data['data']), float(len(priora))/len(data['data'])
		print float(len(priora))/len(data['data'])

	for k,v in sorted(new_counter.iteritems(), key=lambda x: "".join(x[0]), reverse=True):
		#print k,v.items(), v['True']/float(v['False']+ v['True'])
		print v['True']/float(v['False']+ v['True'])

   #return


def nice_result(result):
	for state, value in result:
		print "%s %f"%(state, value)

if __name__ == "__main__":
	#test_meets_rules()
	#main()
	#RunTests()
	GJ, GJFit = mainGJ(sys.argv[1])
	if sys.argv[2] == "GJ":
		algo = GJ
	elif sys.argv[2] == "GJFit":
		algo = GJFit
		
	result =  max_based(algo)
	nice_result(result)
