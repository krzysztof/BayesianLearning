import numpy as np
import sympy as sp
import sys
from math import log, sqrt
from collections import Counter
from GaussJordan import GaussJordanElimination, augment
from fractions import Fraction
from utilities import read_data, match_by_column, get_or_default, binarize, max_based, LEAKDEF

def mainGJ(filename, **kwargs):  # TODO: FIX THE KWARGS !!!
	""" Main execution using GaussJordan elimination"""

	DEBUG = get_or_default(kwargs, 'DEBUG', False)
	data = read_data(filename, ' ')
	c = Counter(data['data'])
	new_counter = match_by_column(c, 4)

	binary_data = binarize(new_counter)

	items = sorted(binary_data.items(), key=lambda x: x[1][2], reverse=True)

	def leak_exponent(k):
		#return (-sum(k)+1,)
		return (1,)
		#return ()

	log_base = 2

	A_vect = [k + leak_exponent(k) for k, v in items if v[0] not in(1.0, 0.0)]
	A = np.array(A_vect) * Fraction(1, 1)

	b_vect = [v[0] for k, v in items if v[0] not in (1.0, 0.0)]
	b_vect = [log(1.0 - b, log_base) for b in b_vect]

	b_cnt = [(v[1], v[2]) for k, v in items if v[0] not in (1.0, 0.0)]

	if DEBUG:
		for i in xrange(A.shape[0]):
			print "b%d"%i, A_vect[i], b_vect[i], b_cnt[i]
	
	b = np.array(sp.symbols('b0:%d' % A.shape[0]))
	subs = dict(zip(b,b_vect))
	subs_cnt = dict(zip(b,b_cnt))
	
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
		return reduce(lambda x, y: x * y, l)

	def _min_fitness(b_val, b_subs_cnt_orig):
		b_subs_cnt = dict((k, v[1]) for k, v in b_subs_cnt_orig.iteritems())
		total = sum(b_subs_cnt.values())
		coeff = [(b.args if b.args else (1, b)) for b in (b_val.args if not type(b_val)==sp.Symbol else [b_val])]
		min_c = min(b_subs_cnt[c[1]] for c in coeff)
		return min_c / float(total)

	def _avg_fitness(b_val, b_subs_cnt_orig):
		b_subs_cnt = dict((k, v[1]) for k, v in b_subs_cnt_orig.iteritems())
		total = sum(b_subs_cnt.values())
		coeff = [(b.args if b.args else (1, b)) for b in (b_val.args if not type(b_val) == sp.Symbol else [b_val])]
		#print coeff
		return sum(b_subs_cnt[s[1]] / float(total) for s in coeff)/ float(sum(abs(s) for s, _ in coeff))
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
		
	result =  max_based(algo, LEAKDEF.LEAK_ONLY)
	nice_result(result)
