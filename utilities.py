from collections import Counter, defaultdict

ZERO_PROB_EPS = 0.01
ONE_PROB_EPS = 0.99
NO_PROB_EPS = 0.5

def read_data(filename, delimiter):
	"""
	rtype: dict
	return: {'header': ['Smoker','BadDiet','Cancer'], 'data':[ ['True','False','Yes'],...], 'domain':{'Smoker':['True','False'],...}}

    ACCEPTS DATA SUCH AS:
    A B C
    a1 b1 c1
    a2 b2 c2
    a3 b3 c3
    ...
	"""

	data = tuple([tuple(line.strip().split(delimiter)) for line in open(filename).readlines() ])
	header, data = data[0], data[1:]
	domain = {}
	for i in range(len(header)):
		domain[header[i]] = tuple(set([ data[j][i] for j in range(len(data)) ]) )

	return {'header': header, 'data': data, 'domain': domain, }


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


def get_or_default(kwargs, key, default):
	return kwargs[key] if key in kwargs else default


def binarize(data_counter):
	binary_data = {}

	def binname(name):
		return {'False': 0, 'True': 1}[name]

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


class LEAKDEF:
	NONE = 0
	HENRION = 1
	DIEZ = 2
	LEAK_ONLY = 3

def max_based(or_cpt, leakdef = LEAKDEF.DIEZ, handle_errors = False):
	leak = or_cpt[-1]
	n = len(or_cpt) - 1
	generator = (list(int(j) for j in bin(i)[2:].zfill(n)) for i in xrange(2**n - 1, -1, -1))

	def fun(g, or_cpt, leak, child_val):
		name = " ".join([("True" if gi == 1 else "False") for gi in g] + [child_val,])

		if leakdef == LEAKDEF.NONE:
			value =  1.0 - reduce(lambda x,y: x*y, [1.0 - a*b for a,b in zip(g, or_cpt)][:n])
		elif leakdef == LEAKDEF.HENRION:
			leak_exp = sum(g) - 1
			value = 1.0 - reduce(lambda x, y: x * y, [1.0 - a*b for a,b in zip(g, or_cpt)][:n])/((1.0 - leak)**leak_exp)
		elif leakdef == LEAKDEF.DIEZ:
			value = 1.0 - reduce(lambda x,y: x*y, [1.0 - a*b for a,b in zip(g, or_cpt)][:n])/(1.0 - leak)
		elif leakdef == LEAKDEF.LEAK_ONLY:
			if sum(g) == 0:
				value = leak
			else:
				value =  1.0 - reduce(lambda x,y: x*y, [1.0 - a*b for a,b in zip(g, or_cpt)][:n])

		if handle_errors:
			if value <= 0.0:
				value = ZERO_PROB_EPS
			elif value >= 1.0:
				value = ONE_PROB_EPS

		return (name, value)
	ret_list = list(fun(g, or_cpt, leak, "True") for g in generator)
	generator = (list(int(j) for j in bin(i)[2:].zfill(n)) for i in xrange(2**n - 1, -1, -1))
	ret_list_false = [(name, 1.0 - val) for name, val in list(fun(g, or_cpt, leak, "False") for g in generator)]
	return ret_list + ret_list_false

def test():
	or_dist = [0.5491949910554562, 0.38953488372093026, 0.0001, 0.06181150550795611, 0.17593313034240754]
	print or_dist
	leak = or_dist[-1]
	print 1.0 - (1.0 - or_dist[2])*(1.0 - or_dist[3])/(1-leak)

	for i in max_based(or_dist, leakdef = LEAKDEF.HENRION, handle_errors=True):
		print i

if __name__ == "__main__":
	test()
