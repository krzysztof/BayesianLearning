from collections import Counter, defaultdict
from itertools import product
from copy import deepcopy

class LEAKDEF:
	NONE = 0
	HENRION = 1
	DIEZ = 2
	LEAK_ONLY = 3

class CPT(object):
	TYPE_NOISY_MAX = 1
	TYPE_CPT = 2
	def __init__(self, matrix, parent_dims, network_type, labels, states):
		self.labels = labels
		self.matrix = matrix
		self.network_type = network_type
		self.parent_dims = parent_dims
		self.states = states

	def __str__(self):
		if self.parent_dims:
			out_matrix = self._split_rows()
		else:
			out_matrix = self.matrix
		return "\n".join(["%s"%(row) for row in out_matrix])

	def _split_rows(self, matrix = None):
		if matrix == None:
			matrix = self.matrix
		new_array = []
		for row in matrix:
			params = row[:-1]
			leak = row[-1]
			sum_begin = 0
			new_row = []
			for pd in self.parent_dims:
				new_row.append(params[sum_begin: sum_begin+pd])
				sum_begin+=pd
			new_row.append([leak,])
			new_array.append(new_row)
		return new_array

	def max_based(self, leakdef = LEAKDEF.DIEZ, handle_errors = False):
		assert self.network_type != self.TYPE_CPT, "Network is already TYPE_CPT"
		cpt_matrix = deepcopy(self.matrix)
		sum_row = [0.0]*len(cpt_matrix[0])
		for r_idx in range(len(cpt_matrix)):
			for c_idx in range(len(cpt_matrix[r_idx])):
				cpt_matrix[r_idx][c_idx] += sum_row[c_idx]
				sum_row[c_idx] = cpt_matrix[r_idx][c_idx]

		splitted_arr = CPT(cpt_matrix, self.parent_dims, self.TYPE_NOISY_MAX, self.labels , self.states)._split_rows()

		# Values array
		max_based_arr = []
		for row in splitted_arr:
			max_based_row = []
			for comb in product(*row):
				max_based_row.append(1.0 - reduce(lambda x,y:x*y, [1.0 - elem for elem in comb]))
			max_based_arr.append(max_based_row)

		# Labels array
		labels_array = []
		for r_idx in range(len(max_based_arr)):
			labels_row = []
			states_lists = self.states[:-1] + [[self.states[-1][r_idx],],]
			for comb in product(*states_lists):
				labels_row.append(" ".join(comb))
				#max_based_row.append(1.0 - reduce(lambda x,y:x*y, [1.0 - elem for elem in comb]))
			labels_array.append(labels_row)

		for i in range(len(max_based_arr)-1, 0, -1):
			max_based_arr[i] = [max_based_arr[i][j] - max_based_arr[i-1][j] for j in range(len(max_based_arr[i]))]

		self.labels_array = labels_array
		new_cpt = CPT(max_based_arr, parent_dims=None, network_type=self.TYPE_CPT, labels=self.labels, states = self.states)
		new_cpt.labels_array = labels_array
		return new_cpt
		
		#return max_based(self, leakdef=leakdef, handle_errors=handle_errors)
	def print_raw(self):
		assert self.network_type == self.TYPE_CPT, "this works only for TYPE_CPT networks"

		ret_str = " ".join(self.labels)+"\n"
		for r_idx in xrange(len(self.matrix)):
			for c_idx in xrange(len(self.matrix[r_idx])):
				ret_str += "%s %s\n" %(self.labels_array[r_idx][c_idx],self.matrix[r_idx][c_idx])
		return ret_str[:-1]
		return "\n".join(" ".join(str(i) for i in row) for row in self.matrix)

def max_based(or_cpt, leakdef = LEAKDEF.DIEZ, handle_errors = False):
	leak = or_cpt[-1]
	params = or_cpt[:-1]
	n = len(params)

	generator = (list(int(j) for j in bin(i)[2:].zfill(n)) for i in xrange(2**n - 1, -1, -1))

	def fun(g, or_params, leak):
		#name = " ".join([("True" if gi == 1 else "False") for gi in g] + [child_val,])
		if leakdef == LEAKDEF.NONE:
			value =  1.0 - reduce(lambda x,y: x*y, [1.0 - a*b[0] for a,b in zip(g, or_params)][:n])
		elif leakdef == LEAKDEF.HENRION:
			leak_exp = sum(g) - 1
			value = 1.0 - reduce(lambda x, y: x * y, [1.0 - a*b[0] for a,b in zip(g, or_params)][:n])/((1.0 - leak)**leak_exp)
		elif leakdef == LEAKDEF.DIEZ:
			value = 1.0 - reduce(lambda x,y: x*y, [1.0 - a*b[0] for a,b in zip(g, or_params)][:n])/(1.0 - leak)
		elif leakdef == LEAKDEF.LEAK_ONLY:
			if sum(g) == 0:
				value = leak
			else:
				value =  1.0 - reduce(lambda x,y: x*y, [1.0 - a*b for a,b in zip(g, or_params)][:n])

		if handle_errors:
			if value <= 0.0:
				value = ZERO_PROB_EPS
			elif value >= 1.0:
				value = ONE_PROB_EPS

		return value
	ret_list = list(fun(g, or_cpt[0][:-1], or_cpt[0][-1]) for g in generator)
	#generator = (list(int(j) for j in bin(i)[2:].zfill(n)) for i in xrange(2**n - 1, -1, -1))
	#ret_list_false = [(name, 1.0 - val) for name, val in list(fun(g, or_cpt, leak, "False") for g in generator)]
	#return ret_list + ret_list_false
	return ret_list

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

def print_CPT(cpt):
	for row_i in range(len(cpt)):
		row = cpt[row_i]
		for param in row[:-1]:
			print "|",
			for state in param:
				print "%.2f"%(state,),
		print "|| %.2f |"%(row[-1],)

def test():
	#or_dist = [
	#	[0.23,0.11,0, 0.25,0,0.15,0,0.04000000000000004,0.13,0,0.005000000000000004],
	#	[0.25,0.14,0,0.55,0,0.6599999999999999,0,0.6299999999999999,0.22,0,0.01500000000000001],
	#	[0.52,0.75,1,0.2,1,0.19,1,0.33,0.65,1, 0.98],
	#	]
	#a = CPT(or_dist, [3, 2, 2, 3])
	#for row in  a._split_rows():
	#	print row
	#b = a.max_based()
	#print b
	#print a
	dist = [[0.61,0, 0.25,0,0.15,0,0.04000000000000004,0,0.01000000000000001],
		[0.39,0,0.75,0,0.85,0,0.96,0,0.99]]
	parent_dims = [2]*4
	a = CPT(dist, parent_dims=parent_dims, network_type=CPT.TYPE_NOISY_MAX, labels=list("ABCDE"), states=[["True", "False"]]*5)
	print a.max_based().print_raw()


	#genie_out = [[0.5311162,0.511579375,0.57507405625,0.4483720000000001,0.4253875,0.5000871250000001,0.3748216000000001,0.3487725,0.433432075,0.2644960000000001,0.23385,0.3334495,0.4580434,0.435461875,0.50885183125,0.3624040000000001,0.3358375,0.422178625,0.2773912,0.2472825,0.3451357749999999,0.149872,0.1144499999999999,0.2295714999999999,0.3910600000000001,0.3656875000000001,0.448148125,0.2836000000000001,0.25375,0.3507625,0.18808,0.15425,0.2641975,0.04480000000000006,0.005000000000000004,0.13435],
	#[0.4624934159999998,0.4690558249999999,0.4123388237499999,0.5179943999999999,0.4726925,0.433664875,0.5932264799999998,0.5544034999999998,0.5036323249999999,0.567336,0.2565499999999999,0.3353105,0.5327396999999999,0.5366081249999999,0.4729936687499999,0.589086,0.5171625,0.482271375,0.6765243,0.6130674999999999,0.5640917249999999,0.607578,0.1505500000000001,0.2926785000000001,0.5966507999999998,0.5970724999999998,0.5276458749999999,0.65172,0.5502499999999999,0.5218375,0.7504739999999999,0.6595499999999999,0.6147724999999999,0.6317999999999999,0.01500000000000001,0.22865],
	#[0.006390384000000001,0.0193648,0.01258712,0.03363360000000001,0.10192,0.06624800000000002,0.03195192,0.09682400000000001,0.06293560000000001,0.168168,0.5096000000000001,0.33124,0.0092169,0.02793,0.0181545,0.04851,0.147,0.09555,0.0460845,0.13965,0.09077250000000001,0.24255,0.735,0.47775,0.0122892,0.03724,0.024206,0.06468,0.196,0.1274,0.06144600000000001,0.1862,0.12103,0.3234,0.98,0.637]]

	#assert len(genie_out) == len(b.matrix), "different sizes"
	#for r_idx in range(len(genie_out)):
	#	assert len(genie_out[r_idx]) == len(b.matrix[r_idx]), "different sizes"
	#	for c_idx in range(len(genie_out[r_idx])):
	#		assert abs(genie_out[r_idx][c_idx] - b.matrix[r_idx][c_idx]) < 10e-5, "different values at %s %s : %s vs %s"%(r_idx, c_idx, genie_out[r_idx][c_idx], b.matrix[r_idx][c_idx])
	


	#print print_CPT(or_dist)
	#print max_based(or_dist)
	#for i in max_based(or_dist, leakdef = LEAKDEF.HENRION, handle_errors=True):
	#	print i

if __name__ == "__main__":
	test()
