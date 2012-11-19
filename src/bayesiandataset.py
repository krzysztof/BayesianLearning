from itertools import combinations
from math import e,log, sqrt
from collections import Counter
import random

import numpy as np
from string import zfill

class ChildNode(object):
	"""
	Structure containing information on ChildNode.
	"""
	def __init__(self, column_index, char_state, parent_indices, parent_char_states):
		"""
		@type column_index: int
		@param column_index: column index of child node

		@type parent_indices: tuple
		@param parent_indices: column indices of direct parent nodes. Tuple of integers.

		@type char_state: int
		@param char_state: characteristic state of child node (used to determine "off" state in NoisyMax)

		@type parent_char_states: tuple
		@param parent_char_states: same as above, but for each of the parents. Tuple of integers.
		"""

		self.column_index = column_index
		self.char_state = char_state
		self.parent_indices = parent_indices
		self.parent_char_states = parent_char_states

	def __str__(self):
		return "Col:%s CharState:%s Parents:%s ParCharStates%s" % (self.column_index, self.char_state, self.parent_indices, self.parent_char_states )


class BayesianDataSet(object):
	"""
	Structure containing information on data generated from bayesian network.
	self.data = data in raw form (2 dimensional tuple)
	"""

	def __init__(self, filename, delimiter=' '):
		"""
		Create BayesianDataSet from a filename.
		Assumes the first line to be column headers (names), and remaining being the records themselves.
		"""

		self.raw_data = tuple([ tuple(line.strip().split(delimiter)) for line in open(filename).readlines() ])
		self.headers = self.raw_data[0]
		self.data = self.raw_data[1:]

		# create a domain
		# e.g. self.domain[0] = ['State0', 'State1', 'State2']
		self.domain = []
		for i in range(len(self.headers)):
			self.domain.append(tuple( set([ self.data[j][i] for j in range(len(self.data)) ]) ))
		self.domain = tuple(self.domain)

		#string to data mapping:
		#self._str_to_num['Column1'] = {'State1':0, 'State2':1}
		self._str_to_num = {}
		for i in range(len(self.headers)):
			self._str_to_num[self.headers[i]] = dict( zip(self.domain[i], range(len(self.domain[i]))) )

		#turn string data into numbers
		number_data = []
		for i in range(len(self.data)):
			number_data.append([ self._str_to_num[self.headers[j]][self.data[i][j]] for j in range(len(self.data[i])) ])

		self.data = number_data

		# dictionary of child nodes information
		# eg. self.children[0] = ChildNode()
		self.children = {}

	def addChildNode(self, column_index, char_state_name, parent_indices, parent_char_states_names):

		char_state = self._str_to_num[self.headers[column_index]][char_state_name]

		parent_char_states = tuple([self._str_to_num[self.headers[parent_indices[i]]][parent_char_states_names[i]] for i in range(len(parent_indices)) ])

		child_node = ChildNode(column_index, char_state, parent_indices, parent_char_states)
		self.children[column_index] = child_node

	def _encode_equation(self, equation, child_node):
		# sizes of each domain
		lengths = [len(self.domain[i]) for i in child_node.parent_indices]

		#each domain values without characteristic states
		sizes = [range(len(self.domain[i])) for i in child_node.parent_indices]
		for i in range(len(sizes)):
			sizes[i].remove(child_node.parent_char_states[i])

		final = []
		for i in range(len(equation)):
			if equation[i] == child_node.parent_char_states[i]:  # append zeros if it's recognized as char state
				final.extend([0]*(lengths[i]-1))
			else:  # encode into number
				final.extend([int(x) for x in list(zfill(bin(2**(sizes[i].index(equation[i])))[2:],lengths[i]-1))])
		return final

	def _decode_equation(self, equation, child_node):
		lengths = [len(self.domain[i]) for i in child_node.parent_indices]
		sizes = [range(len(self.domain[i])) for i in child_node.parent_indices]
		for i in range(len(sizes)):
			sizes[i].remove(child_node.parent_char_states[i])

		tokenized = []
		tmp = 0
		for l in lengths:
			tokenized.append(equation[tmp:tmp+l-1])
			tmp+=l-1
		final = []
		#print equation, tokenized,sizes
		for i in range(len(tokenized)):
			if 1 not in tokenized[i]:
			#if tuple([0]*len(tokenized[i])) == tokenized[i]:
				final.append(child_node.parent_char_states[i])
			else:
				final.append(sizes[i][len(tokenized[i]) - tokenized[i].index(1)-1])
		return tuple(final)

	def _equation_set_valid(self, equation_set):
		"""
		@param equation_set: [((0,0,1,0,1),2341), ...] (of size k)
		"""
		vectors = [eq[0] for eq in equation_set]
		A = np.array(vectors)
		detA = np.linalg.det(A)
		return detA > 10e-7

	def _choose_equations(self, equations, k, i):
		"""
		@param equations: [((0,0,1,0,1),2341), ...]
		@type equations: list
		@param k: number of equations to pick
		@type k: int
		@param i: column in encoded equation (parameter) to choose for
		@type i: int
		"""

		## RANDOM TEST
		#choice = random.sample(equations,k)
		#while(not self._equation_set_valid(choice)):
		#	choice = random.sample(equations, k)
		#return choice
		## RANDOM TEST

		## CHECKING COMBINATIONS
		#t = k
		#chosen = equations[:k]
		#found = False
		#while(t <= len(equations) and not found):
		#	for comb in combinations(equations[:t], k):
		#		if self._equation_set_valid(comb):
		#			chosen = comb
		#			found = True
		#			break;
		#	t+=1
		#return chosen

		## CHECKING COMBINATIONS
		def comparator(a,b):
			if a[0][i] > b[0][i]:
				return -1
			elif a[0][i] < b[0][i]:
				return 1
			elif a[1] > b[1]:
				return -1
			elif a[1] < b[1]:
				return 1
			return 0
		equations2 = list(equations)
		equations2.sort(comparator)
		for eq  in equations2:
			print eq

	def _solve_prod_equation(self, equation_set, leak):
		"""
		@param equation_set: [ ((0,1,0,0,1,0), 0.51), ((0,1,0,0,0,0),0.23), ...]
		@return: solution as a list of floats [0.31, 0.33, 0.21...]
		"""
		a = [equation for equation, value in equation_set]
		b = [log(1.0 - (1.0 - (1.0-value)*(1.0-leak)**(sum(equation)-1))) for equation, value in equation_set]

		a = zip(*a)

		A = np.array(a)
		detA = np.linalg.det(A)
		Y = []
		for i in range(len(a)):
			p = np.array(a[:i] + [b] + a[i+1:])
			detP = np.linalg.det(p)
			Y.append(1.0-e**(detP/detA))
		return Y

	def _solve_prod_equation_single(self, equation_set, leak, param):
		"""
		@param equation_set: [ ((0,1,0,0,1,0), 0.51), ((0,1,0,0,0,0),0.23), ...]
		@return: solution as a list of floats [0.31, 0.33, 0.21...]
		"""
		a = [equation for equation, value in equation_set]
		b = [log(1.0 - (1.0 - (1.0-value)*(1.0-leak)**(sum(equation)-1))) for equation, value in equation_set]

		a = zip(*a)

		A = np.array(a)
		detA = np.linalg.det(A)
		p = np.array(a[:param] + [b] + a[param+1:])
		detP = np.linalg.det(p)
		return 1.0-e**(detP/detA)

	def _get_parameter(self, equation, child_node):
		"""
		@param equation: tuple with trailing one (0,..,1,...0)
		@return ("Node0", "State0")
		"""
		lengths = [len(self.domain[i]) for i in child_node.parent_indices]
		sizes = [range(len(self.domain[i])) for i in child_node.parent_indices]
		for i in range(len(sizes)):
			sizes[i].remove(child_node.parent_char_states[i])

		tokenized = []
		tmp = 0
		for l in lengths:
			tokenized.append(equation[tmp:tmp+l-1])
			tmp+=l-1
		#print equation, tokenized,sizes
		decoded_equation = self._decode_equation(equation, child_node)
		for i in range(len(tokenized)):
			if 1 in tokenized[i]:
				node_name= self.headers[child_node.parent_indices[i]]
				state_name = self.domain[child_node.parent_indices[i]][decoded_equation[i]]
				return (node_name, state_name)

	def countForChild(self, child_column_index):
		assert child_column_index in self.children, "Children of index %d not set!" % ( child_column_index )
		child_node = self.children[child_column_index]

		parent_columns = tuple([tuple([self.data[i][j] for j in child_node.parent_indices]) for i in range(len(self.data))])
		parent_child_columns = tuple([tuple([self.data[i][j] for j in child_node.parent_indices + [child_node.column_index]]) for i in range(len(self.data))])

		parent_counts = Counter(parent_columns)
		parent_child_counts = Counter(parent_child_columns)

		non_char_states = range(len(self.domain[child_node.column_index]))
		non_char_states.remove(child_node.char_state)
		K = sum([len(self.domain[i])-1 for i in child_node.parent_indices])

		final_output = []
		for state in non_char_states:
			for i in range(K):
				param = self._decode_equation([0]*i+[1]+[0]*(K-1-i), child_node)
				nominator = parent_child_counts[param + (state,)]
				denom = parent_counts[param]

				#print counts of each equation
				print [0]*i+[1]+[0]*(K-1-i), denom

				value = float(nominator) / denom
				param_name = self._get_parameter([0]*i+[1]+[0]*(K-1-i), child_node)
				#print ((self.domain[child_node.column_index][state],)+param_name), value,denom
				final_output.append( (((self.domain[child_node.column_index][state],)+param_name), value, ) )
		return final_output

	def solveForChild(self, child_column_index):
		assert child_column_index in self.children, "Children of index %d not set!" % ( child_column_index )
		child_node = self.children[child_column_index]

		# choose only columns from parent indices
		parent_columns = tuple([tuple([self.data[i][j] for j in child_node.parent_indices]) for i in range(len(self.data))])

		# choose only columns from parent indices and child node
		parent_child_columns = tuple([tuple([self.data[i][j] for j in child_node.parent_indices + [child_node.column_index]]) for i in range(len(self.data))])

		# prepare counters for both iterables
		parent_counts = Counter(parent_columns)
		parent_child_counts = Counter(parent_child_columns)

		# counts of parent combinations
		parent_counts_s = sorted(parent_counts.items(), key=lambda i:i[1], reverse=True)

		# calculate LEAKS
		#print child_node.char_state, self.domain[child_node.column_index]
		non_char_states = range(len(self.domain[child_node.column_index]))
		non_char_states.remove(child_node.char_state)

		LEAK = {}
		denom = parent_counts[child_node.parent_char_states]
		for state in non_char_states:
			nominator = parent_child_counts[child_node.parent_char_states + (state,)]
			assert denom > 0, "DENOMINATOR IN COUNTING LEAK IS 0"
			LEAK[state] = float(nominator)/float(denom)

		#number of parameters in data set
		K = sum([len(self.domain[i])-1 for i in child_node.parent_indices])

		# counts of encoded parent combinations minus the LEAK denominator (0,0,0,...0)
		encoded_parent_counts_s = [(self._encode_equation(equation, child_node), count) for equation, count in parent_counts_s if equation!=child_node.parent_char_states]
		#for eq, cnt in encoded_parent_counts_s:
		#	print eq, cnt

		final_output = []
		for state in non_char_states:
			#equations with float values calculated on the right side
			solution = []
			for i in range(K):
				# TODO: turn into SOLVE function (both choosing and solving) because it's easier to solve non KxK matrices (e.g. 2x2 if it's better)
				chosen_equations = self._choose_equations(encoded_parent_counts_s, K, i)
				## PRINT CHOSEN EQUATIONS BELOW
				for eq, cnt in chosen_equations:
					print eq,cnt

				chosen_equations_decoded = []
				for eq, cnt in chosen_equations:
					chosen_equations_decoded.append((self._decode_equation(eq,child_node), cnt,))

				chosen_equations_parametrized = []
				for j in range(len(chosen_equations_decoded)):
					nominator = float(parent_child_counts[tuple(chosen_equations_decoded[j][0])+(state,)])
					denominator = float(chosen_equations_decoded[j][1])
					chosen_equations_parametrized.append((chosen_equations[j][0], nominator/denominator))
				solution.append(self._solve_prod_equation_single(chosen_equations_parametrized, LEAK[state], i))

			for i in range(K):
				param = self._get_parameter([0]*i+[1]+[0]*(K-1-i), child_node)
				final_output.append( (((self.domain[child_node.column_index][state],)+param), solution[i]) )
			#print solution
			#for eq, cnt in chosen_equations_parametrized:
				#print "state%s: %s %s"%(state,eq,cnt)
			return final_output


	def __str__(self):
		return "BDS: columns: %d, rows: %d, header: %s, domain: %s" %(len(self.headers), len(self.data), self.headers, self.domain)
