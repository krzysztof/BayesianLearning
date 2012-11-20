# Learning of Canonical networks
# author: Krzysztof Nowak
from collections import Counter

import sys

from bayesiandataset import BayesianDataSet

def compareNoisyOR():
	"""
	Network CancerOR 10k records
	"""
	real = (
		(('True', 'Smoker', 'True'), 0.61),
		(('True', 'Genetic', 'True'), 0.25),
		(('True', 'CoalWorker', 'True'), 0.15),
		(('True', 'BadDiet', 'True'), 0.04),
	)
	b1 = BayesianDataSet("../data/5n/10k/Network1.txt")
	b1.addChildNode(4,'False', [0,1,2,3], ['False', 'False', 'False', 'False'])
	print "COUNTED"
	for p in b1.countForChild(4):
		print p
	print "SOLVED"
	for p in b1.solveForChild(4):
		print p
	print "REAL"
	for p in real:
		print p

def compareNoisyMAX():
	"""
	Network CancerMax 20k records
	"""
	real = (
		(('Malignant', 'Smoker', 'TwoPacks'), 0.23),
		(('Malignant', 'Smoker', 'OnePack'), 0.11),
		(('Malignant', 'Genetic', 'True'), 0.25),
		(('Malignant', 'CoalWorker', 'True'), 0.15),
		(('Malignant', 'BadDiet', 'Good'), 0.04),
		(('Malignant', 'BadDiet', 'Bad'), 0.13),

		(('Benign', 'Smoker', 'TwoPacks'), 0.25),
		(('Benign', 'Smoker', 'OnePack'), 0.14),
		(('Benign', 'Genetic', 'True'), 0.55),
		(('Benign', 'CoalWorker', 'True'), 0.66),
		(('Benign', 'BadDiet', 'Good'), 0.63),
		(('Benign', 'BadDiet', 'Bad'), 0.22),
	)
	b1 = BayesianDataSet("../data/MAX/Network20k.txt")
	print b1
	b1.addChildNode(4,'No', [0,1,2,3], ['False', 'False', 'False', 'Medium'])
	print "COUNTED"
	for p in b1.countForChild(4):
		print p
	print "SOLVED"
	for p in b1.solveForChild(4):
		print p

	print "REAL"
	for p in real:
		print p

def compareNoisyMaxSimple():
	"""
	Network CancerMax_simple 10k
	"""
	real = (
		(('Malignant', 'Smoker', 'TwoPacks'), 0.23),
		(('Malignant', 'Smoker', 'OnePack'), 0.11),
		(('Malignant', 'Genetic', 'True'), 0.25),
		(('Benign', 'Smoker', 'TwoPacks'), 0.25),
		(('Benign', 'Smoker', 'OnePack'), 0.14),
		(('Benign', 'Genetic', 'True'), 0.55),
	)
	b1 = BayesianDataSet("../data/MAX/CancerMAX_simple10k.txt")
	print b1
	b1.addChildNode(2,'No', [0,1], ['False', 'False'])
	print "COUNTED"
	for p in b1.countForChild(2):
		print p
	print "SOLVED"
	for p in b1.solveForChild(2):
		print p

	print "REAL"
	for p in real:
		print p

def printSomething():
	#b1 = BayesianDataSet("../data/NoisyOR_100k.txt")
	b1 = BayesianDataSet("../data/5n/10k/Network1.txt")
	b1.addChildNode(4,'False', [0,1,2,3], ['False', 'False', 'False', 'False'])
	child_node = b1.children[4]

	parent_columns = tuple([tuple([b1.data[i][j] for j in [0,1,2,3]]) for i in range(len(b1.data))])


	parent_child_columns = tuple([tuple([b1.data[i][j] for j in [0,1,2,3,4]]) for i in range(len(b1.data))])
	parent_child_counts = Counter(parent_child_columns)

	parent_counts_s = sorted(Counter(parent_columns).items(), key=lambda i:i[1], reverse=True)
	encoded_parent_counts_s = [(b1._encode_equation(equation, child_node), count) for equation, count in parent_counts_s if equation!=child_node.parent_char_states]
	for eq, val in encoded_parent_counts_s:
		#print eq,val
		# 0 at the end because state 0 is True
		nom = parent_child_counts[b1._decode_equation(tuple(eq), child_node)+(0,)]
		#denom = b1._decode_equation(tuple(eq))
		denom = val
		#print eq,val, parent_child_counts[nom],parent_child_counts[b1._decode_equation(tuple(eq)+(0,),child_node)], parent_child_counts[b1._decode_equation(tuple(eq)+(1,),child_node)]/ float(val)
		print eq, "n:%s d:%s, val:%s"%(nom ,denom, float(nom)/denom)

def main():
	#printSomething()
	#compareNoisyOR()
	#compareNoisyMaxSimple()
	compareNoisyMAX()


if __name__ == "__main__":
	main()
