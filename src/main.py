# Learning of Canonical networks
# author: Krzysztof Nowak

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
	b1 = BayesianDataSet("../data/MAX/Network20k.txt")
	print b1
	b1.addChildNode(4,'No', [0,1,2,3], ['False', 'False', 'False', 'Medium'])
	print "COUNTED"
	for p in b1.countForChild(4):
		print p
	print "SOLVED"
	for p in b1.solveForChild(4):
		print p

def compareNoisyMaxSimple():
	"""
	Network CancerMax_simple 10k
	"""
	b1 = BayesianDataSet("../data/MAX/CancerMAX_simple10k.txt")
	print b1
	b1.addChildNode(2,'No', [0,1], ['False', 'False'])
	print "COUNTED"
	for p in b1.countForChild(2):
		print p
	print "SOLVED"
	for p in b1.solveForChild(2):
		print p


def main():
	compareNoisyOR()


if __name__ == "__main__":
	main()
