from collections import Counter

from bayesiandataset import BayesianDataSet, ChildNode

def main():
	filename = "../data/Testing100k.txt"
	file_set = open(filename,'r')
	#bn0 = BayesianDataSet(filename)

	#bn0.addChildNode(4,'No', [0,1,2,3], ['False', 'False', 'False', 'Medium'])
	#child_node0 = bn0.children[4]

	#line =
	cnt = Counter([line[:-2] for line in file_set.readlines()])

	#(('TwoPacks', 'OnePack', 'False'), ('False', 'True'),
	# ('True', 'False'), ('Bad', 'Medium', 'Good'), ('Benign', 'Malignant', 'No'))

#	print bn0.domain

	# 0 0 0 0 1 0  -># 5697
	a = cnt["True False False False True"]
	b = cnt["True False False False False"]
	print float(a)/(a+b)
	#a1 = cnt["False False False Good Benign"] # 3650
	#a2 = cnt["False False False Good Malignant"] # 250
	#a3 = cnt["False False False Good No"] # 1797
	#a0 = float(sum([a1,a2,a3]))


	## 1 0 0 0 1 0  -># 1895
	#b1 = cnt["OnePack False False Good Benign"] # 1139
	#b2 = cnt["OnePack False False Good Malignant"] # 295
	#b3 = cnt["OnePack False False Good No"] # 461
	#b0 = float(sum([b1,b2,b3]))

	#print a1,a1/a0, a2,a2/a0, a3,a3/a0, a0
	#print b1,b1/b0, b2,b2/b0, b3,b3/b0, b0


if __name__ == "__main__":
	main()