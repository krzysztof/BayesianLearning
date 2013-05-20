from utilities import read_data, match_by_column, binarize, ZERO_PROB_EPS, ONE_PROB_EPS, NO_PROB_EPS, LEAKDEF, CPT
from collections import Counter, defaultdict
import sys

def main(filename, child_name):
	data = read_data(filename, ' ')
	c = Counter(data['data'])

	n = data['header'].index(child_name)

	new_counter = match_by_column(c, n)

	binary_data = binarize(new_counter)

	items = sorted(binary_data.items(), key=lambda x: x[1][2], reverse=True)
	items_dict = dict(items)

	leak_comb = tuple([0,]*n)
	params = []
	for comb in list(tuple( [0,]*i + [1,] + [0,]*(n-1-i)) for i in xrange(n)) + [leak_comb]:
		#print comb, items_dict[comb]
		if comb not in items_dict:
			params.append(NO_PROB_EPS)
		#elif items_dict[comb][0]==0.0:
		#	params.append(ZERO_PROB_EPS)
		#elif items_dict[comb][0]==1.0:
		#	params.append(ONE_PROB_EPS)
		else:
			params.append(items_dict[comb][0])
	
	leak = params[-1]
	params = reduce( lambda x,y: x+y, [[a,0] for a in params[:-1]]) + [leak,]
	parent_dims = [2]*n
	labels = []
	for h in data['header']:
		labels.append(["True", "False"])
		
	naive_or = CPT([params, [1.0 - p for p in params]], parent_dims, CPT.TYPE_NOISY_MAX, data['header'], labels)
	naive_cpt = naive_or.max_based()

	print naive_cpt.print_raw()

if __name__ == "__main__":
	#main(sys.argv[1],"LungCancer")
	main(sys.argv[1],"C1")
	
