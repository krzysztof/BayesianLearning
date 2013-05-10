from utilities import read_data, match_by_column, binarize, ZERO_PROB_EPS, ONE_PROB_EPS, NO_PROB_EPS, max_based, LEAKDEF
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
				
	for comb, prob in max_based(params, leakdef = LEAKDEF.HENRION, handle_errors = True):
		print comb, prob

if __name__ == "__main__":
	main(sys.argv[1],"LungCancer")
	
