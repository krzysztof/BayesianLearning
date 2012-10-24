import numpy as np
import sys

def generate(data):
	for k,v in data.iteritems():
		AVals, BVals = v
		show(k, AVals, BVals)

def show(name, AVals, BVals):
	import matplotlib.pyplot as plt
	N = len(AVals)
	#AVals = (0.03,0.05,0.07,0.023)
	#AErr = (0.02,0.03,0.05,0.022)

	ind = np.arange(N)  # the x locations for the groups
	width = 0.4       # the width of the bars


	plt.subplot(111)
	rects1 = plt.bar(ind, AVals, width,
		color='g',
		error_kw=dict(elinewidth=2, ecolor='blue'))

	#BVals = (0.03,0.05,0.07,0.023)
	#BErr = (0.02,0.03,0.05,0.022)
	rects2 = plt.bar(ind+width, BVals, width,
		color='r',
		error_kw=dict(elinewidth=2, ecolor='blue'))

	# add some
	plt.ylabel(name)
	plt.title(name)
	plt.xticks(ind+width, ('10k', '5k', '2k', '1k') )

	plt.legend( (rects1[0], rects2[0]), ('New', 'Old') , loc='upper left')

	def autolabel(rects):
		# attach some text labels
		for rect in rects:
			height = rect.get_height()
			plt.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%.3f'%float(height),
				ha='center', va='bottom')

	autolabel(rects1)
	autolabel(rects2)

	#plt.show()
	folder_name = sys.argv[1]
	plt.savefig('graphs/%s/%s.png'%(folder_name, name))
	plt.clf()

if __name__=='__main__':
	data = {}
	#sizes = ['10k', '5k', '2k', '1k']
	sizes = ['10k', '5k', '2k', '1k']
	for size in sizes:
		data[size] = {}
		#file = open(sys.argv[1]) # should take sizes into account
		file = open("summary%s"%size)
		for line in file.readlines():
			#print line
			k,v = line.split(' ')
			data[size][k] = float(v[:-1])
	#for k,v in data['10k'].iteritems():
	#data['10k']['eucl'] = ((1,2,3,4), (4,5,6,7))
	ordered_data = {}
	keys = data['10k'].keys()
	for key in keys:
		ordered_data[key] = [data[size][key] for size in sizes]
	uniq_keys = list(set([key[4:] for key in keys]))
	final_data = {}
	for key in uniq_keys:
		final_data[key] = [ordered_data['new_'+key], ordered_data['old_'+key]]
	generate(final_data)
