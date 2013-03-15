import numpy as np
from math import log,sqrt,e

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
	plt.xticks(ind+width, ('20k', '5k', '2k', '1k') )

	plt.legend( (rects1[0], rects2[0]), ('Counted', 'Solved') , loc='upper left')

	def autolabel(rects):
		# attach some text labels
		for rect in rects:
			height = rect.get_height()
			plt.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%.3f'%float(height),
				ha='center', va='bottom')

	autolabel(rects1)
	autolabel(rects2)

	plt.show()
	#folder_name = sys.argv[1]
	#plt.savefig('graphs/%s/%s.png'%(folder_name, name))
	#plt.clf()

def main():
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

def euclidian_dist(a,b):
	assert len(a) == len(b)
	sum = 0.0;
	for i in range(len(a)):
		sum+= (a[i] - b[i])**2
	return sqrt(sum)

def kl_dist(P, Q):
	assert len(P) == len(Q)
	suma = 0.0;
	for i in range(len(P)):

		q = Q[i]
		if Q[i] < 10e-7:
			q = 10e-7
		tmp =  P[i]*log((P[i]/q),e)
		print tmp,i

		suma += tmp
	return suma

def hellinger_dist(P, Q):
	assert len(P) == len(Q)
	assert all((Q[i]>=0 for i in xrange(len(Q))))
	suma = sum(( (sqrt(P[i]) - sqrt(Q[i]))**2 for i in range(len(P)) ))
	return sqrt(suma) * 0.7071067811865475  # * 1./sqrt(2)


def manual():
	REAL = [0.11,0.23,0.25,0.15,0.04, 0.13]

	C_20k = [0.11490125673249552, 0.2535211267605634,
			 0.273972602739726, 0.18181818181818182,
			 0.04388274530454625, 0.1255203330131284]
	S_20k =[ 0.11692088743870532, 0.22340595108662775,
			 0.22217053415061294, 0.11924477984018245,
			 0.048839931330187047, 0.13005425221980049]

	C_5k = [0.1732283464566929, 0.2153846153846154,
			0.3025210084033613, 0.3333333333333333,
			0.042722664735698766, 0.13003901170351106]
	S_5k = [0.13139681549402915, 0.26807551858501466,
			0.23158693046741452, 0.041916026484861479,
			0.042722664735698745, 0.13003901170351106]

	C_2k = [ 0.09803921568627451, 0.4230769230769231,
			 0.3829787234042553, 0.25,
			 0.04318936877076412,  0.15894039735099338]

	S_2k =[0.13088450292397658, 0.20944622507122512,
		   0.23057259713701439, 0.12905092592592593,
		   0.043189368770764069, 0.15894039735099341]
	C_1k = [ 0.04878048780487805, 0.1875,
			 0.12,  0.0,
			 0.03728813559322034,0.14074074074074075]
	S_1k =[ 0.055697823303457183, 0.33059467918622853,
			0.23293607800650051,  0.22095070422535212,
			0.037288135593220306, 0.14074074074074072]
	print kl_dist(REAL, C_2k)

	sets = [[C_20k, S_20k], [C_5k, S_5k], [C_2k, S_2k], [C_1k, S_1k]]

	eukl_dist = []
	for c,s in sets:
		pair = euclidian_dist(REAL, c), euclidian_dist(REAL, s)
		eukl_dist.append(pair)

	kl_dist_v = []
	for c,s in sets:
		pair = kl_dist(REAL,c), kl_dist(REAL,s)
		kl_dist_v.append(pair)

	hil_dist_v = []
	for c,s in sets:
		pair = hellinger_dist(REAL,c), hellinger_dist(REAL,s)
		hil_dist_v.append(pair)

	a = zip(*eukl_dist)
	b = zip(*kl_dist_v)
	c =  zip(*hil_dist_v)
	print C_20k, S_20k,REAL
	show('Euclidian Distnace', a[0], a[1])
	show('KL Divergence', b[0], b[1])
	show('Hellinger distance', c[0], c[1])

if __name__=='__main__':
	#main() #collects
	manual()

