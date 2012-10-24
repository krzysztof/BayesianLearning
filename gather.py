import sys

datas = []
for file in sys.argv[1:]:
	data = {}
	for line in list(open(file).readlines()[3:]):
		key, value = line.split(' ')
		data[key] = float(value)
	datas.append(data)

N = len(datas)
keys = datas[0].keys()

average_data = {}
for key in keys:
	val = 0.0
	for data in datas:
		try:
			val+=data[key]
		except KeyError:
			print data
	val /= N
	average_data[key]=val

items = sorted(average_data.items(), key=lambda i:i[0])
for k,v in items:
	print k,v
