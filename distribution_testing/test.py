import sys
from collections import Counter, defaultdict

def main():
	dane = open(sys.argv[1],'r').readlines()
	print dane[0]
	dane = dane[1:]

	c = Counter(dane)
	clear_data = sorted(c.items(), key=lambda x:x[0])
	clear_data = [ (k.replace('t','1').replace('f','0')[:-2], v) for k,v in clear_data]
	clear_data_par = [ (" ".join(k.split(" ")[:-1]), v) for k,v in clear_data]
	d = defaultdict(lambda: 0)
	for i in clear_data_par:
		d[i[0]] += i[1]
	idx = 1
	for	k,v in sorted(d.items(), key=lambda x:x[1], reverse=True):
		print "%d\t|\t%s\t|%d\t|%.4f%%" % (idx, k, v, v/float(len(dane)))
		idx+=1
	#for k, v in clear_data:
		#print "::%s::" % k[:-2]
		#print "%s | %4d | %.4f%%" % (k[:-2], v, v/float(len(dane)))
if __name__ == "__main__":
	main()
