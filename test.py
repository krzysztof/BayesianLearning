
import sys, os
for file in sys.argv[1:]:
	name = "/".join(file.split('/')[-3:])
	os.system("python main.py %s > %s" % (file, "output/"+name))
    
#python main.py data/10k/Network$i.txt > output/10k/Network$i.txt
