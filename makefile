#run for single network
run:
	python main.py $(file)
#run for multiple networks
test:
	python test.py data/$(n)n/$(s)/*.txt
#gathering data
gather:
	python gather.py output/$(n)n/$(s)/*.txt
#test and gather:
# usage: make tg n=5 s=10k > summary10k.txt -s
tg:
	python test.py data/$(n)n/$(s)/*.txt
	python gather.py output/$(n)n/$(s)/*.txt
	
create_summary:
	make tg n=5 s=1k > summary1k -s
	make tg n=5 s=2k > summary2k -s
	make tg n=5 s=5k > summary5k -s
	make tg n=5 s=10k > summary10k -s
