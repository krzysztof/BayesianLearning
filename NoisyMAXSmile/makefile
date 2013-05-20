CC = g++
SMILE = smile
SMILEARN = smilearn
CFLAGS = -O3 -ffast-math
#CFLAGS =   

.PHONY: all generator smile_learner naive_learner test clean

all: generator smile_learner naive_learner
	
generator: 
	$(CC) $(CFLAGS) generator.cpp -o generator -l$(SMILE) -l$(SMILEARN)
smile_learner: 
	$(CC) $(CFLAGS) smile_learner.cpp -o smile_learner -l$(SMILE) -l$(SMILEARN)
naive_learner: 
	$(CC) $(CFLAGS) naive_learner.cpp -o naive_learner -l$(SMILE) -l$(SMILEARN)
test:
	$(CC) $(CFLAGS) test.cpp -o test -l$(SMILE) -l$(SMILEARN)
test_run:
	./test 10 ../src/CancerOR.xdsl
clean:
	rm *.o
