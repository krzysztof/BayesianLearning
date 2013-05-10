CC = g++
SMILE = smile
SMILEARN = smilearn
CFLAGS = -O3 -ffast-math

.PHONY: generator learner test clean

generator: 
	$(CC) $(CFLAGS) generator.cpp -o generator -l$(SMILE) -l$(SMILEARN)
learner: 
	$(CC) $(CFLAGS) learner.cpp -o learner -l$(SMILE) -l$(SMILEARN)
test:
	$(CC) $(CFLAGS) test.cpp -o test -l$(SMILE) -l$(SMILEARN)
test_run:
	./test 10 ../src/CancerOR.xdsl
clean:
	rm *.o
