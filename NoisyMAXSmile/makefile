CC = g++
SMILE = smile
SMILEARN = smilearn
CFLAGS = -O3 -ffast-math -Wall
#CFLAGS =   

.PHONY: test Experiment GaussJordan.o Learners.o Generator.o

Experiment: Experiment.cpp Learners.o Generator.o
	$(CC) $(CFLAGS) Experiment.cpp Learners.o GaussJordan.o Generator.o -o Experiment -l$(SMILE) -l$(SMILEARN)

GaussJordan.o: GaussJordan.cpp GaussJordan.h
	$(CC) $(CFLAGS) -c GaussJordan.cpp

Learners.o: GaussJordan.o Learners.h
	$(CC) $(CFLAGS) -c Learners.cpp

Generator.o: Generator.cpp Generator.h
	$(CC) $(CFLAGS) -c Generator.cpp -l$(SMILE) -l$(SMILEARN)

test: GaussJordan.o
	#$(CC) $(CFLAGS) test.cpp -o test -l$(SMILE) -l$(SMILEARN)
	$(CC) $(CFLAGS) testGJ.cpp GaussJordan.o -o testGJ -l$(SMILE) -l$(SMILEARN)
	./testGJ < testGJ.in
clean:
	rm *.o
