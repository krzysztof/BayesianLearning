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
make_experiment:
	./Experiment Test ../src/OR_2.xdsl tmp_results2_coeffs
	./Experiment Test ../src/OR_3.xdsl tmp_results3_coeffs
	./Experiment Test ../src/OR_4.xdsl tmp_results4_coeffs
	./Experiment Test ../src/OR_5.xdsl tmp_results5_coeffs
	./Experiment Test ../src/OR_6.xdsl tmp_results6_coeffs
	./Experiment Test ../src/OR_7.xdsl tmp_results7_coeffs
	./Experiment Test ../src/OR_8.xdsl tmp_results8_coeffs
	./Experiment Test ../src/OR_9.xdsl tmp_results9_coeffs
	./Experiment Test ../src/OR_10.xdsl tmp_results10_coeffs
clean:
	rm *.o
