#include <iostream>
#include <fstream>
#include "smile/smile.h"
#include <stdlib.h>
#include <string>
#include <vector>
#include <time.h>
#include <map>

using namespace std;

//void printCPT(DSL_node *node);

typedef vector<pair<string, double> > vpsd;
typedef vector<string> vs;

double fRand(double fMin, double fMax);


string getRandom(double r, vpsd &items);

typedef map<string, string> keyMap;
typedef map<string, double> valMap;
typedef map< keyMap, valMap > cptMap;
cptMap get_cptmap(DSL_node *node);

void generateDatafile(int n, string networkName, string outFile);
void randomizeNetwork(string networkName, string networkOutputName);

/*
	./generator number_of_records xdls_network_file [seed - optional]
*/
//int main(int argc, char* argv[]) {
//	if(argc > 3) {
//		srand(atoi(argv[3]));
//	} else {
//		srand(time(NULL));
//	}
//	int number_of_records = atoi(argv[1]);
//	string network_model_name = string(argv[2]);
//	generateDatafile(number_of_records, network_model_name);
//
//	return 0;
//}

