#include <iostream>
#include <fstream>
#include "smile/smile.h"
#include <stdlib.h>
#include <string>
#include <vector>
#include <time.h>
#include <map>

using namespace std;

void printCPT(DSL_node *node) {
	DSL_network* net = node->Network(); // node network
	int handle = node->Handle();
	DSL_nodeDefinition *def = node->Definition();
	const DSL_Dmatrix &cpt = *def->GetMatrix();
	const DSL_idArray &outcomes = *def->GetOutcomesNames();
	const DSL_intArray &parents = net->GetParents(handle);
	int parentCount = parents.NumItems();

	DSL_intArray coords;
	for (int elemIdx = 0; elemIdx < cpt.GetSize(); elemIdx ++) {
		cpt.IndexToCoordinates(elemIdx, coords);
		cout << "P(" << node->GetId() << " = " << outcomes[coords[parentCount]] << " | ";
		for (int parentIdx = 0; parentIdx < parentCount; parentIdx ++) {
			 if (parentIdx > 0) cout << ", ";
			 DSL_node *parentNode = net->GetNode(parents[parentIdx]);
			 const DSL_idArray &parentStates = *parentNode->Definition()->GetOutcomesNames();
			 cout << parentNode->GetId() << " = " << parentStates[coords[parentIdx]];
		}
		cout << ") = " << cpt[elemIdx] << endl;
	}
}
typedef vector<pair<string, double> > vpsd;
typedef vector<string> vs;

double fRand(double fMin, double fMax) {
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

string getRandom(double r, vpsd &items) {
	double suma = 0.0;
	for(int i=0 ;i< items.size() ; ++i ) {
		if( (items[i].second + suma) >= r) {
			return items[i].first;
		}
		suma += items[i].second;
	}
	return items[items.size()-1].first;
}


typedef map<string, string> keyMap;
typedef map<string, double> valMap;
typedef map< keyMap, valMap > cptMap;
cptMap get_cptmap(DSL_node *node) {
	DSL_network* net = node->Network(); // node network
	int handle = node->Handle();
	DSL_nodeDefinition *def = node->Definition();
	const DSL_Dmatrix &cpt = *def->GetMatrix();
	const DSL_idArray &outcomes = *def->GetOutcomesNames();
	const DSL_intArray &parents = net->GetParents(handle);
	int parentCount = parents.NumItems();

	DSL_intArray coords;
	cptMap cptmap;
	for (int elemIdx = 0; elemIdx < cpt.GetSize(); elemIdx ++) {
		cpt.IndexToCoordinates(elemIdx, coords);
		keyMap km;
		for (int parentIdx = 0; parentIdx < parentCount; parentIdx ++) {
			DSL_node *parentNode = net->GetNode(parents[parentIdx]);
			const DSL_idArray &parentStates = *parentNode->Definition()->GetOutcomesNames();
			km[string(parentNode->GetId())] = string(parentStates[coords[parentIdx]]);
		}
		if (cptmap.count(km) == 0) {
			valMap vm;
			vm[string(outcomes[coords[parentCount]])] = cpt[elemIdx];
			cptmap[km] = vm;
		} else {
			cptmap[km][string(outcomes[coords[parentCount]])] = cpt[elemIdx];
		}
	}

	return cptmap;
}

void generateDatafile(int n, string inname) {
	DSL_network theNet;
	theNet.ReadFile(inname.c_str());
	int cancer = theNet.FindNode("LungCancer");
	DSL_node *node = theNet.GetNode(cancer);

	DSL_network* net = node->Network(); // node network
	int handle = node->Handle();
	DSL_nodeDefinition *def = node->Definition();
	const DSL_Dmatrix &cpt = *def->GetMatrix();
	const DSL_idArray &outcomes = *def->GetOutcomesNames();
	const DSL_intArray &parents = net->GetParents(handle);
	int parentCount = parents.NumItems();

	DSL_intArray coords;

	double leak = 0.01;

	// PRINT HEADER
		for (int parentIdx = 0; parentIdx < parentCount; parentIdx ++ ) {
			DSL_node *parentNode = net->GetNode(parents[parentIdx]);
		cout << parentNode->GetId() << " ";			 
	}
	cout << node->GetId() << endl;

	cptMap cptmap = get_cptmap(node);

	for(int t=0;t<n;++t) {
		double d;
		vector<bool> line;
		keyMap km;
			for(int parentIdx = 0; parentIdx < parentCount; parentIdx ++ ) {
				DSL_node *parentNode = net->GetNode(parents[parentIdx]);
			const DSL_Dmatrix &parent_cpt = *parentNode->Definition()->GetMatrix();
			d = fRand(0.0, 1.0);
			string value;
				if (d < parent_cpt[0]) {
				value = "True";
				} else {
				value = "False";
			}
			cout << value << " ";
			km[string(parentNode->GetId())] = value;
			}

			d = fRand(0.0, 1.0);
		if(d < cptmap[km]["True"]) {
			cout<< "True" <<endl;
		} else {
			cout<< "False" <<endl;
		}
	}
}


void cancerNet(ifstream ifs){
	DSL_network theNet;
	DSL_stringArray tfnames;
	tfnames.Add("True");
	tfnames.Add("False");

	int smoker = theNet.AddNode(DSL_CPT, "Smoker");
	int genetic = theNet.AddNode(DSL_CPT, "Genetic");
	int coalworker = theNet.AddNode(DSL_CPT, "CoalWorker");
	int baddiet = theNet.AddNode(DSL_CPT, "BadDiet");

	int cancer = theNet.AddNode(DSL_CPT, "LungCancer");

	int parents[] = {smoker, genetic, coalworker, baddiet};

	DSL_doubleArray theProbs;


	theProbs.Flush();
	theProbs.SetSize(2);

	for(int i=0; i<4 ; ++i) {
		
		cin >> theProbs[0];
		theProbs[1] = 1.0 - theProbs[0];
		//cin >> theProbs[1];
		theNet.GetNode(parents[i])->Definition()->SetNumberOfOutcomes(tfnames);
		theNet.GetNode(parents[i])->Definition()->SetDefinition(theProbs);

		theNet.AddArc(parents[i], cancer);
	}

	theProbs.Flush();
	theProbs.SetSize(32);
	for(int i=0;i<32;i=i+2) {
		cin >> theProbs[i];
		theProbs[i+1] = 1.0 - theProbs[i];
		//theProbs[i] = 0.4;
		//theProbs[i+1] = 0.6;
	}

	theNet.GetNode(cancer)->Definition()->SetNumberOfOutcomes(tfnames);
	theNet.GetNode(cancer)->Definition()->SetDefinition(theProbs);
	theNet.GetNode(cancer)->ChangeType(DSL_NOISY_MAX);

	printCPT(theNet.GetNode(cancer));
	theNet.WriteFile("tutorial.dsl");
}

int main(int argc, char* argv[])
{
	//cancerNet();
	//generateDatafile(100, "CancerOR.xdsl", "CancerOR_generatedbyCpp.txt");
	if(argc > 3) {
		srand(atoi(argv[3]));
	} else {
		srand(time(NULL));
	}
	generateDatafile(atoi(argv[1]), argv[2]);
	return 0;
}

