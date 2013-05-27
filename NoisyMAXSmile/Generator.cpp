#include <iostream>
#include <fstream>
#include "smile/smile.h"
#include <stdlib.h>
#include <string>
#include <vector>
#include <time.h>
#include <map>

using namespace std;

//void printCPT(DSL_node *node) {
//	DSL_network* net = node->Network(); // node network
//	int handle = node->Handle();
//	DSL_nodeDefinition *def = node->Definition();
//	const DSL_Dmatrix &cpt = *def->GetMatrix();
//	const DSL_idArray &outcomes = *def->GetOutcomesNames();
//	const DSL_intArray &parents = net->GetParents(handle);
//	int parentCount = parents.NumItems();
//
//	DSL_intArray coords;
//	for (int elemIdx = 0; elemIdx < cpt.GetSize(); elemIdx ++) {
//		cpt.IndexToCoordinates(elemIdx, coords);
//		cout << "P(" << node->GetId() << " = " << outcomes[coords[parentCount]] << " | ";
//		for (int parentIdx = 0; parentIdx < parentCount; parentIdx ++) {
//			 if (parentIdx > 0) cout << ", ";
//			 DSL_node *parentNode = net->GetNode(parents[parentIdx]);
//			 const DSL_idArray &parentStates = *parentNode->Definition()->GetOutcomesNames();
//			 cout << parentNode->GetId() << " = " << parentStates[coords[parentIdx]];
//		}
//		cout << ") = " << cpt[elemIdx] << endl;
//	}
//}

typedef vector<pair<string, double> > vpsd;
typedef vector<string> vs;

double fRand(double fMin, double fMax) {
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

string getRandom(double r, vpsd &items) {
	double suma = 0.0;
	for(unsigned int i=0 ;i< items.size() ; ++i ) {
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

DSL_doubleArray getRandomColumn(int columnSize) {
	double columnSum = 0.0;
	DSL_doubleArray theProbs;
	theProbs.SetSize(columnSize);
	for (int rowIdx = 0 ; rowIdx < columnSize; ++rowIdx) {
		columnSum += (theProbs[rowIdx] = fRand(0.0,1.0));
	}
	double newColumnSum = 0.0;
	for (int rowIdx = 0 ; rowIdx < columnSize - 1; ++rowIdx) {
		newColumnSum += (theProbs[rowIdx] /= columnSum);
	}
	theProbs[columnSize - 1] = 1.0 - newColumnSum;
	return theProbs;
}

void randomizeNetwork( string networkName, string networkOutFile) {
	DSL_network theNet;
	theNet.ReadFile(networkName.c_str());
	int child_node = theNet.FindNode("C1");
	DSL_node *node = theNet.GetNode(child_node);
	DSL_nodeDefinition *def = node->Definition();
	//const DSL_idArray &outcomes = *def->GetOutcomesNames();
	int handle = node->Handle();
	const DSL_intArray &parents = theNet.GetParents(handle);
	int parentCount = parents.NumItems();

	for (int parentIdx = 0; parentIdx < parentCount; parentIdx ++ ) {
		DSL_node *parentNode = theNet.GetNode(parents[parentIdx]);
		DSL_nodeDefinition *parentDef = parentNode ->Definition();
		DSL_Dmatrix &parent_cpt = *parentDef->GetMatrix();
		int columnSize = parent_cpt.GetSize();
		DSL_doubleArray theProbs = getRandomColumn(columnSize);

		//theProbs.SetSize(parent_cpt.GetSize());
		//vector<double>
		//double columnSum = 0.0;
		//for (int rowIdx = 0 ; rowIdx < parent_cpt.GetSize(); ++rowIdx) {
		//	columnSum += (theProbs[rowIdx] = fRand(0.0,1.0));
		//}
		//double newColumnSum = 0.0;
		//for (int rowIdx = 0 ; rowIdx < parent_cpt.GetSize() - 1; ++rowIdx) {
		//	newColumnSum += (theProbs[rowIdx] /= columnSum);
		//}
		//theProbs[parent_cpt.GetSize() - 1] = 1.0 - newColumnSum;

		parentDef->SetDefinition(theProbs);
	}

		DSL_noisyMAX *childMAXDefinition = new DSL_noisyMAX(*(node->Definition()));
		DSL_Dmatrix &childMAXCpt = childMAXDefinition->GetCiWeights();
		int childDimension = def->GetNumberOfOutcomes();
		DSL_doubleArray theProbs;
		theProbs.SetSize(childMAXCpt.GetSize());
		//cout << "childsize"<<childMAXCpt.GetSize() << endl;
		for(int cptIdx = 0 ; cptIdx < childMAXCpt.GetSize(); cptIdx += childDimension) {
			if (childMAXCpt[cptIdx] == 0.0){
				for(int newColIdx = 0 ; newColIdx < childDimension; ++ newColIdx){
					theProbs[cptIdx + newColIdx] = childMAXCpt[cptIdx + newColIdx];
				}
			} else {

				DSL_doubleArray newCol = getRandomColumn(childDimension);
				for(int newColIdx = 0 ; newColIdx < newCol.GetSize(); ++ newColIdx){
					theProbs[cptIdx + newColIdx] = newCol[newColIdx];
				}
			}
			//cout << childMAXCpt[cptIdx] << " ";
		}
		//cout << endl;
		//for(int cptIdx = 0 ; cptIdx < theProbs.GetSize(); ++cptIdx) {
		//	cout << theProbs[cptIdx] << " ";
		//}
		//cout << endl;
		node->Definition()->SetDefinition(theProbs);
		node->Definition()->GetMatrix()->Normalize();


	theNet.WriteFile(networkOutFile.c_str());
}
void generateDatafile(int n, string networkName, string outFile) {
	ofstream outputStream(outFile.c_str());
	DSL_network theNet;
	theNet.ReadFile(networkName.c_str());
	int child_node = theNet.FindNode("C1");
	DSL_node *node = theNet.GetNode(child_node);

	DSL_network* net = node->Network(); // node network
	int handle = node->Handle();
	DSL_nodeDefinition *def = node->Definition();
	//const DSL_Dmatrix &cpt = *def->GetMatrix();
	const DSL_idArray &outcomes = *def->GetOutcomesNames();
	const DSL_intArray &parents = net->GetParents(handle);
	int parentCount = parents.NumItems();

	DSL_intArray coords;

	//double leak = 0.01;

	// PRINT HEADER
	for (int parentIdx = 0; parentIdx < parentCount; parentIdx ++ ) {
		DSL_node *parentNode = net->GetNode(parents[parentIdx]);
		outputStream << parentNode->GetId() << " ";			 
	}
	outputStream << node->GetId() << endl;

	cptMap cptmap = get_cptmap(node);

	for(int t=0;t<n;++t) {
		double d;
		vector<bool> line;
		keyMap km;
		for(int parentIdx = 0; parentIdx < parentCount; parentIdx ++ ) {
			DSL_node *parentNode = net->GetNode(parents[parentIdx]);
			const DSL_Dmatrix &parent_cpt = *parentNode->Definition()->GetMatrix();
			const DSL_idArray &parent_outcomes = *parentNode->Definition()->GetOutcomesNames();
			d = fRand(0.0, 1.0);
			string value;
			if (d < parent_cpt[0]) {
				value = string(parent_outcomes[0]);
			} else {
				value = string(parent_outcomes[1]);
			}
			outputStream << value << " ";
			km[string(parentNode->GetId())] = value;
		}

		d = fRand(0.0, 1.0);
		if(d < cptmap[km][outcomes[0]]) {
			outputStream<< outcomes[0] <<endl;
		} else {
			outputStream<< outcomes[1] <<endl;
		}
	}
	outputStream.close();
}

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

