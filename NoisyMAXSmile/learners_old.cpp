#include <iostream>
#include <fstream>
#include "smile/smile.h"
#include "smile/smilearn.h"
#include <string>
#include <vector>
#include <time.h>
//#include "GaussJordan.h"
#include <algorithm>
#include  "Learners.h"

#define AS_FLOAT(NOM,DENOM,row, col) ((double)(NOM)[(row)][(col)])/((double)(DENOM[(col)]))

using namespace std;

typedef vector<int> vi1D;
typedef vector<vi1D > vi2D;

class LearningInfo {
	public:
		string dataInfile;
		string networkInfile;
		string childName;

		DSL_dataset dataSet;
		DSL_network originalNet;
		vector<DSL_datasetMatch> matches;
		map<int, int> matchNetToData;
		map<int, int> matchDataToNet;

		int childIdx;
		int sumParentDimensions;
		int childDimension;
		DSL_node *childNode;
		DSL_noisyMAX *childMAXDefinition;

		DSL_doubleArray theProbs;
		vi1D parentIndices;
		int numberOfParents;
		map<int, int> distinguishedStates; // map of (nodeIdx, distinguishedStateIdx)
		vi1D parentDimensions;
		vector< DSL_intArray > parentOutcomesStrengths;

		int minimalNumberOfParameters;
		int parametersRowLength;
		vi1D parameterRowOffset;

		int probsSize;

		LearningInfo(string data_infile, string network_infile, string child_name) {

			if (dataSet.ReadFile(data_infile.c_str()) != DSL_OKAY) {
				cout << "Cannot read data file... exiting." << endl;
				exit(1);
			}

			if (originalNet.ReadFile(network_infile.c_str(), DSL_XDSL_FORMAT) != DSL_OKAY) {
				cout << "Cannot read network... exiting." << endl;
				exit(1);
			}

			string err;
			if (dataSet.MatchNetwork(originalNet, matches, err) != DSL_OKAY) {
				cout << "Cannot match network... exiting." << endl;
				exit(1);
			}

			for(unsigned int i=0 ; i < matches.size() ; ++i) {
				matchNetToData[matches[i].node] = matches[i].column;
				matchDataToNet[matches[i].column] = matches[i].node;
			}

			childIdx = originalNet.FindNode(child_name.c_str());
			childNode = originalNet.GetNode(childIdx);

			if (childNode->Definition()->GetType() != (DSL_CHANCE | DSL_DISCRETE | DSL_NOISY_MAX) ){
				cout << "Child should be a NoisyMAX... exiting" << endl;
				// ewentualnie zmienic na noisy-max ręcznie
				exit(1);
			}

			childMAXDefinition = new DSL_noisyMAX(*(childNode->Definition()));

			DSL_intArray &parents = originalNet.GetParents(childNode->Handle());
			numberOfParents = parents.NumItems();
			parentIndices = vector<int>(numberOfParents, 0);
			for(int i=0; i<numberOfParents; ++i)
				parentIndices[i] = parents[i];

			childDimension = childNode->Definition()->GetNumberOfOutcomes();
			parentDimensions = vector<int>(numberOfParents, 0);
			sumParentDimensions = 0;

			parentOutcomesStrengths = vector<DSL_intArray>(numberOfParents);
			minimalNumberOfParameters = 1; // minimal number of unique parameters to calculate (count leak right away)

			for(int parentIdx = 0 ; parentIdx < numberOfParents ; ++parentIdx) {
				DSL_node *parentNode = originalNet.GetNode(parentIndices[parentIdx]);
				sumParentDimensions += (parentDimensions[parentIdx] = parentNode->Definition()->GetNumberOfOutcomes()); //parent dimension is equal to the number of outcomes
				parentOutcomesStrengths[parentIdx] = childMAXDefinition->GetParentOutcomeStrengths(parentIdx);
				//for (int stateIdx=0 ; stateIdx < parentDimensions[parentIdx] ; ++stateIdx)
				//	cout << parentOutcomesStrengths[parentIdx][stateIdx] << " ";
				//cout << endl;
				minimalNumberOfParameters += parentDimensions[parentIdx] - 1; // (each parent dimension reduced by one) because we don't count distinguished states of parents
				distinguishedStates[parentIdx] = parentOutcomesStrengths[parentIdx][parentDimensions[parentIdx] - 1];
			}

			int sumOffset = 0;
			parameterRowOffset = vi1D(numberOfParents + 1, 0); // +1 so we know the offset for LEAK column
			for(int parentIdx = 0; parentIdx < numberOfParents ; ++parentIdx) {
				parameterRowOffset[parentIdx] = sumOffset;
				sumOffset += parentDimensions[parentIdx] - 1;
			}
			parameterRowOffset[numberOfParents] = sumOffset;

			parametersRowLength = minimalNumberOfParameters;
			minimalNumberOfParameters *= (childDimension - 1); // number of unique rows, last row is always 1.0 - sum

		//	DEBUG(minimalNumberOfParameters);
		//	DEBUG(childDimension);
		//	DEBUGV(parentDimensions);
		//	DEBUGV(parameterRowOffset);
			

			//for(int j=0; j< 7 ; ++j) {
			//	DSL_datasetVarInfo vi = ds.GetVariableInfo(j);
			//	cout << "discreete:" << vi.discrete << " id:" << vi.id << endl << " missingInt:" << vi.missingInt << " mF:" << vi.missingFloat << "snames:"<< endl;
			//	for(int i=0;i<vi.stateNames.size(); ++i)
			//		cout << vi.stateNames[i]<< " ";
			//	cout <<endl;
			//}

			//for(int i = 0; i < ds.GetNumberOfRecords(); ++i) { 

			//	vector<int> row(ds.GetNumberOfVariables(), 0);
			//	int sum_ones = 0;

			//	for(int j = 0; j < ds.GetNumberOfVariables(); ++j) {
			//		sum_ones += (row[j] = ds.GetInt(j,i));
			//	}
			//}
			//vector<int> rd = ds.GetIntData(0);
			//cout <<"RDSize:"<<rd.size()<< endl;
			//for(int i=0;i<rd.size();++i) {
			//	cout << vi.stateNames[rd[i]] << endl;
			//}
			//
		}

		void fillCPTWithOR(vf2D &noisyORParameters) { 

			probsSize = childDimension * sumParentDimensions + childDimension; //childDimension[rows] * sumParentDimensions[cols] + leak_column
			theProbs.SetSize(probsSize);

			int probsIdx = 0; //index of cpt vector
			int paramsRowIdx = 0; //index of minimal_parameter array
			int paramsColIdx = 0; //index of minimal_parameter array
			for (int parentIdx = 0; parentIdx < numberOfParents; ++parentIdx) { // for each column
				for(int parentStateIdx = 0; parentStateIdx < parentDimensions[parentIdx] ; ++parentStateIdx) { //for each subset of columns
					double sumColumn = 0.0;
					for(int childStateIdx = 0; childStateIdx < childDimension; ++childStateIdx) { //for each row
						//if (parentStateIdx == distinguishedStates[parentIdx]) { //if it's the distinguished state
						if (parentStateIdx == parentDimensions[parentIdx] - 1) { //if it's the distinguished state
							if (childStateIdx == childDimension - 1) { // if it's both distinguished for child and parent
								theProbs[probsIdx++] = 1;
							} else {
								sumColumn += theProbs[probsIdx++] = 0;
							}
						} else {
							if (childStateIdx == childDimension - 1) { // if it's the child's distinguished state
								theProbs[probsIdx++] = 1.0 - sumColumn;
							} else { // standard parameter
								sumColumn += theProbs[probsIdx++] = noisyORParameters[paramsRowIdx][paramsColIdx];
								++paramsRowIdx;
								if (paramsRowIdx == (childDimension - 1)) { ++paramsColIdx; paramsRowIdx = 0; }
							}
						}
					}
				}
			}

			double sumColumn = 0.0;
			for (int c_idx = 0; c_idx < childDimension; ++c_idx) { // LEAK parameters
				if (c_idx == childDimension - 1) {
					theProbs[probsIdx++] = 1.0 - sumColumn;
				} else {
					sumColumn += theProbs[probsIdx++] = noisyORParameters[paramsRowIdx][paramsColIdx];
					++paramsRowIdx;
					if (paramsRowIdx == (childDimension - 1)) { ++paramsColIdx; paramsRowIdx = 0; }
				}
			}

			childNode->Definition()->SetDefinition(theProbs);
			childNode->Definition()->GetMatrix()->Normalize();
		}
};


void printCPT(DSL_node *node) {
	DSL_network* net = node->Network(); // node network																						 
	int handle = node->Handle();
	DSL_nodeDefinition *def = node->Definition();
	const DSL_Dmatrix &cpt = *def->GetMatrix();
	const DSL_idArray &outcomes = *def->GetOutcomesNames();
	const DSL_intArray &parents = net->GetParents(handle);
	int parentCount = parents.NumItems();

	DSL_intArray coords;

	// for (int elemIdx = 0; elemIdx < cpt.GetSize(); elemIdx ++) { for (int parentIdx = 0; parentIdx < parentCount; parentIdx ++){ } }

	for (int elemIdx = 0; elemIdx < cpt.GetSize(); elemIdx ++) {
		string name = "";
		cpt.IndexToCoordinates(elemIdx, coords);
		//cout << "P(" << node->GetId() << " = " << outcomes[coords[parentCount]] << " | ";
		for (int parentIdx = 0; parentIdx < parentCount; parentIdx ++) {
			DSL_node *parentNode = net->GetNode(parents[parentIdx]);
			if(elemIdx == 0) {
				cout << parentNode->GetId()<< " ";
				if(parentIdx == parentCount-1){ cout<< node->GetId() <<endl;}
			}
			const DSL_idArray &parentStates = *parentNode->Definition()->GetOutcomesNames();
			//cout << parentNode->GetId() << " = " << parentStates[coords[parentIdx]];				 

			name += parentStates[coords[parentIdx]];
			name += " ";
		}
		name += outcomes[coords[parentCount]];
		cout << name << " " << cpt[elemIdx] << endl;
		//cout << ") = " << cpt[elemIdx] << endl;
	}
}

vf1D getWeights(DSL_network &network, string &childName){
	int childIdx = network.FindNode(childName.c_str());
	DSL_node* node = network.GetNode(childIdx);
	int handle = node->Handle();
	DSL_nodeDefinition *def = node->Definition();
	const DSL_Dmatrix &cpt = *def->GetMatrix();
	const DSL_intArray &parents = network.GetParents(handle);
	int parentCount = parents.NumItems();

	DSL_intArray coords;

	unsigned int colSize = def->GetNumberOfOutcomes();
	unsigned int rowSize = cpt.GetSize() / colSize;
	vf1D weights(rowSize, 0.0);

	unsigned int rowIdx = 0;
	for (int elemIdx = 0; elemIdx < cpt.GetSize(); elemIdx += colSize) {
		cpt.IndexToCoordinates(elemIdx, coords);
		double mult = 1.0;
		for (int parentIdx = 0; parentIdx < parentCount; parentIdx ++) {
			DSL_node *parentNode = network.GetNode(parents[parentIdx]);
			const DSL_Dmatrix &parent_cpt = *parentNode->Definition()->GetMatrix();
			mult *= parent_cpt[coords[parentIdx]];
		}
		weights[rowIdx++] = mult;
	}

	return weights;
}

vf2D getCPTArray(DSL_network &network, string &childName){
	//DSL_network* net = node->Network(); // node network																						 
	int childIdx = network.FindNode(childName.c_str());
	DSL_node* node = network.GetNode(childIdx);
	DSL_nodeDefinition *def = node->Definition();
	const DSL_Dmatrix &cpt = *def->GetMatrix();

	DSL_intArray coords;

	unsigned int colSize = def->GetNumberOfOutcomes();
	unsigned int rowSize = cpt.GetSize() / colSize;
	vf2D cpt2(rowSize , vf1D(colSize, 0.0));

	unsigned int colIdx = 0;
	unsigned int rowIdx = 0;
	for (int elemIdx = 0; elemIdx < cpt.GetSize(); elemIdx ++) {
		cpt2[rowIdx][colIdx] = cpt[elemIdx];
		++colIdx;
		if(colIdx == colSize){
			colIdx = 0;
			++rowIdx;
		}
	}

	return cpt2;
}

float EuclidianDistance(vf2D &A, vf2D &B, vf1D &w) {

	double weightsSum = 0.0;
	for(unsigned int i = 0 ; i < w.size(); ++i) {
		weightsSum += w[i];
	}

	if (A.size () != B.size()){
		cout << "Sizes differ! .. exiting" << endl;
		exit(1);
	}

	//for(int i=0;i<A.size(); ++i) { DEBUGV(A[i]); }
	//for(int i=0;i<B.size(); ++i) { DEBUGV(B[i]); }

	double totalSum = 0.0;

	for(unsigned int columnIdx = 0 ; columnIdx < A.size(); ++columnIdx) {
		double colSum = 0.0;
		for(unsigned int rowIdx = 0 ; rowIdx < A[columnIdx].size(); ++rowIdx) {

			//cout << A[columnIdx][rowIdx] << " vs " << B[columnIdx][rowIdx] << endl;
			//cout << A[columnIdx][rowIdx] - B[columnIdx][rowIdx] << endl;
			colSum += pow(A[columnIdx][rowIdx] - B[columnIdx][rowIdx], 2);
		}
		totalSum += colSum * w[columnIdx];
		//totalSum += colSum;
	}
	return totalSum;
	//printCPT(A.GetNode(childIdx));
}

void fillORParam(vector<int> &record, LearningInfo &info, vi2D &noisyORParametersNominator, vi1D &noisyORParametersDenominator) {
	int nonDistinguished = 0;
	int activatedParentIdx = -1;
	int childDataColumn = info.matchNetToData[info.childIdx];
	for(unsigned int i = 0 ; i < record.size() - 1; ++i) {
		if (info.distinguishedStates[i] != record[i]) {
			++nonDistinguished;
			if ( nonDistinguished > 1) 
				return;
			activatedParentIdx = i;
		}
	}

	if (nonDistinguished == 0) { // LEAK

		int paramColIdx = info.parameterRowOffset[info.numberOfParents];
		int childOutcome = record[childDataColumn];
		if (childOutcome < (info.childDimension - 1)) { // if childOutcome is different than distinguished
			noisyORParametersNominator[childOutcome][paramColIdx] += 1;
		}
		noisyORParametersDenominator[paramColIdx] += 1;
	} else if (nonDistinguished == 1) { // Noisy-OR parameter
		int parentState = record[activatedParentIdx];
		int paramColIdx = info.parameterRowOffset[info.matchDataToNet[activatedParentIdx]] + parentState;
		int childOutcome = record[childDataColumn];

		if (childOutcome < (info.childDimension - 1)) { // if childOutcome is different than distinguished
			noisyORParametersNominator[childOutcome][paramColIdx] += 1;
		}
		noisyORParametersDenominator[paramColIdx] += 1;
	}
}

typedef map< vi1D, vi1D > recordMap;

void fillRecordCounter(vi1D &row, LearningInfo &info, recordMap &rm) {

	vi1D childlessRow( row.begin(), row.end());

	int childDataColumn = info.matchNetToData[info.childIdx];
	childlessRow.erase(childlessRow.begin() + childDataColumn);

	if(rm.count(childlessRow) == 0) {
		rm[childlessRow] = vi1D (info.childDimension, 0);
	}

	rm[childlessRow][row[childDataColumn]] += 1;
}

void printRecordMap(recordMap rm) {
	recordMap::iterator iter;
	for(iter = rm.begin(); iter!= rm.end(); ++iter) {
		DEBUGV(iter->first);
		DEBUGV(iter->second);
		cout <<endl;
	}
}

bool recordCounterCompare(pair<vi1D, vi1D> const& p1, pair<vi1D, vi1D> const& p2) {
	unsigned int sum1 = 0;
	unsigned int sum2 = 0;
	for(unsigned int i=0 ; i < p1.second.size(); ++i) { sum1+=p1.second[i]; }
	for(unsigned int i=0 ; i < p2.second.size(); ++i) { sum2+=p2.second[i]; }
	if (sum1 > sum2) { return true; }
	return false;
}

DSL_network OpenNetwork(string network_infile){
	DSL_network originalNet;
	if (originalNet.ReadFile(network_infile.c_str(), DSL_XDSL_FORMAT) != DSL_OKAY) {
		cout << "Cannot read network... exiting." << endl;
		exit(1);
	}

	return originalNet;
}

DSL_network LearnParamsGaussJordan(string data_infile, string network_infile, string child_name) {

	LearningInfo info = LearningInfo( data_infile, network_infile, child_name);

	//assert that we are working on binary case
	if (info.childDimension > 2) {
		cout << "This method works only for binary children..exiting" << endl;
		exit(1);
	}
	for(int parentIdx = 0 ; parentIdx < info.numberOfParents ; ++parentIdx) {
		if (info.parentDimensions[parentIdx] > 2){
			cout << "This method works only for binary children..exiting" << endl;
			exit(1);
		}
	}

	map< vi1D, vi1D > recordCounter;

	for(int i = 0; i < info.dataSet.GetNumberOfRecords(); ++i) { 
		vi1D row(info.dataSet.GetNumberOfVariables(), 0);

		for(int j = 0; j < info.dataSet.GetNumberOfVariables(); ++j) {
			row[j] = info.dataSet.GetInt(j, i);
		}
		fillRecordCounter(row, info, recordCounter);
	}
	//printRecordMap(recordCounter);

	vector< pair < vi1D, vi1D> > integerCoefficients;

	recordMap::iterator iter;
	for(iter = recordCounter.begin(); iter!= recordCounter.end(); ++iter) {

		vi1D newRecord(info.parametersRowLength, 0); //+1 Because of leak
		for(int columnIdx = 0 ; columnIdx < info.numberOfParents; ++ columnIdx) {  // LAZY STATE FLIPPING, look at info.parentOutcomesStrengths
			newRecord[columnIdx] = (iter->first[columnIdx]+1)%2;
		}
		newRecord[info.numberOfParents] = 1; // add leak column
		integerCoefficients.push_back(make_pair(newRecord, iter->second));
	}

	sort(integerCoefficients.begin(), integerCoefficients.end(), recordCounterCompare);
	//for(unsigned int i=0 ; i < integerCoefficients.size(); ++i) {
	//	DEBUGV(integerCoefficients[i].first);
	//	DEBUGV(integerCoefficients[i].second);
	//	cout << endl;
	//}

	vf2D coefficients(integerCoefficients.size(), vf1D(info.parametersRowLength, 0.0)); // +1 Because of leak
	vf1D absoluteTerms(integerCoefficients.size(), 0.0);
	for(unsigned int rowIdx = 0 ; rowIdx < coefficients.size() ; ++rowIdx) {
		for(unsigned int columnIdx = 0 ; columnIdx < coefficients[rowIdx].size() ; ++columnIdx) {
			coefficients[rowIdx][columnIdx] = (double) integerCoefficients[rowIdx].first[columnIdx];
			vi1D &childCounts = integerCoefficients[rowIdx].second;
			absoluteTerms[rowIdx] = ((double)(childCounts[0]))/(childCounts[0] + childCounts[1]);
		}

		//DEBUGV(integerCoefficients[rowIdx].first);
		//DEBUGV(integerCoefficients[rowIdx].second);
		//DEBUGV(coefficients[rowIdx]);
		//DEBUG(absoluteTerms[rowIdx]);
		//cout << endl;
	}

	SystemOfEquations soe (coefficients.size(), info.parametersRowLength);

	soe.setCoefficients(coefficients);
	soe.setAbsoluteTerms(absoluteTerms);

	//cout << soe.print() << endl;

	for(unsigned int rowIdx = 0 ; rowIdx < coefficients.size() ; ++rowIdx) {
		soe.absoluteTerms[rowIdx] = log(1.0 - soe.absoluteTerms[rowIdx]);
	}

	//cout << soe.print() << endl;

	soe.GaussJordanElimination();

	for(unsigned int rowIdx = 0 ; rowIdx < coefficients.size() ; ++rowIdx) {
		soe.absoluteTerms[rowIdx] = 1.0 - exp(soe.absoluteTerms[rowIdx]);
	}
	//cout << soe.print() << endl;

	vector<int> indices = soe.getParameterIndices();

	//vf2D noisyORParameters ( info.parametersRowLength, vf1D(info.childDimension - 1 , 0.0));
	vf2D noisyORParameters (info.childDimension - 1, vf1D(info.parametersRowLength, 0.0));
	for(int paramIdx = 0 ; paramIdx < info.parametersRowLength ; ++paramIdx) {

		t_coeff coefficient = soe.absoluteTerms[indices[paramIdx]];
		if (coefficient <= 0.0) {
			coefficient = 0.0001;
		}
		else if( coefficient >= 1.0) {
			coefficient = 0.9999;
		}

		noisyORParameters[0][paramIdx] = coefficient; // index of 0 because of binary child
	}

	info.fillCPTWithOR(noisyORParameters);

	info.originalNet.WriteFile("debugGJ.dsl");

	return info.originalNet;
}

DSL_network LearnParamsNaive(string data_infile, string network_infile, string child_name) {

	LearningInfo info = LearningInfo( data_infile, network_infile, child_name);

	vi2D noisyORParametersNominator(info.childDimension - 1, vi1D(info.parametersRowLength, 0)); // pairs of (nominator, denominator)
	vi1D noisyORParametersDenominator(info.parametersRowLength, 0);

	for(int i = 0; i < info.dataSet.GetNumberOfRecords(); ++i) { 
		vi1D row(info.dataSet.GetNumberOfVariables(), 0);
		for(int j = 0; j < info.dataSet.GetNumberOfVariables(); ++j)
			row[j] = info.dataSet.GetInt(j, i);
		fillORParam(row, info, noisyORParametersNominator, noisyORParametersDenominator);
	}

	//for(unsigned int i = 0 ; i < noisyORParametersNominator.size() ; ++i) {
	//	DEBUGV(noisyORParametersNominator[i]);
	//}
	//DEBUGV(noisyORParametersDenominator);

	// handle 0.0 and 1.0
	for(int i = 0 ; i < info.parametersRowLength ; ++ i) {
		if ( noisyORParametersDenominator[i] == 0 ) { // if no record for given noisyOR parameter was found => apply uniform distribution
			for(int j = 0 ; j < info.childDimension - 1 ; ++j) {
				noisyORParametersNominator[j][i] = 1;
			}
			noisyORParametersDenominator[i] = info.childDimension;
		} else {
			int sumNominators = 0;
			for(int j = 0 ; j < info.childDimension - 1 ; ++j) {
				sumNominators += noisyORParametersNominator[j][j];
			}

			if (sumNominators == noisyORParametersDenominator[i]) {
				noisyORParametersDenominator[i] += 1;
			}

			for(int j = 0 ; j < info.childDimension - 1 ; ++j) {
				if (noisyORParametersNominator[j][i] == 0) {
					noisyORParametersNominator[j][i] = 1;
					noisyORParametersDenominator[i] += 1;
				}
			}
		}
	}

	//for(unsigned int i = 0 ; i < noisyORParametersNominator.size() ; ++i) {
	//	DEBUGV(noisyORParametersNominator[i]);
	//}
	//DEBUGV(noisyORParametersDenominator);

	vf2D noisyORParameters ( noisyORParametersNominator.size(), vf1D(noisyORParametersNominator[0].size(), 0.0));

	for(unsigned int i = 0 ; i < noisyORParametersNominator.size() ; ++i) {
		for(unsigned int j = 0 ; j < noisyORParametersNominator[i].size() ; ++j) {
			noisyORParameters[i][j] = AS_FLOAT(noisyORParametersNominator, noisyORParametersDenominator, i, j);
		}
	}

	info.fillCPTWithOR(noisyORParameters);

	info.originalNet.WriteFile("debugNAIVE.dsl");

	return info.originalNet;
}

DSL_network LearnParamsEM(string data_infile, string network_infile, string child_name) {
	DSL_dataset ds;
	if (ds.ReadFile(data_infile.c_str()) != DSL_OKAY) {
		cout << "Cannot read data file... exiting." << endl;
		exit(1);
	}

	DSL_network originalNet;
	if (originalNet.ReadFile(network_infile.c_str(), DSL_XDSL_FORMAT) != DSL_OKAY) {
		cout << "Cannot read network... exiting." << endl;
		exit(1);
	}

	int childIdx = originalNet.FindNode(child_name.c_str());
	originalNet.GetNode(childIdx)->ChangeType(DSL_CPT);

	vector<DSL_datasetMatch> matches;
	string err;
	if (ds.MatchNetwork(originalNet, matches, err) != DSL_OKAY) {
		cout << "Cannot match network... exiting." << endl;
		exit(1);
	}

	DSL_em em;
	em.SetUniformizeParameters(true);
	em.SetRandomizeParameters(true);
	em.SetSeed(0);
	em.SetEquivalentSampleSize(1);

	if (em.Learn(ds, originalNet, matches) != DSL_OKAY) {
		cout << "Cannot learn parameters... exiting." << endl;
		exit(1);
	}

	return originalNet;
}

int main(int argc, char* argv[]) {
	ios_base::sync_with_stdio(0);

	string data_infile = string(argv[1]);
	string network_infile = string(argv[2]);
	string child_name = string("C1");

	if (argv[3] == string("EM") || argv[3] == string("Smile")) {
		DSL_network net = LearnParamsEM(data_infile, network_infile, child_name);
		int childIdx = net.FindNode(child_name.c_str());
		DSL_node* childNode = net.GetNode(childIdx);

		if (argv[3] == string("Smile")) {
			childNode->ChangeType(DSL_NOISY_MAX);
		}

		//printCPT(childNode);

	} else if (argv[3] == string("Naive")) {
		DSL_network net = LearnParamsNaive(data_infile, network_infile, child_name);
		//int childIdx = net.FindNode(child_name.c_str());
		//DSL_node* childNode = net.GetNode(childIdx);
		//printCPT(childNode);

	} else if (argv[3] == string("Gauss")) {
		DSL_network net = LearnParamsGaussJordan(data_infile, network_infile, child_name);
		//int childIdx = net.FindNode(child_name.c_str());
		//DSL_node* childNode = net.GetNode(childIdx);
		//printCPT(childNode);
	} else if (argv[3] == string("NvGJ")) { 
		
		DSL_network netGJ = LearnParamsGaussJordan(data_infile, network_infile, child_name);
		DSL_network netNaive = LearnParamsNaive(data_infile, network_infile, child_name);
		DSL_network netEM = LearnParamsEM(data_infile, network_infile, child_name);

		DSL_network netOriginal = OpenNetwork(network_infile);

		vf1D weights = getWeights(netOriginal, child_name);

		vf2D cptGJ = getCPTArray(netGJ, child_name);
		vf2D cptNaive = getCPTArray(netNaive, child_name);
		vf2D cptEM = getCPTArray(netEM, child_name);

		vf2D cptOriginal = getCPTArray(netOriginal, child_name);

		cout << "Gauss:" << EuclidianDistance( cptOriginal, cptGJ, weights) << endl;
		cout << "Naive:" << EuclidianDistance( cptOriginal, cptNaive, weights) << endl;
		cout << "EM:" << EuclidianDistance( cptOriginal, cptEM, weights) << endl;

		//cout << "Euclidian distance (orig vs GJ):    " << EuclidianDistance(original, netGJ, child_name) << endl;
		//cout << "Euclidian distance (orig vs Naive): " << EuclidianDistance(original, netNaive, child_name) << endl;
		//int childIdx = net.FindNode(child_name.c_str());

	}

	return 0;
}

