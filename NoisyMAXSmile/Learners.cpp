#include "Learners.h"

LearningInfo::LearningInfo(string data_infile, string network_infile, string child_name) {

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

LearningInfo::~LearningInfo(){
	delete childMAXDefinition;
}

void LearningInfo::fillCPTWithOR(vf2D &noisyORParameters) { 

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

bool recordCounterCompare::operator()(const pair<vi1D, vi1D>& p1, const pair<vi1D, vi1D>& p2) {
	unsigned int sum1 = 0;
	unsigned int sum2 = 0;
	for(unsigned int i=0 ; i < p1.second.size(); ++i) { sum1+=p1.second[i]; }
	for(unsigned int i=0 ; i < p2.second.size(); ++i) { sum2+=p2.second[i]; }
	return sum1 > sum2;
}

bool recordCounterCompare2::operator()(const pair<vi1D, vi1D>& p1, const pair<vi1D, vi1D>& p2) {
	unsigned int sum1 = 0;
	unsigned int coefficient1 = 0;
	unsigned int sum2 = 0;
	unsigned int coefficient2 = 0;
	for(unsigned int i=0 ; i < p1.first.size(); ++i) { coefficient1+=abs(p1.first[i]); }
	for(unsigned int i=0 ; i < p2.first.size(); ++i) { coefficient2+=abs(p2.first[i]); }

	//cout << "Cef1:" << coefficient1 <<endl;
	//cout << "Cef2:" << coefficient2 <<endl;

	for(unsigned int i=0 ; i < p1.second.size(); ++i) { sum1+=p1.second[i]; }
	for(unsigned int i=0 ; i < p2.second.size(); ++i) { sum2+=p2.second[i]; }
	sum1 /= (double)(coefficient1);
	sum2 /= (double)(coefficient2);
	return sum1 > sum2;
}

//struct recordCounterCompare2 : public std::binary_function<pair<vi1D, vi1D> ,pair<vi1D, vi1D>,bool>
//{
//	inline bool operator()(const pair<vi1D, vi1D>& p1, const pair<vi1D, vi1D>& p2) {
//		unsigned int sum1 = 0;
//		unsigned int sum2 = 0;
//		for(unsigned int i=0 ; i < p1.second.size(); ++i) { sum1+=p1.second[i]; }
//		for(unsigned int i=0 ; i < p2.second.size(); ++i) { sum2+=p2.second[i]; }
//		return sum1 > sum2;
//	}
//};



//bool recordCounterCompare(pair<vi1D, vi1D> const& p1, pair<vi1D, vi1D> const& p2) {
//	unsigned int sum1 = 0;
//	unsigned int sum2 = 0;
//	for(unsigned int i=0 ; i < p1.second.size(); ++i) { sum1+=p1.second[i]; }
//	for(unsigned int i=0 ; i < p2.second.size(); ++i) { sum2+=p2.second[i]; }
//	if (sum1 > sum2) { return true; }
//	return false;
//}

//bool recordCounterCompare2(pair<vi1D, vi1D> const& p1, pair<vi1D, vi1D> const& p2) {
//	unsigned int sum1 = 0;
//	unsigned int coefficient1 = 0;
//	unsigned int sum2 = 0;
//	unsigned int coefficient2 = 0;
//	for(unsigned int i=0 ; i < p1.first.size(); ++i) { coefficient1+=abs(p1.first[i]); }
//	for(unsigned int i=0 ; i < p2.first.size(); ++i) { coefficient2+=abs(p2.first[i]); }
//
//	//cout << "Cef1:" << coefficient1 <<endl;
//	//cout << "Cef2:" << coefficient2 <<endl;
//
//	for(unsigned int i=0 ; i < p1.second.size(); ++i) { sum1+=p1.second[i]; }
//	for(unsigned int i=0 ; i < p2.second.size(); ++i) { sum2+=p2.second[i]; }
//	sum1 /= (double)(coefficient1);
//	sum2 /= (double)(coefficient2);
//	if (sum1 > sum2) { return true; }
//	return false;
//}

DSL_network OpenNetwork(string network_infile){
	DSL_network originalNet;
	if (originalNet.ReadFile(network_infile.c_str(), DSL_XDSL_FORMAT) != DSL_OKAY) {
		cout << "Cannot read network... exiting." << endl;
		exit(1);
	}

	return originalNet;
}

DSL_network LearnParamsGaussJordan(string data_infile, string network_infile, string child_name, int comparator, bool DEBUG_ON) {

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
		int sumOnes = 0;
		for(int columnIdx = 0 ; columnIdx < info.numberOfParents; ++ columnIdx) {  // LAZY STATE FLIPPING, look at info.parentOutcomesStrengths
			newRecord[columnIdx] = (iter->first[columnIdx]+1)%2;
			sumOnes += newRecord[columnIdx];
		}
		//newRecord[info.numberOfParents] = 1; // add leak column
		newRecord[info.numberOfParents] = -(sumOnes-1); // add leak column
		integerCoefficients.push_back(make_pair(newRecord, iter->second));
	}

	if (DEBUG_ON){
		cout << "Before sorting:" << endl;
		for(unsigned int i=0 ; i < integerCoefficients.size(); ++i) {
			for(unsigned int j=0 ; j < integerCoefficients[i].first.size(); ++j) {
				cout << integerCoefficients[i].first[j] << " ";
			}

			int denom = 0;
			for(unsigned int j=0 ; j < integerCoefficients[i].second.size(); ++j) {
				denom += integerCoefficients[i].second[j];
			}
			cout << "|" << integerCoefficients[i].second[0] << " / " << denom << " = " << (double)(integerCoefficients[i].second[0])/denom << endl;

			//DEBUGV(integerCoefficients[i].first);
			//DEBUGV(integerCoefficients[i].second);
			//cout << endl;
		}
	}

	if (comparator == 0){
		sort(integerCoefficients.begin(), integerCoefficients.end(), recordCounterCompare());
	} else if (comparator == 1){
		sort(integerCoefficients.begin(), integerCoefficients.end(), recordCounterCompare2());
	}
	
	if (DEBUG_ON){
		cout << "After sorting:" << endl;
		for(unsigned int i=0 ; i < integerCoefficients.size(); ++i) {
			for(unsigned int j=0 ; j < integerCoefficients[i].first.size(); ++j) {
				cout << integerCoefficients[i].first[j] << " ";
			}

			int denom = 0;
			for(unsigned int j=0 ; j < integerCoefficients[i].second.size(); ++j) {
				denom += integerCoefficients[i].second[j];
			}
			cout << "|" << integerCoefficients[i].second[0] << " / " << denom << " = " << (double)(integerCoefficients[i].second[0])/denom << endl;

			//DEBUGV(integerCoefficients[i].first);
			//DEBUGV(integerCoefficients[i].second);
			//cout << endl;
		}
	}

	//if (DEBUG_ON){
	//	for(unsigned int i=0 ; i < integerCoefficients.size(); ++i) {
	//		DEBUGV(integerCoefficients[i].first);
	//		DEBUGV(integerCoefficients[i].second);
	//		cout << endl;
	//	}
	//}

	vf2D coefficients(integerCoefficients.size(), vf1D(info.parametersRowLength, 0.0)); // +1 Because of leak
	vf1D absoluteTerms(integerCoefficients.size(), 0.0);
	vi2D termCounts(integerCoefficients.size(), vi1D(info.childDimension, 0));
	for(unsigned int rowIdx = 0 ; rowIdx < coefficients.size() ; ++rowIdx) {
		for(unsigned int columnIdx = 0 ; columnIdx < coefficients[rowIdx].size() ; ++columnIdx) {
			coefficients[rowIdx][columnIdx] = (double) integerCoefficients[rowIdx].first[columnIdx];
		}
		vi1D &childCounts = integerCoefficients[rowIdx].second;
		absoluteTerms[rowIdx] = ((double)(childCounts[0]))/(childCounts[0] + childCounts[1]);
		termCounts[rowIdx] = childCounts;
		//for(unsigned int termIdx = 0 ; termIdx < coefficients[rowIdx].second.size() ; ++termIdx) {
		//	termCounts[rowIdx][termIdx] = coefficients[rowIdx].second[termIdx];
		//}

		//DEBUGV(integerCoefficients[rowIdx].first);
		//DEBUGV(integerCoefficients[rowIdx].second);
		//DEBUGV(coefficients[rowIdx]);
		//DEBUG(absoluteTerms[rowIdx]);
		//cout << endl;
	}

	SystemOfEquations soe (coefficients.size(), info.parametersRowLength);

	soe.setCoefficients(coefficients);
	soe.setAbsoluteTerms(absoluteTerms);
	soe.setTermCounts(termCounts);

	//if (DEBUG_ON) {
	//	cout << soe.print() << endl;
	//}

	for(unsigned int rowIdx = 0 ; rowIdx < coefficients.size() ; ++rowIdx) {
		soe.absoluteTerms[rowIdx] = log(1.0 - soe.absoluteTerms[rowIdx]);
	}

	//if (DEBUG_ON) {
	//	cout << soe.print() << endl;
	//}

	int returnCode = soe.GaussJordanElimination();
	//if (returnCode != GJ_OK){
	//	cout << "GaussJordanElimination returned code" << returnCode << " exiting .."; 
	//	exit(1);
	//}

	for(unsigned int rowIdx = 0 ; rowIdx < coefficients.size() ; ++rowIdx) {
		soe.absoluteTerms[rowIdx] = 1.0 - exp(soe.absoluteTerms[rowIdx]);
	}

	if (DEBUG_ON) {
		cout << soe.print(true) << endl;
	}

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

DSL_network LearnParamsNaive(string data_infile, string network_infile, string child_name, bool DEBUG_ON) {

	LearningInfo info = LearningInfo( data_infile, network_infile, child_name);

	vi2D noisyORParametersNominator(info.childDimension - 1, vi1D(info.parametersRowLength, 0)); // pairs of (nominator, denominator)
	vi1D noisyORParametersDenominator(info.parametersRowLength, 0);

	for(int i = 0; i < info.dataSet.GetNumberOfRecords(); ++i) { 
		vi1D row(info.dataSet.GetNumberOfVariables(), 0);
		for(int j = 0; j < info.dataSet.GetNumberOfVariables(); ++j)
			row[j] = info.dataSet.GetInt(j, i);
		fillORParam(row, info, noisyORParametersNominator, noisyORParametersDenominator);
	}

	if (DEBUG_ON){
		for(unsigned int i = 0 ; i < noisyORParametersNominator.size() ; ++i) {
			DEBUGV(noisyORParametersNominator[i]);
		}
		DEBUGV(noisyORParametersDenominator);
	}

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

	if (DEBUG_ON) {
		for(unsigned int i = 0 ; i < noisyORParametersNominator.size() ; ++i) {
			DEBUGV(noisyORParametersNominator[i]);
		}
		DEBUGV(noisyORParametersDenominator);
	}

	vf2D noisyORParameters ( noisyORParametersNominator.size(), vf1D(noisyORParametersNominator[0].size(), 0.0));

	for(unsigned int i = 0 ; i < noisyORParametersNominator.size() ; ++i) {
		for(unsigned int j = 0 ; j < noisyORParametersNominator[i].size() ; ++j) {
			noisyORParameters[i][j] = AS_FLOAT(noisyORParametersNominator, noisyORParametersDenominator, i, j);
		}
	}

	if (DEBUG_ON) {
		DEBUGV(noisyORParameters[0]);
	}

	info.fillCPTWithOR(noisyORParameters);

	info.originalNet.WriteFile("debugNAIVE.dsl");

	return info.originalNet;
}

DSL_network LearnParamsEM(string data_infile, string network_infile, string child_name, bool DEBUG_ON) {
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

DSL_network LearnParamsSmile(string data_infile, string network_infile, string child_name, bool DEBUG_ON) {
	DSL_network em_net = LearnParamsEM(data_infile, network_infile, child_name);

	int childIdx = em_net.FindNode(child_name.c_str());
	DSL_node* childNode = em_net.GetNode(childIdx);
	childNode->ChangeType(DSL_NOISY_MAX);
	return em_net;
}


