#include <iostream>
#include <fstream>
#include "smile/smile.h"
#include "smile/smilearn.h"
#include <string>
#include <vector>
#include <time.h>
#include "GaussJordan.h"
#include <algorithm>

#define AS_FLOAT(NOM,DENOM,row, col) ((double)(NOM)[(row)][(col)])/((double)(DENOM[(col)]))

using namespace std;

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

		LearningInfo(string data_infile, string network_infile, string child_name);
		~LearningInfo();

		void fillCPTWithOR(vf2D &noisyORParameters);
};

void printCPT(DSL_node *node);
//vf1D getWeights(DSL_network &network, string &childName); 
//vf2D getCPTArray(DSL_network &network, string &childName);

//float EuclidianDistance(vf2D &A, vf2D &B, vf1D &w);

void fillORParam(vector<int> &record, LearningInfo &info, vi2D &noisyORParametersNominator, vi1D &noisyORParametersDenominator);

typedef map< vi1D, vi1D > recordMap;

void fillRecordCounter(vi1D &row, LearningInfo &info, recordMap &rm);

void printRecordMap(recordMap rm) ;

struct recordCounterCompare : public std::binary_function<pair<vi1D, vi1D> ,pair<vi1D, vi1D>,bool>
{
	inline bool operator()(const pair<vi1D, vi1D>& p1, const pair<vi1D, vi1D>& p2);
};

struct recordCounterCompare2 : public std::binary_function<pair<vi1D, vi1D> ,pair<vi1D, vi1D>,bool>
{
	inline bool operator()(const pair<vi1D, vi1D>& p1, const pair<vi1D, vi1D>& p2);
};


//struct recordCounterCompare : public std::binary_function<pair<vi1D, vi1D> ,pair<vi1D, vi1D>,bool>;
//bool recordCounterCompare(pair<vi1D, vi1D> const& p1, pair<vi1D, vi1D> const& p2);
//bool recordCounterCompare2(pair<vi1D, vi1D> const& p1, pair<vi1D, vi1D> const& p2);

DSL_network OpenNetwork(string network_infile);

DSL_network LearnParamsGaussJordan(string data_infile, string network_infile, string child_name, int comparator = 0 , bool DEBUG_ON = false);
DSL_network LearnParamsNaive(string data_infile, string network_infile, string child_name, bool DEBUG_ON = false);
DSL_network LearnParamsEM(string data_infile, string network_infile, string child_name, bool DEBUG_ON = false);
DSL_network LearnParamsSmile(string data_infile, string network_infile, string child_name, bool DEBUG_ON = false);
