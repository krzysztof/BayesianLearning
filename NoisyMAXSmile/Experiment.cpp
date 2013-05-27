#include <iostream>
#include <fstream>
#include "smile/smile.h"
#include "smile/smilearn.h"
#include <string>
#include <vector>
#include <time.h>
//#include "GaussJordan.h"
#include <algorithm>
#include "Learners.h"
#include "Generator.h"

#define AS_FLOAT(NOM,DENOM,row, col) ((double)(NOM)[(row)][(col)])/((double)(DENOM[(col)]))

using namespace std;

typedef vector<int> vi1D;
typedef vector<vi1D > vi2D;

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

double OldEuclidianDistance(vf2D &A, vf2D &B) {
	double totalSum = 0.0;
	for(unsigned int columnIdx = 0 ; columnIdx < A.size(); ++columnIdx) {
		totalSum += pow(A[columnIdx][0] - B[columnIdx][0],2);
	}
	return sqrt(totalSum);
}

double OldHellingerDistance(vf2D &A, vf2D &B) {
	double totalSum = 0.0;
	for(unsigned int columnIdx = 0 ; columnIdx < A.size(); ++columnIdx) {
		totalSum += pow(sqrt(A[columnIdx][0]) - sqrt(B[columnIdx][0]), 2);
		//totalSum += colSum;
	}
	return sqrt(totalSum) * 0.7071067811865475;// * w[columnIdx];
	//return totalSum;
}

double HellingerDistance(vf2D &A, vf2D &B, vf1D &w) {
	double totalSum = 0.0;
	for(unsigned int columnIdx = 0 ; columnIdx < A.size(); ++columnIdx) {
		double partialSum = 0.0;
		for(unsigned int rowIdx = 0 ; rowIdx < A[columnIdx].size(); ++rowIdx) {

			partialSum += pow(sqrt(A[columnIdx][rowIdx]) - sqrt(B[columnIdx][rowIdx]), 2);
		}
		totalSum += sqrt(partialSum) * 0.7071067811865475 * w[columnIdx];
	}
	return totalSum;
}

double EuclidianDistance(vf2D &A, vf2D &B, vf1D &w) {

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



int main(int argc, char* argv[]) {
	ios_base::sync_with_stdio(0);
	srand(time(NULL));

	string child_name = string("C1");

	if (argv[1] == string("Help")){
		cout << "Usage: " << endl;
		cout << argv[0] << " EM network_file data_file" << endl;
		cout << argv[0] << " Smile network_file data_file" << endl;
		cout << argv[0] << " Naive network_file data_file" << endl;
		cout << argv[0] << " Gauss network_file data_file" << endl;
		cout << argv[0] << " Test network_file output_file" << endl;
		cout << argv[0] << " Generate network_file N output_file" << endl;
		cout << argv[0] << " GenerateRandomNet network_file output_file" << endl;
		exit(1);
	}

	
	if (argv[1] == string("EM") || argv[1] == string("Smile")) {

		string network_infile = string(argv[2]);
		string data_infile = string(argv[3]);
		DSL_network net = LearnParamsEM(data_infile, network_infile, child_name);
		int childIdx = net.FindNode(child_name.c_str());
		DSL_node* childNode = net.GetNode(childIdx);

		if (argv[1] == string("Smile")) {
			childNode->ChangeType(DSL_NOISY_MAX);
		}

	} else if (argv[1] == string("Naive")) {

		string network_infile = string(argv[2]);
		string data_infile = string(argv[3]);
		DSL_network net = LearnParamsNaive(data_infile, network_infile, child_name, true);
		DSL_network netOriginal = OpenNetwork(network_infile);

		vf2D cptNaive = getCPTArray(net, child_name);
		vf2D cptOriginal = getCPTArray(netOriginal, child_name);
		vf1D weights = getWeights(netOriginal, child_name);

		cout << "OldEucl/Naive:" << OldEuclidianDistance( cptOriginal, cptNaive) << endl;
		cout << "Eucl/Naive:" << EuclidianDistance( cptOriginal, cptNaive, weights) << endl;
		//int childIdx = net.FindNode(child_name.c_str());
		//DSL_node* childNode = net.GetNode(childIdx);
		//printCPT(childNode);

	} else if (argv[1] == string("Gauss")) {

		string network_infile = string(argv[2]);
		string data_infile = string(argv[3]);
		DSL_network net = LearnParamsGaussJordan(data_infile, network_infile, child_name, 0, true);
		DSL_network netOriginal = OpenNetwork(network_infile);

		vf2D cptGJ = getCPTArray(net, child_name);
		vf2D cptOriginal = getCPTArray(netOriginal, child_name);
		vf1D weights = getWeights(netOriginal, child_name);
		cout << "OldEucl/Gauss:" << OldEuclidianDistance( cptOriginal, cptGJ) << endl;
		cout << "Eucl/Gauss:" << EuclidianDistance( cptOriginal, cptGJ, weights) << endl;

	} else if (argv[1] == string("Gauss2")) {
		string network_infile = string(argv[2]);
		string data_infile = string(argv[3]);
		DSL_network net = LearnParamsGaussJordan(data_infile, network_infile, child_name, 1, true);
		DSL_network netOriginal = OpenNetwork(network_infile);

		vf2D cptGJ = getCPTArray(net, child_name);
		vf2D cptOriginal = getCPTArray(netOriginal, child_name);
		vf1D weights = getWeights(netOriginal, child_name);
		cout << "OldEucl/Gauss:" << OldEuclidianDistance( cptOriginal, cptGJ) << endl;
		cout << "Eucl/Gauss:" << EuclidianDistance( cptOriginal, cptGJ, weights) << endl;
	} else if (argv[1] == string("GenerateRandomNet")) {
		string network_infile = string(argv[2]);
		string generatorOutfile = string(argv[3]);
		randomizeNetwork(network_infile, generatorOutfile);
		//generateDatafile(n, network_infile, generatorOutfile);
	} else if (argv[1] == string("Generate")) {
		string network_infile = string(argv[2]);
		int n = atoi(argv[3]);
		string generatorOutfile = string(argv[4]);
		randomizeNetwork(network_infile, string("randomized_net.xdsl"));
		generateDatafile(n, network_infile, generatorOutfile);
	} else if (argv[1] == string("MassGenerate")) {
		//int n = atoi(argv[4]);
		//int reps = atoi(argv[5]);
		//for(int r=0; r < reps ; ++r){
		//	//string generatorOutfile = string("./data/OR_") + n +
		//	//generateDatafile(n, network_infile, generatorOutfile);
		//}
		
	} else if (argv[1] == string("Test")) { 
		string network_infile = string(argv[2]);

		ofstream outputStream(argv[3]);
		//int arr[] = {100,1000,10000,100000};
		int arr[] = {500,1000,1500,2000};
		
		int arrSize = 4;

		int reps = 100;

		for(int ni=0; ni < arrSize; ++ni) {
			int n = arr[ni];

			outputStream <<"CASE " << n << endl;
			cout << network_infile << " doing n=" << n << endl;
			for(int r=0; r < reps; r++) {
				string tmp_datafile = string("tmp_DATAFILE.txt");
				string tmp_network = string("tmp_NETWORK.xdsl");

				randomizeNetwork(network_infile, tmp_network);
				//tmp_network = network_infile;

				generateDatafile(n, tmp_network, tmp_datafile);

				string data_infile = tmp_datafile;

				DSL_network netOriginal = OpenNetwork(tmp_network);
				vf2D cptOriginal = getCPTArray(netOriginal, child_name);
				vf1D weights = getWeights(netOriginal, child_name);


				DSL_network netGJ = LearnParamsGaussJordan(data_infile, tmp_network, child_name);
				vf2D cptGauss = getCPTArray(netGJ, child_name);

				DSL_network netGJ2 = LearnParamsGaussJordan(data_infile, tmp_network, child_name, 1, false);
				vf2D cptGauss2 = getCPTArray(netGJ2, child_name);

				DSL_network netNaive = LearnParamsNaive(data_infile, tmp_network, child_name);
				vf2D cptNaive = getCPTArray(netNaive, child_name);

				//outputStream << EuclidianDistance( cptOriginal, cptGauss, weights) << " ";
				//outputStream << EuclidianDistance( cptOriginal, cptGauss2, weights) << " ";
				//outputStream << EuclidianDistance( cptOriginal, cptNaive, weights) << endl;
				
				outputStream << OldEuclidianDistance(cptOriginal, cptGauss) << " ";
				outputStream << OldEuclidianDistance(cptOriginal, cptGauss2) << " ";
				outputStream << OldEuclidianDistance(cptOriginal, cptNaive) << endl;

				//outputStream << OldHellingerDistance(cptOriginal, cptGauss) << " ";
				//outputStream << OldHellingerDistance(cptOriginal, cptGauss2) << " ";
				//outputStream << OldHellingerDistance(cptOriginal, cptNaive) << endl;



				//DSL_network netSmile = LearnParamsSmile(data_infile, tmp_network, child_name);
				//vf2D cptSmile = getCPTArray(netSmile, child_name);

				////outputStream << EuclidianDistance( cptOriginal, cptSmile, weights) << endl;
				//outputStream << OldEuclidianDistance(cptOriginal, cptSmile) << endl;

				//outputStream << "OldEucl/Gauss:" << OldEuclidianDistance( cptOriginal, cptGauss) << endl;
				//outputStream << "OldEucl/Naive:" << OldEuclidianDistance( cptOriginal, cptNaive) << endl;
				//outputStream << "OldEucl/Smile:" << OldEuclidianDistance( cptOriginal, cptSmile) << endl;
			}
			//outputStream << "N:" << n << endl;
			//outputStream << "Gauss:" << EuclidianDistance( cptOriginal, cptGJ, weights) << endl;
			//outputStream << "Naive:" << EuclidianDistance( cptOriginal, cptNaive, weights) << endl;
			//outputStream << "EM:" << EuclidianDistance( cptOriginal, cptEM, weights) << endl;
			//WYPISAC OBIE MAPY
		}
		outputStream <<"CASE 0" << endl;
		outputStream.close();
	}

	return 0;
}

