#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include "GaussJordan.h"

using namespace std;

void failTest(string testNumber, string error) {
	cout << "Failed test '" << testNumber << "' :" << error << endl;
}

void passTest(string testNumber) {
	cout << "Passed test '" << testNumber << "'" <<  endl;
}

void readSystemOfEquations(vf2D &coefficients, vf1D &absoluteTerms) {
		for(unsigned int r = 0 ; r < coefficients.size() ; ++r) {
			for(unsigned int c = 0 ; c< coefficients[0].size() ; ++c) {
				cin >> coefficients[r][c]; 
			}
			cin >> absoluteTerms[r]; 
		}
}

int main(int argc, char* argv[]) {
	string test_name;
	while(true) {
		cin >> test_name;
		if (test_name == string("EXIT"))
			break;
		int cols, rows;
		cin >> rows >> cols;
		vf2D coefficients( rows, vf1D(cols, 0.0));
		vf1D absoluteTerms( rows, 0.0);

		readSystemOfEquations(coefficients, absoluteTerms);

		SystemOfEquations soe(rows, cols);
		soe.setCoefficients(coefficients);
		soe.setAbsoluteTerms(absoluteTerms);

		//cout << soe.print();
		//cout << "After" << endl;

		unsigned int result = soe.GaussJordanElimination();

		unsigned int out_result;
		cin >> out_result;

		if (result != out_result) {
			ostringstream error;
			error << "Return codes differ, got: " << result << ", expected: " << out_result;
			failTest(test_name, error.str());
			continue;
		}


		readSystemOfEquations(coefficients, absoluteTerms);

		SystemOfEquations out_soe(rows, cols);
		out_soe.setCoefficients(coefficients);
		out_soe.setAbsoluteTerms(absoluteTerms);
		
		if ( soe != out_soe ){
			ostringstream error;
			error << "Systems of equations differ.";
			cout << "Result: " << soe.print() << endl;
			cout << "Expected: " << out_soe.print() << endl;
			failTest(test_name, error.str());
			continue;
		}
		passTest(test_name);
	}

	return 0;
}
