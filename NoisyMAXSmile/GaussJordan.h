#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>

#define DEBUG(x) cout << #x << ": " << (x) << endl
#define DEBUGV(x) cout << #x << ":"; for(unsigned int ___idx = 0; ___idx < (x).size() ; ++ ___idx){ cout << " " << (x)[___idx];} cout << endl

using namespace std;

typedef double t_coeff;
typedef vector<t_coeff> vf1D;
typedef vector< vf1D > vf2D;

typedef vector<int> vi1D;
typedef vector<vi1D > vi2D;

const t_coeff EPSILON = 10e-7;

const unsigned int GJ_OK = 0;
const unsigned int GJ_PIVOT_ERROR = 1;

struct Equation {
	vf1D coefficients;
	vf1D linearCombination;
	double absoluteTerm;
	vi1D termCounts;
	double fitness;
};

class SystemOfEquations {
	public:
		unsigned int numberOfRows;
		unsigned int numberOfColumns;
		vf2D coefficients;
		vf2D linearCombination;
		vf1D absoluteTerms;
		vi2D termCounts;
		vf1D termFitness;

		vector<Equation > equations;
		
		SystemOfEquations(unsigned int rows, unsigned int cols);	
		void setCoefficients(vf2D &originalCoefficients);
		void setAbsoluteTerms(vf1D &originalAbsoluteTerms);
		void setTermCounts(vi2D &termCounts);

		void computeFitness();
		void sortEquations();
		string print(bool printLinearCombination = false);
		int GaussJordanElimination(unsigned int rowLimiter = 0);
		bool operator==(const SystemOfEquations& rhs);
		bool operator!=(const SystemOfEquations& rhs);
		vector<int> getParameterIndices();
};

