#ifndef GAUSSJORDAN_H
#define GAUSSJORDAN_H
#include "GaussJordan.h"
#endif

SystemOfEquations::SystemOfEquations(unsigned int numberOfRows, unsigned int numberOfColumns) {
	this->numberOfRows = numberOfRows;
	this->numberOfColumns = numberOfColumns;
	this->coefficients = vf2D(numberOfRows, vf1D(numberOfColumns, 0.0));
	this->absoluteTerms = vf1D(numberOfRows, 0.0);
}

void SystemOfEquations::setCoefficients(vf2D &originalCoefficients) {
	for(unsigned int r = 0 ; r < originalCoefficients.size() ; ++r) {
		for(unsigned int c = 0 ; c < originalCoefficients[r].size(); ++c) {
			this->coefficients[r][c] = originalCoefficients[r][c];
		}
	}
}

void SystemOfEquations::setAbsoluteTerms(vf1D &originalAbsoluteTerms) {
	for(unsigned int r = 0 ; r < originalAbsoluteTerms.size() ; ++r) {
		this->absoluteTerms[r] = originalAbsoluteTerms[r];
	}
}

string SystemOfEquations::print(bool printLinearCombination) {
	ostringstream output;
	output << endl;
	for(unsigned int rowIdx = 0 ; rowIdx < this->numberOfRows; ++rowIdx) {
		for(unsigned int columnIdx = 0; columnIdx < this->numberOfColumns; ++columnIdx) {
			output << this->coefficients[rowIdx][columnIdx] << " ";
		}
		output << "| " << this->absoluteTerms[rowIdx];
		if (printLinearCombination) {
			output << " : ";
			for(unsigned int combIdx = 0; combIdx < this->linearCombination[rowIdx].size(); ++combIdx) {
				output << linearCombination[rowIdx][combIdx] << " ";
			}
			output << " {";
			for(int termIdx = 0 ; termIdx < this->termCounts[rowIdx].size() ; ++termIdx) {
				output << termCounts[rowIdx][termIdx] << ", ";
			}
			output << "} (" << termFitness[rowIdx] << ")";
		}
		output << endl;
	}
	return output.str();
}

bool SystemOfEquations::operator==(const SystemOfEquations& rhs) {

	for (unsigned int rowIdx = 0 ; rowIdx < this->numberOfRows ; ++ rowIdx) {
		for (unsigned int columnIdx = 0 ; columnIdx < this->numberOfColumns ; ++ columnIdx) {
			if(abs(this->coefficients[rowIdx][columnIdx] - rhs.coefficients[rowIdx][columnIdx]) > EPSILON) {
				return false;
			}
		}

		if(abs(this->absoluteTerms[rowIdx] - rhs.absoluteTerms[rowIdx]) > EPSILON) {
			return false;
		}

	}
	return true;
}

bool SystemOfEquations::operator!=(const SystemOfEquations& rhs) {
	return !this->operator==(rhs);
}

void SystemOfEquations::setTermCounts(vi2D &termCounts) {
	this->termCounts = termCounts;
	this->computeFitness();
}

void SystemOfEquations::computeFitness() {
	this->termFitness = vf1D(termCounts.size(), 0.0);
	for(int rowIdx = 0 ; rowIdx < termCounts.size() ; ++rowIdx) {
		double sumCoeffs = 0;
		for(int colIdx = 0 ; colIdx < this->coefficients[rowIdx].size() ; ++colIdx) {
			sumCoeffs += abs(coefficients[rowIdx][colIdx]);
		}
		double termSum = 0;
		for(int termIdx = 0 ; termIdx < this->termCounts[rowIdx].size() ; ++termIdx) {
			termSum += termCounts[rowIdx][termIdx];
		}
		termFitness[rowIdx] = termSum/(sumCoeffs*sumCoeffs);
		//double x =(double)(termCounts[i][0]);
		//double n =(double)(termCounts[i][0] + termCounts[i][1]);
		//termFitness[i] = sqrt( (0.3 - ( (x / n) * (1.0 - x/n) ) ) / n);
	}
}

void SystemOfEquations::sortEquations(){
	vi1D indices(this->coefficients.size(), 0);
	for(int i=0 ; i < this->coefficients.size() ; ++i) {
		indices[i]=i;
	}
}

int SystemOfEquations::GaussJordanElimination(unsigned int rowLimiter) { 
	if (rowLimiter == 0) { 
		rowLimiter = this->numberOfRows; //number of rows to consider (it replaces numberOfRows in further computations)
	}
	const unsigned int INF = 9999999;

	unsigned int returnCode = GJ_OK;

	this->linearCombination = vf2D(rowLimiter, vf1D (rowLimiter, 0.0));
	for(unsigned int diagIdx = 0; diagIdx < rowLimiter; ++diagIdx) {
		linearCombination[diagIdx][diagIdx] = 1.0;
	}
	vector<unsigned int> occupiedRows(rowLimiter, INF); // rows occupied by pivots
	for(unsigned int columnIdx = 0; columnIdx < this->numberOfColumns; ++columnIdx) {

		unsigned int pivotIdx = INF;

		for(unsigned int rowIdx = 0; rowIdx < rowLimiter; ++rowIdx){
			if ((occupiedRows[rowIdx] == INF) && (abs(this->coefficients[rowIdx][columnIdx]) > EPSILON)){
				pivotIdx = rowIdx;
				occupiedRows[pivotIdx] = columnIdx;
				break;
			}
		}

		if (pivotIdx == INF) {
			returnCode |= GJ_PIVOT_ERROR;
			continue;
		}

		t_coeff pivot = this->coefficients[pivotIdx][columnIdx];
		this->absoluteTerms[pivotIdx] /= pivot;
		for(unsigned int combIdx = 0; combIdx < rowLimiter; ++combIdx){
			linearCombination[pivotIdx][combIdx] /= pivot;
		}
		for(unsigned int c = 0 ; c < this->numberOfColumns ; ++c) {
			this->coefficients[pivotIdx][c] /= pivot;
		}

		for(unsigned int r = 0 ; r < rowLimiter; ++r) {
			if (r == pivotIdx) { continue; }

			if (abs(this->coefficients[r][columnIdx]) > EPSILON) {
				t_coeff tmpCoeff = this->coefficients[r][columnIdx];

				for(unsigned int r2 = 0; r2 < rowLimiter; ++r2){
					linearCombination[r][r2] -= linearCombination[pivotIdx][r2] * tmpCoeff;
				}

				for(unsigned int c = 0 ; c < this->numberOfColumns ; ++c) {
					this->coefficients[r][c] -= this->coefficients[pivotIdx][c] * tmpCoeff;
				}
				this->absoluteTerms[r] -= this->absoluteTerms[pivotIdx] * tmpCoeff;
			}
		}
	}

	return returnCode;
}

vector<int> SystemOfEquations::getParameterIndices() {
	vector<int> indices(numberOfColumns, -1);
	for(unsigned int columnIdx = 0 ; columnIdx < this->numberOfColumns ; ++columnIdx) {
		for(unsigned int rowIdx = 0 ; rowIdx < this->numberOfRows ; ++rowIdx) {
			if (abs(coefficients[rowIdx][columnIdx]) > EPSILON) {
				indices[columnIdx] = rowIdx;
				break;
			}
		}
	}
	return indices;
}
