#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <vector>
#include <iostream>
#include <fstream>


using namespace std;

const int D = 600;

typedef complex<double> cdouble;

const cdouble I(0.0,1.0);

const double mass = 1.0;
const double omega = 1.0;
const double alpha = 2.0;
const double hbar = 1.0;//6.582119 * pow(10,-16); 
const double xstart = -1.0;
const double x0 = -4.0;
const double xD = 4.0;
const double deltaX = (xD - x0)/D;
const double T = 2.0 * M_PI;
const double deltaT = T/256.0;
const double A = 1;
const double B = 2;
const double C = 4;

vector<cdouble> vecMatMul(vector<cdouble> v, vector<cdouble> m);
void printCVec(vector<cdouble> v, string fileName);
vector<cdouble> vecNorm(vector<cdouble> v);
double getAvgX(vector<cdouble> v);
double getV(vector<cdouble> v);
double getK(vector<cdouble> v);
double V(double x);

ofstream oFile;

int main() {
	//Init discrete psi vector
	vector<cdouble> psiDisc;
	cdouble tempPsiVal;
	double x;
	for(int i=0; i<D; i++) {
		x = x0 + i*deltaX;
		tempPsiVal = pow((pow(8,.5)/M_PI),2)*exp(-1.0*pow(2,.5)*pow((x-xstart),2.0));
		psiDisc.push_back(tempPsiVal); // * conj(tempPsiVal));
	}
	
	cdouble A = pow((2*M_PI*hbar*I*deltaT/mass),.5);
	vector<cdouble> transitionMatrix;
	double xi;
	double xj;
	for(int i=0; i<D; i++) {
		for(int j=0; j<D; j++) {
			xi = x0 + i*deltaX;
			xj = x0 + j*deltaX;
			transitionMatrix.push_back((1.0/A) * exp( (I*deltaT/hbar) * ( (.5*mass*pow(xj-xi,2)/pow(deltaT,2)) - V((xi+xj)/2) )));
		}
	}
	vector<cdouble> psi = psiDisc;
	for(int n=0; n<=350; n++) {
		psi = vecNorm(psi);
		string file = "out" + to_string(n) + ".txt";
		printCVec(psi, file);
		
		for(int i=0; i<2; i++) {
			psi = vecMatMul(psi,transitionMatrix);
			for(int i=0; i<D; i++) {
				psi.at(i) = psi.at(i)*deltaX;
			}
			/*
			oFile.open("out0.txt", ios::app);
			oFile << n << "	" << getV(psi) << endl;
			oFile.close();
			oFile.open("out1.txt", ios::app);
			oFile << n << "	" << getK(psi) << endl;
			oFile.close();
			oFile.open("out2.txt", ios::app);
			oFile << n << "	" << getV(psi) + getK(psi) << endl;
			oFile.close();
			cout << n << "	" << getAvgX(psi) << endl;
			*/
			//psi = vecNorm(psi);
		}
		

	}
	
	//printCVec(psi);
	//printCVec(transitionMatrix);
	
	
}

double V(double x) {
	//return(pow(x,4) - 2.0*pow(x,2) + 1);
	return( pow( (1.0/pow(A,2)*pow(x,2) - pow(B,.5)) , 2) );
}

vector<cdouble> vecMatMul(vector<cdouble> v, vector<cdouble> m) {
	vector<cdouble> newVec;
	for(int i=0; i<D; i++) {
		cdouble sum = (0.0,0.0);
		for(int j=0; j<D; j++) {
			sum+= m.at(i*D + j)*v.at(j);
		}
		newVec.push_back(sum);
	}
	return newVec;
}

void printCVec(vector<cdouble> v, string fileName) {
	double i = 0;
	double x;
	for(cdouble val : v) {
		x = x0 + i*deltaX;
		i++;
		oFile.open(fileName, ios::app);
		oFile << x << "	" << real(val*conj(val)) << endl;
		oFile.close();
	}
}

vector<cdouble> vecNorm(vector<cdouble> v) {
	double sum = 0.0;
	for(val : v) {
		sum += norm(val)*
		deltaX;
	}
	double mag = pow(sum,.5);
	for(int i=0; i<D; i++) {
		v.at(i) = v.at(i) / mag;
	}
	return v;
}

double getAvgX(vector<cdouble> v) {
	int currMaxProb = 0;
	for(int i=0; i<D; i++) {
		//cout << real(v.at(i)*conj(v.at(i))) << "	" << real(v.at(currMaxProb)*conj(v.at(currMaxProb))) << endl;
		if( real(v.at(i)*conj(v.at(i))) > real(v.at(currMaxProb)*conj(v.at(currMaxProb))) ) {
			currMaxProb = i;
		}
	}
	return(x0 + currMaxProb*deltaX);
}

double getV(vector<cdouble> v){
	double sum = 0.0;
	for(int i=0;i<D;i++) {
		double x = xstart + i*deltaX;
		sum += V(x) * real(v.at(i) * conj(v.at(i))) * deltaX;
	}
	return sum; //* pow(deltaX, 3);

}
double getK(vector<cdouble> v){
	double sum = 0.0;
	for(int i=1;i<D-1;i++) {
		cdouble minus;
		cdouble plus;
		
		minus = v.at(i-1);
		
		plus = v.at(i+1);

		
		sum += (real(conj(v.at(i)) * ( (minus + plus - (v.at(i) * 2.0)) / deltaX ) ) );

	}
	return sum * -1.0 * pow(hbar,2) / (2*mass);


}

