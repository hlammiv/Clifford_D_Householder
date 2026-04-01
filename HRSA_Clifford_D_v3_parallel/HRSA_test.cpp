#include<iostream>
#include "householder_search.h"
#include<vector>
#include<array>
#include<cmath>
#include<complex>
#include<numeric>
#include<atomic>
#include<csignal>
#include "cyclotomic_int9.h"
#include "Z9chi.h"
//#include "z3++.h"


std::atomic<bool> interrupted(false); // For clean exiting with Ctrl+C if execution takes too long

using namespace std;

int main(int argc, char* argv[]){
	
	signal(SIGINT, handleCtrlC);
	
	if(argc < 4 || argc > 5){
		cout << "INVALID INPUT: Do ./ESA_test 'theta' 'epsilon' 'maximum f'  [optional] 'contraction factor'." << endl;
	}
	
	double theta = atof(argv[1]);
	double epsilon = atof(argv[2]);
	int max_f = atoi(argv[3]);
	double c = 1.0; // if no correction factor specified, correction factor is taken to be just 1.
	
	if(argc == 5){
		c = atof(argv[4]);
		
		if (c <= 0.0 || c > 1.0){
			cout << "Invalid choice of c. c must lie in interval (0,1]." << endl;
			exit(0); 
		}
	}
	
	
	

	array<ringZ9chi,3> output2 = HRSA(-1*theta, epsilon, max_f, c);
	
	
	
	cout << "x_1: ";
	output2[0].print(); 
	cout << " " << output2[0].abs_val() << "\n";
	
	cout << "x_2: ";
	output2[1].print(); 
	cout << " " << output2[1].abs_val() << "\n";
	
	cout << "x_3: ";
	output2[2].print(); 
	cout << " " << output2[2].abs_val() << "\n";
		
		
	cout << "Target x_1 = e^(i theta/2) : " << cos(theta/2.0) << " + i" << sin(theta/2) << endl; 
	cout << "Found x_1 (should be close to Target x_1): " << real(output2[0].toComplexDouble()) << " + i" << imag(output2[0].toComplexDouble()) << endl;
		
	cout << "Found x_2 (should be close to -1): " << real(output2[1].toComplexDouble()) << " + i" << imag(output2[1].toComplexDouble()) << endl;
	
	cout << "Found x_3 (should be close to 0): " << real(output2[2].toComplexDouble()) << " + i" << imag(output2[2].toComplexDouble()) << endl;
	
	cout 	<< "Squared Frobenius norm of vector (should be exactly 2): " 
			<< output2[0].abs_val_sq() + output2[1].abs_val_sq() + output2[2].abs_val_sq() << endl;

	// Full matrix check: compute ||X_{(0,1)} (I - u u^dagger) - R_{(0,1)}^Z(theta)||_F
	// and verify it is less than epsilon. This is the actual guarantee the algorithm
	// is designed to provide, as opposed to the proxy vector condition in epsTest.
	matrixFrobeniusCheck(output2, theta, epsilon);

	return 0;
}
