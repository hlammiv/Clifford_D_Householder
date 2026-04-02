#include<iostream>
#include "householder_search.h"
#include "decompose.h"
#include<vector>
#include<array>
#include<cmath>
#include<complex>
#include<numeric>
#include<atomic>
#include<csignal>
#include<cstring>
#include "cyclotomic_int9.h"
#include "Z9chi.h"


std::atomic<bool> interrupted(false); // For clean exiting with Ctrl+C if execution takes too long

using namespace std;

int main(int argc, char* argv[]){
	
	signal(SIGINT, handleCtrlC);
	
	if(argc < 4){
		cout << "Usage: ./HRSA_tester theta epsilon max_f [c] [--max-solns N]" << endl;
		cout << "  theta          : target rotation angle" << endl;
		cout << "  epsilon        : desired precision" << endl;
		cout << "  max_f          : maximum f level to search" << endl;
		cout << "  c              : contraction factor in (0,1], default 1.0" << endl;
		cout << "  --max-solns N  : collect N candidates at winning f, pick lowest D-count" << endl;
		cout << "                   (default: 1 = original behavior, first solution)" << endl;
		return 1;
	}
	
	double theta = atof(argv[1]);
	double epsilon = atof(argv[2]);
	int max_f = atoi(argv[3]);
	double c = 1.0;
	int max_solns = 1;  // default: original behavior
	
	// Parse optional positional arg (c) and --max-solns flag
	for(int i = 4; i < argc; ++i){
		if(strcmp(argv[i], "--max-solns") == 0 && i+1 < argc){
			max_solns = atoi(argv[++i]);
			if(max_solns < 1) max_solns = 1;
		} else if(argv[i][0] != '-'){
			// positional: must be c
			c = atof(argv[i]);
			if (c <= 0.0 || c > 1.0){
				cout << "Invalid choice of c. c must lie in interval (0,1]." << endl;
				return 1;
			}
		}
	}

	array<ringZ9chi,3> output2;
	if(max_solns > 1){
		output2 = HRSA_bestD(-1*theta, epsilon, max_f, c, max_solns);
	} else {
		output2 = HRSA(-1*theta, epsilon, max_f, c);
	}
	
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
	
	cout << "Squared Frobenius norm of vector (should be exactly 2): " 
	     << output2[0].abs_val_sq() + output2[1].abs_val_sq() + output2[2].abs_val_sq() << endl;

	matrixFrobeniusCheck(output2, theta, epsilon);

	// --- D-gate decomposition ---
	cout << "\n=== Gate Decomposition ===" << endl;
	Mat3 V = buildUnitary(output2);
	
	cout << "sde_chi of V[0][0]: " << sdeChiFull(V.m[0][0]) << endl;
	
	DecompResult dr = decompose(V);
	if (dr.success) {
		cout << "Decomposition successful!" << endl;
		cout << "  sde_chi of input: " << dr.sde_chi << endl;
		cout << "  Total D gates:    " << dr.D_count << endl;
		cout << "  Peeling steps:    " << dr.steps.size() << endl;
		for (size_t i = 0; i < dr.steps.size(); ++i) {
			const GateStep& gs = dr.steps[i];
			cout << "    Step " << i+1 << ": H D("
			     << gs.a0 << "," << gs.a1 << "," << gs.a2
			     << ") Dgate^" << gs.eps << " X^" << gs.delta << endl;
		}
	} else {
		cout << "Decomposition failed at sde_chi=" << dr.sde_chi << endl;
	}

	return 0;
}
