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
		cout << "Usage: ./HRSA_tester theta epsilon max_f [c] [flags]" << endl;
		cout << "  theta            : target rotation angle" << endl;
		cout << "  epsilon          : desired precision" << endl;
		cout << "  max_f            : maximum f level to search" << endl;
		cout << "  c                : contraction factor in (0,1], default 1.0" << endl;
		cout << "  --max-solns N    : collect N HRSA candidates, pick lowest D-count" << endl;
		cout << "  --max-direct K   : direct search up to K D-gates (default 2, 0=Clifford only)" << endl;
		cout << "  --no-direct      : skip direct search, go straight to HRSA" << endl;
		return 1;
	}
	
	double theta = atof(argv[1]);
	double epsilon = atof(argv[2]);
	int max_f = atoi(argv[3]);
	double c = 1.0;
	int max_solns = 1;  // default: original behavior
	int max_direct = 2; // default: direct search up to k=2
	
	// Parse optional positional arg (c) and flags
	for(int i = 4; i < argc; ++i){
		if(strcmp(argv[i], "--max-solns") == 0 && i+1 < argc){
			max_solns = atoi(argv[++i]);
			if(max_solns < 1) max_solns = 1;
		} else if(strcmp(argv[i], "--max-direct") == 0 && i+1 < argc){
			max_direct = atoi(argv[++i]);
			if(max_direct < 0) max_direct = 0;
		} else if(strcmp(argv[i], "--no-direct") == 0){
			max_direct = -1;  // skip direct search entirely
		} else if(argv[i][0] != '-'){
			// positional: must be c
			c = atof(argv[i]);
			if (c <= 0.0 || c > 1.0){
				cout << "Invalid choice of c. c must lie in interval (0,1]." << endl;
				return 1;
			}
		}
	}

	// Track results across all phases for the summary block
	string method = "none";
	int result_D = -1;
	double result_frob = -1.0;
	int result_f = -1;
	bool found = false;

	// --- Phase 1: Direct gate-sequence search for low D-count ---
	if (max_direct >= 0 && !found) {
		cout << "\n=== Phase 1: Direct search (D-count 0-" << max_direct << ") ===" << endl;
		DirectSearchResult dsr = directSearch(theta, epsilon, max_direct);
	
		if (dsr.success) {
			method = "Direct";
			result_D = dsr.D_count;
			result_frob = dsr.frob_dist;
			result_f = 0;
			found = true;
		} else {
			cout << "\nNo low-D solution found via direct search." << endl;
		}
	}

	// --- Phase 2: Diagonal synthesis ---
	if (!found) {
		cout << "\n=== Phase 2: Diagonal synthesis ===" << endl;
		DiagSearchResult dsr = diagSearch(theta, epsilon, max_f);
		if (dsr.success) {
			// Decompose to get D-count
			Mat3 V;
			V.m[0][0] = dsr.diag_entries[0];
			V.m[1][1] = dsr.diag_entries[1];
			V.m[2][2] = dsr.diag_entries[2];
			DecompResult dr = decompose(V, true);

			method = "Diagonal(f=" + to_string(dsr.f_level) + ")";
			result_D = dr.success ? dr.D_count : -1;
			result_frob = dsr.frob_dist;
			result_f = dsr.f_level;
			found = true;
		} else {
			cout << "\nDiagonal synthesis did not find solution. Falling back to HRSA..." << endl;
		}
	}

	// --- Phase 3: HRSA Householder search ---
	if (!found) {
		cout << "\n=== Phase 3: HRSA search ===" << endl;
		array<ringZ9chi,3> output2;
		if(max_solns > 1){
			output2 = HRSA_bestD(-1*theta, epsilon, max_f, c, max_solns);
		} else {
			output2 = HRSA(-1*theta, epsilon, max_f, c);
		}

		// Check if HRSA found something (non-zero vector)
		bool hrsa_found = !(output2[0].isZero() && output2[1].isZero() && output2[2].isZero());

		if (hrsa_found) {
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

			Mat3 V = buildUnitary(output2);
			DecompResult dr = decompose(V, true);

			// Compute Frobenius distance for summary

			// Get f level from the exponent
			result_f = output2[0].getExp();

			method = "HRSA(f=" + to_string(result_f) + ")";
			result_D = dr.success ? dr.D_count : -1;
			// Use Matrix Frobenius distance (computed by matrixFrobeniusCheck)
			// Recompute it here for the summary
			{
				using cd = complex<double>;
				array<cd,3> uc = {output2[0].toComplexDouble(), output2[1].toComplexDouble(), output2[2].toComplexDouble()};
				cd H[3][3];
				for(int i=0;i<3;++i) for(int j=0;j<3;++j)
					H[i][j] = (i==j ? cd(1,0) : cd(0,0)) - uc[i]*conj(uc[j]);
				cd XH[3][3];
				for(int j=0;j<3;++j){ XH[0][j]=H[1][j]; XH[1][j]=H[0][j]; XH[2][j]=H[2][j]; }
				cd R[3][3]={};
				R[0][0]=polar(1.0,-theta/2.0); R[1][1]=polar(1.0,theta/2.0); R[2][2]=cd(1,0);
				double frob_sq=0;
				for(int i=0;i<3;++i) for(int j=0;j<3;++j){
					cd d=XH[i][j]-R[i][j]; frob_sq+=d.real()*d.real()+d.imag()*d.imag();
				}
				result_frob = sqrt(frob_sq);
			}
			found = true;
		}
	}

	// =========================================================
	//  PARSEABLE SUMMARY — parse_results.py reads these lines
	// =========================================================
	cout << "\n=== Summary ===" << endl;
	if (found) {
		cout << "Epsilon Diff. Val.: " << result_frob << " Eps. Cond.: " << epsilon*epsilon/(8.0*c*c) << endl;
		cout << "Matrix Frobenius distance ||X_(0,1) H - R_(0,1)^Z(theta)||_F = "
		     << result_frob << "   (target epsilon = " << epsilon << ", passes: "
		     << (result_frob < epsilon ? "YES" : "NO") << ")" << endl;
		cout << "Total D gates:    " << result_D << endl;
		cout << "Method: " << method << endl;
		cout << "Success!" << endl;
	} else {
		cout << "Failure." << endl;
	}

	return 0;
}
