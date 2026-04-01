#include "householder_search.h"
#include <algorithm>
#include <mutex>
#include <omp.h>

using namespace std;

// FILTER_REL_TOL: relative tolerance added to the candidate pre-filters in
// entryEnumeration to avoid discarding borderline candidates due to floating-point
// rounding error in the filter comparisons.
//
// The three filters each compare a computed float quantity against a threshold T:
//   x1, x3 filters:  computed_val < eps_cond          (T = eps_cond)
//   x2 filter:       computed_val < f_pow * eps_cond   (T = f_pow * eps_cond)
//
// A candidate whose true value is just barely below T might be computed as just
// barely above T due to float rounding, and incorrectly excluded. The filter buffer
// must therefore scale with T, not be a fixed constant. We widen each threshold by
// a relative factor (1 + FILTER_REL_TOL) so the buffer shrinks with eps_cond.
//
// At eps=0.1:   eps_cond~1.25e-3, buffer~1.25e-9  (was 1e-4 with old fixed FPRL)
// At eps=0.01:  eps_cond~1.25e-5, buffer~1.25e-11 (old FPRL was 8x eps_cond!)
// At eps=1e-10: eps_cond~1.25e-21,buffer~1.25e-27 (old FPRL was 8e16x eps_cond!)
//
// 1e-6 provides ~1 million times the actual float rounding error (~2e-16 relative),
// a comfortable margin with no risk of letting genuinely-invalid candidates through.
//
// FPRL must NOT appear in epsTest — that is the strict mathematical acceptance gate.
// The pre-filters are intentionally conservative (they may admit extra candidates),
// and epsTest is what enforces the actual epsilon condition.
const double FILTER_REL_TOL = 1e-6;
extern std::atomic<bool> interrupted;

// Global mutex protecting all cout output. Without this, concurrent epsTest calls
// from the parallel x1 loop produce interleaved/garbled output.
static mutex cout_mutex;

array<ringZ9chi,3> HRSA(double theta, double epsilon, int max_f, double c){
	int f = 0;
	ringZ9chi zero;
	array<ringZ9chi,3> answer;
	vector<ringZ9> x1_cands, x2_cands;
	// OPT: lookup is now a hash map keyed by quad() value instead of a flat vector.
	// This turns the inner x_3 scan from O(|lookup|) to O(1) average, eliminating
	// the dominant cost of the original triple-nested loop.
	unordered_multimap<int,ringZ9> lookup;
	
	complex<double> angle_dir(cos(theta/2.0), sin(-theta/2.0));
	
	answer[0] = zero; answer[1] = zero; answer[2] = zero;
	
	while( f <= max_f){
		// OPT: use three_power() for exact integer arithmetic instead of floating-point pow().
		// static_cast<int>(pow(3,2*f)) can silently return the wrong value for larger f
		// due to floating-point rounding (e.g., 3^10 = 59049.9999... -> 59048).
		int f_pow_sq = three_power(f) * three_power(f);
		ringZ9 f_pow_sq_doub(2*f_pow_sq);
		
		cout << "Enumerating entry candidates for f = " << f << "." << endl;
		entryEnumeration(x1_cands, x2_cands, lookup, theta, epsilon, f, c);
		
		// Parallelise over x1_cands when f >= 2: each x1 is fully independent.
		// lookup is read-only at this point, so concurrent equal_range calls are safe.
		// 'found' signals the first thread to pass epsTest; others check and exit early.
		// Below f=2 the candidate lists are small and thread overhead exceeds the gain,
		// so we fall through to a simple serial loop for f=0 and f=1.
		atomic<bool> found(false);
		mutex answer_mutex;
		array<ringZ9chi,3> local_best;

		#pragma omp parallel for schedule(dynamic) if(f >= 2) \
		    shared(found, answer_mutex, local_best)
		for(int i = 0; i < (int)x1_cands.size(); ++i){

			if(found || interrupted) continue;

			const ringZ9& x_1 = x1_cands[i];
			int q1 = x_1.quad();

			for(const ringZ9& x_2 : x2_cands){

				if(found || interrupted) break;

				int q2 = x_2.quad();
				if(q1 + q2 > 2 * f_pow_sq) continue;

				// The hash map gives us all x3 candidates with quad(x3) == target_norm.
				// quad(x3) is only the constant term of complexConj(x3)*x3; that product
				// is not always a rational integer. We must verify the full ring equality
				// complexConj(x1)*x1 + complexConj(x2)*x2 + complexConj(x3)*x3 == 2*3^{2f}
				// (as a ring element, not just its constant term) to guarantee |u|^2 = 2
				// exactly. The hash map is a fast pre-filter; this exact check culls false
				// positives that would produce non-unitary H matrices.
				int target_norm = 2*f_pow_sq - q1 - q2;
				ringZ9 x3sq = f_pow_sq_doub - x_1.complexConj()*x_1 - x_2.complexConj()*x_2;

				auto range = lookup.equal_range(target_norm);
				for(auto it = range.first; it != range.second; ++it){

					if(found) break;

					const ringZ9& x_3 = it->second;

					// Exact ring check: reject x3 candidates where the ring product
					// has non-zero higher-order terms (i.e. is not a rational integer).
					if(!(x_3.complexConj()*x_3 == x3sq)) continue;

					array<ringZ9chi,3> candidate = {
						ringZ9chi(x_1,f),
						ringZ9chi(x_2,f),
						ringZ9chi(x_3,f)
					};

					// epsTest is lock-free: no cout inside, pure arithmetic.
					if(epsTest(candidate[0], candidate[1], candidate[2], f, angle_dir, epsilon, c)){
						lock_guard<mutex> lk(answer_mutex);
						if(!found){
							found = true;
							local_best = candidate;
						}
					}
				}
			}
		}

		if(found){
			answer = local_best;
			// Log the winning eps_diff now that we're back in the serial region.
			double eps_diff = pow(abs(answer[0].toComplexDouble() - angle_dir),2)
			                + (answer[1] + ringZ9chi(ringZ9(1),0)).abs_val_sq()
			                + answer[2].abs_val_sq();
			cout << "Epsilon Diff. Val.: " << eps_diff
			     << " Eps. Cond.: " << epsilon*epsilon/(8.0*c*c) << endl;
			cout << "Success!" << endl;
			return answer;
		}
		if(interrupted) return answer;
		
		f++;
		x1_cands.clear();
		x2_cands.clear();
		lookup.clear();
	}
	cout << "Failure." << endl;
	return answer;
}

void entryEnumeration(	vector<ringZ9>& x1_cands, 
						vector<ringZ9>& x2_cands, 
						unordered_multimap<int,ringZ9>& lookup, 
						double theta, double epsilon, int f, double c){
	complex<double> target_ang(cos(theta/2.0),sin(-theta/2.0));
	int f_pow = three_power(f);
	
	int A = 4 * f_pow*f_pow;
	
	double eps_cond = epsilon*epsilon/(8.0*c*c);
	
	// OPT: precompute 1/3^f once for converting ringZ9 coefficients to complex doubles
	// without constructing and normalizing a full ringZ9chi object each time.
	double inv_3f = 1.0;
	for (int i = 0; i < f; ++i) inv_3f /= 3.0;

	// Precomputed trig tables matching ringZ9chi::real_part/imag_part static arrays.
	static const double cos_vals[6] = {
		1.0, 0.766044443118978, 0.17364817766693041,
		-0.5, -0.9396926207859083, -0.9396926207859084,
	};
	static const double sin_vals[6] = {
		0.0, 0.6427876096865393, 0.984807753012208,
		0.8660254037844387, 0.3420201433256689, -0.34202014332566866,
	};

	ringZ9 test;	
	int max_a3 = static_cast<int>(ceil(sqrt(A / 3.0)));

	// Parallelise the outermost loop over a3. Each a3 value is fully independent.
	// Threads accumulate into thread-local vectors/maps to avoid contention, then
	// merge into the shared outputs after the parallel region.
	// dynamic scheduling balances load: work per a3 varies widely since inner budgets
	// collapse quickly for large |a3|.
	vector<vector<ringZ9>> tl_x1(omp_get_max_threads());
	vector<vector<ringZ9>> tl_x2(omp_get_max_threads());
	vector<unordered_multimap<int,ringZ9>> tl_lookup(omp_get_max_threads());

	#pragma omp parallel for schedule(dynamic) if(f >= 2) shared(tl_x1, tl_x2, tl_lookup)
	for (int a3 = -max_a3; a3 <= max_a3; ++a3){

		if(interrupted) continue;

		int tid = omp_get_thread_num();
		auto& loc_x1     = tl_x1[tid];
		auto& loc_x2     = tl_x2[tid];
		auto& loc_lookup = tl_lookup[tid];
		
		int budget_a3 = A - 3*a3*a3;
		if (budget_a3 < 0) continue;
		
		int max_a4 = static_cast<int>(ceil(sqrt(budget_a3)));
		for (int a4 = -max_a4; a4 <= max_a4; ++a4){
			int budget_a4 = budget_a3 - 3*a4*a4;
			if (budget_a4 < 0) continue;
			
			int max_a5 = static_cast<int>(ceil(sqrt(budget_a4)));
			for (int a5 = -max_a5; a5 <= max_a5; ++a5) {
				int budget_a5 = budget_a4 - 3*a5*a5;
				if (budget_a5 < 0) continue;

				int max_b = static_cast<int>(ceil(sqrt(budget_a5)));
				for (int b0 = -max_b; b0 <= max_b; ++b0) {
					int budget_b0 = budget_a5 - b0*b0;
					if ((b0 + a3) % 2 != 0 || budget_b0 < 0) continue;
					int a0 = (b0 + a3) / 2;

					int max_b1 = static_cast<int>(ceil(sqrt(budget_b0)));
					for (int b1 = -max_b1; b1 <= max_b1; ++b1) {
						
						if(interrupted) goto next_a3;
						
						int budget_b1 = budget_b0 - b1*b1;
						if ((b1 + a4) % 2 != 0 || budget_b1 < 0) continue;
						int a1 = (b1 + a4) / 2;

						int max_b2 = static_cast<int>(ceil(sqrt(budget_b1)));
						for (int b2 = -max_b2; b2 <= max_b2; ++b2) {
							if ((b2 + a5) % 2 != 0) continue;
							int a2 = (b2 + a5) / 2;

							int arr[9] = {a0, a1, a2, a3, a4, a5, 0, 0, 0};
							ringZ9 local_test(arr);

							double re_raw = 0.0, im_raw = 0.0;
							for (int k = 0; k < 6; ++k) {
								re_raw += cos_vals[k] * arr[k];
								im_raw += sin_vals[k] * arr[k];
							}
							double re = re_raw * inv_3f;
							double im = im_raw * inv_3f;
							double abs_sq = re*re + im*im;

							if( pow(abs(complex<double>(re,im) - target_ang), 2) < eps_cond * (1.0 + FILTER_REL_TOL)){
								loc_x1.push_back(local_test);
							}
							{
								ringZ9chi x2_chi(local_test, f);
								auto x2_shifted = x2_chi + ringZ9chi(ringZ9(1), 0);
								if( x2_shifted.abs_val_sq() < eps_cond * (1.0 + FILTER_REL_TOL) ){
									loc_x2.push_back(local_test);
								}
							}
							if( abs_sq < eps_cond * (1.0 + FILTER_REL_TOL)){
								loc_lookup.emplace(local_test.quad(), local_test);
							}
						} // b2
					} // b1
				} // b0
			} // a5
		} // a4
		next_a3:; // label for interrupted early-exit from b1 loop
	} // a3 (parallel)

	// Merge thread-local results into the shared output containers.
	for(int t = 0; t < omp_get_max_threads(); ++t){
		x1_cands.insert(x1_cands.end(), tl_x1[t].begin(), tl_x1[t].end());
		x2_cands.insert(x2_cands.end(), tl_x2[t].begin(), tl_x2[t].end());
		lookup.insert(tl_lookup[t].begin(), tl_lookup[t].end());
	}
	
	// Sort x1_cands by distance to e^(i*theta/2) and x2_cands by distance to -1.
	// This causes the algorithm to try the most promising candidates first, which
	// tends to find a passing triple earlier — often at a lower f-level than the
	// unsorted original.
	//
	// Both sorts are safe with respect to correctness: they reorder candidates but
	// never remove any. For a given f, every valid (x1,x2,x3) triple is still
	// reachable. The sort cannot cause the algorithm to advance to a higher f than
	// necessary. It may return a DIFFERENT valid triple at the same f (potentially
	// with a larger eps_diff than the original), but that triple still passes epsTest.
	sort(x1_cands.begin(), x1_cands.end(), [&target_ang](const ringZ9& a, const ringZ9& b) {
		return abs(a.toComplexDouble() - target_ang) < abs(b.toComplexDouble() - target_ang);
	});
	complex<double> neg_one(-1.0, 0.0);
	sort(x2_cands.begin(), x2_cands.end(), [&neg_one](const ringZ9& a, const ringZ9& b) {
		return abs(a.toComplexDouble() - neg_one) < abs(b.toComplexDouble() - neg_one);
	});
	// lookup is a hash map; order is irrelevant for O(1) lookups by key.
	
	cout << "Finished enumerating entry possibilities for f = " << f << "." << endl;
}


bool epsTest(ringZ9chi u_1, ringZ9chi u_2, ringZ9chi u_3, int f, complex<double> angle_dir, double epsilon, double c){
	
	double eps_diff = pow(abs(u_1.toComplexDouble() - angle_dir),2) 
	                + (u_2 + ringZ9chi(ringZ9(1),0)).abs_val_sq() 
	                + u_3.abs_val_sq();
	double eps_cond = epsilon*epsilon/(8.0*c*c);

	// No cout here: this function is called from parallel threads and must be
	// lock-free. The caller logs the result once a solution is confirmed.
	return eps_diff < eps_cond;
}


int three_power(int n){
	
	int prod = 1;
	
	if(n < 0){
		cout << "Negative argument passed to three_power(int n)" << endl;
		return 0;
	}
	
	for(int i = 0; i < n; i++){
		prod = prod*3;
	}
	
	return prod;
	
}

void handleCtrlC(int sig) {
    interrupted = true;
}

void matrixFrobeniusCheck(const std::array<ringZ9chi,3>& u, double theta, double epsilon){
	// Build H = I_3 - u * u^dagger as a 3x3 complex<double> matrix.
	// u^dagger means conjugate-transpose; since u is a column vector,
	// H_{ij} = delta_{ij} - u_i * conj(u_j).
	using cd = complex<double>;
	array<cd,3> uc = {
		u[0].toComplexDouble(),
		u[1].toComplexDouble(),
		u[2].toComplexDouble()
	};

	// H[i][j] = delta_ij - uc[i]*conj(uc[j])
	cd H[3][3];
	for(int i = 0; i < 3; ++i)
		for(int j = 0; j < 3; ++j)
			H[i][j] = (i==j ? cd(1,0) : cd(0,0)) - uc[i] * conj(uc[j]);

	// X_{(0,1)} = [[0,1,0],[1,0,0],[0,0,1]]  (SWAP of rows 0 and 1).
	// This is the fixed correction unitary such that X_{(0,1)} * H ~ R_{(0,1)}^Z(theta)
	// when u ~ [e^{itheta/2}, -1, 0]. It is theta-independent.
	// XH[i][j] = sum_k X[i][k] * H[k][j]; with X = SWAP this exchanges rows 0 and 1:
	cd XH[3][3];
	for(int j = 0; j < 3; ++j){
		XH[0][j] = H[1][j];  // row 0 of X*H = row 1 of H
		XH[1][j] = H[0][j];  // row 1 of X*H = row 0 of H
		XH[2][j] = H[2][j];  // row 2 unchanged
	}

	// R_{(0,1)}^Z(theta) = Diag(e^{-itheta/2}, e^{itheta/2}, 1).
	// Takes the EXTERNAL theta (before negation by the test driver).
	cd R[3][3] = {};
	R[0][0] = polar(1.0, -theta/2.0);
	R[1][1] = polar(1.0,  theta/2.0);
	R[2][2] = cd(1.0, 0.0);

	// Frobenius norm: ||XH - R||_F = sqrt(sum_{ij} |XH_{ij} - R_{ij}|^2)
	double frob_sq = 0.0;
	for(int i = 0; i < 3; ++i)
		for(int j = 0; j < 3; ++j){
			cd diff = XH[i][j] - R[i][j];
			frob_sq += diff.real()*diff.real() + diff.imag()*diff.imag();
		}
	double frob = sqrt(frob_sq);

	cout << "Matrix Frobenius distance ||X_(0,1) H - R_(0,1)^Z(theta)||_F = "
	     << frob << "   (target epsilon = " << epsilon << ", passes: "
	     << (frob < epsilon ? "YES" : "NO") << ")" << endl;
}
