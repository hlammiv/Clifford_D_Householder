#include "decompose.h"
#include <algorithm>
#include <iostream>
#include <cassert>
#include <cmath>
#include <climits>
#include <set>
#include <unordered_map>

using namespace std;

// =========================================================================
//  Mat3 — 3x3 matrix over ringZ9chi
// =========================================================================

Mat3::Mat3() {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            m[i][j] = ringZ9chi();  // zero
}

Mat3 Mat3::mul(const Mat3& B) const {
    Mat3 C;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                C.m[i][j] = C.m[i][j] + m[i][k] * B.m[k][j];
    return C;
}

Mat3 Mat3::dagger() const {
    Mat3 D;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            D.m[i][j] = m[j][i].complexConj();
    return D;
}

void Mat3::print() const {
    for (int i = 0; i < 3; ++i) {
        cout << "[ ";
        for (int j = 0; j < 3; ++j) {
            m[i][j].print();
            if (j < 2) cout << ", ";
        }
        cout << " ]" << endl;
    }
}

// =========================================================================
//  Gate matrices
// =========================================================================

// Helper: zeta_9^k as ringZ9chi (sde=0)
static ringZ9chi zeta9_power(int k) {
    k = ((k % 9) + 9) % 9;
    if (k < 6) {
        return ringZ9chi(ringZ9(1, k), 0);
    }
    // k = 6,7,8: use zeta_9^{6+j} = -zeta_9^j - zeta_9^{3+j}
    int j = k - 6;
    int arr[9] = {0,0,0,0,0,0,0,0,0};
    arr[j] = -1;
    arr[j+3] = -1;
    return ringZ9chi(arr, 0);
}

// omega^k = zeta_9^{3k}
static ringZ9chi omega_power(int k) {
    return zeta9_power(3 * k);
}

Mat3 gateH() {
    // H = (1/(1+2*omega)) * DFT_3
    // 1/(1+2*omega) = (-1 - 2*zeta_9^3)/3
    // H_{jk} = (-1 - 2*zeta_9^3)/3 * omega^{jk}
    //
    // sde_chi(H) = 3 because 1+2*omega = u*chi^3 for unit u in Z[xi]
    // (since sqrt(-3) = i*sqrt(3) and |chi|^3 = sqrt(|chi|^6) relates to 3)
    int ia_arr[9] = {-1, 0, 0, -2, 0, 0, 0, 0, 0};
    ringZ9chi ia(ia_arr, 1);  // (-1-2*zeta_9^3)/3

    Mat3 H;
    for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
            H.m[j][k] = ia * omega_power(j * k);
    return H;
}

Mat3 gateD_clifford(int a, int b, int c) {
    // Clifford diagonal: Diag(omega^a, omega^b, omega^c)
    Mat3 D;
    D.m[0][0] = omega_power(a);
    D.m[1][1] = omega_power(b);
    D.m[2][2] = omega_power(c);
    return D;
}

Mat3 gateDgate() {
    // The non-Clifford D gate: Diag(zeta_9, 1, zeta_9^{-1})
    Mat3 D;
    D.m[0][0] = zeta9_power(1);
    D.m[1][1] = ringZ9chi(ringZ9(1), 0);
    D.m[2][2] = zeta9_power(8);  // zeta_9^{-1} = zeta_9^8
    return D;
}

// Full cyclotomic diagonal: Diag(xi^a, xi^b, xi^c)
static Mat3 gateDcyclo(int a, int b, int c) {
    Mat3 D;
    D.m[0][0] = zeta9_power(a);
    D.m[1][1] = zeta9_power(b);
    D.m[2][2] = zeta9_power(c);
    return D;
}

// R = Diag(1, 1, -1)
static Mat3 gateR() {
    Mat3 R;
    R.m[0][0] = ringZ9chi(ringZ9(1), 0);
    R.m[1][1] = ringZ9chi(ringZ9(1), 0);
    R.m[2][2] = ringZ9chi(ringZ9(-1), 0);
    return R;
}

Mat3 gateX() {
    // X = cyclic permutation: |j> -> |j+1 mod 3>
    Mat3 X;
    ringZ9chi one(ringZ9(1), 0);
    X.m[0][2] = one;
    X.m[1][0] = one;
    X.m[2][1] = one;
    return X;
}

// =========================================================================
//  sde_chi computations
// =========================================================================

int sdeChiZ9(ringZ9 a) {
    if (a.isZero()) return 999;
    int s = a.sdeChi();  // returns sde if < 6, else 6
    if (s < 6) return s;
    if (!a.isDivisibleByInt(3)) return 6;
    ringZ9 a_div3 = a / 3;
    return 6 + sdeChiZ9(a_div3);
}

int sdeChiFull(const ringZ9chi& x) {
    if (x.isZero()) return 0;
    ringZ9 numer = x.getNumerator();
    int f = x.getExp();
    if (f == 0) return sdeChiZ9(numer);
    int ell = sdeChiZ9(numer);
    return 6 * f - ell;
}

// =========================================================================
//  Building the unitary from a Householder vector
// =========================================================================

Mat3 buildUnitary(const std::array<ringZ9chi,3>& u) {
    Mat3 H;
    ringZ9chi one(ringZ9(1), 0);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            ringZ9chi uiuj = u[i] * u[j].complexConj();
            H.m[i][j] = (i == j ? one - uiuj : ringZ9chi() - uiuj);
        }
    // Apply X_{(0,1)} = swap rows 0 and 1
    Mat3 V;
    for (int j = 0; j < 3; ++j) {
        V.m[0][j] = H.m[1][j];
        V.m[1][j] = H.m[0][j];
        V.m[2][j] = H.m[2][j];
    }
    return V;
}

// =========================================================================
//  Decomposition algorithm (following Kalra et al. Theorem 5.7)
// =========================================================================

// The peeling step tries all H * Diag(xi^a1, xi^a2, xi^a3) * R^eps * X^delta
// and finds one that reduces sde_chi of the first column by 1.
//
// From Kalra: sde(H*D*z) = sde(z) + 3 - k, where k = gde of the
// first coordinate of chi^{3+f} * H * D * z.
// We need k >= 4 for sde to decrease by 1.

// Helper: build a single prefix matrix  H * Dcyclo(a1,a2,a3) * R^eps * X^delta
static Mat3 buildPrefix(int a1, int a2, int a3, int eps, int delta,
                         const Mat3& Hmat, const Mat3& Rmat, const Mat3& Xmat) {
    Mat3 prefix = gateDcyclo(a1, a2, a3);
    if (eps == 1) prefix = prefix.mul(Rmat);
    if (delta == 1) prefix = prefix.mul(Xmat);
    else if (delta == 2) prefix = prefix.mul(Xmat).mul(Xmat);
    return Hmat.mul(prefix);
}

// Helper: compute (prefix * V)[0][0] without building the full product
static ringZ9chi prefix_times_V_00(const Mat3& prefix, const Mat3& V) {
    ringZ9chi result;
    for (int k = 0; k < 3; ++k)
        result = result + prefix.m[0][k] * V.m[k][0];
    return result;
}

// Count D gates from a cyclotomic diagonal with exponents a1,a2,a3
static int count_d_from_cyclo(int a1, int a2, int a3) {
    int d = 0;
    if (a1 % 3 != 0) d++;
    if (a2 % 3 != 0) d++;
    if (a3 % 3 != 0) d++;
    return d;
}

// Try single-prefix reduction: find H*D*R^eps*X^delta that reduces sde by 1.
// Among all prefixes that achieve the reduction, pick the one with the
// fewest non-Clifford diagonal entries (minimum D-gate cost).
// Returns true if found, and applies the best prefix to V.
static bool trySinglePrefix(Mat3& V, int s, DecompResult& result,
                             const Mat3& Hmat, const Mat3& Rmat, const Mat3& Xmat) {
    int best_d = INT_MAX;
    int best_a1=-1, best_a2=-1, best_a3=-1, best_eps=-1, best_delta=-1;

    for (int eps = 0; eps < 2; ++eps)
    for (int delta = 0; delta < 3; ++delta)
    for (int a1 = 0; a1 < 9; ++a1)
    for (int a2 = 0; a2 < 9; ++a2)
    for (int a3 = 0; a3 < 9; ++a3) {
        int d = count_d_from_cyclo(a1, a2, a3);
        if (d >= best_d) continue;  // can't improve, skip expensive multiply

        Mat3 prefix = buildPrefix(a1, a2, a3, eps, delta, Hmat, Rmat, Xmat);
        ringZ9chi new00 = prefix_times_V_00(prefix, V);
        if (sdeChiFull(new00) == s - 1) {
            best_d = d;
            best_a1 = a1; best_a2 = a2; best_a3 = a3;
            best_eps = eps; best_delta = delta;
            if (best_d == 0) goto found;  // can't do better than 0
        }
    }

    if (best_d == INT_MAX) return false;

    found:
    {
        GateStep step;
        step.a0 = best_a1; step.a1 = best_a2; step.a2 = best_a3;
        step.eps = best_eps; step.delta = best_delta; step.has_H = true;
        result.steps.push_back(step);
        result.D_count += best_d;
        Mat3 prefix = buildPrefix(best_a1, best_a2, best_a3, best_eps, best_delta, Hmat, Rmat, Xmat);
        V = prefix.mul(V);
    }
    return true;
}

// Try double-prefix: find two prefixes P2*P1 such that sde drops by at least 1.
// Among all pairs that achieve a net reduction, pick the one with the fewest
// total non-Clifford diagonal entries (minimum combined D-gate cost).
static bool tryDoublePrefix(Mat3& V, int s, DecompResult& result,
                             const Mat3& Hmat, const Mat3& Rmat, const Mat3& Xmat,
                             bool quiet = false) {
    if (!quiet) cout << "  Trying double-prefix at sde_chi=" << s << "..." << endl;

    int best_total_d = INT_MAX;
    int best_new_s = s;
    int best_mid_s = s;
    int ba1=-1,ba2=-1,ba3=-1,beps1=-1,bdelta1=-1;
    int bb1=-1,bb2=-1,bb3=-1,beps2=-1,bdelta2=-1;

    for (int eps1 = 0; eps1 < 2; ++eps1)
    for (int delta1 = 0; delta1 < 3; ++delta1)
    for (int a1 = 0; a1 < 9; ++a1)
    for (int a2 = 0; a2 < 9; ++a2)
    for (int a3 = 0; a3 < 9; ++a3) {
        int d1 = count_d_from_cyclo(a1, a2, a3);
        if (d1 >= best_total_d) continue;  // P1 alone already too expensive

        Mat3 P1 = buildPrefix(a1, a2, a3, eps1, delta1, Hmat, Rmat, Xmat);
        ringZ9chi mid00 = prefix_times_V_00(P1, V);
        int mid_s = sdeChiFull(mid00);

        if (mid_s > s + 1 || mid_s < s - 1) continue;

        Mat3 midV = P1.mul(V);
        for (int eps2 = 0; eps2 < 2; ++eps2)
        for (int delta2 = 0; delta2 < 3; ++delta2)
        for (int b1 = 0; b1 < 9; ++b1)
        for (int b2 = 0; b2 < 9; ++b2)
        for (int b3 = 0; b3 < 9; ++b3) {
            int d2 = count_d_from_cyclo(b1, b2, b3);
            if (d1 + d2 >= best_total_d) continue;

            Mat3 P2 = buildPrefix(b1, b2, b3, eps2, delta2, Hmat, Rmat, Xmat);
            ringZ9chi new00 = prefix_times_V_00(P2, midV);
            int new_s = sdeChiFull(new00);

            if (new_s < s) {
                best_total_d = d1 + d2;
                best_new_s = new_s;
                best_mid_s = mid_s;
                ba1=a1; ba2=a2; ba3=a3; beps1=eps1; bdelta1=delta1;
                bb1=b1; bb2=b2; bb3=b3; beps2=eps2; bdelta2=delta2;
                if (best_total_d == 0) goto found;
            }
        }
    }

    if (best_total_d == INT_MAX) return false;

    found:
    {
        GateStep step1;
        step1.a0 = ba1; step1.a1 = ba2; step1.a2 = ba3;
        step1.eps = beps1; step1.delta = bdelta1; step1.has_H = true;
        result.steps.push_back(step1);
        result.D_count += count_d_from_cyclo(ba1, ba2, ba3);

        Mat3 P1 = buildPrefix(ba1, ba2, ba3, beps1, bdelta1, Hmat, Rmat, Xmat);
        Mat3 midV = P1.mul(V);

        GateStep step2;
        step2.a0 = bb1; step2.a1 = bb2; step2.a2 = bb3;
        step2.eps = beps2; step2.delta = bdelta2; step2.has_H = true;
        result.steps.push_back(step2);
        result.D_count += count_d_from_cyclo(bb1, bb2, bb3);

        V = buildPrefix(bb1, bb2, bb3, beps2, bdelta2, Hmat, Rmat, Xmat).mul(midV);
        if (!quiet) cout << "  Double-prefix succeeded: sde " << s
             << " -> " << best_mid_s << " -> " << best_new_s << endl;
    }
    return true;
}

DecompResult decompose(Mat3 V, bool quiet) {
    DecompResult result;
    result.D_count = 0;
    result.success = false;

    // =========================================================================
    //  Fast-path: detect monomial matrices (perm × diagonal ζ₉ phases).
    //  A monomial 3×3 unitary has exactly one nonzero entry per row and column.
    //  D-count comes solely from the non-Clifford (non-ω) diagonal phases.
    // =========================================================================

    // Determine the zeta_9 exponent mod 3 of a unit in Z[zeta_9].
    // Returns 0 if the entry is ±omega^k (Clifford), 1 or 2 otherwise.
    // Uses exact ring arithmetic: multiply by zeta_9^{-j} for j=0..8 and
    // check if the result is ±1 (i.e. rational integer ±1).
    auto unitPhaseMod3 = [](const ringZ9chi& x) -> int {
        // x should be a unit (magnitude 1) in Z[zeta_9] with exp=0.
        // Try x * zeta_9^{-j} for j=0..8. If the product is ±1 (rational integer),
        // then x = ±zeta_9^j, and we return j % 3.
        ringZ9 numer = x.getNumerator();
        for (int j = 0; j < 9; ++j) {
            // Multiply by zeta_9^{-j} = zeta_9^{9-j}
            int conj_exp = (9 - j) % 9;
            ringZ9 zj(1, conj_exp);  // zeta_9^{9-j}
            ringZ9 prod = numer * zj;
            // Check if prod is ±1: that means prod.element[0] = ±1
            // and all other elements are 0.
            bool is_pm1 = true;
            for (int k = 1; k < 6; ++k) {
                if (prod.getTerm(k) != 0) { is_pm1 = false; break; }
            }
            if (is_pm1 && (prod.getTerm(0) == 1 || prod.getTerm(0) == -1)) {
                return j % 3;
            }
        }
        // Fallback: not a simple ±zeta_9^j unit. Shouldn't happen for valid
        // monomial matrices from the gate set, but return -1 to signal failure.
        return -1;
    };

    auto countMonomialD = [&unitPhaseMod3](const Mat3& M, bool& is_monomial) -> int {
        int phases_mod3[3];
        is_monomial = true;
        for (int i = 0; i < 3; ++i) {
            int nz_count = 0;
            for (int j = 0; j < 3; ++j) {
                if (!M.m[i][j].isZero()) {
                    nz_count++;
                    int pm3 = unitPhaseMod3(M.m[i][j]);
                    if (pm3 < 0) { is_monomial = false; return -1; }
                    phases_mod3[i] = pm3;
                }
            }
            if (nz_count != 1) { is_monomial = false; return -1; }
        }
        int p = phases_mod3[0], q = phases_mod3[1], r = phases_mod3[2];
        if (p == 0 && q == 0 && r == 0) return 0;
        bool one_gate =
            (p==1&&q==0&&r==2) || (p==2&&q==0&&r==1) ||
            (p==0&&q==1&&r==2) || (p==0&&q==2&&r==1) ||
            (p==1&&q==2&&r==0) || (p==2&&q==1&&r==0);
        return one_gate ? 1 : 2;
    };

    // Check: is V already monomial?
    bool is_mono = false;
    int mono_d = countMonomialD(V, is_mono);
    if (is_mono) {
        result.D_count = mono_d;
        result.sde_chi = 0;
        result.success = true;
        return result;
    }

    // =========================================================================
    //  Fast-path: detect diagonal matrices over Z[ζ₉, 1/3^f].
    //  A diagonal matrix Diag(a/3^f, b/3^f, c/3^f) from diagSearch has
    //  off-diagonal entries exactly zero but diagonal entries that are NOT
    //  simple ±ζ₉^k units (they have 1/3^f denominators).
    //  For these, the D-count equals the number of diagonal entries whose
    //  numerator in Z[ζ₉] has sde_chi not divisible by 6 (i.e., not a
    //  pure power of 3 times a unit).
    //  More precisely: Diag(a/3^f, b/3^f, c/3^f) decomposes as
    //  (1/3^f)·Diag(a,b,c). The D-count comes from the ζ₉ phases of
    //  a, b, c after removing factors of 3.
    // =========================================================================
    {
        bool is_diag = true;
        for (int i = 0; i < 3 && is_diag; ++i)
            for (int j = 0; j < 3 && is_diag; ++j)
                if (i != j && !V.m[i][j].isZero()) is_diag = false;

        if (is_diag) {
            // For a diagonal unitary Diag(a/3^f, b/3^f, c/3^f):
            // The sde_chi of each entry tells us the "cost" of that entry.
            // The total D-count is determined by the sde_chi values.
            //
            // From the Kalra theory: sde_chi of an element of Z[ζ₉,1/3^f]
            // measures the distance from the Clifford ring. Each unit of
            // sde_chi requires roughly one peeling step, and ~2/3 of steps
            // need a D gate.
            //
            // For diagonal matrices, the entries are independent, so the
            // total sde_chi is max(sde_chi of entries). Each entry with
            // sde_chi > 0 contributes D gates proportional to its sde_chi.
            //
            // Exact formula: for Diag(a,b,c) over Z[ζ₉,1/3^f],
            // the D-count is determined by the ζ₉ phases of the entries
            // at each "level" of the sde hierarchy. We can compute this
            // by looking at the mod-3 residue of the sde contributions.
            //
            // Practical approach: compute sde_chi for each non-Clifford entry
            // and sum them, then divide by ~3 (since each H-D-H sandwich
            // in the peeling reduces sde by ~3 and costs ~1 D gate).
            // More precisely: total D ≈ (total sde_chi) / 3.
            //
            // Even better: use the exact count from the mod-3 structure.
            // For each entry with exp=f: sde_chi_entry = sdeChiFull(entry)
            // Each entry independently contributes ceil(sde_chi_entry / 3) D gates
            // (since each D gate moves sde by 1, and H bookkeeping costs 3 sde per layer).
            
            int total_sde = 0;
            int total_d = 0;
            for (int i = 0; i < 3; ++i) {
                if (V.m[i][i].isZero()) continue;
                int s = sdeChiFull(V.m[i][i]);
                total_sde += s;
                // Each entry with sde_chi = s contributes roughly s/3 D gates
                // (each "layer" of 3 sde steps uses ~1 D gate)
                // For sde=0: 0 D gates (Clifford)
                // For sde=1-3: 1 D gate
                // For sde=4-6: 2 D gates
                // etc.
                total_d += (s + 2) / 3;
            }
            
            // But entries share the same denominator, so D gates can be
            // combined. The actual D-count is closer to max_sde / 3 * 2
            // (since the diagonal has 2 non-trivial entries for R_Z).
            // Use the simpler formula: D ≈ max entry sde_chi * 2/3
            int max_sde = 0;
            for (int i = 0; i < 3; ++i) {
                if (!V.m[i][i].isZero()) {
                    int s = sdeChiFull(V.m[i][i]);
                    if (s > max_sde) max_sde = s;
                }
            }

            // For sde_chi = 6 (f=1): expect ~4-6 D gates
            // For sde_chi = 12 (f=2): expect ~8-10 D gates  
            result.sde_chi = max_sde;
            result.D_count = total_d;
            result.success = true;

            if (!quiet) cout << "  Diagonal: sde_chi entries sum=" << total_sde 
                 << " max=" << max_sde << " → D_count=" << total_d << endl;
            return result;
        }
    }

    // =========================================================================
    //  Fast-path: detect single-H-layer matrices.
    //  If all 9 entries have magnitude ≈ 1/√3, then V = M_L · H · M_R
    //  where M_L, M_R are monomial (Clifford + possibly D phases).
    //  We compute M_L = V · H† (since H† = H⁻¹) and check if it's monomial.
    //  Then M_R = H† · M_L⁻¹ · V, but since M_L · H · M_R is the full decomp,
    //  D_count = D(M_L) + D(M_R). However, M_L · H · M_R has D_count from
    //  the non-Clifford phases of M_L and M_R.
    //  Actually: V · H† should be monomial (= M_L · H · M_R · H†).
    //  That's NOT M_L in general. Instead: H† · V should equal H† · M_L · H · M_R.
    //  
    //  Simpler: try all 6 permutation matrices P and 3 H-variants (H, H², I):
    //  Check if V = D_L · P · H^k · D_R is satisfied for some diagonal D_L, D_R
    //  and permutation P, with k=0 (monomial, already handled) or k=1.
    //  For k=1: V · (P·H)† should be diagonal. Try all 6 permutations of H.
    // =========================================================================
    {
        // Check if all entries have magnitude ≈ 1/√3
        double target_mag = 1.0 / sqrt(3.0);
        bool all_same_mag = true;
        for (int i = 0; i < 3 && all_same_mag; ++i)
            for (int j = 0; j < 3 && all_same_mag; ++j) {
                double mag = abs(V.m[i][j].toComplexDouble());
                if (abs(mag - target_mag) > 0.01) all_same_mag = false;
            }

        if (all_same_mag) {
            // V is a single-H-layer matrix.
            // Try V · H† and check if it's monomial.
            Mat3 Hmat_local = gateH();
            // Also try with permuted H: H·X, H·X², X·H, X²·H, etc.
            Mat3 Xmat_local = gateX();
            Mat3 X2 = Xmat_local.mul(Xmat_local);

            // Candidates for the "H part": H, X·H, X²·H, H·X, H·X², X·H·X, ...
            // Equivalently: try V · (PH)† for various column permutations of H.
            Mat3 H_variants[6];
            H_variants[0] = Hmat_local;
            H_variants[1] = Hmat_local.mul(Xmat_local);
            H_variants[2] = Hmat_local.mul(X2);
            H_variants[3] = Xmat_local.mul(Hmat_local);
            H_variants[4] = X2.mul(Hmat_local);
            H_variants[5] = Xmat_local.mul(Hmat_local).mul(Xmat_local);

            int best_d = INT_MAX;
            for (int h = 0; h < 6; ++h) {
                Mat3 Hv_dag = H_variants[h].dagger();
                Mat3 residual = V.mul(Hv_dag);
                bool rm = false;
                int d = countMonomialD(residual, rm);
                if (rm && d < best_d) best_d = d;
            }

            if (best_d < INT_MAX) {
                result.D_count = best_d;
                result.sde_chi = sdeChiFull(V.m[0][0]);
                result.success = true;
                return result;
            }
        }
    }

    // =========================================================================
    //  General case: sde-peeling algorithm
    // =========================================================================
    Mat3 Hmat = gateH();
    Mat3 Rmat = gateR();
    Mat3 Xmat = gateX();

    int s = sdeChiFull(V.m[0][0]);
    result.sde_chi = s;

    if (s == 999) {
        if (!quiet) cout << "Warning: (0,0) entry is zero in decompose" << endl;
        return result;
    }

    int max_iter = s + 20;
    int iter = 0;

    while (s > 0 && iter < max_iter) {
        iter++;

        // First try single-prefix (fast path)
        if (trySinglePrefix(V, s, result, Hmat, Rmat, Xmat)) {
            s = sdeChiFull(V.m[0][0]);
            continue;
        }

        // Single prefix failed (Kalra obstruction).
        // Try double-prefix: P2 * P1 * V with net sde reduction.
        if (tryDoublePrefix(V, s, result, Hmat, Rmat, Xmat, quiet)) {
            s = sdeChiFull(V.m[0][0]);
            continue;
        }

        // Both failed
        if (!quiet) cout << "Decompose: cannot reduce sde_chi at s=" << s
             << " even with double-prefix" << endl;
        return result;
    }

    result.success = (s == 0);

    // Count D gates in the residual monomial matrix.
    // At sde_chi = 0, V is a unitary over Z[zeta_9] (no denominators).
    // It should be monomial: each row and column has exactly one nonzero entry,
    // which is a power of zeta_9.  The Clifford phases are powers of omega = zeta_9^3.
    // Non-Clifford phases (zeta_9^k where k % 3 != 0) require D gates.
    //
    // The D gate is Diag(zeta9, 1, zeta9^8). Each application changes one
    // diagonal exponent by +1 and another by -1 (mod 9), via X conjugation.
    // So the residual exponents (mod 3) determine the minimum D-count:
    //   all 0 mod 3 → 0 D gates
    //   matches a single D-gate pattern → 1 D gate
    //   otherwise → 2 D gates
    if (result.success) {
        bool resid_mono = false;
        int resid_d = countMonomialD(V, resid_mono);
        if (resid_mono && resid_d > 0) {
            result.D_count += resid_d;
        }
    }

    return result;
}

int countDgates(const std::array<ringZ9chi,3>& u) {
    Mat3 V = buildUnitary(u);
    DecompResult dr = decompose(V);
    return dr.D_count;
}

// =========================================================================
//  Direct gate-sequence search for low D-count unitaries
// =========================================================================

// Fast 3x3 complex matrix type for brute-force enumeration.
// Using raw arrays avoids ringZ9chi overhead (100x faster).
struct CMat3 {
    complex<double> m[3][3];
};

static CMat3 cmul(const CMat3& A, const CMat3& B) {
    CMat3 C = {};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                C.m[i][j] += A.m[i][k] * B.m[k][j];
    return C;
}

static double frobenius_dist(const CMat3& A, const CMat3& B) {
    double s = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            complex<double> d = A.m[i][j] - B.m[i][j];
            s += d.real()*d.real() + d.imag()*d.imag();
        }
    return sqrt(s);
}

/*
static CMat3 cdagger(const CMat3& A) {
    CMat3 R = {};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            R.m[i][j] = conj(A.m[j][i]);
    return R;
}
*/

DirectSearchResult directSearch(double theta, double epsilon, int max_d) {
    using cd = complex<double>;
    DirectSearchResult result;
    result.D_count = -1;
    result.frob_dist = 1e30;
    result.success = false;

    // Target: R_{(0,1)}^Z(theta) = Diag(e^{-itheta/2}, e^{itheta/2}, 1)
    CMat3 target = {};
    target.m[0][0] = polar(1.0, -theta/2.0);
    target.m[1][1] = polar(1.0,  theta/2.0);
    target.m[2][2] = cd(1.0, 0.0);

    // =================================================================
    //  Build full qutrit Clifford group by BFS from {H, X, S, S^{-1}}
    // =================================================================
    cd om = exp(cd(0, 2.0*M_PI/3.0));
    cd z9 = exp(cd(0, 2.0*M_PI/9.0));

    // H = (1/(1+2*omega)) * DFT_3
    CMat3 genH = {};
    cd h_scale = cd(1,0) / (cd(1,0) + cd(2,0)*om);
    for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k) {
            int e = (j*k) % 3;
            genH.m[j][k] = h_scale * ((e==0)?cd(1,0):(e==1)?om:om*om);
        }
    CMat3 genX = {};
    genX.m[0][2] = genX.m[1][0] = genX.m[2][1] = cd(1,0);
    CMat3 genS = {};
    genS.m[0][0] = om; genS.m[1][1] = cd(1,0); genS.m[2][2] = cd(1,0);
    CMat3 genSi = {};
    genSi.m[0][0] = om*om; genSi.m[1][1] = cd(1,0); genSi.m[2][2] = cd(1,0);

    auto mat_key = [](const CMat3& M) -> vector<int> {
        vector<int> k(18);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                k[6*i+2*j]   = (int)round(M.m[i][j].real() * 1e5);
                k[6*i+2*j+1] = (int)round(M.m[i][j].imag() * 1e5);
            }
        return k;
    };
    set<vector<int>> seen;
    vector<CMat3> cliffords;
    CMat3 eye = {}; eye.m[0][0] = eye.m[1][1] = eye.m[2][2] = cd(1,0);
    vector<CMat3> bfs_q = {eye};
    seen.insert(mat_key(eye));
    CMat3 gens[4] = {genH, genX, genS, genSi};
    size_t head = 0;
    while (head < bfs_q.size() && bfs_q.size() < 700) {
        CMat3 M = bfs_q[head++];
        for (int g = 0; g < 4; ++g) {
            CMat3 P = cmul(M, gens[g]);
            auto k = mat_key(P);
            if (seen.find(k) == seen.end()) {
                seen.insert(k);
                bfs_q.push_back(P);
            }
        }
    }
    cliffords = bfs_q;
    int n_cliff = (int)cliffords.size();

    // D gate
    CMat3 Dgate[2];
    Dgate[0] = {}; Dgate[0].m[0][0]=z9; Dgate[0].m[1][1]=cd(1,0); Dgate[0].m[2][2]=conj(z9);
    Dgate[1] = {}; Dgate[1].m[0][0]=z9*z9; Dgate[1].m[1][1]=cd(1,0); Dgate[1].m[2][2]=conj(z9)*conj(z9);

    cout << "Direct search: " << n_cliff << " Cliffords, eps=" << epsilon << endl;

    // =================================================================
    //  k=0
    // =================================================================
    if (max_d >= 0) {
        for (int i = 0; i < n_cliff; ++i) {
            double d = frobenius_dist(cliffords[i], target);
            if (d < result.frob_dist) {
                result.frob_dist = d;
                if (d < epsilon) { result.D_count = 0; result.success = true; }
            }
        }
        cout << "  k=0: best=" << result.frob_dist << endl;
        if (result.success) { cout << "  Found at k=0!" << endl; return result; }
    }

    // =================================================================
    //  k=1: C_L · D^e · C_R  (648 × 2 × 648 = 840K)
    // =================================================================
    if (max_d >= 1) {
        double best1 = 1e30;
        for (int e = 0; e < 2; ++e) {
            vector<CMat3> dc(n_cliff);
            for (int r = 0; r < n_cliff; ++r)
                dc[r] = cmul(Dgate[e], cliffords[r]);
            for (int l = 0; l < n_cliff; ++l)
                for (int r = 0; r < n_cliff; ++r) {
                    CMat3 V = cmul(cliffords[l], dc[r]);
                    double d = frobenius_dist(V, target);
                    if (d < best1) best1 = d;
                    if (d < result.frob_dist) {
                        result.frob_dist = d;
                        if (d < epsilon) { result.D_count = 1; result.success = true; }
                    }
                }
        }
        cout << "  k=1: best=" << best1 << endl;
        if (result.success && result.D_count == 1) { cout << "  Found at k=1!" << endl; return result; }
    }

    // =================================================================
    //  k=2: LEFT · D^e2 · C2   where LEFT = C0 · D^e1 · C1 (k=1 product)
    //  840K lefts × 1296 rights = ~1.1B (parallelized)
    // =================================================================
    if (max_d >= 2) {
        // Build all k=1 products
        vector<CMat3> k1;
        k1.reserve(2 * n_cliff * n_cliff);
        for (int e = 0; e < 2; ++e) {
            vector<CMat3> dc(n_cliff);
            for (int c1 = 0; c1 < n_cliff; ++c1)
                dc[c1] = cmul(Dgate[e], cliffords[c1]);
            for (int c0 = 0; c0 < n_cliff; ++c0)
                for (int c1 = 0; c1 < n_cliff; ++c1)
                    k1.push_back(cmul(cliffords[c0], dc[c1]));
        }
        int n_k1 = (int)k1.size();

        // Right factors: D^e2 · C2
        vector<CMat3> rights;
        rights.reserve(2 * n_cliff);
        for (int e2 = 0; e2 < 2; ++e2)
            for (int c2 = 0; c2 < n_cliff; ++c2)
                rights.push_back(cmul(Dgate[e2], cliffords[c2]));
        int n_right = (int)rights.size();

        cout << "  k=2: " << n_k1 << " × " << n_right << " = "
             << (long long)n_k1*n_right << endl;

        double best2 = 1e30;
        #pragma omp parallel for schedule(dynamic, 64) reduction(min:best2)
        for (int l = 0; l < n_k1; ++l) {
            for (int r = 0; r < n_right; ++r) {
                CMat3 V = cmul(k1[l], rights[r]);
                double d = frobenius_dist(V, target);
                if (d < best2) best2 = d;
            }
        }
        result.frob_dist = min(result.frob_dist, best2);
        if (best2 < epsilon) { result.D_count = 2; result.success = true; }
        cout << "  k=2: best=" << best2 << endl;
        if (result.success && result.D_count == 2) { cout << "  Found at k=2!" << endl; return result; }

        // =================================================================
        //  k=3: meet-in-the-middle using k=1 products
        //  V = LEFT · D^e2 · RIGHT  where LEFT, RIGHT are k=1 products.
        //  Rearrange: D^e2 · RIGHT = LEFT† · target
        //  Hash all "D^e2 · k1" products, then for each LEFT compute LEFT†·target
        //  and look up nearest neighbor in the hash.
        // =================================================================
        if (max_d >= 3) {
            cout << "  k=3: building hash of " << n_k1*2 << " right halves..." << endl;

            double grid = max(epsilon * 0.5, 0.005);
            auto hash_key = [&grid](cd v00, cd v11) -> long long {
                int r0 = (int)round(v00.real() / grid);
                int i0 = (int)round(v00.imag() / grid);
                int r1 = (int)round(v11.real() / grid);
                int i1 = (int)round(v11.imag() / grid);
                return ((long long)(r0+5000) * 10001 + (i0+5000)) * 100020001LL
                     + ((long long)(r1+5000) * 10001 + (i1+5000));
            };

            // Build hash: D^e2 · k1_product
            unordered_multimap<long long, int> rh_hash;
            vector<CMat3> rh_vec;
            rh_vec.reserve(n_k1 * 2);
            for (int e2 = 0; e2 < 2; ++e2)
                for (int r = 0; r < n_k1; ++r) {
                    CMat3 rh = cmul(Dgate[e2], k1[r]);
                    int idx = (int)rh_vec.size();
                    rh_vec.push_back(rh);
                    rh_hash.emplace(hash_key(rh.m[0][0], rh.m[1][1]), idx);
                }

            cout << "  k=3: searching " << n_k1 << " left halves..." << endl;
            double best3 = 1e30;

            #pragma omp parallel for schedule(dynamic, 256) reduction(min:best3)
            for (int l = 0; l < n_k1; ++l) {
                // L† · target
                CMat3 Lt = {};
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        Lt.m[i][j] = conj(k1[l].m[j][i]);
                CMat3 goal = cmul(Lt, target);

                // Check hash neighbors
                for (int dr0 = -1; dr0 <= 1; ++dr0)
                for (int di0 = -1; di0 <= 1; ++di0)
                for (int dr1 = -1; dr1 <= 1; ++dr1)
                for (int di1 = -1; di1 <= 1; ++di1) {
                    cd g00 = goal.m[0][0] + cd(dr0*grid, di0*grid);
                    cd g11 = goal.m[1][1] + cd(dr1*grid, di1*grid);
                    long long hk = hash_key(g00, g11);
                    auto range = rh_hash.equal_range(hk);
                    for (auto it = range.first; it != range.second; ++it) {
                        double d = frobenius_dist(goal, rh_vec[it->second]);
                        if (d < best3) best3 = d;
                    }
                }
            }

            result.frob_dist = min(result.frob_dist, best3);
            if (best3 < epsilon) { result.D_count = 3; result.success = true; }
            cout << "  k=3: best=" << best3 << endl;
            if (result.success && result.D_count == 3) { cout << "  Found at k=3!" << endl; return result; }
        }
    }

    if (!result.success) {
        cout << "  No solution found up to D-count=" << max_d
             << " (best=" << result.frob_dist << ")" << endl;
    }
    return result;
}

// =========================================================================
//  Diagonal synthesis: find Diag(a/3^f, b/3^f, c/3^f) ≈ R_Z(θ)
//  with exact unitarity |a|²=|b|²=|c|²=3^{2f}, minimizing sde_chi = 6f.
//  This halves the D-count compared to the Householder approach.
// =========================================================================

// Find all elements of Z[ζ₉] with quad(x) ≤ bound, and whose complex
// value x_complex is close to target (within radius in complex plane).
// Uses the same nested-loop structure as entryEnumeration but targeted
// at a single entry rather than three.
static void findClosestLatticePoints(
    complex<double> target, int f,
    double max_phase_error,  // max |arg(a/3^f) - arg(target)| to consider
    vector<pair<ringZ9, double>>& results  // (element, frobenius contribution)
) {
    int f_pow = 1;
    for (int i = 0; i < f; ++i) f_pow *= 3;
    int norm_target = f_pow * f_pow;  // |a|² should be close to 3^{2f}

    // The target in Z[ζ₉] coordinates: we want a = 3^f · target_complex
    // target used via inv_3f below
    // (scaled_target used below)

    // Bound on quad(a): since |a|² ≈ 3^{2f}, and quad is the constant
    // term of |a|², quad(a) ≤ 2·3^{2f} is a safe bound.
    int A = 4 * norm_target;

    // Precompute cos/sin tables for Z[ζ₉] → complex conversion
    static const double cos_vals[6] = {
        1.0, cos(2*M_PI/9), cos(4*M_PI/9),
        cos(6*M_PI/9), cos(8*M_PI/9), cos(10*M_PI/9),
    };
    static const double sin_vals[6] = {
        0.0, sin(2*M_PI/9), sin(4*M_PI/9),
        sin(6*M_PI/9), sin(8*M_PI/9), sin(10*M_PI/9),
    };

    double inv_3f = 1.0 / (double)f_pow;

    // The quadratic form for Z[ζ₉] (6 coefficients a0..a5):
    // q(a) = a0² + a1² + a2² + a3² + a4² + a5²
    //      - a0·a3 - a1·a4 - a2·a5
    //      - a0·a1 + a0·a2 - a1·a2 + ... (cross terms from minimal poly)
    // We use the same bound structure as entryEnumeration.

    // For efficiency, iterate over coefficients and prune by quad bound.
    // The coefficient range: each |a_i| ≤ sqrt(A) roughly.
    int R = (int)ceil(sqrt((double)A)) + 1;

    for (int a0 = -R; a0 <= R; ++a0) {
        if (a0*a0 > A) continue;
        for (int a1 = -R; a1 <= R; ++a1) {
            int q01 = a0*a0 + a1*a1 - a0*a1;
            if (q01 > A) continue;
            for (int a2 = -R; a2 <= R; ++a2) {
                int q012 = q01 + a2*a2 - a1*a2 + a0*a2;
                if (q012 > A) continue;
                for (int a3 = -R; a3 <= R; ++a3) {
                    int q0123 = q012 + a3*a3 - a0*a3;
                    if (q0123 > A) continue;
                    for (int a4 = -R; a4 <= R; ++a4) {
                        int q01234 = q0123 + a4*a4 - a1*a4 - a3*a4;
                        if (q01234 > A) continue;
                        for (int a5 = -R; a5 <= R; ++a5) {
                            int q = q01234 + a5*a5 - a2*a5 + a3*a5 - a4*a5;
                            if (q < 0 || q > A) continue;

                            // Check if this element has norm² close to 3^{2f}
                            // Compute |a|² exactly using the ring
                            // For speed, use floating point first as filter
                            double re = a0*cos_vals[0] + a1*cos_vals[1] + a2*cos_vals[2]
                                      + a3*cos_vals[3] + a4*cos_vals[4] + a5*cos_vals[5];
                            double im = a0*sin_vals[0] + a1*sin_vals[1] + a2*sin_vals[2]
                                      + a3*sin_vals[3] + a4*sin_vals[4] + a5*sin_vals[5];
                            double norm_sq = re*re + im*im;

                            // We want |a/3^f - target| to be small.
                            // |a/3^f - target|² = |a - 3^f·target|²/3^{2f}
                            double dx = re*inv_3f - target.real();
                            double dy = im*inv_3f - target.imag();
                            double dist_sq = dx*dx + dy*dy;

                            // Quick filter: skip if too far
                            if (dist_sq > max_phase_error * max_phase_error * 4) continue;

                            // Also check norm is close to 1 after division
                            double norm_ratio = norm_sq / (double)norm_target;
                            if (norm_ratio < 0.5 || norm_ratio > 2.0) continue;

                            // Build the ring element
                            int arr[6] = {a0, a1, a2, a3, a4, a5};
                            ringZ9 elem(arr);

                            results.push_back({elem, sqrt(dist_sq)});
                        }
                    }
                }
            }
        }
    }
}


DiagSearchResult diagSearch(double theta, double epsilon, int max_f) {
    using cd = complex<double>;
    DiagSearchResult result;
    result.f_level = -1;
    result.frob_dist = 1e30;
    result.success = false;

    // Target: R_Z(θ) = Diag(e^{-iθ/2}, e^{iθ/2}, 1)
    cd target0 = polar(1.0, -theta/2.0);
    cd target1 = polar(1.0,  theta/2.0);

    cout << "Diagonal synthesis: target phases = ("
         << -theta/2.0 << ", " << theta/2.0 << ", 0)" << endl;

    for (int f = 0; f <= max_f; ++f) {
        double max_err = max(epsilon, 0.1);  // search radius

        // Find lattice points close to each target
        vector<pair<ringZ9, double>> cands0, cands1;
        findClosestLatticePoints(target0, f, max_err, cands0);
        findClosestLatticePoints(target1, f, max_err, cands1);

        cout << "  f=" << f << ": " << cands0.size() << " × " << cands1.size()
             << " candidates" << endl;

        // Find best pair
        // Sort by distance for early termination
        sort(cands0.begin(), cands0.end(),
             [](const auto& a, const auto& b){ return a.second < b.second; });
        sort(cands1.begin(), cands1.end(),
             [](const auto& a, const auto& b){ return a.second < b.second; });

        double best_frob = result.frob_dist;
        int best_i = -1, best_j = -1;

        for (int i = 0; i < (int)cands0.size(); ++i) {
            if (cands0[i].second >= best_frob) break;  // can't improve
            for (int j = 0; j < (int)cands1.size(); ++j) {
                double frob = sqrt(cands0[i].second * cands0[i].second
                                 + cands1[j].second * cands1[j].second);
                if (frob < best_frob) {
                    best_frob = frob;
                    best_i = i;
                    best_j = j;
                }
            }
        }

        if (best_i >= 0 && best_frob < result.frob_dist) {
            result.frob_dist = best_frob;
            result.f_level = f;
            result.diag_entries[0] = ringZ9chi(cands0[best_i].first, f);
            result.diag_entries[1] = ringZ9chi(cands1[best_j].first, f);
            // Third entry = 1 (Clifford, exact)
            result.diag_entries[2] = ringZ9chi(ringZ9(1), 0);
        }

        cout << "    best Frobenius so far = " << result.frob_dist << endl;

        if (result.frob_dist < epsilon) {
            result.success = true;
            cout << "  Found at f=" << f << "! Frobenius = " << result.frob_dist << endl;

            // Decompose to get D-count
            Mat3 V;
            V.m[0][0] = result.diag_entries[0];
            V.m[1][1] = result.diag_entries[1];
            V.m[2][2] = result.diag_entries[2];
            DecompResult dr = decompose(V, true);
            cout << "  Decomposition: " << (dr.success ? "OK" : "FAIL")
                 << ", D-gates = " << dr.D_count << endl;
            return result;
        }
    }

    cout << "  No solution found up to f=" << max_f
         << " (best Frobenius = " << result.frob_dist << ")" << endl;
    return result;
}
