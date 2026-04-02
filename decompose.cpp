#include "decompose.h"
#include <iostream>
#include <cassert>
#include <cmath>

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
// Returns true if found, and applies the prefix to V.
static bool trySinglePrefix(Mat3& V, int s, DecompResult& result,
                             const Mat3& Hmat, const Mat3& Rmat, const Mat3& Xmat) {
    for (int eps = 0; eps < 2; ++eps)
    for (int delta = 0; delta < 3; ++delta)
    for (int a1 = 0; a1 < 9; ++a1)
    for (int a2 = 0; a2 < 9; ++a2)
    for (int a3 = 0; a3 < 9; ++a3) {
        Mat3 prefix = buildPrefix(a1, a2, a3, eps, delta, Hmat, Rmat, Xmat);
        ringZ9chi new00 = prefix_times_V_00(prefix, V);
        if (sdeChiFull(new00) == s - 1) {
            GateStep step;
            step.a0 = a1; step.a1 = a2; step.a2 = a3;
            step.eps = eps; step.delta = delta; step.has_H = true;
            result.steps.push_back(step);
            result.D_count += count_d_from_cyclo(a1, a2, a3);
            V = prefix.mul(V);
            return true;
        }
    }
    return false;
}

// Try double-prefix: find two prefixes P2*P1 such that sde drops by at least 1.
// This is the fallback when the Kalra obstruction blocks single-prefix.
// We try all P1 that keep sde the same or raise it by at most 1,
// then for each such P1*V, try all P2 that bring sde below the original.
static bool tryDoublePrefix(Mat3& V, int s, DecompResult& result,
                             const Mat3& Hmat, const Mat3& Rmat, const Mat3& Xmat) {
    cout << "  Trying double-prefix at sde_chi=" << s << "..." << endl;

    // First pass: collect all P1 that give sde in {s-1, s, s+1}
    // (s-1 would have been caught by single, so really {s, s+1})
    // For efficiency, we only try a subset: fix eps1/delta1 and vary a's
    for (int eps1 = 0; eps1 < 2; ++eps1)
    for (int delta1 = 0; delta1 < 3; ++delta1)
    for (int a1 = 0; a1 < 9; ++a1)
    for (int a2 = 0; a2 < 9; ++a2)
    for (int a3 = 0; a3 < 9; ++a3) {
        Mat3 P1 = buildPrefix(a1, a2, a3, eps1, delta1, Hmat, Rmat, Xmat);
        ringZ9chi mid00 = prefix_times_V_00(P1, V);
        int mid_s = sdeChiFull(mid00);

        // Only consider P1 that keep sde manageable (same or +1)
        if (mid_s > s + 1 || mid_s < s - 1) continue;

        // Now try all P2 on top of P1*V
        Mat3 midV = P1.mul(V);
        for (int eps2 = 0; eps2 < 2; ++eps2)
        for (int delta2 = 0; delta2 < 3; ++delta2)
        for (int b1 = 0; b1 < 9; ++b1)
        for (int b2 = 0; b2 < 9; ++b2)
        for (int b3 = 0; b3 < 9; ++b3) {
            Mat3 P2 = buildPrefix(b1, b2, b3, eps2, delta2, Hmat, Rmat, Xmat);
            ringZ9chi new00 = prefix_times_V_00(P2, midV);
            int new_s = sdeChiFull(new00);

            // Accept any net reduction from original s
            if (new_s < s) {
                // Record step 1
                GateStep step1;
                step1.a0 = a1; step1.a1 = a2; step1.a2 = a3;
                step1.eps = eps1; step1.delta = delta1; step1.has_H = true;
                result.steps.push_back(step1);
                result.D_count += count_d_from_cyclo(a1, a2, a3);

                // Record step 2
                GateStep step2;
                step2.a0 = b1; step2.a1 = b2; step2.a2 = b3;
                step2.eps = eps2; step2.delta = delta2; step2.has_H = true;
                result.steps.push_back(step2);
                result.D_count += count_d_from_cyclo(b1, b2, b3);

                V = P2.mul(midV);
                cout << "  Double-prefix succeeded: sde " << s
                     << " -> " << mid_s << " -> " << new_s << endl;
                return true;
            }
        }
    }
    return false;
}

DecompResult decompose(Mat3 V) {
    DecompResult result;
    result.D_count = 0;
    result.success = false;

    Mat3 Hmat = gateH();
    Mat3 Rmat = gateR();
    Mat3 Xmat = gateX();

    int s = sdeChiFull(V.m[0][0]);
    result.sde_chi = s;

    if (s == 999) {
        cout << "Warning: (0,0) entry is zero in decompose" << endl;
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
        if (tryDoublePrefix(V, s, result, Hmat, Rmat, Xmat)) {
            s = sdeChiFull(V.m[0][0]);
            continue;
        }

        // Both failed
        cout << "Decompose: cannot reduce sde_chi at s=" << s
             << " even with double-prefix" << endl;
        return result;
    }

    result.success = (s == 0);
    return result;
}

int countDgates(const std::array<ringZ9chi,3>& u) {
    Mat3 V = buildUnitary(u);
    DecompResult dr = decompose(V);
    return dr.D_count;
}
