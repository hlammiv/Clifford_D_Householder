#ifndef DECOMPOSE_H
#define DECOMPOSE_H

#include <array>
#include <vector>
#include <string>
#include "cyclotomic_int9.h"
#include "Z9chi.h"

/**
 * @file decompose.h
 * @brief Exact synthesis of a unitary in U(3, Z[zeta_9, 1/3]) into a product
 *        of Clifford+D generators, following the sde-reduction approach of
 *        Kalra et al. and the C+R decomposition of Gustafson et al.
 *
 * The algorithm iteratively peels off generators  H · D(a0,a1,a2) · Dgate^eps · X^delta
 * that reduce sde_chi of the (0,0) entry by 1, until sde_chi = 0.
 * The remaining Clifford element is identified by lookup.
 *
 * @author Claude (with Chris H.)
 * @date April 2026
 */

/// A 3x3 matrix over ringZ9chi.
struct Mat3 {
    ringZ9chi m[3][3];

    Mat3();  ///< Zero matrix
    Mat3 mul(const Mat3& B) const;      ///< Matrix multiplication
    Mat3 dagger() const;                ///< Conjugate-transpose
    void print() const;
};

/// One step in the gate sequence: H D(a0,a1,a2) Dgate^eps X^delta
struct GateStep {
    int a0, a1, a2;  ///< parameters for D(a0,a1,a2) = Diag(omega^a0, omega^a1, omega^a2)
    int eps;          ///< exponent on the D gate (0 or 1 or 2)
    int delta;        ///< exponent on X gate (0, 1, or 2)
    bool has_H;       ///< whether the Hadamard is included
};

/// Result of the decomposition
struct DecompResult {
    std::vector<GateStep> steps;
    int D_count;             ///< Total number of D gates used
    int sde_chi;             ///< sde_chi of the input unitary
    bool success;
};

// ---- Fixed gate matrices ----

/// Qutrit Hadamard (DFT matrix) over Z[zeta_9, 1/3]
Mat3 gateH();

/// Diagonal Clifford gate D(a,b,c) = Diag(omega^a, omega^b, omega^c)
/// where omega = zeta_9^3 = e^(2*pi*i/3)
Mat3 gateD_clifford(int a, int b, int c);

/// The non-Clifford D gate: Diag(zeta_9, 1, zeta_9^{-1}) = Diag(zeta_9, 1, zeta_9^8)
Mat3 gateDgate();

/// The X (cyclic shift) gate
Mat3 gateX();

// ---- sde_chi for ringZ9chi ----

/// Compute sde_chi for an element of Z[zeta_9, 1/3].
/// Extends sdeChi() (which handles sde < 6) to arbitrary sde
/// by using: sde_chi(a/3^f) = sde_chi(a) + 6*f  when a is in reduced form.
/// Since the existing sdeChi() returns min(sde_chi(a), 6) for a in Z[zeta_9],
/// we handle the case sde >= 6 by dividing by 3 and recursing.
int sdeChiFull(const ringZ9chi& x);

/// Compute sde_chi of numerator in Z[zeta_9], handling sde >= 6
/// by dividing out factors of 3 (each = chi^6 * unit).
int sdeChiZ9(ringZ9 a);

// ---- Decomposition ----

/**
 * @brief Decompose a 3x3 unitary over Z[zeta_9, 1/3] into Clifford+D generators.
 *
 * Uses the iterative sde-reduction algorithm: at each step, try all possible
 * prefixes  H · D(a0,a1,a2) · Dgate^eps · X^delta  and find one that reduces
 * sde_chi of the (0,0) entry by 1.  Each step with eps>0 contributes to the
 * D-gate count.
 *
 * @param V  The unitary to decompose (constructed from the Householder vector).
 * @return DecompResult containing the gate sequence and D-count.
 */
DecompResult decompose(Mat3 V, bool quiet = false);

/**
 * @brief Build the unitary V = X_{(0,1)} (I - u u^dagger) from a Householder vector u.
 */
Mat3 buildUnitary(const std::array<ringZ9chi,3>& u);

/**
 * @brief Convenience: run HRSA, build the unitary, decompose it, return D-count.
 */
int countDgates(const std::array<ringZ9chi,3>& u);

/**
 * @brief Direct gate-sequence search for low D-count unitaries.
 *
 * Enumerates all products  C₀·D^{e₁}·C₁·D^{e₂}·...·D^{eₖ}·Cₖ
 * where each Cᵢ is a monomial Clifford (permutation × ω-diagonal)
 * and eᵢ ∈ {1,2}, searching for one within epsilon of R_Z(θ) in
 * Frobenius distance.
 *
 * Searches k=0,1,...,max_d in order, returning the first hit (minimum D-count).
 * k=0: 162 candidates (instant)
 * k=1: ~52K candidates (instant)
 * k=2: ~17M candidates (~2 seconds)
 * k=3: meet-in-the-middle (~4 seconds)
 *
 * @param theta     Target rotation angle for R_{(0,1)}^Z(θ).
 * @param epsilon   Desired Frobenius distance tolerance.
 * @param max_d     Maximum D-count to search (recommend 2 or 3).
 * @param best_frob [out] The Frobenius distance of the best match found.
 * @return The best unitary found, or zero matrix if none within epsilon.
 */
struct DirectSearchResult {
    Mat3 V;
    int D_count;
    double frob_dist;
    bool success;
};

DirectSearchResult directSearch(double theta, double epsilon, int max_d);

/**
 * @brief Diagonal synthesis: find Diag(a/3^f, b/3^f, 1) ≈ R_Z(θ) directly.
 *
 * Instead of the Householder construction (which builds a full-rank projector
 * at sde_chi = 12f), this searches for diagonal unitaries at sde_chi = 6f.
 * This halves the D-gate cost: ~6 D-gates at ε≈0.003 vs ~20+ for Householder.
 *
 * For each f = 0,1,2,...,max_f:
 *   Enumerates Z[ζ₉] elements a,b with a/3^f ≈ e^{-iθ/2}, b/3^f ≈ e^{iθ/2}
 *   and checks if the resulting diagonal is within epsilon.
 */
struct DiagSearchResult {
    ringZ9chi diag_entries[3];
    int f_level;
    double frob_dist;
    bool success;
};

DiagSearchResult diagSearch(double theta, double epsilon, int max_f);

#endif
