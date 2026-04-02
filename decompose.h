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
DecompResult decompose(Mat3 V);

/**
 * @brief Build the unitary V = X_{(0,1)} (I - u u^dagger) from a Householder vector u.
 */
Mat3 buildUnitary(const std::array<ringZ9chi,3>& u);

/**
 * @brief Convenience: run HRSA, build the unitary, decompose it, return D-count.
 */
int countDgates(const std::array<ringZ9chi,3>& u);

#endif
