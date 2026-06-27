#ifndef CANONICAL_H
#define CANONICAL_H

#include "params.h"

/*
 * Construct the three projective-path bases U, V, and W.
 *
 * Storage:
 *   M contains n consecutive row-major n x n slices.
 *   U, V, W are row-major n x n matrices.
 *   Their columns are u_i, v_i, and w_i, respectively.
 *
 * Return:
 *    0  success
 *   -1  invalid input or path construction failure
 *
 * Preconditions:
 *   1 <= n && n <= MEDS_n
 *   M contains n^3 canonical field elements
 *   u1 contains n canonical field elements
 *   U, V, W each provide n^2 writable field elements
 *   U, V, W, M, and u1 do not overlap
 *
 * The function is variable-time.
 *
 * On failure, U, V, and W are left unchanged.
 * The representative u1 is preserved as the first column of U;
 * it is not projectively normalized here.
 */
int canonical_build_uvw_vartime(
    Fq *U,
    Fq *V,
    Fq *W,
    const Fq *M,
    const Fq *u1,
    int n);

#endif
