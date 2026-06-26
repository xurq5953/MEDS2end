#ifndef MATRIXELIM_H
#define MATRIXELIM_H

#include "matrixmod.h"

int pmod_mat_rank_vartime(
    const pmod_mat_t *A,
    int n);

int pmod_mat_is_invertible_vartime(
    const pmod_mat_t *A,
    int n);

int pmod_mat_inv_vartime(
    pmod_mat_t *A_inv,
    const pmod_mat_t *A,
    int n);

int pmod_mat_right_kernel_corank1_vartime(
    GFq_t *kernel,
    const pmod_mat_t *A,
    int n);

int pmod_mat_left_kernel_corank1_vartime(
    GFq_t *kernel,
    const pmod_mat_t *A,
    int n);

int pmod_vec_in_span_vartime(
    const GFq_t *v,
    const GFq_t *basis,
    int basis_count,
    int n);

int pmod_vec_extends_independent_set_vartime(
    const GFq_t *v,
    const GFq_t *basis,
    int basis_count,
    int n);

#endif
