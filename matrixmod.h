#ifndef MATRIXMOD_H
#define MATRIXMOD_H

#include <stdio.h>
#include <stdint.h>

#include "field.h"

#define pmod_mat_entry(M, M_r, M_c, r, c) M[(M_c)*(r)+(c)]

#define pmod_mat_set_entry(M, M_r, M_c, r, c, v) (M[(M_c)*(r)+(c)] = v)

#define pmod_mat_t GFq_t
#define Fq GFq_t

void pmod_mat_print(const pmod_mat_t *M, int M_r, int M_c);
void pmod_mat_fprint(FILE *stream, const pmod_mat_t *M, int M_r, int M_c);

void pmod_mat_zero(pmod_mat_t *A, int n);
void pmod_mat_copy(pmod_mat_t *dst, const pmod_mat_t *src, int n);
void pmod_mat_identity(pmod_mat_t *A, int n);

void pmod_mat_transpose(
    pmod_mat_t *AT,
    const pmod_mat_t *A,
    int n);

void pmod_mat_mul(
    pmod_mat_t *C,
    const pmod_mat_t *A,
    const pmod_mat_t *B,
    int n);

void pmod_mat_mul_rect(
    Fq *C,
    int C_r,
    int C_c,
    const Fq *A,
    int A_r,
    int A_c,
    const Fq *B,
    int B_r,
    int B_c);

void pmod_mat_vec_mul(
    GFq_t *out,
    const pmod_mat_t *A,
    const GFq_t *v,
    int n);

void pmod_mat_transpose_vec_mul(
    GFq_t *out,
    const pmod_mat_t *A,
    const GFq_t *v,
    int n);

void pmod_mat_set_col(
    pmod_mat_t *A,
    int col,
    const GFq_t *v,
    int n);

void pmod_mat_get_col(
    GFq_t *v,
    const pmod_mat_t *A,
    int col,
    int n);

void pmod_mat_linear_combination(
    pmod_mat_t *out,
    const pmod_mat_t *matrices,
    const GFq_t *coeffs,
    int count,
    int n);

void pmod_mat_diag_scale(
    pmod_mat_t *out,
    const GFq_t *left_diag,
    const pmod_mat_t *A,
    const GFq_t *right_diag,
    int n);

/*
 * Legacy MEDS systematic-form helpers.
 * Do not use for rank or kernel computation.
 */
int pmod_mat_syst_ct(Fq *M, int M_r, int M_c);
int pmod_mat_row_echelon_ct(Fq *M, int M_r, int M_c);
int pmod_mat_back_substitution_ct(Fq *M, int M_r, int M_c);

int pmod_mat_inv(Fq *B, const Fq *A, int A_r, int A_c);

#endif
