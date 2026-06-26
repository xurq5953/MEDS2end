#ifndef MATRIXMOD_H
#define MATRIXMOD_H

#include <stdio.h>
#include <stdint.h>

#include "params.h"

#define pmod_mat_entry(M, M_r, M_c, r, c) M[(M_c)*(r)+(c)]

#define pmod_mat_set_entry(M, M_r, M_c, r, c, v) (M[(M_c)*(r)+(c)] = v)

#define Fq GFq_t

void pmod_mat_print(Fq *M, int M_r, int M_c);
void pmod_mat_fprint(FILE *stream, Fq *M, int M_r, int M_c);

void pmod_mat_mul(Fq *C, int C_r, int C_c, Fq *A, int A_r, int A_c, Fq *B, int B_r, int B_c);

int pmod_mat_syst_ct(Fq *M, int M_r, int M_c);
int pmod_mat_row_echelon_ct(Fq *M, int M_r, int M_c);
int pmod_mat_back_substitution_ct(Fq *M, int M_r, int M_c);

int pmod_mat_inv(Fq *B, Fq *A, int A_r, int A_c);

GFq_t GF_inv(GFq_t val);

#endif

