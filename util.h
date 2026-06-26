#ifndef UTILS_H
#define UTILS_H

#include "params.h"
#include "fips202.h"
#include "matrixmod.h"

void XOF(uint8_t **buf, size_t *length, const uint8_t *seed, size_t seed_len, int num);

Fq rnd_GF(keccak_state *shake);

void rnd_sys_mat(Fq *M, int M_r, int M_c, const uint8_t *seed, size_t seed_len);

void rnd_inv_matrix(Fq *M, int M_r, int M_c, uint8_t *seed, size_t seed_len);

int parse_hash(const uint8_t *digest, int digest_len, uint8_t *h, int len_h);

int solve(Fq *A, Fq *B_inv, Fq *G0prime, Fq Amm);

void pi(Fq *Gout, Fq *A, Fq *B, Fq *G);
void phi(Fq *Gout, Fq *A, Fq *B, Fq *C, Fq *G);

#endif
