#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "log.h"

#include "fips202.h"

#include "params.h"

#include "api.h"
#include "randombytes.h"

#include "meds.h"

#include "util.h"
#include "bitstream.h"

#include "matrixmod.h"

static void encode_mat_to_bs(bitstream_t *bs, const Fq *M, int elem_count)
{
  for (int j = 0; j < elem_count; j++)
    bs_write(bs, M[j], GFq_bits);

  bs_finalize(bs);
}

static void decode_mat_from_bs(bitstream_t *bs, Fq *M, int elem_count)
{
  for (int j = 0; j < elem_count; j++)
    M[j] = bs_read(bs, GFq_bits);

  bs_finalize(bs);
}

static void derive_round_matrices(
    Fq *A,
    Fq *B,
    Fq *C,
    const uint8_t round_seed[MEDS_round_seed_bytes],
    const uint8_t alpha[MEDS_salt_bytes],
    uint32_t round_index
  )
{
  uint8_t round_input[MEDS_salt_bytes + MEDS_round_seed_bytes + 4];
  uint8_t sigma_A[MEDS_round_seed_bytes];
  uint8_t sigma_B[MEDS_round_seed_bytes];
  uint8_t sigma_C[MEDS_round_seed_bytes];

  memcpy(round_input, alpha, MEDS_salt_bytes);
  memcpy(round_input + MEDS_salt_bytes, round_seed, MEDS_round_seed_bytes);
  round_input[MEDS_salt_bytes + MEDS_round_seed_bytes + 0] = (uint8_t)(round_index & 0xff);
  round_input[MEDS_salt_bytes + MEDS_round_seed_bytes + 1] = (uint8_t)((round_index >> 8) & 0xff);
  round_input[MEDS_salt_bytes + MEDS_round_seed_bytes + 2] = (uint8_t)((round_index >> 16) & 0xff);
  round_input[MEDS_salt_bytes + MEDS_round_seed_bytes + 3] = (uint8_t)((round_index >> 24) & 0xff);

  XOF((uint8_t*[]){sigma_A, sigma_B, sigma_C},
      (size_t[]){MEDS_round_seed_bytes, MEDS_round_seed_bytes, MEDS_round_seed_bytes},
      round_input, sizeof(round_input),
      3);

  rnd_inv_matrix(A, MEDS_m, MEDS_m, sigma_A, MEDS_round_seed_bytes);
  rnd_inv_matrix(B, MEDS_n, MEDS_n, sigma_B, MEDS_round_seed_bytes);
  rnd_inv_matrix(C, MEDS_k, MEDS_k, sigma_C, MEDS_round_seed_bytes);
}

static void encode_G_full(uint8_t out[MEDS_G_BYTES], const Fq *G)
{
  memset(out, 0, MEDS_G_BYTES);

  bitstream_t bs;

  bs_init(&bs, out, MEDS_G_BYTES);

  for (int j = 0; j < MEDS_k * MEDS_m * MEDS_n; j++)
    bs_write(&bs, G[j], GFq_bits);

  bs_finalize(&bs);
}

static void load_public_key_matrices(Fq *G[MEDS_s], const unsigned char *pk)
{
  rnd_sys_mat(G[0], MEDS_k, MEDS_m*MEDS_n, pk, MEDS_pub_seed_bytes);

  bitstream_t bs;

  bs_init(&bs, (uint8_t*)pk + MEDS_pub_seed_bytes, MEDS_PK_BYTES - MEDS_pub_seed_bytes);

  for (int si = 1; si < MEDS_s; si++)
  {
    for (int j = 0; j < MEDS_k * MEDS_m * MEDS_n; j++)
      G[si][j] = bs_read(&bs, GFq_bits);

    bs_finalize(&bs);
  }
}

static void load_secret_key_matrices(
    Fq *A_inv[MEDS_s],
    Fq *B_inv[MEDS_s],
    Fq *C_inv[MEDS_s],
    const unsigned char *sk_body
  )
{
  bitstream_t bs;

  bs_init(&bs, (uint8_t*)sk_body, MEDS_SK_BYTES - MEDS_sec_seed_bytes - MEDS_pub_seed_bytes);

  for (int si = 1; si < MEDS_s; si++)
  {
    for (int j = 0; j < MEDS_m*MEDS_m; j++)
      A_inv[si][j] = bs_read(&bs, GFq_bits);

    bs_finalize(&bs);
  }

  for (int si = 1; si < MEDS_s; si++)
  {
    for (int j = 0; j < MEDS_n*MEDS_n; j++)
      B_inv[si][j] = bs_read(&bs, GFq_bits);

    bs_finalize(&bs);
  }

  for (int si = 1; si < MEDS_s; si++)
  {
    for (int j = 0; j < MEDS_k*MEDS_k; j++)
      C_inv[si][j] = bs_read(&bs, GFq_bits);

    bs_finalize(&bs);
  }
}

int crypto_sign_keypair(
    unsigned char *pk,
    unsigned char *sk
  )
{
  uint8_t delta[MEDS_sec_seed_bytes];

  randombytes(delta, MEDS_sec_seed_bytes);


  Fq G_data[MEDS_k * MEDS_m * MEDS_n * MEDS_s];
  Fq *G[MEDS_s];

  for (int i = 0; i < MEDS_s; i++)
    G[i] = &G_data[i * MEDS_k * MEDS_m * MEDS_n];

  uint8_t sigma_G0[MEDS_pub_seed_bytes];
  uint8_t sigma[MEDS_sec_seed_bytes];

  XOF((uint8_t*[]){sigma_G0, sigma},
       (size_t[]){MEDS_pub_seed_bytes, MEDS_sec_seed_bytes},
       delta, MEDS_sec_seed_bytes,
       2);

  LOG_VEC(sigma, MEDS_sec_seed_bytes);
  LOG_VEC_FMT(sigma_G0, MEDS_pub_seed_bytes, "sigma_G0");


  rnd_sys_mat(G[0], MEDS_k, MEDS_m*MEDS_n, sigma_G0, MEDS_pub_seed_bytes);

  LOG_MAT(G[0], MEDS_k, MEDS_m*MEDS_n);

  Fq A_inv_data[MEDS_s * MEDS_m * MEDS_m];
  Fq B_inv_data[MEDS_s * MEDS_n * MEDS_n];
  Fq C_inv_data[MEDS_s * MEDS_k * MEDS_k];

  Fq *A_inv[MEDS_s];
  Fq *B_inv[MEDS_s];
  Fq *C_inv[MEDS_s];

  for (int i = 0; i < MEDS_s; i++)
  {
    A_inv[i] = &A_inv_data[i * MEDS_m * MEDS_m];
    B_inv[i] = &B_inv_data[i * MEDS_n * MEDS_n];
    C_inv[i] = &C_inv_data[i * MEDS_k * MEDS_k];
  }

  for (int i = 1; i < MEDS_s; i++)
  {
    Fq A[MEDS_m * MEDS_m] = {0};
    Fq B[MEDS_n * MEDS_n] = {0};
    Fq C[MEDS_k * MEDS_k] = {0};

    while (1 == 1) // redo generation for this index until success
    {
      uint8_t sigma_Ai[MEDS_sec_seed_bytes];
      uint8_t sigma_Bi[MEDS_sec_seed_bytes];
      uint8_t sigma_Ci[MEDS_sec_seed_bytes];

      XOF((uint8_t*[]){sigma_Ai, sigma_Bi, sigma_Ci, sigma},
          (size_t[]){MEDS_sec_seed_bytes, MEDS_sec_seed_bytes, MEDS_sec_seed_bytes, MEDS_sec_seed_bytes},
          sigma, MEDS_sec_seed_bytes,
          4);

      rnd_inv_matrix(A, MEDS_m, MEDS_m, sigma_Ai, MEDS_sec_seed_bytes);
      rnd_inv_matrix(B, MEDS_n, MEDS_n, sigma_Bi, MEDS_sec_seed_bytes);
      rnd_inv_matrix(C, MEDS_k, MEDS_k, sigma_Ci, MEDS_sec_seed_bytes);

      if (pmod_mat_inv(A_inv[i], A, MEDS_m, MEDS_m) < 0)
      {
        LOG("no inv A_inv");
        continue;
      }

      if (pmod_mat_inv(B_inv[i], B, MEDS_n, MEDS_n) < 0)
      {
        LOG("no inv B_inv");
        continue;
      }

      if (pmod_mat_inv(C_inv[i], C, MEDS_k, MEDS_k) < 0)
      {
        LOG("no inv C_inv");
        continue;
      }

      LOG_MAT_FMT(A, MEDS_m, MEDS_m, "A[%i]", i);
      LOG_MAT_FMT(A_inv[i], MEDS_m, MEDS_m, "A_inv[%i]", i);
      LOG_MAT_FMT(B, MEDS_n, MEDS_n, "B[%i]", i);
      LOG_MAT_FMT(B_inv[i], MEDS_n, MEDS_n, "B_inv[%i]", i);
      LOG_MAT_FMT(C, MEDS_k, MEDS_k, "C[%i]", i);
      LOG_MAT_FMT(C_inv[i], MEDS_k, MEDS_k, "C_inv[%i]", i);


      phi(G[i], A, B, C, G[0]);

      LOG_MAT_FMT(G[i], MEDS_k, MEDS_m*MEDS_n, "G[%i]", i);

      // successfull generated G[s]; break out of while loop
      break;
    }
  }


  // copy pk data
  {
    uint8_t *tmp_pk = pk;

    memcpy(tmp_pk, sigma_G0, MEDS_pub_seed_bytes);
    LOG_VEC(tmp_pk, MEDS_pub_seed_bytes, "sigma_G0 (pk)");
    tmp_pk += MEDS_pub_seed_bytes;

    bitstream_t bs;

    bs_init(&bs, tmp_pk, MEDS_PK_BYTES - MEDS_pub_seed_bytes);

    for (int si = 1; si < MEDS_s; si++)
    {
      for (int j = 0; j < MEDS_k * MEDS_m * MEDS_n; j++)
        bs_write(&bs, G[si][j], GFq_bits);

      bs_finalize(&bs);
    }

    LOG_VEC(tmp_pk, MEDS_PK_BYTES - MEDS_pub_seed_bytes, "G[1:] (pk)");
    tmp_pk += MEDS_PK_BYTES - MEDS_pub_seed_bytes;

    LOG_HEX(pk, MEDS_PK_BYTES);

    if (MEDS_PK_BYTES != MEDS_pub_seed_bytes + bs.byte_pos + (bs.bit_pos > 0 ? 1 : 0))
    {
      fprintf(stderr, "ERROR: MEDS_PK_BYTES and actual pk size do not match! %i vs %i\n", MEDS_PK_BYTES, MEDS_pub_seed_bytes + bs.byte_pos+(bs.bit_pos > 0 ? 1 : 0));
      fprintf(stderr, "%i %i\n", MEDS_pub_seed_bytes + bs.byte_pos, MEDS_pub_seed_bytes + bs.byte_pos + (bs.bit_pos > 0 ? 1 : 0));
      return -1;
    }
  }

  // copy sk data
  {
    memcpy(sk, delta, MEDS_sec_seed_bytes);
    memcpy(sk + MEDS_sec_seed_bytes, sigma_G0, MEDS_pub_seed_bytes);

    bitstream_t bs;

    bs_init(&bs, sk + MEDS_sec_seed_bytes + MEDS_pub_seed_bytes, MEDS_SK_BYTES - MEDS_sec_seed_bytes - MEDS_pub_seed_bytes);

    for (int si = 1; si < MEDS_s; si++)
    {
      for (int j = 0; j < MEDS_m*MEDS_m; j++)
        bs_write(&bs, A_inv[si][j], GFq_bits);

      bs_finalize(&bs);
    }

    for (int si = 1; si < MEDS_s; si++)
    {
      for (int j = 0; j < MEDS_n*MEDS_n; j++)
        bs_write(&bs, B_inv[si][j], GFq_bits);

      bs_finalize(&bs);
    }

    for (int si = 1; si < MEDS_s; si++)
    {
      for (int j = 0; j < MEDS_k*MEDS_k; j++)
        bs_write(&bs, C_inv[si][j], GFq_bits);

      bs_finalize(&bs);
    }

    if (MEDS_SK_BYTES != MEDS_sec_seed_bytes + MEDS_pub_seed_bytes + bs.byte_pos + (bs.bit_pos > 0 ? 1 : 0))
    {
      fprintf(stderr, "ERROR: MEDS_SK_BYTES and actual sk size do not match! %i vs %i\n", MEDS_SK_BYTES, MEDS_sec_seed_bytes + MEDS_pub_seed_bytes + bs.byte_pos+(bs.bit_pos > 0 ? 1 : 0));
      fprintf(stderr, "%i %i\n", MEDS_sec_seed_bytes + MEDS_pub_seed_bytes + bs.byte_pos, MEDS_sec_seed_bytes + MEDS_pub_seed_bytes + bs.byte_pos + (bs.bit_pos > 0 ? 1 : 0));
      return -1;
    }

    LOG_HEX(sk, MEDS_SK_BYTES);
  }

  return 0;
}

int crypto_sign(
    unsigned char *sm, unsigned long long *smlen,
    const unsigned char *m, unsigned long long mlen,
    const unsigned char *sk
  )
{
  uint8_t delta[MEDS_sec_seed_bytes];

  randombytes(delta, MEDS_sec_seed_bytes);


  // skip secret seed
  sk += MEDS_sec_seed_bytes;

  Fq G_0[MEDS_k * MEDS_m * MEDS_n];


  rnd_sys_mat(G_0, MEDS_k, MEDS_m*MEDS_n, sk, MEDS_pub_seed_bytes);

  sk += MEDS_pub_seed_bytes;


  Fq A_inv_data[MEDS_s * MEDS_m * MEDS_m];
  Fq B_inv_data[MEDS_s * MEDS_n * MEDS_n];
  Fq C_inv_data[MEDS_s * MEDS_k * MEDS_k];

  Fq *A_inv[MEDS_s];
  Fq *B_inv[MEDS_s];
  Fq *C_inv[MEDS_s];

  for (int i = 0; i < MEDS_s; i++)
  {
    A_inv[i] = &A_inv_data[i * MEDS_m * MEDS_m];
    B_inv[i] = &B_inv_data[i * MEDS_n * MEDS_n];
    C_inv[i] = &C_inv_data[i * MEDS_k * MEDS_k];
  }

  load_secret_key_matrices(A_inv, B_inv, C_inv, sk);


  for (int i = 1; i < MEDS_s; i++)
    LOG_MAT(A_inv[i], MEDS_m, MEDS_m);

  for (int i = 1; i < MEDS_s; i++)
    LOG_MAT(B_inv[i], MEDS_n, MEDS_n);

  for (int i = 1; i < MEDS_s; i++)
    LOG_MAT(C_inv[i], MEDS_k, MEDS_k);

  LOG_MAT(G_0, MEDS_k, MEDS_m*MEDS_m);

  LOG_VEC(delta, MEDS_sec_seed_bytes);


  uint8_t alpha[MEDS_salt_bytes];
  uint8_t round_seeds[MEDS_t][MEDS_round_seed_bytes];

  XOF((uint8_t*[]){alpha, (uint8_t*)round_seeds},
      (size_t[]){MEDS_salt_bytes, MEDS_t * MEDS_round_seed_bytes},
      delta, MEDS_sec_seed_bytes,
      2);


  Fq A_tilde_data[MEDS_t * MEDS_m * MEDS_m];
  Fq B_tilde_data[MEDS_t * MEDS_n * MEDS_n];
  Fq C_tilde_data[MEDS_t * MEDS_k * MEDS_k];

  Fq *A_tilde[MEDS_t];
  Fq *B_tilde[MEDS_t];
  Fq *C_tilde[MEDS_t];

  for (int i = 0; i < MEDS_t; i++)
  {
    A_tilde[i] = &A_tilde_data[i * MEDS_m * MEDS_m];
    B_tilde[i] = &B_tilde_data[i * MEDS_n * MEDS_n];
    C_tilde[i] = &C_tilde_data[i * MEDS_k * MEDS_k];
  }

  keccak_state h_shake;
  shake256_init(&h_shake);

  for (int i = 0; i < MEDS_t; i++)
  {
    Fq G_tilde_ti[MEDS_k * MEDS_m * MEDS_n];

    derive_round_matrices(A_tilde[i], B_tilde[i], C_tilde[i], round_seeds[i], alpha, (uint32_t)i);

    LOG_MAT_FMT(A_tilde[i], MEDS_m, MEDS_m, "A_tilde[%i]", i);
    LOG_MAT_FMT(B_tilde[i], MEDS_n, MEDS_n, "B_tilde[%i]", i);
    LOG_MAT_FMT(C_tilde[i], MEDS_k, MEDS_k, "C_tilde[%i]", i);


    phi(G_tilde_ti, A_tilde[i], B_tilde[i], C_tilde[i], G_0);


    LOG_MAT_FMT(G_tilde_ti, MEDS_k, MEDS_m*MEDS_n, "G_tilde[%i]", i);

    uint8_t bs_buf[MEDS_G_BYTES];

    encode_G_full(bs_buf, G_tilde_ti);
    shake256_absorb(&h_shake, bs_buf, MEDS_G_BYTES);
  }

  shake256_absorb(&h_shake, (uint8_t*)m, mlen);

  shake256_finalize(&h_shake);

  uint8_t digest[MEDS_digest_bytes];

  shake256_squeeze(digest, MEDS_digest_bytes, &h_shake);

  LOG_VEC(digest, MEDS_digest_bytes);


  uint8_t h[MEDS_t];

  parse_hash(digest, MEDS_digest_bytes, h, MEDS_t);

  LOG_VEC(h, MEDS_t);


  memset(sm, 0, MEDS_SIG_BYTES);

  bitstream_t bs;

  bs_init(&bs, sm, MEDS_RESPONSE_BYTES);

  size_t response_count = 0;

  for (int i = 0; i < MEDS_t; i++)
  {
    if (h[i] > 0)
    {
      Fq mu[MEDS_m*MEDS_m];
      Fq nu[MEDS_n*MEDS_n];
      Fq eta[MEDS_k*MEDS_k];

      pmod_mat_mul(mu, MEDS_m, MEDS_m, A_tilde[i], MEDS_m, MEDS_m, A_inv[h[i]], MEDS_m, MEDS_m);
      pmod_mat_mul(nu, MEDS_n, MEDS_n, B_inv[h[i]], MEDS_n, MEDS_n, B_tilde[i], MEDS_n, MEDS_n);
      pmod_mat_mul(eta, MEDS_k, MEDS_k, C_tilde[i], MEDS_k, MEDS_k, C_inv[h[i]], MEDS_k, MEDS_k);

      LOG_MAT(mu, MEDS_m, MEDS_m);
      LOG_MAT(nu, MEDS_n, MEDS_n);
      LOG_MAT(eta, MEDS_k, MEDS_k);

      encode_mat_to_bs(&bs, mu, MEDS_m*MEDS_m);
      encode_mat_to_bs(&bs, nu, MEDS_n*MEDS_n);
      encode_mat_to_bs(&bs, eta, MEDS_k*MEDS_k);

      response_count++;
    }
  }

  if (response_count != MEDS_w || bs.byte_pos != MEDS_RESPONSE_BYTES || bs.bit_pos != 0)
    return -1;

  uint8_t *seed_out = sm + MEDS_ZERO_SEED_OFFSET;
  size_t zero_seed_count = 0;

  for (int i = 0; i < MEDS_t; i++)
  {
    if (h[i] == 0)
    {
      memcpy(seed_out, round_seeds[i], MEDS_round_seed_bytes);
      seed_out += MEDS_round_seed_bytes;
      zero_seed_count++;
    }
  }

  if (zero_seed_count != MEDS_ZERO_SEED_COUNT || seed_out != sm + MEDS_DIGEST_OFFSET)
    return -1;

  memcpy(sm + MEDS_DIGEST_OFFSET, digest, MEDS_digest_bytes);
  memcpy(sm + MEDS_SALT_OFFSET, alpha, MEDS_salt_bytes);
  memcpy(sm + MEDS_SIG_BYTES, m, mlen);

  *smlen = MEDS_SIG_BYTES + mlen;

  LOG_HEX(sm, MEDS_SIG_BYTES + mlen);

  return 0;
}

int crypto_sign_open(
    unsigned char *m, unsigned long long *mlen,
    const unsigned char *sm, unsigned long long smlen,
    const unsigned char *pk
  )
{
  LOG_HEX(sm, smlen);

  if (smlen < MEDS_SIG_BYTES)
    return -1;

  Fq G_data[MEDS_k*MEDS_m*MEDS_n * MEDS_s];
  Fq *G[MEDS_s];

  for (int i = 0; i < MEDS_s; i++)
    G[i] = &G_data[i * MEDS_k * MEDS_m * MEDS_n];


  load_public_key_matrices(G, pk);

  for (int i = 0; i < MEDS_s; i++)
    LOG_MAT_FMT(G[i], MEDS_k, MEDS_m*MEDS_n, "G[%i]", i);

  const uint8_t *digest = sm + MEDS_DIGEST_OFFSET;

  const uint8_t *alpha = sm + MEDS_SALT_OFFSET;

  LOG_VEC(digest, MEDS_digest_bytes);

  uint8_t h[MEDS_t];

  parse_hash(digest, MEDS_digest_bytes, h, MEDS_t);


  bitstream_t bs;

  bs_init(&bs, (uint8_t*)sm, MEDS_RESPONSE_BYTES);

  const uint8_t *seed_in = sm + MEDS_ZERO_SEED_OFFSET;
  size_t response_count = 0;
  size_t zero_seed_count = 0;

  Fq G_hat_i[MEDS_k*MEDS_m*MEDS_n];

  Fq mu[MEDS_m*MEDS_m];
  Fq nu[MEDS_n*MEDS_n];
  Fq eta[MEDS_k*MEDS_k];

  keccak_state shake;
  shake256_init(&shake);

  for (int i = 0; i < MEDS_t; i++)
  {
    if (h[i] > 0)
    {
      decode_mat_from_bs(&bs, mu, MEDS_m*MEDS_m);
      decode_mat_from_bs(&bs, nu, MEDS_n*MEDS_n);
      decode_mat_from_bs(&bs, eta, MEDS_k*MEDS_k);


      LOG_MAT_FMT(mu, MEDS_m, MEDS_m, "mu[%i]", i);
      LOG_MAT_FMT(nu, MEDS_n, MEDS_n, "nu[%i]", i);
      LOG_MAT_FMT(eta, MEDS_k, MEDS_k, "eta[%i]", i);


      Fq inv_mu[MEDS_m*MEDS_m];
      Fq inv_nu[MEDS_n*MEDS_n];
      Fq inv_eta[MEDS_k*MEDS_k];

      if (pmod_mat_inv(inv_mu, mu, MEDS_m, MEDS_m) < 0)
        return -1;

      if (pmod_mat_inv(inv_nu, nu, MEDS_n, MEDS_n) < 0)
        return -1;

      if (pmod_mat_inv(inv_eta, eta, MEDS_k, MEDS_k) < 0)
        return -1;


      phi(G_hat_i, mu, nu, eta, G[h[i]]);

      LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);

      response_count++;
    }
    else
    {
      LOG_VEC_FMT(seed_in, MEDS_round_seed_bytes, "seeds[%i]", i);

      Fq A_tilde[MEDS_m*MEDS_m];
      Fq B_tilde[MEDS_n*MEDS_n];
      Fq C_tilde[MEDS_k*MEDS_k];

      derive_round_matrices(A_tilde, B_tilde, C_tilde, seed_in, alpha, (uint32_t)i);

      LOG_MAT_FMT(A_tilde, MEDS_m, MEDS_m, "A_tilde[%i]", i);
      LOG_MAT_FMT(B_tilde, MEDS_n, MEDS_n, "B_tilde[%i]", i);
      LOG_MAT_FMT(C_tilde, MEDS_k, MEDS_k, "C_tilde[%i]", i);


      phi(G_hat_i, A_tilde, B_tilde, C_tilde, G[0]);

      LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);

      seed_in += MEDS_round_seed_bytes;
      zero_seed_count++;
    }


    uint8_t bs_buf[MEDS_G_BYTES];

    encode_G_full(bs_buf, G_hat_i);
 
    shake256_absorb(&shake, bs_buf, MEDS_G_BYTES);
  }

  if (response_count != MEDS_w ||
      zero_seed_count != MEDS_ZERO_SEED_COUNT ||
      seed_in != sm + MEDS_DIGEST_OFFSET ||
      bs.byte_pos != MEDS_RESPONSE_BYTES ||
      bs.bit_pos != 0)
    return -1;

  shake256_absorb(&shake, (uint8_t*)(sm + MEDS_SIG_BYTES), smlen - MEDS_SIG_BYTES);

  shake256_finalize(&shake);

  uint8_t check[MEDS_digest_bytes];

  shake256_squeeze(check, MEDS_digest_bytes, &shake);

  if (memcmp(digest, check, MEDS_digest_bytes) != 0)
  {
    fprintf(stderr, "Signature verification failed!\n");

    return -1;
  }

  memcpy(m, (uint8_t*)(sm + MEDS_SIG_BYTES), smlen - MEDS_SIG_BYTES);
  *mlen = smlen - MEDS_SIG_BYTES;

  return 0;
}
