#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "randombytes.h"

#include "params.h"
#include "api.h"
#include "meds.h"
#include "util.h"

static void mat_identity(pmod_mat_t *M, int dim)
{
  for (int r = 0; r < dim; r++)
    for (int c = 0; c < dim; c++)
      M[r * dim + c] = (r == c) ? 1 : 0;
}

static void mat_identity_with_swap01(pmod_mat_t *M, int dim)
{
  mat_identity(M, dim);

  M[0 * dim + 0] = 0;
  M[0 * dim + 1] = 1;
  M[1 * dim + 0] = 1;
  M[1 * dim + 1] = 0;
}

static int mat_equal(const pmod_mat_t *A, const pmod_mat_t *B, int len)
{
  for (int i = 0; i < len; i++)
    if (A[i] != B[i])
      return 0;

  return 1;
}

static void fill_test_G(pmod_mat_t *G)
{
  for (int i = 0; i < MEDS_k * MEDS_m * MEDS_n; i++)
    G[i] = (pmod_mat_t)((i + 1) % MEDS_p);
}

static int test_phi_identity_C(void)
{
  uint8_t seedA[MEDS_sec_seed_bytes] = {1};
  uint8_t seedB[MEDS_sec_seed_bytes] = {2};
  uint8_t seedG[MEDS_pub_seed_bytes] = {3};

  pmod_mat_t A[MEDS_m * MEDS_m];
  pmod_mat_t B[MEDS_n * MEDS_n];
  pmod_mat_t C[MEDS_k * MEDS_k];
  pmod_mat_t G[MEDS_k * MEDS_m * MEDS_n];
  pmod_mat_t G_pi[MEDS_k * MEDS_m * MEDS_n];
  pmod_mat_t G_phi[MEDS_k * MEDS_m * MEDS_n];

  rnd_inv_matrix(A, MEDS_m, MEDS_m, seedA, sizeof(seedA));
  rnd_inv_matrix(B, MEDS_n, MEDS_n, seedB, sizeof(seedB));
  rnd_sys_mat(G, MEDS_k, MEDS_m * MEDS_n, seedG, sizeof(seedG));
  mat_identity(C, MEDS_k);

  pi(G_pi, A, B, G);
  phi(G_phi, A, B, C, G);

  return mat_equal(G_pi, G_phi, MEDS_k * MEDS_m * MEDS_n) ? 0 : -1;
}

static int test_phi_swap_C(void)
{
  pmod_mat_t A[MEDS_m * MEDS_m];
  pmod_mat_t B[MEDS_n * MEDS_n];
  pmod_mat_t C[MEDS_k * MEDS_k];
  pmod_mat_t G[MEDS_k * MEDS_m * MEDS_n];
  pmod_mat_t G_phi[MEDS_k * MEDS_m * MEDS_n];

  mat_identity(A, MEDS_m);
  mat_identity(B, MEDS_n);
  mat_identity_with_swap01(C, MEDS_k);
  fill_test_G(G);

  phi(G_phi, A, B, C, G);

  const int block_len = MEDS_m * MEDS_n;

  if (!mat_equal(&G_phi[0 * block_len], &G[1 * block_len], block_len))
    return -1;

  if (!mat_equal(&G_phi[1 * block_len], &G[0 * block_len], block_len))
    return -1;

  for (int i = 2; i < MEDS_k; i++)
    if (!mat_equal(&G_phi[i * block_len], &G[i * block_len], block_len))
      return -1;

  return 0;
}

static int test_phi_composition(void)
{
  uint8_t seedA[MEDS_sec_seed_bytes] = {10};
  uint8_t seedB[MEDS_sec_seed_bytes] = {11};
  uint8_t seedC[MEDS_sec_seed_bytes] = {12};
  uint8_t seedAt[MEDS_sec_seed_bytes] = {13};
  uint8_t seedBt[MEDS_sec_seed_bytes] = {14};
  uint8_t seedCt[MEDS_sec_seed_bytes] = {15};
  uint8_t seedG[MEDS_pub_seed_bytes] = {16};

  pmod_mat_t A[MEDS_m * MEDS_m];
  pmod_mat_t B[MEDS_n * MEDS_n];
  pmod_mat_t C[MEDS_k * MEDS_k];
  pmod_mat_t At[MEDS_m * MEDS_m];
  pmod_mat_t Bt[MEDS_n * MEDS_n];
  pmod_mat_t Ct[MEDS_k * MEDS_k];
  pmod_mat_t A_inv[MEDS_m * MEDS_m];
  pmod_mat_t B_inv[MEDS_n * MEDS_n];
  pmod_mat_t C_inv[MEDS_k * MEDS_k];
  pmod_mat_t mu[MEDS_m * MEDS_m];
  pmod_mat_t nu[MEDS_n * MEDS_n];
  pmod_mat_t eta[MEDS_k * MEDS_k];
  pmod_mat_t G[MEDS_k * MEDS_m * MEDS_n];
  pmod_mat_t G_public[MEDS_k * MEDS_m * MEDS_n];
  pmod_mat_t G_left[MEDS_k * MEDS_m * MEDS_n];
  pmod_mat_t G_right[MEDS_k * MEDS_m * MEDS_n];

  rnd_inv_matrix(A, MEDS_m, MEDS_m, seedA, sizeof(seedA));
  rnd_inv_matrix(B, MEDS_n, MEDS_n, seedB, sizeof(seedB));
  rnd_inv_matrix(C, MEDS_k, MEDS_k, seedC, sizeof(seedC));
  rnd_inv_matrix(At, MEDS_m, MEDS_m, seedAt, sizeof(seedAt));
  rnd_inv_matrix(Bt, MEDS_n, MEDS_n, seedBt, sizeof(seedBt));
  rnd_inv_matrix(Ct, MEDS_k, MEDS_k, seedCt, sizeof(seedCt));
  rnd_sys_mat(G, MEDS_k, MEDS_m * MEDS_n, seedG, sizeof(seedG));

  if (pmod_mat_inv(A_inv, A, MEDS_m, MEDS_m) < 0)
    return -1;

  if (pmod_mat_inv(B_inv, B, MEDS_n, MEDS_n) < 0)
    return -1;

  if (pmod_mat_inv(C_inv, C, MEDS_k, MEDS_k) < 0)
    return -1;

  pmod_mat_mul(mu, MEDS_m, MEDS_m, At, MEDS_m, MEDS_m, A_inv, MEDS_m, MEDS_m);
  pmod_mat_mul(nu, MEDS_n, MEDS_n, B_inv, MEDS_n, MEDS_n, Bt, MEDS_n, MEDS_n);
  pmod_mat_mul(eta, MEDS_k, MEDS_k, Ct, MEDS_k, MEDS_k, C_inv, MEDS_k, MEDS_k);

  phi(G_public, A, B, C, G);
  phi(G_left, mu, nu, eta, G_public);
  phi(G_right, At, Bt, Ct, G);

  return mat_equal(G_left, G_right, MEDS_k * MEDS_m * MEDS_n) ? 0 : -1;
}

static int run_phi_tests(void)
{
  if (test_phi_identity_C() != 0)
  {
    fprintf(stderr, "test_phi_identity_C failed\n");
    return -1;
  }

  if (test_phi_swap_C() != 0)
  {
    fprintf(stderr, "test_phi_swap_C failed\n");
    return -1;
  }

  if (test_phi_composition() != 0)
  {
    fprintf(stderr, "test_phi_composition failed\n");
    return -1;
  }

  printf("phi tests passed\n");
  return 0;
}

double osfreq(void);

long long cpucycles(void)
{
  unsigned long long result;
  asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
      : "=a" (result) ::  "%rdx");
  return result;
}

int main(int argc, char *argv[])
{
  printf("paramter set: %s\n\n", MEDS_name);

  long long time = 0;
  long long keygen_time = 0xfffffffffffffff;
  long long sign_time = 0xfffffffffffffff;
  long long verify_time = 0xfffffffffffffff;

  int rounds = 1;

  if (argc > 1)
    rounds = atoi(argv[1]);

  unsigned char entropy_input[48] = {0};

  randombytes_init(entropy_input, NULL, 256);

  if (run_phi_tests() != 0)
    return 1;

  char msg[4] = "Test";

  printf("pk:  %i bytes\n", CRYPTO_PUBLICKEYBYTES);
  printf("sk:  %i bytes\n", CRYPTO_SECRETKEYBYTES);
  printf("sig: %i bytes\n", CRYPTO_BYTES);
  printf("\n");

  for (int round = 0; round < rounds; round++)
  {
    uint8_t sk[CRYPTO_SECRETKEYBYTES] = {0};
    uint8_t pk[CRYPTO_PUBLICKEYBYTES] = {0};

    time = -cpucycles();
    crypto_sign_keypair(pk, sk);
    time += cpucycles();

    if (time < keygen_time) keygen_time = time;

    uint8_t sig[CRYPTO_BYTES + sizeof(msg)] = {0};
    unsigned long long sig_len = sizeof(sig);

    time = -cpucycles();
    crypto_sign(sig, &sig_len, (const unsigned char *)msg, sizeof(msg), sk);
    time += cpucycles();

    if (time < sign_time) sign_time = time;

    unsigned char msg_out[4];
    unsigned long long msg_out_len = sizeof(msg_out);


    time = -cpucycles();
    int ret = crypto_sign_open(msg_out, &msg_out_len, sig, sizeof(sig), pk);
    time += cpucycles();

    if (time < verify_time) verify_time = time;

    if (ret == 0)
      printf("success\n");
    else
      printf("!!! FAILED !!!\n");
  }

  double freq = osfreq();

  printf("\n");
  printf("Time (min of %i runs):\n", rounds);
  printf("keygen: %f   (%llu cycles)\n", keygen_time / freq, keygen_time);
  printf("sign:   %f   (%llu cycles)\n", sign_time / freq, sign_time);
  printf("verify: %f   (%llu cycles)\n", verify_time / freq, verify_time);

  return 0;
}
