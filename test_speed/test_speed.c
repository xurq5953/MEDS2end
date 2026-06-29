#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../api.h"
#include "../canonical.h"
#include "../field.h"
#include "../matrixelim.h"
#include "../matrixmod.h"
#include "../params.h"
#include "../randombytes.h"
#include "../triform.h"
#include "cpucycles.h"
#include "speed_print.h"

#define DEFAULT_SPEED_ROUNDS 10
#define DEFAULT_PROTOCOL_ROUNDS 16

static uint64_t *cycles = NULL;
static int speed_rounds = DEFAULT_SPEED_ROUNDS;
static int protocol_rounds = DEFAULT_PROTOCOL_ROUNDS;

typedef int (*speed_fn)(void *ctx);

typedef struct {
  uint32_t state;
} speed_rng_t;

typedef struct {
  Fq a;
  Fq b;
  Fq c;
  Fq vec[TRINE_n];
  Fq vec2[TRINE_n];
  Fq mat[TRINE_n * TRINE_n];
  Fq mat2[TRINE_n * TRINE_n];
  Fq mat3[TRINE_n * TRINE_n];
  Fq form[TRINE_n * TRINE_n * TRINE_n];
  Fq form2[TRINE_n * TRINE_n * TRINE_n];
  Fq U[TRINE_n * TRINE_n];
  Fq V[TRINE_n * TRINE_n];
  Fq W[TRINE_n * TRINE_n];
  unsigned char pk[CRYPTO_PUBLICKEYBYTES];
  unsigned char sk[CRYPTO_SECRETKEYBYTES];
  unsigned char sm[CRYPTO_BYTES + 32];
  unsigned char msg[32];
  unsigned char opened[32];
  unsigned long long smlen;
  unsigned long long opened_len;
} speed_ctx_t;

static uint32_t rng_next(speed_rng_t *rng)
{
  uint32_t x = rng->state;

  if (x == 0)
    x = 0x9e3779b9u;

  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  rng->state = x;

  return x;
}

static Fq rng_fq(speed_rng_t *rng)
{
  return (Fq)(rng_next(rng) % TRINE_q);
}

static void fill_random_fq(Fq *out, size_t count, speed_rng_t *rng)
{
  for (size_t i = 0; i < count; i++)
    out[i] = rng_fq(rng);
}

static void fill_basis_vector(Fq *v, int index, int n)
{
  memset(v, 0, (size_t)n * sizeof(*v));
  v[index] = 1;
}

static size_t triform_ref_index(int slice, int row, int col, int n)
{
  return (size_t)slice * (size_t)n * (size_t)n
       + (size_t)row * (size_t)n
       + (size_t)col;
}

static void force_coordinate_path_constraints(Fq *M, int n)
{
  for (int i = 0; i < n; i++)
  {
    for (int slice = 0; slice < n; slice++)
      M[triform_ref_index(slice, i, i, n)] = 0;

    for (int row = 0; row < n; row++)
      M[triform_ref_index(i, row, i, n)] = 0;

    if (i < n - 1)
      for (int col = 0; col < n; col++)
        M[triform_ref_index(i, i + 1, col, n)] = 0;
  }
}

static void make_builduvw_fixture(Fq *M, Fq *u1, int n)
{
  Fq U[TRINE_n * TRINE_n];
  Fq V[TRINE_n * TRINE_n];
  Fq W[TRINE_n * TRINE_n];

  for (uint32_t attempt = 0; attempt < 1024u; attempt++)
  {
    speed_rng_t rng = {0x5eed1234u + attempt};

    fill_random_fq(M, (size_t)n * (size_t)n * (size_t)n, &rng);
    force_coordinate_path_constraints(M, n);
    fill_basis_vector(u1, 0, n);

    if (canonical_build_uvw_vartime(U, V, W, M, u1, n) == 0)
      return;
  }

  fprintf(stderr, "failed to construct BuildUVW speed fixture\n");
  exit(1);
}

static void make_invertible_matrix(Fq *A, int n)
{
  pmod_mat_identity(A, n);

  for (int row = 0; row < n; row++)
    for (int col = 0; col < n; col++)
      if (row != col)
        A[row * n + col] = (Fq)((row + 2 * col + 1) % 7);
}

static void make_corank_one_matrix(Fq *A, int n)
{
  memset(A, 0, (size_t)n * (size_t)n * sizeof(*A));

  for (int i = 0; i < n - 1; i++)
    A[i * n + i] = 1;
}

static int time_function(const char *label, speed_fn fn, void *ctx, int rounds)
{
  if (rounds < 2)
    rounds = 2;

  for (int i = 0; i < rounds; i++)
  {
    cycles[i] = (uint64_t)cpucycles();
    if (fn(ctx) != 0)
    {
      fprintf(stderr, "%s failed at benchmark round %d\n", label, i);
      return -1;
    }
  }

  cycles[rounds] = (uint64_t)cpucycles();
  print_results(label, cycles, (size_t)rounds + 1u);
  return 0;
}

static int bench_keypair(void *arg)
{
  speed_ctx_t *ctx = arg;
  if (crypto_sign_keypair(ctx->pk, ctx->sk) != 0)
    return -1;

  return 0;
}

static int bench_sign(void *arg)
{
  speed_ctx_t *ctx = arg;
  ctx->smlen = 0;

  if (crypto_sign(ctx->sm, &ctx->smlen, ctx->msg, sizeof(ctx->msg), ctx->sk) != 0)
    return -1;

  if (ctx->smlen !=
      (unsigned long long)CRYPTO_BYTES +
      (unsigned long long)sizeof(ctx->msg))
    return -1;

  return 0;
}

static int bench_open(void *arg)
{
  speed_ctx_t *ctx = arg;
  memset(ctx->opened, 0xa5, sizeof(ctx->opened));
  ctx->opened_len = 0;

  if (crypto_sign_open(ctx->opened, &ctx->opened_len, ctx->sm, ctx->smlen, ctx->pk) != 0)
    return -1;

  if (ctx->opened_len != (unsigned long long)sizeof(ctx->msg))
    return -1;

  if (memcmp(ctx->opened, ctx->msg, sizeof(ctx->msg)) != 0)
    return -1;

  return 0;
}

static int bench_gf_add(void *arg)
{
  speed_ctx_t *ctx = arg;
  ctx->c = GF_add(ctx->a, ctx->b);
  return 0;
}

static int bench_gf_mul(void *arg)
{
  speed_ctx_t *ctx = arg;
  ctx->c = GF_mul(ctx->a, ctx->b);
  return 0;
}

static int bench_gf_inv(void *arg)
{
  speed_ctx_t *ctx = arg;
  ctx->c = GF_inv(ctx->a);
  return 0;
}

static int bench_mat_mul(void *arg)
{
  speed_ctx_t *ctx = arg;
  pmod_mat_mul(ctx->mat3, ctx->mat, ctx->mat2, TRINE_n);
  return 0;
}

static int bench_mat_vec_mul(void *arg)
{
  speed_ctx_t *ctx = arg;
  pmod_mat_vec_mul(ctx->vec2, ctx->mat, ctx->vec, TRINE_n);
  return 0;
}

static int bench_mat_rank(void *arg)
{
  speed_ctx_t *ctx = arg;
  ctx->c = (Fq)pmod_mat_rank_vartime(ctx->mat, TRINE_n);
  return 0;
}

static int bench_mat_inv(void *arg)
{
  speed_ctx_t *ctx = arg;
  (void)pmod_mat_inv_vartime(ctx->mat3, ctx->mat, TRINE_n);
  return 0;
}

static int bench_right_kernel(void *arg)
{
  speed_ctx_t *ctx = arg;
  (void)pmod_mat_right_kernel_corank1_vartime(ctx->vec2, ctx->mat2, TRINE_n);
  return 0;
}

static int bench_triform_matrix_at_w(void *arg)
{
  speed_ctx_t *ctx = arg;
  triform_matrix_at_w(ctx->mat3, ctx->form, ctx->vec, TRINE_n);
  return 0;
}

static int bench_triform_phi_u(void *arg)
{
  speed_ctx_t *ctx = arg;
  triform_phi_u(ctx->mat3, ctx->form, ctx->vec, TRINE_n);
  return 0;
}

static int bench_triform_phi_v(void *arg)
{
  speed_ctx_t *ctx = arg;
  triform_phi_v(ctx->mat3, ctx->form, ctx->vec, TRINE_n);
  return 0;
}

static int bench_triform_eval(void *arg)
{
  speed_ctx_t *ctx = arg;
  ctx->c = triform_eval(ctx->form, ctx->vec, ctx->vec2, ctx->vec, TRINE_n);
  return 0;
}

static int bench_triform_action_pullback(void *arg)
{
  speed_ctx_t *ctx = arg;
  triform_action_pullback(ctx->form2, ctx->form, ctx->mat, ctx->mat2, ctx->mat3, TRINE_n);
  return 0;
}

static int bench_builduvw(void *arg)
{
  speed_ctx_t *ctx = arg;
  (void)canonical_build_uvw_vartime(ctx->U, ctx->V, ctx->W, ctx->form, ctx->vec, TRINE_n);
  return 0;
}

static void init_context(speed_ctx_t *ctx)
{
  speed_rng_t rng = {0x12345678u};

  memset(ctx, 0, sizeof(*ctx));
  ctx->a = 7;
  ctx->b = 11;

  fill_random_fq(ctx->vec, TRINE_n, &rng);
  fill_random_fq(ctx->vec2, TRINE_n, &rng);
  fill_random_fq(ctx->form, (size_t)TRINE_n * TRINE_n * TRINE_n, &rng);
  make_invertible_matrix(ctx->mat, TRINE_n);
  make_invertible_matrix(ctx->mat2, TRINE_n);
  make_invertible_matrix(ctx->mat3, TRINE_n);
  make_corank_one_matrix(ctx->mat2, TRINE_n);
  make_builduvw_fixture(ctx->form, ctx->vec, TRINE_n);

  for (size_t i = 0; i < sizeof(ctx->msg); i++)
    ctx->msg[i] = (unsigned char)i;

  if (crypto_sign_keypair(ctx->pk, ctx->sk) != 0)
  {
    fprintf(stderr, "crypto_sign_keypair failed during speed setup\n");
    exit(1);
  }

  ctx->smlen = 0;
  if (crypto_sign(ctx->sm, &ctx->smlen, ctx->msg, sizeof(ctx->msg), ctx->sk) != 0)
  {
    fprintf(stderr, "crypto_sign failed during speed setup\n");
    exit(1);
  }

  if (ctx->smlen !=
      (unsigned long long)CRYPTO_BYTES +
      (unsigned long long)sizeof(ctx->msg))
  {
    fprintf(stderr, "crypto_sign produced an unexpected signed-message length\n");
    exit(1);
  }

  memset(ctx->opened, 0xa5, sizeof(ctx->opened));
  ctx->opened_len = 0;
  if (crypto_sign_open(ctx->opened, &ctx->opened_len, ctx->sm, ctx->smlen, ctx->pk) != 0)
  {
    fprintf(stderr, "crypto_sign_open failed during speed setup\n");
    exit(1);
  }

  if (ctx->opened_len != (unsigned long long)sizeof(ctx->msg) ||
      memcmp(ctx->opened, ctx->msg, sizeof(ctx->msg)) != 0)
  {
    fprintf(stderr, "crypto_sign_open recovered the wrong message during speed setup\n");
    exit(1);
  }
}

int main(int argc, char **argv)
{
  speed_ctx_t *ctx = NULL;

  if (argc > 1)
    speed_rounds = atoi(argv[1]);

  if (argc > 2)
    protocol_rounds = atoi(argv[2]);

  if (speed_rounds < 2)
    speed_rounds = 2;

  if (protocol_rounds < 2)
    protocol_rounds = 2;

  cycles = calloc((size_t)speed_rounds + 1u, sizeof(*cycles));
  if (protocol_rounds > speed_rounds)
  {
    free(cycles);
    cycles = calloc((size_t)protocol_rounds + 1u, sizeof(*cycles));
  }

  if (cycles == NULL)
  {
    fprintf(stderr, "failed to allocate cycle buffer\n");
    return 1;
  }

  ctx = calloc(1u, sizeof(*ctx));
  if (ctx == NULL)
  {
    fprintf(stderr, "failed to allocate speed context\n");
    free(cycles);
    return EXIT_FAILURE;
  }

  init_context(ctx);

  printf("name: %s\n", CRYPTO_ALGNAME);
  printf(
      "parameters (n, q, r, K, X): (%d, %d, %d, %d, %d)\n",
      TRINE_n,
      TRINE_q,
      TRINE_r,
      TRINE_K,
      TRINE_X);
  printf("pk: %d bytes\n", TRINE_PK_BYTES);
  printf("sk: %d bytes\n", TRINE_SK_BYTES);
  printf("sig: %d bytes\n", TRINE_SIG_BYTES);
  printf("rounds: core=%d protocol=%d\n\n", speed_rounds, protocol_rounds);

  if (time_function("crypto_sign_keypair", bench_keypair, ctx, protocol_rounds) != 0)
    goto failure;
  if (time_function("crypto_sign", bench_sign, ctx, protocol_rounds) != 0)
    goto failure;
  if (time_function("crypto_sign_open", bench_open, ctx, protocol_rounds) != 0)
    goto failure;

  if (time_function("GF_add", bench_gf_add, ctx, speed_rounds) != 0)
    goto failure;
  if (time_function("GF_mul", bench_gf_mul, ctx, speed_rounds) != 0)
    goto failure;
  if (time_function("GF_inv", bench_gf_inv, ctx, speed_rounds) != 0)
    goto failure;
  if (time_function("pmod_mat_mul", bench_mat_mul, ctx, speed_rounds) != 0)
    goto failure;
  if (time_function("pmod_mat_vec_mul", bench_mat_vec_mul, ctx, speed_rounds) != 0)
    goto failure;
  if (time_function("pmod_mat_rank_vartime", bench_mat_rank, ctx, speed_rounds) != 0)
    goto failure;
  if (time_function("pmod_mat_inv_vartime", bench_mat_inv, ctx, speed_rounds) != 0)
    goto failure;
  if (time_function("pmod_mat_right_kernel_corank1_vartime", bench_right_kernel, ctx, speed_rounds) != 0)
    goto failure;
  if (time_function("triform_matrix_at_w", bench_triform_matrix_at_w, ctx, speed_rounds) != 0)
    goto failure;
  if (time_function("triform_phi_u", bench_triform_phi_u, ctx, speed_rounds) != 0)
    goto failure;
  if (time_function("triform_phi_v", bench_triform_phi_v, ctx, speed_rounds) != 0)
    goto failure;
  if (time_function("triform_eval", bench_triform_eval, ctx, speed_rounds) != 0)
    goto failure;
  if (time_function("triform_action_pullback", bench_triform_action_pullback, ctx, speed_rounds) != 0)
    goto failure;
  if (time_function("canonical_build_uvw_vartime", bench_builduvw, ctx, speed_rounds) != 0)
    goto failure;

  free(ctx);
  free(cycles);
  return EXIT_SUCCESS;

failure:
  free(ctx);
  free(cycles);
  return EXIT_FAILURE;
}
