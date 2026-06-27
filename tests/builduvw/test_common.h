#ifndef TEST_BUILDUVW_COMMON_H
#define TEST_BUILDUVW_COMMON_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "canonical.h"
#include "field.h"
#include "matrixelim.h"
#include "matrixmod.h"
#include "triform.h"

#define CHECK(cond) do { \
  if (!(cond)) { \
    fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #cond); \
    exit(1); \
  } \
} while (0)

#define CHECK_CTX(cond, seed, n, round, label) do { \
  if (!(cond)) { \
    fprintf(stderr, \
        "%s:%d: check failed: %s\n" \
        "  label=%s seed=0x%08x n=%d round=%d\n", \
        __FILE__, __LINE__, #cond, \
        (label), (unsigned)(seed), (n), (round)); \
    exit(1); \
  } \
} while (0)

typedef struct {
  uint32_t state;
} test_rng_t;

static inline void test_rng_init(test_rng_t *rng, uint32_t seed)
{
  rng->state = (seed == 0) ? 0x9e3779b9u : seed;
}

static inline uint32_t test_rng_next(test_rng_t *rng)
{
  CHECK(rng->state != 0);

  uint32_t x = rng->state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  rng->state = x;
  return x;
}

static inline Fq test_rng_fq(test_rng_t *rng)
{
  return (Fq)(test_rng_next(rng) % MEDS_p);
}

static inline size_t ref_index(int slice, int row, int col, int n)
{
  return (size_t)slice * (size_t)n * (size_t)n
       + (size_t)row * (size_t)n
       + (size_t)col;
}

static inline size_t mat_count(int n)
{
  return (size_t)n * (size_t)n;
}

static inline size_t form_count(int n)
{
  return (size_t)n * (size_t)n * (size_t)n;
}

static inline void fill_random(Fq *out, size_t count, test_rng_t *rng)
{
  for (size_t i = 0; i < count; i++)
    out[i] = test_rng_fq(rng);
}

static inline void fill_constant(Fq *out, size_t count, Fq value)
{
  for (size_t i = 0; i < count; i++)
    out[i] = value;
}

static inline void fill_zero(Fq *out, size_t count)
{
  memset(out, 0, count * sizeof(*out));
}

static inline void fill_basis_vector(Fq *v, int index, int n)
{
  fill_zero(v, (size_t)n);
  v[index] = 1;
}

static inline void vec_scale(Fq *out, Fq scalar, const Fq *v, int n)
{
  for (int i = 0; i < n; i++)
    out[i] = GF_mul(scalar, v[i]);
}

static inline int vec_nonzero(const Fq *v, int n)
{
  for (int i = 0; i < n; i++)
    if (v[i] != 0)
      return 1;

  return 0;
}

static inline void expect_vec_equal(const Fq *a, const Fq *b, int n)
{
  for (int i = 0; i < n; i++)
    CHECK(a[i] == b[i]);
}

static inline void expect_mat_equal(const Fq *a, const Fq *b, int n)
{
  for (int i = 0; i < n * n; i++)
    CHECK(a[i] == b[i]);
}

static inline void expect_mat_equal_ctx(
    const Fq *a,
    const Fq *b,
    int n,
    uint32_t seed,
    int round,
    const char *label)
{
  for (int i = 0; i < n * n; i++)
    CHECK_CTX(a[i] == b[i], seed, n, round, label);
}

static inline int vectors_projectively_equal(
    const Fq *a,
    const Fq *b,
    int n)
{
  int pivot = -1;

  if (!vec_nonzero(a, n) || !vec_nonzero(b, n))
    return 0;

  for (int i = 0; i < n; i++)
  {
    if (a[i] != 0 && b[i] != 0)
    {
      pivot = i;
      break;
    }
  }

  if (pivot < 0)
    return 0;

  for (int i = 0; i < n; i++)
  {
    Fq left = GF_mul(a[i], b[pivot]);
    Fq right = GF_mul(b[i], a[pivot]);

    if (left != right)
      return 0;
  }

  return 1;
}

static inline void expect_projective_columns_equal(
    const Fq *A,
    const Fq *B,
    int n)
{
  Fq a_col[MEDS_n];
  Fq b_col[MEDS_n];

  for (int col = 0; col < n; col++)
  {
    pmod_mat_get_col(a_col, A, col, n);
    pmod_mat_get_col(b_col, B, col, n);
    CHECK(vectors_projectively_equal(a_col, b_col, n));
  }
}

static inline int outputs_unchanged(
    const Fq *U,
    const Fq *V,
    const Fq *W,
    Fq sentinel,
    int n)
{
  for (int i = 0; i < n * n; i++)
    if (U[i] != sentinel || V[i] != sentinel || W[i] != sentinel)
      return 0;

  return 1;
}

static inline void force_coordinate_path_constraints(Fq *M, int n)
{
  for (int i = 0; i < n; i++)
  {
    for (int slice = 0; slice < n; slice++)
      M[ref_index(slice, i, i, n)] = 0;

    for (int row = 0; row < n; row++)
      M[ref_index(i, row, i, n)] = 0;

    if (i < n - 1)
      for (int col = 0; col < n; col++)
        M[ref_index(i, i + 1, col, n)] = 0;
  }
}

static inline void make_coordinate_path_candidate(Fq *M, Fq *u1, int n, uint32_t seed)
{
  test_rng_t rng;

  test_rng_init(&rng, seed);
  fill_random(M, form_count(n), &rng);
  force_coordinate_path_constraints(M, n);
  fill_basis_vector(u1, 0, n);
}

static inline void make_initial_phi_u_full_rank(Fq *M, Fq *u1, int n)
{
  fill_zero(M, form_count(n));
  fill_basis_vector(u1, 0, n);

  for (int row = 0; row < n; row++)
    M[ref_index(row, 0, row, n)] = 1;
}

static inline void make_initial_phi_u_corank_two(Fq *M, Fq *u1, int n)
{
  fill_zero(M, form_count(n));
  fill_basis_vector(u1, 0, n);

  for (int row = 0; row < n - 2; row++)
    M[ref_index(row, 0, row, n)] = 1;
}

static inline void make_u_dependency_fixture(Fq *M, Fq *u1, int n)
{
  fill_zero(M, form_count(n));
  fill_basis_vector(u1, 0, n);

  for (int a = 0; a < n; a++)
    for (int b = 0; b < n; b++)
      if (b != a)
        M[ref_index((a + b) % n, a, b, n)] = 1;
}

#endif
