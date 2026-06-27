#ifndef TEST_TRIFORM_COMMON_H
#define TEST_TRIFORM_COMMON_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "field.h"
#include "matrixmod.h"
#include "triform.h"

#define CHECK(cond) do { \
  if (!(cond)) { \
    fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #cond); \
    exit(1); \
  } \
} while (0)

typedef struct {
  uint32_t state;
} test_rng_t;

static inline uint32_t test_rng_next(test_rng_t *rng)
{
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

static inline void fill_random(Fq *out, size_t count, test_rng_t *rng)
{
  for (size_t i = 0; i < count; i++)
    out[i] = test_rng_fq(rng);
}

static inline void fill_zero(Fq *out, size_t count)
{
  memset(out, 0, count * sizeof(*out));
}

static inline void fill_identity(Fq *out, int n)
{
  fill_zero(out, triform_slice_size(n));

  for (int i = 0; i < n; i++)
    out[i * n + i] = 1;
}

static inline int vec_is_nonzero(const Fq *v, int n)
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

static inline void expect_triform_equal(const Fq *a, const Fq *b, int n)
{
  for (size_t i = 0; i < triform_element_count(n); i++)
    CHECK(a[i] == b[i]);
}

static inline void ref_mat_vec_mul(Fq *out, const Fq *A, const Fq *v, int n)
{
  for (int row = 0; row < n; row++)
  {
    Fq acc = 0;

    for (int col = 0; col < n; col++)
      acc = GF_add(acc, GF_mul(A[row * n + col], v[col]));

    out[row] = acc;
  }
}

static inline void ref_mat_transpose_vec_mul(Fq *out, const Fq *A, const Fq *v, int n)
{
  for (int col = 0; col < n; col++)
  {
    Fq acc = 0;

    for (int row = 0; row < n; row++)
      acc = GF_add(acc, GF_mul(A[row * n + col], v[row]));

    out[col] = acc;
  }
}

static inline void ref_mat_mul(Fq *out, const Fq *A, const Fq *B, int n)
{
  for (int row = 0; row < n; row++)
    for (int col = 0; col < n; col++)
    {
      Fq acc = 0;

      for (int k = 0; k < n; k++)
        acc = GF_add(acc, GF_mul(A[row * n + k], B[k * n + col]));

      out[row * n + col] = acc;
    }
}

static inline void ref_matrix_at_w(Fq *out, const Fq *M, const Fq *w, int n)
{
  for (int row = 0; row < n; row++)
    for (int col = 0; col < n; col++)
    {
      Fq acc = 0;

      for (int slice = 0; slice < n; slice++)
      {
        size_t idx = triform_index(slice, row, col, n);
        acc = GF_add(acc, GF_mul(w[slice], M[idx]));
      }

      out[row * n + col] = acc;
    }
}

static inline void ref_phi_u(Fq *out, const Fq *M, const Fq *u, int n)
{
  for (int slice = 0; slice < n; slice++)
    for (int row = 0; row < n; row++)
    {
      Fq acc = 0;

      for (int col = 0; col < n; col++)
      {
        size_t idx = triform_index(slice, col, row, n);
        acc = GF_add(acc, GF_mul(M[idx], u[col]));
      }

      out[row * n + slice] = acc;
    }
}

static inline void ref_phi_v(Fq *out, const Fq *M, const Fq *v, int n)
{
  for (int slice = 0; slice < n; slice++)
    for (int row = 0; row < n; row++)
    {
      Fq acc = 0;

      for (int col = 0; col < n; col++)
      {
        size_t idx = triform_index(slice, row, col, n);
        acc = GF_add(acc, GF_mul(M[idx], v[col]));
      }

      out[row * n + slice] = acc;
    }
}

static inline Fq ref_eval(const Fq *M, const Fq *u, const Fq *v, const Fq *w, int n)
{
  Fq acc = 0;

  for (int slice = 0; slice < n; slice++)
    for (int row = 0; row < n; row++)
      for (int col = 0; col < n; col++)
      {
        Fq term = GF_mul(u[row], M[triform_index(slice, row, col, n)]);
        term = GF_mul(term, v[col]);
        term = GF_mul(term, w[slice]);
        acc = GF_add(acc, term);
      }

  return acc;
}

static inline Fq ref_bilinear_eval(const Fq *A, const Fq *u, const Fq *v, int n)
{
  Fq acc = 0;

  for (int row = 0; row < n; row++)
    for (int col = 0; col < n; col++)
    {
      Fq term = GF_mul(u[row], A[row * n + col]);
      term = GF_mul(term, v[col]);
      acc = GF_add(acc, term);
    }

  return acc;
}

static inline Fq ref_mat_form_eval(const Fq *A, const Fq *left, const Fq *right, int n)
{
  return ref_bilinear_eval(A, left, right, n);
}

static inline Fq ref_right_eval(const Fq *A, const Fq *left, const Fq *right, int n)
{
  return ref_bilinear_eval(A, left, right, n);
}

static inline void ref_action_pullback(
    Fq *out,
    const Fq *M,
    const Fq *A,
    const Fq *B,
    const Fq *C,
    int n)
{
  for (int out_slice = 0; out_slice < n; out_slice++)
    for (int row = 0; row < n; row++)
      for (int col = 0; col < n; col++)
      {
        Fq acc = 0;

        for (int in_slice = 0; in_slice < n; in_slice++)
          for (int a_row = 0; a_row < n; a_row++)
            for (int b_col = 0; b_col < n; b_col++)
            {
              Fq term = GF_mul(A[a_row * n + row], C[in_slice * n + out_slice]);
              term = GF_mul(term, M[triform_index(in_slice, a_row, b_col, n)]);
              term = GF_mul(term, B[b_col * n + col]);
              acc = GF_add(acc, term);
            }

        out[triform_index(out_slice, row, col, n)] = acc;
      }
}

static inline void vec_add(Fq *out, const Fq *a, const Fq *b, int n)
{
  for (int i = 0; i < n; i++)
    out[i] = GF_add(a[i], b[i]);
}

static inline void vec_scale(Fq *out, Fq scalar, const Fq *v, int n)
{
  for (int i = 0; i < n; i++)
    out[i] = GF_mul(scalar, v[i]);
}

#endif
