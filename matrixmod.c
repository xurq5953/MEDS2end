#include <stdio.h>
#include <string.h>

#include "matrixmod.h"

void pmod_mat_print(const pmod_mat_t *M, int M_r, int M_c)
{
  pmod_mat_fprint(stdout, M, M_r, M_c);
}

void pmod_mat_fprint(FILE *stream, const pmod_mat_t *M, int M_r, int M_c)
{
  for (int r = 0; r < M_r; r++)
  {
    fprintf(stream, "[");
    for (int c = 0; c < M_c-1; c++)
      fprintf(stream, "%4i ", pmod_mat_entry(M, M_r, M_c, r, c));
    fprintf(stream, "%4i", pmod_mat_entry(M, M_r, M_c, r, M_c-1));
    fprintf(stream, "]\n");
  }
}

void pmod_mat_zero(pmod_mat_t *A, int n)
{
  memset(A, 0, (size_t)n * n * sizeof(*A));
}

void pmod_mat_copy(pmod_mat_t *dst, const pmod_mat_t *src, int n)
{
  memcpy(dst, src, (size_t)n * n * sizeof(*dst));
}

void pmod_mat_identity(pmod_mat_t *A, int n)
{
  pmod_mat_zero(A, n);

  for (int i = 0; i < n; i++)
    A[i * n + i] = 1;
}

void pmod_mat_transpose(
    pmod_mat_t *AT,
    const pmod_mat_t *A,
    int n)
{
  if (AT == A)
  {
    for (int r = 0; r < n; r++)
      for (int c = r + 1; c < n; c++)
      {
        GFq_t tmp = AT[r * n + c];
        AT[r * n + c] = AT[c * n + r];
        AT[c * n + r] = tmp;
      }

    return;
  }

  for (int r = 0; r < n; r++)
    for (int c = 0; c < n; c++)
      AT[c * n + r] = A[r * n + c];
}

void pmod_mat_mul(
    pmod_mat_t *C,
    const pmod_mat_t *A,
    const pmod_mat_t *B,
    int n)
{
  pmod_mat_t tmp[n * n];

  for (int r = 0; r < n; r++)
    for (int c = 0; c < n; c++)
    {
      uint64_t acc = 0;

      for (int k = 0; k < n; k++)
        acc += (uint64_t)A[r * n + k] * B[k * n + c];

      tmp[r * n + c] = (GFq_t)(acc % MEDS_p);
    }

  memcpy(C, tmp, sizeof(tmp));
}

void pmod_mat_mul_rect(
    Fq *C,
    int C_r,
    int C_c,
    const Fq *A,
    int A_r,
    int A_c,
    const Fq *B,
    int B_r,
    int B_c)
{
  GFq_t tmp[C_r*C_c];

  for (int c = 0; c < C_c; c++)
    for (int r = 0; r < C_r; r++)
    {
      uint64_t val = 0;

      for (int i = 0; i < A_c; i++)
        val = (val + (uint64_t)pmod_mat_entry(A, A_r, A_c, r, i) * (uint64_t)pmod_mat_entry(B, B_r, B_c, i, c));

      tmp[r*C_c + c] = val % MEDS_p;
    }

  for (int c = 0; c < C_c; c++)
    for (int r = 0; r < C_r; r++)
      pmod_mat_set_entry(C, C_r, C_c, r, c, tmp[r*C_c + c]);
}

void pmod_mat_vec_mul(
    GFq_t *out,
    const pmod_mat_t *A,
    const GFq_t *v,
    int n)
{
  GFq_t tmp[n];

  for (int r = 0; r < n; r++)
  {
    uint64_t acc = 0;

    for (int c = 0; c < n; c++)
      acc += (uint64_t)A[r * n + c] * v[c];

    tmp[r] = (GFq_t)(acc % MEDS_p);
  }

  memcpy(out, tmp, sizeof(tmp));
}

void pmod_mat_transpose_vec_mul(
    GFq_t *out,
    const pmod_mat_t *A,
    const GFq_t *v,
    int n)
{
  GFq_t tmp[n];

  for (int c = 0; c < n; c++)
  {
    uint64_t acc = 0;

    for (int r = 0; r < n; r++)
      acc += (uint64_t)A[r * n + c] * v[r];

    tmp[c] = (GFq_t)(acc % MEDS_p);
  }

  memcpy(out, tmp, sizeof(tmp));
}

void pmod_mat_set_col(
    pmod_mat_t *A,
    int col,
    const GFq_t *v,
    int n)
{
  for (int r = 0; r < n; r++)
    A[r * n + col] = v[r];
}

void pmod_mat_get_col(
    GFq_t *v,
    const pmod_mat_t *A,
    int col,
    int n)
{
  for (int r = 0; r < n; r++)
    v[r] = A[r * n + col];
}

void pmod_mat_linear_combination(
    pmod_mat_t *out,
    const pmod_mat_t *matrices,
    const GFq_t *coeffs,
    int count,
    int n)
{
  for (int pos = 0; pos < n * n; pos++)
  {
    uint64_t acc = 0;

    for (int i = 0; i < count; i++)
      acc += (uint64_t)coeffs[i] * matrices[i * n * n + pos];

    out[pos] = (GFq_t)(acc % MEDS_p);
  }
}

void pmod_mat_diag_scale(
    pmod_mat_t *out,
    const GFq_t *left_diag,
    const pmod_mat_t *A,
    const GFq_t *right_diag,
    int n)
{
  pmod_mat_t tmp[n * n];

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    {
      GFq_t x = GF_mul(left_diag[i], A[i * n + j]);
      tmp[i * n + j] = GF_mul(x, right_diag[j]);
    }

  memcpy(out, tmp, sizeof(tmp));
}

int pmod_mat_syst_ct(Fq *M, int M_r, int M_c)
{
  if (pmod_mat_row_echelon_ct(M, M_r, M_c) < 0)
    return -1;

  return pmod_mat_back_substitution_ct(M, M_r, M_c);
}

int pmod_mat_row_echelon_ct(Fq *M, int M_r, int M_c)
{
  for (int r = 0; r < M_r; r++)
  {
    // Add lower rows to the pivot row until the pivot becomes nonzero.
    for (int r2 = r + 1; r2 < M_r; r2++)
    {
      uint64_t pivot_entry = pmod_mat_entry(M, M_r, M_c, r, r);

      for (int c = r; c < M_c; c++)
      {
        uint64_t candidate_entry = pmod_mat_entry(M, M_r, M_c, r2, c);
        uint64_t current_entry = pmod_mat_entry(M, M_r, M_c, r, c);
        uint64_t updated_entry =
            (current_entry + candidate_entry * (pivot_entry == 0)) % MEDS_p;

        pmod_mat_set_entry(M, M_r, M_c, r, c, updated_entry);
      }
    }

    uint64_t pivot = pmod_mat_entry(M, M_r, M_c, r, r);

    if (pivot == 0)
      return -1;

    uint64_t pivot_inv = GF_inv((GFq_t)pivot);

    // Normalize the pivot row.
    for (int c = r; c < M_c; c++)
    {
      uint64_t row_entry = pmod_mat_entry(M, M_r, M_c, r, c);
      uint64_t normalized_entry = (row_entry * pivot_inv) % MEDS_p;

      pmod_mat_set_entry(M, M_r, M_c, r, c, normalized_entry);
    }

    // Eliminate entries below the pivot.
    for (int r2 = r + 1; r2 < M_r; r2++)
    {
      uint64_t factor = pmod_mat_entry(M, M_r, M_c, r2, r);

      for (int c = r; c < M_c; c++)
      {
        uint64_t pivot_row_entry =
            pmod_mat_entry(M, M_r, M_c, r, c);
        uint64_t target_row_entry =
            pmod_mat_entry(M, M_r, M_c, r2, c);

        int64_t reduced_entry =
            (int64_t)target_row_entry -
            (int64_t)((pivot_row_entry * factor) % MEDS_p);

        reduced_entry += MEDS_p * (reduced_entry < 0);

        pmod_mat_set_entry(
            M, M_r, M_c, r2, c, (GFq_t)reduced_entry);
      }
    }
  }

  return 0;
}

int pmod_mat_back_substitution_ct(Fq *M, int M_r, int M_c)
{
  for (int r = M_r - 1; r >= 0; r--)
  {
    for (int r2 = 0; r2 < r; r2++)
    {
      uint64_t factor = pmod_mat_entry(M, M_r, M_c, r2, r);

      uint64_t pivot_entry =
          pmod_mat_entry(M, M_r, M_c, r, r);
      uint64_t target_pivot_entry =
          pmod_mat_entry(M, M_r, M_c, r2, r);

      int64_t reduced_pivot_entry =
          (int64_t)target_pivot_entry -
          (int64_t)((pivot_entry * factor) % MEDS_p);

      reduced_pivot_entry += MEDS_p * (reduced_pivot_entry < 0);

      pmod_mat_set_entry(
          M, M_r, M_c, r2, r, (GFq_t)reduced_pivot_entry);

      for (int c = M_r; c < M_c; c++)
      {
        uint64_t pivot_row_entry =
            pmod_mat_entry(M, M_r, M_c, r, c);
        uint64_t target_row_entry =
            pmod_mat_entry(M, M_r, M_c, r2, c);

        int64_t reduced_entry =
            (int64_t)target_row_entry -
            (int64_t)((pivot_row_entry * factor) % MEDS_p);

        reduced_entry += MEDS_p * (reduced_entry < 0);

        pmod_mat_set_entry(
            M, M_r, M_c, r2, c, (GFq_t)reduced_entry);
      }
    }
  }

  return 0;
}