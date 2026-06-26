#include <string.h>

#include "matrixelim.h"

typedef struct {
  int rank;
  int pivot_cols[MEDS_n];
} pmod_elim_info_t;

static void swap_rows(pmod_mat_t *A, int r1, int r2, int n)
{
  if (r1 == r2)
    return;

  for (int c = 0; c < n; c++)
  {
    GFq_t tmp = A[r1 * n + c];
    A[r1 * n + c] = A[r2 * n + c];
    A[r2 * n + c] = tmp;
  }
}

static void scale_row(pmod_mat_t *A, int row, GFq_t scale, int n)
{
  for (int c = 0; c < n; c++)
    A[row * n + c] = GF_mul(A[row * n + c], scale);
}

static void row_sub_mul(pmod_mat_t *A, int dst, int src, GFq_t factor, int n)
{
  for (int c = 0; c < n; c++)
  {
    GFq_t t = GF_mul(factor, A[src * n + c]);
    A[dst * n + c] = GF_sub(A[dst * n + c], t);
  }
}

static pmod_elim_info_t pmod_mat_rref_core_vartime(
    pmod_mat_t *A,
    pmod_mat_t *rhs,
    int n)
{
  pmod_elim_info_t info = {0};
  int pivot_row = 0;

  for (int col = 0; col < n && pivot_row < n; col++)
  {
    int pivot = -1;

    for (int r = pivot_row; r < n; r++)
    {
      if (A[r * n + col] != 0)
      {
        pivot = r;
        break;
      }
    }

    if (pivot < 0)
      continue;

    swap_rows(A, pivot_row, pivot, n);
    if (rhs != NULL)
      swap_rows(rhs, pivot_row, pivot, n);

    GFq_t inv = GF_inv(A[pivot_row * n + col]);
    scale_row(A, pivot_row, inv, n);
    if (rhs != NULL)
      scale_row(rhs, pivot_row, inv, n);

    for (int r = 0; r < n; r++)
    {
      if (r == pivot_row)
        continue;

      GFq_t factor = A[r * n + col];
      if (factor == 0)
        continue;

      row_sub_mul(A, r, pivot_row, factor, n);
      if (rhs != NULL)
        row_sub_mul(rhs, r, pivot_row, factor, n);
    }

    info.pivot_cols[pivot_row] = col;
    pivot_row++;
  }

  info.rank = pivot_row;
  return info;
}

int pmod_mat_rank_vartime(
    const pmod_mat_t *A,
    int n)
{
  pmod_mat_t tmp[n * n];
  pmod_mat_copy(tmp, A, n);

  pmod_elim_info_t info = pmod_mat_rref_core_vartime(tmp, NULL, n);

  return info.rank;
}

int pmod_mat_is_invertible_vartime(
    const pmod_mat_t *A,
    int n)
{
  return pmod_mat_rank_vartime(A, n) == n;
}

int pmod_mat_inv_vartime(
    pmod_mat_t *A_inv,
    const pmod_mat_t *A,
    int n)
{
  pmod_mat_t left[n * n];
  pmod_mat_t right[n * n];

  pmod_mat_copy(left, A, n);
  pmod_mat_identity(right, n);

  pmod_elim_info_t info = pmod_mat_rref_core_vartime(left, right, n);

  if (info.rank != n)
    return -1;

  if (A_inv != NULL)
    pmod_mat_copy(A_inv, right, n);

  return 0;
}

int pmod_mat_inv(Fq *B, const Fq *A, int A_r, int A_c)
{
  if (A_r != A_c)
    return -1;

  return pmod_mat_inv_vartime(B, A, A_r);
}

int pmod_mat_right_kernel_corank1_vartime(
    GFq_t *kernel,
    const pmod_mat_t *A,
    int n)
{
  pmod_mat_t R[n * n];
  pmod_mat_copy(R, A, n);

  pmod_elim_info_t info = pmod_mat_rref_core_vartime(R, NULL, n);

  if (info.rank != n - 1)
    return -1;

  int is_pivot[MEDS_n] = {0};

  for (int i = 0; i < info.rank; i++)
    is_pivot[info.pivot_cols[i]] = 1;

  int free_col = -1;

  for (int c = 0; c < n; c++)
  {
    if (!is_pivot[c])
    {
      free_col = c;
      break;
    }
  }

  if (free_col < 0)
    return -1;

  memset(kernel, 0, (size_t)n * sizeof(*kernel));
  kernel[free_col] = 1;

  for (int i = 0; i < info.rank; i++)
  {
    int pivot_col = info.pivot_cols[i];
    kernel[pivot_col] = GF_neg(R[i * n + free_col]);
  }

  return 0;
}

int pmod_mat_left_kernel_corank1_vartime(
    GFq_t *kernel,
    const pmod_mat_t *A,
    int n)
{
  pmod_mat_t AT[n * n];
  pmod_mat_transpose(AT, A, n);

  return pmod_mat_right_kernel_corank1_vartime(kernel, AT, n);
}

int pmod_vec_in_span_vartime(
    const GFq_t *v,
    const GFq_t *basis,
    int basis_count,
    int n)
{
  if (basis_count > n)
    return 1;

  pmod_mat_t B[n * n];
  pmod_mat_zero(B, n);

  for (int i = 0; i < basis_count; i++)
    pmod_mat_set_col(B, i, basis + i * n, n);

  int rank_before = pmod_mat_rank_vartime(B, n);

  if (basis_count >= n)
    return 1;

  pmod_mat_set_col(B, basis_count, v, n);

  int rank_after = pmod_mat_rank_vartime(B, n);

  return rank_after == rank_before;
}

int pmod_vec_extends_independent_set_vartime(
    const GFq_t *v,
    const GFq_t *basis,
    int basis_count,
    int n)
{
  if (basis_count >= n)
    return 0;

  pmod_mat_t B[n * n];
  pmod_mat_zero(B, n);

  for (int i = 0; i < basis_count; i++)
    pmod_mat_set_col(B, i, basis + i * n, n);

  pmod_mat_set_col(B, basis_count, v, n);

  return pmod_mat_rank_vartime(B, n) == basis_count + 1;
}
