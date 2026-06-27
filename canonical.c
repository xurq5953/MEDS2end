#include "canonical.h"

#include <string.h>

#include "matrixelim.h"
#include "matrixmod.h"
#include "triform.h"

static int vector_is_nonzero(
    const Fq *v,
    int n)
{
  for (int i = 0; i < n; i++)
    if (v[i] != 0)
      return 1;

  return 0;
}

static void vector_list_append(
    Fq *list,
    int index,
    const Fq *v,
    int n)
{
  memcpy(
      list + (size_t)index * (size_t)n,
      v,
      (size_t)n * sizeof(*v));
}

static void vector_list_to_matrix(
    Fq *A,
    const Fq *list,
    int n)
{
  pmod_mat_zero(A, n);

  for (int col = 0; col < n; col++)
    pmod_mat_set_col(
        A,
        col,
        list + (size_t)col * (size_t)n,
        n);
}

int canonical_build_uvw_vartime(
    Fq *U,
    Fq *V,
    Fq *W,
    const Fq *M,
    const Fq *u1,
    int n)
{
  if (U == NULL ||
      V == NULL ||
      W == NULL ||
      M == NULL ||
      u1 == NULL)
    return -1;

  if (n < 1 || n > MEDS_n)
    return -1;

  if (!vector_is_nonzero(u1, n))
    return -1;

  Fq LU[n * n];
  Fq LV[n * n];
  Fq LW[n * n];
  Fq current_v[n];
  Fq current_w[n];
  Fq next_u[n];
  Fq next_v[n];
  Fq matrix_tmp[n * n];

  int u_count = 1;
  int v_count = 0;
  int w_count = 0;

  vector_list_append(LU, 0, u1, n);

  if (triform_phi_u_lker_corank1_vartime(
          current_v,
          M,
          u1,
          n) != 0)
    return -1;

  for (int i = 0; i < n; i++)
  {
    if (!pmod_vec_extends_independent_set_vartime(
            current_v,
            LV,
            v_count,
            n))
      return -1;

    vector_list_append(LV, v_count, current_v, n);
    v_count++;

    if (triform_phi_v_rker_corank1_vartime(
            current_w,
            M,
            current_v,
            n) != 0)
      return -1;

    if (!pmod_vec_extends_independent_set_vartime(
            current_w,
            LW,
            w_count,
            n))
      return -1;

    vector_list_append(LW, w_count, current_w, n);
    w_count++;

    if (i == n - 1)
      break;

    if (triform_matrix_at_w_lker_corank1_vartime(
            next_u,
            M,
            current_w,
            n) != 0)
      return -1;

    /*
     * This both verifies Phi_U(next_u) has corank one and caches
     * v_{i+1}, avoiding a second elimination in the next iteration.
     */
    if (triform_phi_u_lker_corank1_vartime(
            next_v,
            M,
            next_u,
            n) != 0)
      return -1;

    if (!pmod_vec_extends_independent_set_vartime(
            next_u,
            LU,
            u_count,
            n))
      return -1;

    vector_list_append(LU, u_count, next_u, n);
    u_count++;

    memcpy(current_v, next_v, (size_t)n * sizeof(*current_v));
  }

  if (u_count != n || v_count != n || w_count != n)
    return -1;

  vector_list_to_matrix(matrix_tmp, LU, n);
  if (!pmod_mat_is_invertible_vartime(matrix_tmp, n))
    return -1;

  vector_list_to_matrix(matrix_tmp, LV, n);
  if (!pmod_mat_is_invertible_vartime(matrix_tmp, n))
    return -1;

  vector_list_to_matrix(matrix_tmp, LW, n);
  if (!pmod_mat_is_invertible_vartime(matrix_tmp, n))
    return -1;

  vector_list_to_matrix(U, LU, n);
  vector_list_to_matrix(V, LV, n);
  vector_list_to_matrix(W, LW, n);

  return 0;
}
