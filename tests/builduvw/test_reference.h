#ifndef TEST_BUILDUVW_REFERENCE_H
#define TEST_BUILDUVW_REFERENCE_H

#include "test_common.h"

typedef enum {
  REF_BUILD_OK = 0,
  REF_BUILD_BAD_INPUT,
  REF_BUILD_ZERO_U1,
  REF_BUILD_PHI_U_NOT_CORANK1,
  REF_BUILD_V_DEPENDENT,
  REF_BUILD_PHI_V_NOT_CORANK1,
  REF_BUILD_W_DEPENDENT,
  REF_BUILD_MW_NOT_CORANK1,
  REF_BUILD_NEXT_PHI_U_NOT_CORANK1,
  REF_BUILD_U_DEPENDENT,
  REF_BUILD_FINAL_SINGULAR
} ref_build_status_t;

static inline const char *ref_build_status_name(ref_build_status_t status)
{
  switch (status)
  {
    case REF_BUILD_OK:
      return "OK";
    case REF_BUILD_BAD_INPUT:
      return "BAD_INPUT";
    case REF_BUILD_ZERO_U1:
      return "ZERO_U1";
    case REF_BUILD_PHI_U_NOT_CORANK1:
      return "PHI_U_NOT_CORANK1";
    case REF_BUILD_V_DEPENDENT:
      return "V_DEPENDENT";
    case REF_BUILD_PHI_V_NOT_CORANK1:
      return "PHI_V_NOT_CORANK1";
    case REF_BUILD_W_DEPENDENT:
      return "W_DEPENDENT";
    case REF_BUILD_MW_NOT_CORANK1:
      return "MW_NOT_CORANK1";
    case REF_BUILD_NEXT_PHI_U_NOT_CORANK1:
      return "NEXT_PHI_U_NOT_CORANK1";
    case REF_BUILD_U_DEPENDENT:
      return "U_DEPENDENT";
    case REF_BUILD_FINAL_SINGULAR:
      return "FINAL_SINGULAR";
  }

  return "UNKNOWN";
}

static inline void ref_vector_list_append(
    Fq *list,
    int index,
    const Fq *v,
    int n)
{
  memcpy(list + (size_t)index * (size_t)n, v, (size_t)n * sizeof(*v));
}

static inline void ref_vector_list_to_matrix(Fq *A, const Fq *list, int n)
{
  pmod_mat_zero(A, n);

  for (int col = 0; col < n; col++)
    pmod_mat_set_col(A, col, list + (size_t)col * (size_t)n, n);
}

static inline int ref_extends_independent_set(
    const Fq *v,
    const Fq *packed_basis,
    int count,
    int n)
{
  Fq B[MEDS_n * MEDS_n];

  if (count >= n)
    return 0;

  pmod_mat_zero(B, n);

  for (int i = 0; i < count; i++)
    pmod_mat_set_col(B, i, packed_basis + (size_t)i * (size_t)n, n);

  pmod_mat_set_col(B, count, v, n);

  return pmod_mat_rank_vartime(B, n) == count + 1;
}

static inline ref_build_status_t ref_phi_u_lker(
    Fq *kernel,
    const Fq *M,
    const Fq *u,
    int n,
    ref_build_status_t failure_status)
{
  Fq phi_u[MEDS_n * MEDS_n];

  triform_phi_u(phi_u, M, u, n);

  if (pmod_mat_rank_vartime(phi_u, n) != n - 1)
    return failure_status;

  if (pmod_mat_left_kernel_corank1_vartime(kernel, phi_u, n) != 0)
    return failure_status;

  return REF_BUILD_OK;
}

static inline ref_build_status_t ref_build_uvw(
    Fq *U,
    Fq *V,
    Fq *W,
    const Fq *M,
    const Fq *u1,
    int n)
{
  if (U == NULL || V == NULL || W == NULL || M == NULL || u1 == NULL)
    return REF_BUILD_BAD_INPUT;

  if (n < 1 || n > MEDS_n)
    return REF_BUILD_BAD_INPUT;

  if (!vec_nonzero(u1, n))
    return REF_BUILD_ZERO_U1;

  Fq LU[MEDS_n * MEDS_n];
  Fq LV[MEDS_n * MEDS_n];
  Fq LW[MEDS_n * MEDS_n];
  Fq current_v[MEDS_n];
  Fq current_w[MEDS_n];
  Fq next_u[MEDS_n];
  Fq next_v[MEDS_n];
  Fq phi_v[MEDS_n * MEDS_n];
  Fq matrix_at_w[MEDS_n * MEDS_n];
  Fq matrix_tmp[MEDS_n * MEDS_n];
  ref_build_status_t status;
  int u_count = 1;
  int v_count = 0;
  int w_count = 0;

  ref_vector_list_append(LU, 0, u1, n);

  status = ref_phi_u_lker(
      current_v,
      M,
      u1,
      n,
      REF_BUILD_PHI_U_NOT_CORANK1);

  if (status != REF_BUILD_OK)
    return status;

  for (int i = 0; i < n; i++)
  {
    if (!ref_extends_independent_set(current_v, LV, v_count, n))
      return REF_BUILD_V_DEPENDENT;

    ref_vector_list_append(LV, v_count, current_v, n);
    v_count++;

    triform_phi_v(phi_v, M, current_v, n);
    if (pmod_mat_rank_vartime(phi_v, n) != n - 1)
      return REF_BUILD_PHI_V_NOT_CORANK1;

    if (pmod_mat_right_kernel_corank1_vartime(current_w, phi_v, n) != 0)
      return REF_BUILD_PHI_V_NOT_CORANK1;

    if (!ref_extends_independent_set(current_w, LW, w_count, n))
      return REF_BUILD_W_DEPENDENT;

    ref_vector_list_append(LW, w_count, current_w, n);
    w_count++;

    if (i == n - 1)
      break;

    triform_matrix_at_w(matrix_at_w, M, current_w, n);
    if (pmod_mat_rank_vartime(matrix_at_w, n) != n - 1)
      return REF_BUILD_MW_NOT_CORANK1;

    if (pmod_mat_left_kernel_corank1_vartime(next_u, matrix_at_w, n) != 0)
      return REF_BUILD_MW_NOT_CORANK1;

    status = ref_phi_u_lker(
        next_v,
        M,
        next_u,
        n,
        REF_BUILD_NEXT_PHI_U_NOT_CORANK1);

    if (status != REF_BUILD_OK)
      return status;

    if (!ref_extends_independent_set(next_u, LU, u_count, n))
      return REF_BUILD_U_DEPENDENT;

    ref_vector_list_append(LU, u_count, next_u, n);
    u_count++;
    memcpy(current_v, next_v, (size_t)n * sizeof(*current_v));
  }

  if (u_count != n || v_count != n || w_count != n)
    return REF_BUILD_FINAL_SINGULAR;

  ref_vector_list_to_matrix(matrix_tmp, LU, n);
  if (!pmod_mat_is_invertible_vartime(matrix_tmp, n))
    return REF_BUILD_FINAL_SINGULAR;

  ref_vector_list_to_matrix(matrix_tmp, LV, n);
  if (!pmod_mat_is_invertible_vartime(matrix_tmp, n))
    return REF_BUILD_FINAL_SINGULAR;

  ref_vector_list_to_matrix(matrix_tmp, LW, n);
  if (!pmod_mat_is_invertible_vartime(matrix_tmp, n))
    return REF_BUILD_FINAL_SINGULAR;

  ref_vector_list_to_matrix(U, LU, n);
  ref_vector_list_to_matrix(V, LV, n);
  ref_vector_list_to_matrix(W, LW, n);

  return REF_BUILD_OK;
}

static inline void require_success_fixture(
    Fq *M,
    Fq *u1,
    int n,
    uint32_t base_seed)
{
  Fq U[MEDS_n * MEDS_n];
  Fq V[MEDS_n * MEDS_n];
  Fq W[MEDS_n * MEDS_n];
  ref_build_status_t last_status = REF_BUILD_BAD_INPUT;

  if (n == 1)
  {
    fill_zero(M, form_count(n));
    u1[0] = 1;
    return;
  }

  for (uint32_t attempt = 0; attempt < 1024u; attempt++)
  {
    uint32_t seed = base_seed + attempt;

    make_coordinate_path_candidate(M, u1, n, seed);
    last_status = ref_build_uvw(U, V, W, M, u1, n);

    if (last_status == REF_BUILD_OK)
      return;
  }

  fprintf(
      stderr,
      "failed to find success fixture: n=%d seed=0x%08x last=%s\n",
      n,
      (unsigned)base_seed,
      ref_build_status_name(last_status));
  exit(1);
}

static inline void verify_success_properties(
    const Fq *M,
    const Fq *u1,
    const Fq *U,
    const Fq *V,
    const Fq *W,
    int n)
{
  Fq u_col[MEDS_n];
  Fq v_col[MEDS_n];
  Fq w_col[MEDS_n];
  Fq next_u[MEDS_n];
  Fq phi_u[MEDS_n * MEDS_n];
  Fq phi_v[MEDS_n * MEDS_n];
  Fq matrix_at_w[MEDS_n * MEDS_n];
  Fq left_value[MEDS_n];
  Fq right_value[MEDS_n];

  pmod_mat_get_col(u_col, U, 0, n);
  expect_vec_equal(u_col, u1, n);

  CHECK(pmod_mat_is_invertible_vartime(U, n));
  CHECK(pmod_mat_is_invertible_vartime(V, n));
  CHECK(pmod_mat_is_invertible_vartime(W, n));

  for (int i = 0; i < n; i++)
  {
    pmod_mat_get_col(u_col, U, i, n);
    pmod_mat_get_col(v_col, V, i, n);
    pmod_mat_get_col(w_col, W, i, n);

    CHECK(vec_nonzero(u_col, n));
    CHECK(vec_nonzero(v_col, n));
    CHECK(vec_nonzero(w_col, n));

    triform_phi_u(phi_u, M, u_col, n);
    CHECK(pmod_mat_rank_vartime(phi_u, n) == n - 1);
    pmod_mat_transpose_vec_mul(left_value, phi_u, v_col, n);
    for (int j = 0; j < n; j++)
      CHECK(left_value[j] == 0);

    triform_phi_v(phi_v, M, v_col, n);
    CHECK(pmod_mat_rank_vartime(phi_v, n) == n - 1);
    pmod_mat_vec_mul(right_value, phi_v, w_col, n);
    for (int j = 0; j < n; j++)
      CHECK(right_value[j] == 0);

    if (i < n - 1)
    {
      pmod_mat_get_col(next_u, U, i + 1, n);
      triform_matrix_at_w(matrix_at_w, M, w_col, n);
      CHECK(pmod_mat_rank_vartime(matrix_at_w, n) == n - 1);
      pmod_mat_transpose_vec_mul(left_value, matrix_at_w, next_u, n);
      for (int j = 0; j < n; j++)
        CHECK(left_value[j] == 0);
    }
  }
}

static inline void compare_reference_and_production(
    const Fq *M,
    const Fq *u1,
    int n)
{
  Fq U[MEDS_n * MEDS_n];
  Fq V[MEDS_n * MEDS_n];
  Fq W[MEDS_n * MEDS_n];
  Fq ref_U[MEDS_n * MEDS_n];
  Fq ref_V[MEDS_n * MEDS_n];
  Fq ref_W[MEDS_n * MEDS_n];
  ref_build_status_t ref_status = ref_build_uvw(ref_U, ref_V, ref_W, M, u1, n);
  int prod_status = canonical_build_uvw_vartime(U, V, W, M, u1, n);

  CHECK((ref_status == REF_BUILD_OK) == (prod_status == 0));

  if (ref_status == REF_BUILD_OK)
  {
    expect_projective_columns_equal(U, ref_U, n);
    expect_projective_columns_equal(V, ref_V, n);
    expect_projective_columns_equal(W, ref_W, n);
    pmod_mat_get_col(ref_U, U, 0, n);
    expect_vec_equal(ref_U, u1, n);
  }
}

#endif
