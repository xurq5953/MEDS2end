#include "test_common.h"

static void make_corank_one_form(Fq *M, int n)
{
  fill_zero(M, triform_element_count(n));

  for (int i = 0; i < n - 1; i++)
  {
    M[triform_index(0, i, i, n)] = 1;
    M[triform_index(i, 0, i, n)] = 1;
    M[triform_index(i, i, 0, n)] = 1;
  }
}

static void make_full_rank_phi_u_form(Fq *M, int n)
{
  fill_zero(M, triform_element_count(n));

  for (int i = 0; i < n; i++)
    M[triform_index(i, 0, i, n)] = 1;
}

static void make_full_rank_phi_v_form(Fq *M, int n)
{
  fill_zero(M, triform_element_count(n));

  for (int i = 0; i < n; i++)
    M[triform_index(i, i, 0, n)] = 1;
}

static void make_full_rank_matrix_w_form(Fq *M, int n)
{
  fill_zero(M, triform_element_count(n));

  for (int i = 0; i < n; i++)
    M[triform_index(0, i, i, n)] = 1;
}

static void make_corank_two_form(Fq *M, int n)
{
  fill_zero(M, triform_element_count(n));

  for (int i = 0; i < n - 2; i++)
  {
    M[triform_index(0, i, i, n)] = 1;
    M[triform_index(i, 0, i, n)] = 1;
    M[triform_index(i, i, 0, n)] = 1;
  }
}

static void expect_left_kernel(const Fq *kernel, const Fq *A, int n)
{
  Fq value[MEDS_n];
  ref_mat_transpose_vec_mul(value, A, kernel, n);

  for (int i = 0; i < n; i++)
    CHECK(value[i] == 0);
}

static void expect_right_kernel(const Fq *kernel, const Fq *A, int n)
{
  Fq value[MEDS_n];
  ref_mat_vec_mul(value, A, kernel, n);

  for (int i = 0; i < n; i++)
    CHECK(value[i] == 0);
}

static void test_corank_for_n(int n)
{
  Fq M[MEDS_n * MEDS_n * MEDS_n];
  Fq point[MEDS_n];
  Fq kernel[MEDS_n];
  Fq derived[MEDS_n * MEDS_n];

  fill_zero(point, (size_t)n);
  point[0] = 1;

  make_corank_one_form(M, n);

  CHECK(triform_phi_u_lker_corank1_vartime(kernel, M, point, n) == 0);
  CHECK(vec_is_nonzero(kernel, n));
  triform_phi_u(derived, M, point, n);
  expect_left_kernel(kernel, derived, n);

  CHECK(triform_phi_v_rker_corank1_vartime(kernel, M, point, n) == 0);
  CHECK(vec_is_nonzero(kernel, n));
  triform_phi_v(derived, M, point, n);
  expect_right_kernel(kernel, derived, n);

  CHECK(triform_matrix_at_w_lker_corank1_vartime(kernel, M, point, n) == 0);
  CHECK(vec_is_nonzero(kernel, n));
  triform_matrix_at_w(derived, M, point, n);
  expect_left_kernel(kernel, derived, n);

  make_full_rank_phi_u_form(M, n);
  CHECK(triform_phi_u_lker_corank1_vartime(kernel, M, point, n) == -1);

  make_full_rank_phi_v_form(M, n);
  CHECK(triform_phi_v_rker_corank1_vartime(kernel, M, point, n) == -1);

  make_full_rank_matrix_w_form(M, n);
  CHECK(triform_matrix_at_w_lker_corank1_vartime(kernel, M, point, n) == -1);

  if (n > 1)
  {
    fill_zero(M, triform_element_count(n));
    CHECK(triform_phi_u_lker_corank1_vartime(kernel, M, point, n) == -1);
    CHECK(triform_phi_v_rker_corank1_vartime(kernel, M, point, n) == -1);
    CHECK(triform_matrix_at_w_lker_corank1_vartime(kernel, M, point, n) == -1);
  }

  if (n > 2)
  {
    make_corank_two_form(M, n);
    CHECK(triform_phi_u_lker_corank1_vartime(kernel, M, point, n) == -1);
    CHECK(triform_phi_v_rker_corank1_vartime(kernel, M, point, n) == -1);
    CHECK(triform_matrix_at_w_lker_corank1_vartime(kernel, M, point, n) == -1);
  }
}

int main(void)
{
  test_corank_for_n(1);
  test_corank_for_n(2);
  test_corank_for_n(3);
  test_corank_for_n(MEDS_n);
  return 0;
}
