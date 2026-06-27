#include "test_common.h"

static void test_zero_form_for_n(int n)
{
  test_rng_t rng;
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq out_form[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq A[MEDS_n * MEDS_n] = {0};
  Fq B[MEDS_n * MEDS_n] = {0};
  Fq C[MEDS_n * MEDS_n] = {0};
  Fq u[MEDS_n] = {0};
  Fq v[MEDS_n] = {0};
  Fq w[MEDS_n] = {0};
  Fq matrix_at_w[MEDS_n * MEDS_n] = {0};
  Fq phi_u[MEDS_n * MEDS_n] = {0};
  Fq phi_v[MEDS_n * MEDS_n] = {0};

  test_rng_init(&rng, (uint32_t)(0xed000000u + (uint32_t)n));
  fill_random(A, triform_slice_size(n), &rng);
  fill_random(B, triform_slice_size(n), &rng);
  fill_random(C, triform_slice_size(n), &rng);
  fill_random(u, (size_t)n, &rng);
  fill_random(v, (size_t)n, &rng);
  fill_random(w, (size_t)n, &rng);

  triform_matrix_at_w(matrix_at_w, M, w, n);
  triform_phi_u(phi_u, M, u, n);
  triform_phi_v(phi_v, M, v, n);
  triform_action_pullback(out_form, M, A, B, C, n);

  expect_zero_mat(matrix_at_w, n);
  expect_zero_mat(phi_u, n);
  expect_zero_mat(phi_v, n);
  expect_zero_triform(out_form, n);
  CHECK(triform_eval(M, u, v, w, n) == 0);
}

static void test_zero_vectors_for_n(int n)
{
  test_rng_t rng;
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq u[MEDS_n] = {0};
  Fq v[MEDS_n] = {0};
  Fq w[MEDS_n] = {0};
  Fq zero[MEDS_n] = {0};
  Fq matrix_at_w[MEDS_n * MEDS_n] = {0};
  Fq phi_u[MEDS_n * MEDS_n] = {0};
  Fq phi_v[MEDS_n * MEDS_n] = {0};

  test_rng_init(&rng, (uint32_t)(0xed100000u + (uint32_t)n));
  fill_random(M, triform_element_count(n), &rng);
  fill_random(u, (size_t)n, &rng);
  fill_random(v, (size_t)n, &rng);
  fill_random(w, (size_t)n, &rng);

  triform_phi_u(phi_u, M, zero, n);
  expect_zero_mat(phi_u, n);
  CHECK(triform_eval(M, zero, v, w, n) == 0);

  triform_phi_v(phi_v, M, zero, n);
  expect_zero_mat(phi_v, n);
  CHECK(triform_eval(M, u, zero, w, n) == 0);

  triform_matrix_at_w(matrix_at_w, M, zero, n);
  expect_zero_mat(matrix_at_w, n);
  CHECK(triform_eval(M, u, v, zero, n) == 0);
}

static void test_single_nonzero_slice_for_n(int n)
{
  test_rng_t rng;
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq w[MEDS_n] = {0};
  Fq matrix_at_w[MEDS_n * MEDS_n] = {0};
  int selected = n - 1;

  test_rng_init(&rng, (uint32_t)(0xed200000u + (uint32_t)n));
  fill_random(
      M + (size_t)selected * triform_slice_size(n),
      triform_slice_size(n),
      &rng);

  fill_basis_vector(w, selected, n);
  triform_matrix_at_w(matrix_at_w, M, w, n);
  expect_mat_equal(
      matrix_at_w,
      M + (size_t)selected * triform_slice_size(n),
      n);

  for (int slice = 0; slice < n; slice++)
  {
    if (slice == selected)
      continue;

    fill_basis_vector(w, slice, n);
    triform_matrix_at_w(matrix_at_w, M, w, n);
    expect_zero_mat(matrix_at_w, n);
  }
}

static void test_single_nonzero_tensor_element_for_n(int n)
{
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq u[MEDS_n] = {0};
  Fq v[MEDS_n] = {0};
  Fq w[MEDS_n] = {0};
  Fq value = (Fq)(MEDS_p - 2);
  int slice = n - 1;
  int row = (n > 1) ? 1 : 0;
  int col = (n > 2) ? 2 : 0;

  M[ref_triform_index(slice, row, col, n)] = value;

  fill_basis_vector(u, row, n);
  fill_basis_vector(v, col, n);
  fill_basis_vector(w, slice, n);
  CHECK(triform_eval(M, u, v, w, n) == value);

  if (n > 1)
  {
    fill_basis_vector(u, (row + 1) % n, n);
    CHECK(triform_eval(M, u, v, w, n) == 0);
  }

  fill_basis_vector(u, row, n);
  if (n > 1)
  {
    fill_basis_vector(v, (col + 1) % n, n);
    CHECK(triform_eval(M, u, v, w, n) == 0);
  }

  fill_basis_vector(v, col, n);
  if (n > 1)
  {
    fill_basis_vector(w, (slice + 1) % n, n);
    CHECK(triform_eval(M, u, v, w, n) == 0);
  }
}

static void test_max_canonical_for_n(int n)
{
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq u[MEDS_n] = {0};
  Fq v[MEDS_n] = {0};
  Fq w[MEDS_n] = {0};
  Fq matrix_at_w[MEDS_n * MEDS_n] = {0};
  Fq phi_u[MEDS_n * MEDS_n] = {0};
  Fq phi_v[MEDS_n * MEDS_n] = {0};
  Fq ref_matrix[MEDS_n * MEDS_n] = {0};
  Fq ref_u[MEDS_n * MEDS_n] = {0};
  Fq ref_v[MEDS_n * MEDS_n] = {0};

  fill_constant(M, triform_element_count(n), (Fq)(MEDS_p - 1));
  fill_constant(u, (size_t)n, (Fq)(MEDS_p - 1));
  fill_constant(v, (size_t)n, (Fq)(MEDS_p - 1));
  fill_constant(w, (size_t)n, (Fq)(MEDS_p - 1));

  triform_matrix_at_w(matrix_at_w, M, w, n);
  triform_phi_u(phi_u, M, u, n);
  triform_phi_v(phi_v, M, v, n);
  ref_matrix_at_w(ref_matrix, M, w, n);
  ref_phi_u(ref_u, M, u, n);
  ref_phi_v(ref_v, M, v, n);

  expect_mat_equal(matrix_at_w, ref_matrix, n);
  expect_mat_equal(phi_u, ref_u, n);
  expect_mat_equal(phi_v, ref_v, n);
  CHECK(triform_eval(M, u, v, w, n) == ref_eval(M, u, v, w, n));

  for (int i = 0; i < n * n; i++)
  {
    CHECK(matrix_at_w[i] < MEDS_p);
    CHECK(phi_u[i] < MEDS_p);
    CHECK(phi_v[i] < MEDS_p);
  }
}

static void test_all_basis_vectors_for_n(int n)
{
  test_rng_t rng;
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq u[MEDS_n] = {0};
  Fq v[MEDS_n] = {0};
  Fq w[MEDS_n] = {0};

  test_rng_init(&rng, (uint32_t)(0xed300000u + (uint32_t)n));
  fill_random(M, triform_element_count(n), &rng);

  for (int slice = 0; slice < n; slice++)
    for (int row = 0; row < n; row++)
      for (int col = 0; col < n; col++)
      {
        fill_basis_vector(u, row, n);
        fill_basis_vector(v, col, n);
        fill_basis_vector(w, slice, n);
        CHECK(triform_eval(M, u, v, w, n) ==
              M[ref_triform_index(slice, row, col, n)]);
      }
}

int main(void)
{
  int dims[] = {1, 2, 3, 4, 5, MEDS_n};
  int basis_dims[] = {2, 3, 4};

  for (size_t i = 0; i < sizeof(dims) / sizeof(dims[0]); i++)
  {
    test_zero_form_for_n(dims[i]);
    test_zero_vectors_for_n(dims[i]);
    test_single_nonzero_slice_for_n(dims[i]);
    test_single_nonzero_tensor_element_for_n(dims[i]);
    test_max_canonical_for_n(dims[i]);
  }

  for (size_t i = 0; i < sizeof(basis_dims) / sizeof(basis_dims[0]); i++)
    test_all_basis_vectors_for_n(basis_dims[i]);

  return 0;
}
