#include "test_common.h"

static void test_action_for_n(int n)
{
  test_rng_t rng;
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq out[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq ref[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq A[MEDS_n * MEDS_n] = {0};
  Fq B[MEDS_n * MEDS_n] = {0};
  Fq C[MEDS_n * MEDS_n] = {0};
  Fq x[MEDS_n] = {0};
  Fq y[MEDS_n] = {0};
  Fq z[MEDS_n] = {0};
  Fq Ax[MEDS_n] = {0};
  Fq By[MEDS_n] = {0};
  Fq Cz[MEDS_n] = {0};

  test_rng_init(&rng, (uint32_t)(0xabcdef01u + (uint32_t)n));
  fill_random(M, triform_element_count(n), &rng);
  fill_random(A, triform_slice_size(n), &rng);
  fill_random(B, triform_slice_size(n), &rng);
  fill_random(C, triform_slice_size(n), &rng);
  fill_random(x, (size_t)n, &rng);
  fill_random(y, (size_t)n, &rng);
  fill_random(z, (size_t)n, &rng);

  triform_action_pullback(out, M, A, B, C, n);
  if (n <= 5)
  {
    ref_action_pullback(ref, M, A, B, C, n);
    expect_triform_equal(out, ref, n);
  }

  ref_mat_vec_mul(Ax, A, x, n);
  ref_mat_vec_mul(By, B, y, n);
  ref_mat_vec_mul(Cz, C, z, n);
  CHECK(triform_eval(out, x, y, z, n) == triform_eval(M, Ax, By, Cz, n));

  fill_identity(A, n);
  fill_identity(B, n);
  fill_identity(C, n);
  triform_action_pullback(out, M, A, B, C, n);
  expect_triform_equal(out, M, n);
}

static void test_composition_for_n(int n)
{
  test_rng_t rng;
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq first[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq second[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq direct[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq A[MEDS_n * MEDS_n] = {0};
  Fq B[MEDS_n * MEDS_n] = {0};
  Fq C[MEDS_n * MEDS_n] = {0};
  Fq D[MEDS_n * MEDS_n] = {0};
  Fq E[MEDS_n * MEDS_n] = {0};
  Fq F[MEDS_n * MEDS_n] = {0};
  Fq AD[MEDS_n * MEDS_n] = {0};
  Fq BE[MEDS_n * MEDS_n] = {0};
  Fq CF[MEDS_n * MEDS_n] = {0};

  test_rng_init(&rng, (uint32_t)(0x10203040u + (uint32_t)n));
  fill_random(M, triform_element_count(n), &rng);
  fill_random(A, triform_slice_size(n), &rng);
  fill_random(B, triform_slice_size(n), &rng);
  fill_random(C, triform_slice_size(n), &rng);
  fill_random(D, triform_slice_size(n), &rng);
  fill_random(E, triform_slice_size(n), &rng);
  fill_random(F, triform_slice_size(n), &rng);

  triform_action_pullback(first, M, A, B, C, n);
  triform_action_pullback(second, first, D, E, F, n);

  ref_mat_mul(AD, A, D, n);
  ref_mat_mul(BE, B, E, n);
  ref_mat_mul(CF, C, F, n);
  triform_action_pullback(direct, M, AD, BE, CF, n);

  expect_triform_equal(second, direct, n);
}

static void test_action_a_only_for_n(int n)
{
  test_rng_t rng;
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq out[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq expected[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq A[MEDS_n * MEDS_n] = {0};
  Fq A_transpose[MEDS_n * MEDS_n] = {0};
  Fq B[MEDS_n * MEDS_n] = {0};
  Fq C[MEDS_n * MEDS_n] = {0};

  test_rng_init(&rng, (uint32_t)(0x0a001000u + (uint32_t)n));
  fill_random(M, triform_element_count(n), &rng);
  fill_nonsymmetric_matrix(A, n);
  fill_identity(B, n);
  fill_identity(C, n);

  for (int row = 0; row < n; row++)
    for (int col = 0; col < n; col++)
      A_transpose[row * n + col] = A[col * n + row];

  triform_action_pullback(out, M, A, B, C, n);

  for (int slice = 0; slice < n; slice++)
    ref_mat_mul(
        expected + (size_t)slice * triform_slice_size(n),
        A_transpose,
        M + (size_t)slice * triform_slice_size(n),
        n);

  expect_triform_equal(out, expected, n);
}

static void test_action_b_only_for_n(int n)
{
  test_rng_t rng;
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq out[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq expected[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq A[MEDS_n * MEDS_n] = {0};
  Fq B[MEDS_n * MEDS_n] = {0};
  Fq C[MEDS_n * MEDS_n] = {0};

  test_rng_init(&rng, (uint32_t)(0x0b001000u + (uint32_t)n));
  fill_random(M, triform_element_count(n), &rng);
  fill_identity(A, n);
  fill_nonsymmetric_matrix(B, n);
  fill_identity(C, n);

  triform_action_pullback(out, M, A, B, C, n);

  for (int slice = 0; slice < n; slice++)
    ref_mat_mul(
        expected + (size_t)slice * triform_slice_size(n),
        M + (size_t)slice * triform_slice_size(n),
        B,
        n);

  expect_triform_equal(out, expected, n);
}

static void test_action_c_only_for_n(int n)
{
  test_rng_t rng;
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq out[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq expected[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq coefficients[MEDS_n] = {0};
  Fq A[MEDS_n * MEDS_n] = {0};
  Fq B[MEDS_n * MEDS_n] = {0};
  Fq C[MEDS_n * MEDS_n] = {0};

  test_rng_init(&rng, (uint32_t)(0x0c001000u + (uint32_t)n));
  fill_random(M, triform_element_count(n), &rng);
  fill_identity(A, n);
  fill_identity(B, n);
  fill_nonsymmetric_matrix(C, n);

  if (n >= 2)
    CHECK(C[1] != C[n]);

  triform_action_pullback(out, M, A, B, C, n);

  for (int out_slice = 0; out_slice < n; out_slice++)
  {
    for (int in_slice = 0; in_slice < n; in_slice++)
      coefficients[in_slice] = C[in_slice * n + out_slice];

    ref_matrix_at_w(
        expected + (size_t)out_slice * triform_slice_size(n),
        M,
        coefficients,
        n);
  }

  expect_triform_equal(out, expected, n);
}

static void test_action_c_column_selection_for_n(int n)
{
  test_rng_t rng;
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq out[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq expected[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq A[MEDS_n * MEDS_n] = {0};
  Fq B[MEDS_n * MEDS_n] = {0};
  Fq C[MEDS_n * MEDS_n] = {0};
  int in_slice = n - 1;
  int out_slice = (n > 1) ? 1 : 0;

  test_rng_init(&rng, (uint32_t)(0xc0010000u + (uint32_t)n));
  fill_random(M, triform_element_count(n), &rng);
  fill_identity(A, n);
  fill_identity(B, n);
  C[in_slice * n + out_slice] = 1;

  triform_action_pullback(out, M, A, B, C, n);
  memcpy(
      expected + (size_t)out_slice * triform_slice_size(n),
      M + (size_t)in_slice * triform_slice_size(n),
      triform_slice_size(n) * sizeof(*expected));

  expect_triform_equal(out, expected, n);
}

static void test_action_asymmetric_ac_for_n(int n)
{
  test_rng_t rng;
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq out[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq expected[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq A[MEDS_n * MEDS_n] = {0};
  Fq B[MEDS_n * MEDS_n] = {0};
  Fq C[MEDS_n * MEDS_n] = {0};

  test_rng_init(&rng, (uint32_t)(0xac001000u + (uint32_t)n));
  fill_random(M, triform_element_count(n), &rng);
  fill_nonsymmetric_matrix(A, n);
  fill_identity(B, n);
  fill_nonsymmetric_matrix(C, n);

  triform_action_pullback(out, M, A, B, C, n);
  ref_action_pullback(expected, M, A, B, C, n);
  expect_triform_equal(out, expected, n);
}

int main(void)
{
  int dims[] = {2, 3, 5, MEDS_n};

  for (size_t i = 0; i < sizeof(dims) / sizeof(dims[0]); i++)
  {
    test_action_for_n(dims[i]);
    test_composition_for_n(dims[i]);
    test_action_a_only_for_n(dims[i]);
    test_action_b_only_for_n(dims[i]);
    test_action_c_only_for_n(dims[i]);
    test_action_c_column_selection_for_n(dims[i]);
  }

  test_action_asymmetric_ac_for_n(2);
  test_action_asymmetric_ac_for_n(3);
  test_action_asymmetric_ac_for_n(5);

  return 0;
}
