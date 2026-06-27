#include "test_common.h"

static Fq dot(const Fq *a, const Fq *b, int n)
{
  Fq acc = 0;

  for (int i = 0; i < n; i++)
    acc = GF_add(acc, GF_mul(a[i], b[i]));

  return acc;
}

static void test_derived_for_n(int n)
{
  test_rng_t rng;
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
  Fq temp[MEDS_n] = {0};

  test_rng_init(&rng, (uint32_t)(0x12345678u + (uint32_t)n));
  fill_random(M, triform_element_count(n), &rng);
  fill_random(u, (size_t)n, &rng);
  fill_random(v, (size_t)n, &rng);
  fill_random(w, (size_t)n, &rng);

  triform_matrix_at_w(matrix_at_w, M, w, n);
  ref_matrix_at_w(ref_matrix, M, w, n);
  expect_mat_equal(matrix_at_w, ref_matrix, n);

  triform_phi_u(phi_u, M, u, n);
  ref_phi_u(ref_u, M, u, n);
  expect_mat_equal(phi_u, ref_u, n);

  triform_phi_v(phi_v, M, v, n);
  ref_phi_v(ref_v, M, v, n);
  expect_mat_equal(phi_v, ref_v, n);

  CHECK(triform_eval(M, u, v, w, n) == ref_eval(M, u, v, w, n));
  CHECK(triform_eval(M, u, v, w, n) == ref_mat_form_eval(matrix_at_w, u, v, n));

  ref_mat_vec_mul(temp, phi_u, w, n);
  CHECK(triform_eval(M, u, v, w, n) == dot(v, temp, n));

  ref_mat_vec_mul(temp, phi_v, w, n);
  CHECK(triform_eval(M, u, v, w, n) == dot(u, temp, n));
}

static void test_trilinearity_for_n(int n)
{
  test_rng_t rng;
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq u0[MEDS_n] = {0};
  Fq u1[MEDS_n] = {0};
  Fq v0[MEDS_n] = {0};
  Fq v1[MEDS_n] = {0};
  Fq w0[MEDS_n] = {0};
  Fq w1[MEDS_n] = {0};
  Fq sum[MEDS_n] = {0};
  Fq scaled[MEDS_n] = {0};
  Fq scalar = 7;

  test_rng_init(&rng, (uint32_t)(0x87654321u + (uint32_t)n));
  fill_random(M, triform_element_count(n), &rng);
  fill_random(u0, (size_t)n, &rng);
  fill_random(u1, (size_t)n, &rng);
  fill_random(v0, (size_t)n, &rng);
  fill_random(v1, (size_t)n, &rng);
  fill_random(w0, (size_t)n, &rng);
  fill_random(w1, (size_t)n, &rng);

  vec_add(sum, u0, u1, n);
  CHECK(triform_eval(M, sum, v0, w0, n) ==
        GF_add(triform_eval(M, u0, v0, w0, n), triform_eval(M, u1, v0, w0, n)));

  vec_add(sum, v0, v1, n);
  CHECK(triform_eval(M, u0, sum, w0, n) ==
        GF_add(triform_eval(M, u0, v0, w0, n), triform_eval(M, u0, v1, w0, n)));

  vec_add(sum, w0, w1, n);
  CHECK(triform_eval(M, u0, v0, sum, n) ==
        GF_add(triform_eval(M, u0, v0, w0, n), triform_eval(M, u0, v0, w1, n)));

  vec_scale(scaled, scalar, u0, n);
  CHECK(triform_eval(M, scaled, v0, w0, n) == GF_mul(scalar, triform_eval(M, u0, v0, w0, n)));

  vec_scale(scaled, scalar, v0, n);
  CHECK(triform_eval(M, u0, scaled, w0, n) == GF_mul(scalar, triform_eval(M, u0, v0, w0, n)));

  vec_scale(scaled, scalar, w0, n);
  CHECK(triform_eval(M, u0, v0, scaled, n) == GF_mul(scalar, triform_eval(M, u0, v0, w0, n)));
}

static void check_scaled_matrix_at_w(
    const Fq *M,
    const Fq *w,
    Fq scalar,
    int n)
{
  Fq scaled_w[MEDS_n] = {0};
  Fq lhs[MEDS_n * MEDS_n] = {0};
  Fq rhs0[MEDS_n * MEDS_n] = {0};
  Fq rhs[MEDS_n * MEDS_n] = {0};

  vec_scale(scaled_w, scalar, w, n);
  triform_matrix_at_w(lhs, M, scaled_w, n);
  triform_matrix_at_w(rhs0, M, w, n);
  mat_scale(rhs, scalar, rhs0, n);
  expect_mat_equal(lhs, rhs, n);
}

static void check_scaled_phi_u(
    const Fq *M,
    const Fq *u,
    Fq scalar,
    int n)
{
  Fq scaled_u[MEDS_n] = {0};
  Fq lhs[MEDS_n * MEDS_n] = {0};
  Fq rhs0[MEDS_n * MEDS_n] = {0};
  Fq rhs[MEDS_n * MEDS_n] = {0};

  vec_scale(scaled_u, scalar, u, n);
  triform_phi_u(lhs, M, scaled_u, n);
  triform_phi_u(rhs0, M, u, n);
  mat_scale(rhs, scalar, rhs0, n);
  expect_mat_equal(lhs, rhs, n);
}

static void check_scaled_phi_v(
    const Fq *M,
    const Fq *v,
    Fq scalar,
    int n)
{
  Fq scaled_v[MEDS_n] = {0};
  Fq lhs[MEDS_n * MEDS_n] = {0};
  Fq rhs0[MEDS_n * MEDS_n] = {0};
  Fq rhs[MEDS_n * MEDS_n] = {0};

  vec_scale(scaled_v, scalar, v, n);
  triform_phi_v(lhs, M, scaled_v, n);
  triform_phi_v(rhs0, M, v, n);
  mat_scale(rhs, scalar, rhs0, n);
  expect_mat_equal(lhs, rhs, n);
}

static void test_derived_linearity_for_n(int n)
{
  test_rng_t rng;
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq w0[MEDS_n] = {0};
  Fq w1[MEDS_n] = {0};
  Fq u0[MEDS_n] = {0};
  Fq u1[MEDS_n] = {0};
  Fq v0[MEDS_n] = {0};
  Fq v1[MEDS_n] = {0};
  Fq sum[MEDS_n] = {0};
  Fq lhs[MEDS_n * MEDS_n] = {0};
  Fq rhs0[MEDS_n * MEDS_n] = {0};
  Fq rhs1[MEDS_n * MEDS_n] = {0};
  Fq rhs[MEDS_n * MEDS_n] = {0};
  Fq scalars[] = {0, 1, 7, (Fq)(MEDS_p - 1)};

  test_rng_init(&rng, (uint32_t)(0x31415926u + (uint32_t)n));
  fill_random(M, triform_element_count(n), &rng);
  fill_random(w0, (size_t)n, &rng);
  fill_random(w1, (size_t)n, &rng);
  fill_random(u0, (size_t)n, &rng);
  fill_random(u1, (size_t)n, &rng);
  fill_random(v0, (size_t)n, &rng);
  fill_random(v1, (size_t)n, &rng);

  vec_add(sum, w0, w1, n);
  triform_matrix_at_w(lhs, M, sum, n);
  triform_matrix_at_w(rhs0, M, w0, n);
  triform_matrix_at_w(rhs1, M, w1, n);
  mat_add(rhs, rhs0, rhs1, n);
  expect_mat_equal(lhs, rhs, n);

  vec_add(sum, u0, u1, n);
  triform_phi_u(lhs, M, sum, n);
  triform_phi_u(rhs0, M, u0, n);
  triform_phi_u(rhs1, M, u1, n);
  mat_add(rhs, rhs0, rhs1, n);
  expect_mat_equal(lhs, rhs, n);

  vec_add(sum, v0, v1, n);
  triform_phi_v(lhs, M, sum, n);
  triform_phi_v(rhs0, M, v0, n);
  triform_phi_v(rhs1, M, v1, n);
  mat_add(rhs, rhs0, rhs1, n);
  expect_mat_equal(lhs, rhs, n);

  for (size_t i = 0; i < sizeof(scalars) / sizeof(scalars[0]); i++)
  {
    check_scaled_matrix_at_w(M, w0, scalars[i], n);
    check_scaled_phi_u(M, u0, scalars[i], n);
    check_scaled_phi_v(M, v0, scalars[i], n);
  }
}

int main(void)
{
  int dims[] = {1, 2, 3, 4, 5, MEDS_n};

  for (size_t i = 0; i < sizeof(dims) / sizeof(dims[0]); i++)
  {
    test_derived_for_n(dims[i]);
    test_trilinearity_for_n(dims[i]);
    test_derived_linearity_for_n(dims[i]);
  }

  return 0;
}
