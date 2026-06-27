#include "test_common.h"

static Fq dot(const Fq *a, const Fq *b, int n)
{
  Fq acc = 0;

  for (int i = 0; i < n; i++)
    acc = GF_add(acc, GF_mul(a[i], b[i]));

  return acc;
}

static void random_derived_round(uint32_t seed, int n, int round, test_rng_t *rng)
{
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq u[MEDS_n] = {0};
  Fq v[MEDS_n] = {0};
  Fq w[MEDS_n] = {0};
  Fq other[MEDS_n] = {0};
  Fq sum[MEDS_n] = {0};
  Fq matrix_at_w[MEDS_n * MEDS_n] = {0};
  Fq phi_u[MEDS_n * MEDS_n] = {0};
  Fq phi_v[MEDS_n * MEDS_n] = {0};
  Fq ref_matrix[MEDS_n * MEDS_n] = {0};
  Fq ref_u[MEDS_n * MEDS_n] = {0};
  Fq ref_v[MEDS_n * MEDS_n] = {0};
  Fq rhs0[MEDS_n * MEDS_n] = {0};
  Fq rhs1[MEDS_n * MEDS_n] = {0};
  Fq rhs[MEDS_n * MEDS_n] = {0};
  Fq temp[MEDS_n] = {0};
  Fq eval_value = 0;

  fill_random(M, triform_element_count(n), rng);
  fill_random(u, (size_t)n, rng);
  fill_random(v, (size_t)n, rng);
  fill_random(w, (size_t)n, rng);

  triform_matrix_at_w(matrix_at_w, M, w, n);
  ref_matrix_at_w(ref_matrix, M, w, n);
  expect_mat_equal_ctx(matrix_at_w, ref_matrix, n, seed, round, "M(w)-reference");

  triform_phi_u(phi_u, M, u, n);
  ref_phi_u(ref_u, M, u, n);
  expect_mat_equal_ctx(phi_u, ref_u, n, seed, round, "Phi_U-reference");

  triform_phi_v(phi_v, M, v, n);
  ref_phi_v(ref_v, M, v, n);
  expect_mat_equal_ctx(phi_v, ref_v, n, seed, round, "Phi_V-reference");

  eval_value = triform_eval(M, u, v, w, n);
  CHECK_CTX(eval_value == ref_eval(M, u, v, w, n), seed, n, round, "eval-reference");
  CHECK_CTX(eval_value == ref_mat_form_eval(matrix_at_w, u, v, n), seed, n, round, "eval-M(w)");

  ref_mat_vec_mul(temp, phi_u, w, n);
  CHECK_CTX(eval_value == dot(v, temp, n), seed, n, round, "eval-Phi_U");

  ref_mat_vec_mul(temp, phi_v, w, n);
  CHECK_CTX(eval_value == dot(u, temp, n), seed, n, round, "eval-Phi_V");

  fill_random(other, (size_t)n, rng);
  vec_add(sum, w, other, n);
  triform_matrix_at_w(matrix_at_w, M, sum, n);
  triform_matrix_at_w(rhs0, M, w, n);
  triform_matrix_at_w(rhs1, M, other, n);
  mat_add(rhs, rhs0, rhs1, n);
  expect_mat_equal_ctx(matrix_at_w, rhs, n, seed, round, "M(w)-additivity");

  fill_random(other, (size_t)n, rng);
  vec_add(sum, u, other, n);
  triform_phi_u(phi_u, M, sum, n);
  triform_phi_u(rhs0, M, u, n);
  triform_phi_u(rhs1, M, other, n);
  mat_add(rhs, rhs0, rhs1, n);
  expect_mat_equal_ctx(phi_u, rhs, n, seed, round, "Phi_U-additivity");

  fill_random(other, (size_t)n, rng);
  vec_add(sum, v, other, n);
  triform_phi_v(phi_v, M, sum, n);
  triform_phi_v(rhs0, M, v, n);
  triform_phi_v(rhs1, M, other, n);
  mat_add(rhs, rhs0, rhs1, n);
  expect_mat_equal_ctx(phi_v, rhs, n, seed, round, "Phi_V-additivity");
}

static void run_random_derived(int n, int rounds, uint32_t group)
{
  uint32_t seed = 0x55667788u ^ ((uint32_t)n << 16) ^ group;
  test_rng_t rng;

  test_rng_init(&rng, seed);

  for (int round = 0; round < rounds; round++)
    random_derived_round(seed, n, round, &rng);
}

static void random_action_round(
    uint32_t seed,
    int n,
    int round,
    int check_six_loop_reference,
    test_rng_t *rng)
{
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

  fill_random(M, triform_element_count(n), rng);
  fill_random(A, triform_slice_size(n), rng);
  fill_random(B, triform_slice_size(n), rng);
  fill_random(C, triform_slice_size(n), rng);
  fill_random(x, (size_t)n, rng);
  fill_random(y, (size_t)n, rng);
  fill_random(z, (size_t)n, rng);

  triform_action_pullback(out, M, A, B, C, n);

  ref_mat_vec_mul(Ax, A, x, n);
  ref_mat_vec_mul(By, B, y, n);
  ref_mat_vec_mul(Cz, C, z, n);
  CHECK_CTX(
      triform_eval(out, x, y, z, n) == triform_eval(M, Ax, By, Cz, n),
      seed,
      n,
      round,
      "pullback-eval-identity");

  if (check_six_loop_reference)
  {
    ref_action_pullback(ref, M, A, B, C, n);
    expect_triform_equal_ctx(out, ref, n, seed, round, "pullback-six-loop");
  }
}

static void run_random_action(
    int n,
    int rounds,
    int check_six_loop_reference,
    uint32_t group)
{
  uint32_t seed = 0x55667788u ^ ((uint32_t)n << 16) ^ group;
  test_rng_t rng;

  test_rng_init(&rng, seed);

  for (int round = 0; round < rounds; round++)
    random_action_round(seed, n, round, check_six_loop_reference, &rng);
}

int main(void)
{
  int small_derived_dims[] = {1, 2, 3, 4, 5};
  int small_action_dims[] = {1, 2, 3, 4};

  for (size_t i = 0; i < sizeof(small_derived_dims) / sizeof(small_derived_dims[0]); i++)
    run_random_derived(small_derived_dims[i], 1000, 0x00000001u);

  run_random_derived(MEDS_n, 100, 0x00000002u);

  for (size_t i = 0; i < sizeof(small_action_dims) / sizeof(small_action_dims[0]); i++)
    run_random_action(small_action_dims[i], 100, 1, 0x00000003u);

  run_random_action(MEDS_n, 10, 0, 0x00000004u);

  return 0;
}
