#include "test_common.h"

static void test_action_for_n(int n)
{
  test_rng_t rng = {(uint32_t)(0xabcdef01u + (uint32_t)n)};
  Fq M[MEDS_n * MEDS_n * MEDS_n];
  Fq out[MEDS_n * MEDS_n * MEDS_n];
  Fq ref[MEDS_n * MEDS_n * MEDS_n];
  Fq A[MEDS_n * MEDS_n];
  Fq B[MEDS_n * MEDS_n];
  Fq C[MEDS_n * MEDS_n];
  Fq x[MEDS_n];
  Fq y[MEDS_n];
  Fq z[MEDS_n];
  Fq Ax[MEDS_n];
  Fq By[MEDS_n];
  Fq Cz[MEDS_n];

  fill_random(M, triform_element_count(n), &rng);
  fill_random(A, triform_slice_size(n), &rng);
  fill_random(B, triform_slice_size(n), &rng);
  fill_random(C, triform_slice_size(n), &rng);
  fill_random(x, (size_t)n, &rng);
  fill_random(y, (size_t)n, &rng);
  fill_random(z, (size_t)n, &rng);

  triform_action_pullback(out, M, A, B, C, n);
  ref_action_pullback(ref, M, A, B, C, n);
  expect_triform_equal(out, ref, n);

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
  test_rng_t rng = {(uint32_t)(0x10203040u + (uint32_t)n)};
  Fq M[MEDS_n * MEDS_n * MEDS_n];
  Fq first[MEDS_n * MEDS_n * MEDS_n];
  Fq second[MEDS_n * MEDS_n * MEDS_n];
  Fq direct[MEDS_n * MEDS_n * MEDS_n];
  Fq A[MEDS_n * MEDS_n];
  Fq B[MEDS_n * MEDS_n];
  Fq C[MEDS_n * MEDS_n];
  Fq D[MEDS_n * MEDS_n];
  Fq E[MEDS_n * MEDS_n];
  Fq F[MEDS_n * MEDS_n];
  Fq AD[MEDS_n * MEDS_n];
  Fq BE[MEDS_n * MEDS_n];
  Fq CF[MEDS_n * MEDS_n];

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

int main(void)
{
  int dims[] = {1, 2, 3, MEDS_n};

  for (size_t i = 0; i < sizeof(dims) / sizeof(dims[0]); i++)
  {
    test_action_for_n(dims[i]);
    test_composition_for_n(dims[i]);
  }

  return 0;
}
