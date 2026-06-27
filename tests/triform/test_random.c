#include "test_common.h"

static void test_random_for_n(int n)
{
  test_rng_t rng = {(uint32_t)(0x55667788u + (uint32_t)n)};
  Fq M[MEDS_n * MEDS_n * MEDS_n];
  Fq u[MEDS_n];
  Fq v[MEDS_n];
  Fq w[MEDS_n];
  Fq matrix_at_w[MEDS_n * MEDS_n];
  Fq phi_u[MEDS_n * MEDS_n];
  Fq phi_v[MEDS_n * MEDS_n];
  Fq ref_matrix[MEDS_n * MEDS_n];
  Fq ref_u[MEDS_n * MEDS_n];
  Fq ref_v[MEDS_n * MEDS_n];

  for (int round = 0; round < 20; round++)
  {
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
  }
}

int main(void)
{
  test_random_for_n(1);
  test_random_for_n(2);
  test_random_for_n(3);
  test_random_for_n(MEDS_n);
  return 0;
}
