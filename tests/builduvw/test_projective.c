#include "test_reference.h"

static void test_repeatability_for_n(int n)
{
  Fq U0[MEDS_n * MEDS_n];
  Fq V0[MEDS_n * MEDS_n];
  Fq W0[MEDS_n * MEDS_n];
  Fq U1[MEDS_n * MEDS_n];
  Fq V1[MEDS_n * MEDS_n];
  Fq W1[MEDS_n * MEDS_n];
  Fq M[MEDS_n * MEDS_n * MEDS_n];
  Fq u1[MEDS_n];

  require_success_fixture(M, u1, n, 0x33330000u + (uint32_t)n);

  CHECK(canonical_build_uvw_vartime(U0, V0, W0, M, u1, n) == 0);
  CHECK(canonical_build_uvw_vartime(U1, V1, W1, M, u1, n) == 0);
  expect_mat_equal(U0, U1, n);
  expect_mat_equal(V0, V1, n);
  expect_mat_equal(W0, W1, n);
}

static void test_projective_scaling_for_n(int n)
{
  Fq U0[MEDS_n * MEDS_n];
  Fq V0[MEDS_n * MEDS_n];
  Fq W0[MEDS_n * MEDS_n];
  Fq U1[MEDS_n * MEDS_n];
  Fq V1[MEDS_n * MEDS_n];
  Fq W1[MEDS_n * MEDS_n];
  Fq M[MEDS_n * MEDS_n * MEDS_n];
  Fq u1[MEDS_n];
  Fq scaled_u1[MEDS_n];
  Fq first_col[MEDS_n];
  Fq alpha = 7;

  require_success_fixture(M, u1, n, 0x44440000u + (uint32_t)n);
  vec_scale(scaled_u1, alpha, u1, n);

  CHECK(canonical_build_uvw_vartime(U0, V0, W0, M, u1, n) == 0);
  CHECK(canonical_build_uvw_vartime(U1, V1, W1, M, scaled_u1, n) == 0);
  expect_projective_columns_equal(U0, U1, n);
  expect_projective_columns_equal(V0, V1, n);
  expect_projective_columns_equal(W0, W1, n);

  pmod_mat_get_col(first_col, U1, 0, n);
  expect_vec_equal(first_col, scaled_u1, n);
}

int main(void)
{
  int dims[] = {1, 2, 3, 5, MEDS_n};

  for (size_t i = 0; i < sizeof(dims) / sizeof(dims[0]); i++)
  {
    test_repeatability_for_n(dims[i]);
    test_projective_scaling_for_n(dims[i]);
  }

  return 0;
}
