#include "test_reference.h"

static void test_n1_zero_form_success(void)
{
  Fq U[1];
  Fq V[1];
  Fq W[1];
  Fq M[1] = {0};
  Fq u1[1] = {1};

  CHECK(canonical_build_uvw_vartime(U, V, W, M, u1, 1) == 0);
  CHECK(U[0] == 1);
  CHECK(V[0] != 0);
  CHECK(W[0] != 0);
  CHECK(pmod_mat_is_invertible_vartime(U, 1));
  CHECK(pmod_mat_is_invertible_vartime(V, 1));
  CHECK(pmod_mat_is_invertible_vartime(W, 1));

  u1[0] = 0;
  CHECK(canonical_build_uvw_vartime(U, V, W, M, u1, 1) == -1);
}

static void test_success_for_n(int n)
{
  Fq U[MEDS_n * MEDS_n];
  Fq V[MEDS_n * MEDS_n];
  Fq W[MEDS_n * MEDS_n];
  Fq ref_U[MEDS_n * MEDS_n];
  Fq ref_V[MEDS_n * MEDS_n];
  Fq ref_W[MEDS_n * MEDS_n];
  Fq M[MEDS_n * MEDS_n * MEDS_n];
  Fq u1[MEDS_n];

  require_success_fixture(M, u1, n, 0x22220000u + (uint32_t)n);

  CHECK(ref_build_uvw(ref_U, ref_V, ref_W, M, u1, n) == REF_BUILD_OK);
  CHECK(canonical_build_uvw_vartime(U, V, W, M, u1, n) == 0);
  verify_success_properties(M, u1, U, V, W, n);
  expect_projective_columns_equal(U, ref_U, n);
  expect_projective_columns_equal(V, ref_V, n);
  expect_projective_columns_equal(W, ref_W, n);
}

int main(void)
{
  int dims[] = {1, 2, 3, 4, 5, MEDS_n};

  test_n1_zero_form_success();

  for (size_t i = 0; i < sizeof(dims) / sizeof(dims[0]); i++)
    test_success_for_n(dims[i]);

  return 0;
}
