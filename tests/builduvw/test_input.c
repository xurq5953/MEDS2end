#include "test_reference.h"

static void test_bad_dimensions(void)
{
  Fq U[MEDS_n * MEDS_n] = {0};
  Fq V[MEDS_n * MEDS_n] = {0};
  Fq W[MEDS_n * MEDS_n] = {0};
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq u1[MEDS_n] = {0};

  CHECK(canonical_build_uvw_vartime(U, V, W, M, u1, 0) == -1);
  CHECK(canonical_build_uvw_vartime(U, V, W, M, u1, -1) == -1);
  CHECK(canonical_build_uvw_vartime(U, V, W, M, u1, MEDS_n + 1) == -1);
}

static void test_null_pointers(void)
{
  Fq U[MEDS_n * MEDS_n] = {0};
  Fq V[MEDS_n * MEDS_n] = {0};
  Fq W[MEDS_n * MEDS_n] = {0};
  Fq M[MEDS_n * MEDS_n * MEDS_n] = {0};
  Fq u1[MEDS_n] = {0};

  CHECK(canonical_build_uvw_vartime(NULL, V, W, M, u1, 1) == -1);
  CHECK(canonical_build_uvw_vartime(U, NULL, W, M, u1, 1) == -1);
  CHECK(canonical_build_uvw_vartime(U, V, NULL, M, u1, 1) == -1);
  CHECK(canonical_build_uvw_vartime(U, V, W, NULL, u1, 1) == -1);
  CHECK(canonical_build_uvw_vartime(U, V, W, M, NULL, 1) == -1);
}

static void test_zero_representative_preserves_outputs(void)
{
  Fq U[MEDS_n * MEDS_n];
  Fq V[MEDS_n * MEDS_n];
  Fq W[MEDS_n * MEDS_n];
  Fq M[MEDS_n * MEDS_n * MEDS_n];
  Fq u1[MEDS_n];
  int n = 3;

  fill_constant(U, mat_count(n), 123);
  fill_constant(V, mat_count(n), 123);
  fill_constant(W, mat_count(n), 123);
  fill_random(M, form_count(n), &(test_rng_t){1});
  fill_zero(u1, (size_t)n);

  CHECK(canonical_build_uvw_vartime(U, V, W, M, u1, n) == -1);
  CHECK(outputs_unchanged(U, V, W, 123, n));
}

static void test_inputs_not_modified(void)
{
  Fq U[MEDS_n * MEDS_n];
  Fq V[MEDS_n * MEDS_n];
  Fq W[MEDS_n * MEDS_n];
  Fq M[MEDS_n * MEDS_n * MEDS_n];
  Fq M_copy[MEDS_n * MEDS_n * MEDS_n];
  Fq u1[MEDS_n];
  Fq u1_copy[MEDS_n];
  int n = 4;

  require_success_fixture(M, u1, n, 0x11110000u);
  memcpy(M_copy, M, form_count(n) * sizeof(*M));
  memcpy(u1_copy, u1, (size_t)n * sizeof(*u1));

  CHECK(canonical_build_uvw_vartime(U, V, W, M, u1, n) == 0);
  CHECK(memcmp(M, M_copy, form_count(n) * sizeof(*M)) == 0);
  CHECK(memcmp(u1, u1_copy, (size_t)n * sizeof(*u1)) == 0);

  make_initial_phi_u_full_rank(M, u1, n);
  memcpy(M_copy, M, form_count(n) * sizeof(*M));
  memcpy(u1_copy, u1, (size_t)n * sizeof(*u1));
  fill_constant(U, mat_count(n), 77);
  fill_constant(V, mat_count(n), 77);
  fill_constant(W, mat_count(n), 77);

  CHECK(canonical_build_uvw_vartime(U, V, W, M, u1, n) == -1);
  CHECK(memcmp(M, M_copy, form_count(n) * sizeof(*M)) == 0);
  CHECK(memcmp(u1, u1_copy, (size_t)n * sizeof(*u1)) == 0);
  CHECK(outputs_unchanged(U, V, W, 77, n));
}

int main(void)
{
  test_bad_dimensions();
  test_null_pointers();
  test_zero_representative_preserves_outputs();
  test_inputs_not_modified();
  return 0;
}
