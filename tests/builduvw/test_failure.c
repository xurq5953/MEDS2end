#include "test_reference.h"

static void expect_failure_preserves_outputs(const Fq *M, const Fq *u1, int n)
{
  Fq U[MEDS_n * MEDS_n];
  Fq V[MEDS_n * MEDS_n];
  Fq W[MEDS_n * MEDS_n];

  fill_constant(U, mat_count(n), 91);
  fill_constant(V, mat_count(n), 91);
  fill_constant(W, mat_count(n), 91);

  CHECK(canonical_build_uvw_vartime(U, V, W, M, u1, n) == -1);
  CHECK(outputs_unchanged(U, V, W, 91, n));
}

static void test_initial_phi_u_full_rank(void)
{
  Fq M[MEDS_n * MEDS_n * MEDS_n];
  Fq u1[MEDS_n];
  int n = 4;

  make_initial_phi_u_full_rank(M, u1, n);
  CHECK(ref_build_uvw((Fq[MEDS_n * MEDS_n]){0}, (Fq[MEDS_n * MEDS_n]){0}, (Fq[MEDS_n * MEDS_n]){0}, M, u1, n) ==
        REF_BUILD_PHI_U_NOT_CORANK1);
  expect_failure_preserves_outputs(M, u1, n);
}

static void test_initial_phi_u_corank_two(void)
{
  Fq M[MEDS_n * MEDS_n * MEDS_n];
  Fq u1[MEDS_n];
  int n = 4;

  make_initial_phi_u_corank_two(M, u1, n);
  CHECK(ref_build_uvw((Fq[MEDS_n * MEDS_n]){0}, (Fq[MEDS_n * MEDS_n]){0}, (Fq[MEDS_n * MEDS_n]){0}, M, u1, n) ==
        REF_BUILD_PHI_U_NOT_CORANK1);
  expect_failure_preserves_outputs(M, u1, n);
}

static void test_initial_phi_u_zero(void)
{
  Fq M[MEDS_n * MEDS_n * MEDS_n];
  Fq u1[MEDS_n];
  int n = 3;

  fill_zero(M, form_count(n));
  fill_basis_vector(u1, 0, n);
  CHECK(ref_build_uvw((Fq[MEDS_n * MEDS_n]){0}, (Fq[MEDS_n * MEDS_n]){0}, (Fq[MEDS_n * MEDS_n]){0}, M, u1, n) ==
        REF_BUILD_PHI_U_NOT_CORANK1);
  expect_failure_preserves_outputs(M, u1, n);
}

static void test_u_dependency_fixture(void)
{
  Fq M[MEDS_n * MEDS_n * MEDS_n];
  Fq u1[MEDS_n];
  int dims[] = {3, 5};

  for (size_t i = 0; i < sizeof(dims) / sizeof(dims[0]); i++)
  {
    int n = dims[i];
    make_u_dependency_fixture(M, u1, n);
    CHECK(ref_build_uvw((Fq[MEDS_n * MEDS_n]){0}, (Fq[MEDS_n * MEDS_n]){0}, (Fq[MEDS_n * MEDS_n]){0}, M, u1, n) ==
          REF_BUILD_U_DEPENDENT);
    expect_failure_preserves_outputs(M, u1, n);
  }
}

int main(void)
{
  test_initial_phi_u_full_rank();
  test_initial_phi_u_corank_two();
  test_initial_phi_u_zero();
  test_u_dependency_fixture();
  return 0;
}
