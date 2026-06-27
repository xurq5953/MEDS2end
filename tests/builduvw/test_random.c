#include "test_reference.h"

static void make_random_case(Fq *M, Fq *u1, int n, test_rng_t *rng, int category)
{
  fill_random(M, form_count(n), rng);
  fill_basis_vector(u1, 0, n);

  switch (category)
  {
    case 0:
      break;
    case 1:
      make_initial_phi_u_full_rank(M, u1, n);
      break;
    case 2:
      force_coordinate_path_constraints(M, n);
      break;
    case 3:
      if (n >= 2)
        make_initial_phi_u_corank_two(M, u1, n);
      break;
    case 4:
      for (size_t i = 0; i < form_count(n); i++)
        if ((i % 5u) != 0u)
          M[i] = 0;
      break;
    default:
      fill_constant(M, form_count(n), (Fq)(MEDS_p - 1));
      force_coordinate_path_constraints(M, n);
      break;
  }
}

static void diff_round(int n, uint32_t seed, int round, test_rng_t *rng)
{
  Fq U[MEDS_n * MEDS_n];
  Fq V[MEDS_n * MEDS_n];
  Fq W[MEDS_n * MEDS_n];
  Fq ref_U[MEDS_n * MEDS_n];
  Fq ref_V[MEDS_n * MEDS_n];
  Fq ref_W[MEDS_n * MEDS_n];
  Fq M[MEDS_n * MEDS_n * MEDS_n];
  Fq u1[MEDS_n];
  int category = (int)(test_rng_next(rng) % 6u);
  ref_build_status_t ref_status;
  int prod_status;

  make_random_case(M, u1, n, rng, category);
  ref_status = ref_build_uvw(ref_U, ref_V, ref_W, M, u1, n);
  prod_status = canonical_build_uvw_vartime(U, V, W, M, u1, n);

  if ((ref_status == REF_BUILD_OK) != (prod_status == 0))
  {
    fprintf(
        stderr,
        "status mismatch seed=0x%08x n=%d round=%d category=%d ref=%s prod=%d\n",
        (unsigned)seed,
        n,
        round,
        category,
        ref_build_status_name(ref_status),
        prod_status);
    exit(1);
  }

  if (ref_status == REF_BUILD_OK)
  {
    expect_projective_columns_equal(U, ref_U, n);
    expect_projective_columns_equal(V, ref_V, n);
    expect_projective_columns_equal(W, ref_W, n);
    verify_success_properties(M, u1, U, V, W, n);
  }
}

static void run_diff_rounds(int n, int rounds, uint32_t group)
{
  uint32_t seed = 0x55550000u ^ ((uint32_t)n << 8) ^ group;
  test_rng_t rng;

  test_rng_init(&rng, seed);

  for (int round = 0; round < rounds; round++)
    diff_round(n, seed, round, &rng);
}

static void run_success_rounds(int n, int rounds, uint32_t group)
{
  for (int round = 0; round < rounds; round++)
  {
    uint32_t seed = 0x66660000u ^ ((uint32_t)n << 8) ^ group ^ (uint32_t)round;
    Fq U[MEDS_n * MEDS_n];
    Fq V[MEDS_n * MEDS_n];
    Fq W[MEDS_n * MEDS_n];
    Fq M[MEDS_n * MEDS_n * MEDS_n];
    Fq u1[MEDS_n];

    require_success_fixture(M, u1, n, seed);
    CHECK_CTX(
        canonical_build_uvw_vartime(U, V, W, M, u1, n) == 0,
        seed,
        n,
        round,
        "success-fixture");
    verify_success_properties(M, u1, U, V, W, n);
  }
}

int main(void)
{
  int small_dims[] = {1, 2, 3, 4, 5};

  for (size_t i = 0; i < sizeof(small_dims) / sizeof(small_dims[0]); i++)
  {
    run_diff_rounds(small_dims[i], 120, 0x01u);
    run_success_rounds(small_dims[i], 20, 0x02u);
  }

  run_diff_rounds(MEDS_n, 20, 0x03u);
  run_success_rounds(MEDS_n, 3, 0x04u);

  return 0;
}
