#include "trine_expand.h"

#include <stdint.h>
#include <string.h>

#include "fips202.h"
#include "matrixelim.h"
#include "triform.h"
#include "util.h"

enum
{
  TRINE_PURPOSE_PUBLIC_SEED = 1,
  TRINE_PURPOSE_BASE_FORM = 2,
  TRINE_PURPOSE_MATRIX_A = 3,
  TRINE_PURPOSE_MATRIX_B = 4,
  TRINE_PURPOSE_MATRIX_C = 5
};

static const uint8_t TRINE_EXPAND_PREFIX[] =
    "MEDS2END-TRINE-EXPAND-v1";

static int role_to_purpose(
    trine_matrix_role_t role,
    uint8_t *purpose)
{
  switch (role)
  {
    case TRINE_MATRIX_A:
      *purpose = TRINE_PURPOSE_MATRIX_A;
      return 0;
    case TRINE_MATRIX_B:
      *purpose = TRINE_PURPOSE_MATRIX_B;
      return 0;
    case TRINE_MATRIX_C:
      *purpose = TRINE_PURPOSE_MATRIX_C;
      return 0;
    default:
      return -1;
  }
}

static void store_u32_le(
    uint8_t out[4],
    uint32_t value)
{
  out[0] = (uint8_t)(value & 0xffu);
  out[1] = (uint8_t)((value >> 8) & 0xffu);
  out[2] = (uint8_t)((value >> 16) & 0xffu);
  out[3] = (uint8_t)((value >> 24) & 0xffu);
}

static void init_domain_shake(
    keccak_state *shake,
    uint8_t purpose,
    uint32_t index,
    const uint8_t *seed,
    size_t seed_len)
{
  uint8_t index_le[4];

  store_u32_le(index_le, index);

  shake256_init(shake);
  shake256_absorb(
      shake,
      TRINE_EXPAND_PREFIX,
      sizeof(TRINE_EXPAND_PREFIX) - 1u);
  shake256_absorb(shake, &purpose, 1u);
  shake256_absorb(shake, index_le, sizeof(index_le));
  shake256_absorb(shake, seed, seed_len);
  shake256_finalize(shake);
}

static int ranges_overlap(
    const Fq *a,
    const Fq *b,
    size_t count)
{
  if (a == NULL || b == NULL)
    return 0;

  const uintptr_t a_begin = (uintptr_t)a;
  const uintptr_t a_end = a_begin + count * sizeof(*a);
  const uintptr_t b_begin = (uintptr_t)b;
  const uintptr_t b_end = b_begin + count * sizeof(*b);

  return a_begin < b_end && b_begin < a_end;
}

int trine_expand_public_seed(
    uint8_t out_public_seed[MEDS_pub_seed_bytes],
    const uint8_t secret_seed[MEDS_sec_seed_bytes])
{
  if (out_public_seed == NULL || secret_seed == NULL)
    return -1;

  uint8_t tmp[MEDS_pub_seed_bytes];
  keccak_state shake;

  init_domain_shake(
      &shake,
      TRINE_PURPOSE_PUBLIC_SEED,
      0,
      secret_seed,
      MEDS_sec_seed_bytes);
  shake256_squeeze(tmp, sizeof(tmp), &shake);

  memcpy(out_public_seed, tmp, sizeof(tmp));
  return 0;
}

int trine_expand_base_form(
    Fq *out_base_form,
    const uint8_t public_seed[MEDS_pub_seed_bytes],
    int n)
{
  if (out_base_form == NULL || public_seed == NULL)
    return -1;

  if (n < 1 || n > MEDS_n)
    return -1;

  Fq tmp[triform_element_count(n)];
  keccak_state shake;

  init_domain_shake(
      &shake,
      TRINE_PURPOSE_BASE_FORM,
      0,
      public_seed,
      MEDS_pub_seed_bytes);

  for (size_t i = 0; i < triform_element_count(n); i++)
    tmp[i] = rnd_GF(&shake);

  memcpy(out_base_form, tmp, sizeof(tmp));
  return 0;
}

int trine_expand_secret_matrix_pair_vartime(
    Fq *out_matrix,
    Fq *out_inverse,
    const uint8_t secret_seed[MEDS_sec_seed_bytes],
    trine_matrix_role_t role,
    uint32_t index,
    int n)
{
  if (secret_seed == NULL)
    return -1;

  if (out_matrix == NULL && out_inverse == NULL)
    return -1;

  if (n < 1 || n > MEDS_n)
    return -1;

  if (index >= (uint32_t)MEDS_X)
    return -1;

  const size_t matrix_elements = (size_t)n * (size_t)n;
  if (ranges_overlap(out_matrix, out_inverse, matrix_elements))
    return -1;

  uint8_t purpose;
  if (role_to_purpose(role, &purpose) != 0)
    return -1;

  uint8_t matrix_seed[MEDS_sec_seed_bytes];
  keccak_state shake;
  Fq matrix_tmp[matrix_elements];
  Fq inverse_tmp[matrix_elements];

  init_domain_shake(&shake, purpose, index, secret_seed, MEDS_sec_seed_bytes);
  shake256_squeeze(matrix_seed, sizeof(matrix_seed), &shake);

  rnd_inv_matrix(matrix_tmp, n, n, matrix_seed, sizeof(matrix_seed));

  if (out_inverse != NULL &&
      pmod_mat_inv_vartime(inverse_tmp, matrix_tmp, n) != 0)
    return -1;

  if (out_matrix != NULL)
    memcpy(out_matrix, matrix_tmp, sizeof(matrix_tmp));

  if (out_inverse != NULL)
    memcpy(out_inverse, inverse_tmp, sizeof(inverse_tmp));

  return 0;
}
