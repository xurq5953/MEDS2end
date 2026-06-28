#include <stdint.h>
#include <string.h>

#include "fips202.h"

#include "log.h"

#include "util.h"
#include "meds.h"
#include "matrixmod.h"
#include "matrixelim.h"


void XOF(uint8_t **buf, size_t *length, const uint8_t *seed, size_t seed_len, int num)
{
  keccak_state shake;

  shake256_absorb_once(&shake, seed, seed_len);

  for (int i = 0; i < num; i++)
    shake256_squeeze(buf[i], length[i], &shake);
}

Fq rnd_GF(keccak_state *shake)
{
  Fq val = MEDS_p;

  while (val >= MEDS_p)
  {
    uint8_t data[sizeof(Fq)];
    
    shake256_squeeze(data, sizeof(Fq), shake);

    val = 0;
    for (int i = 0; i < sizeof(Fq); i++)
      val |= data[i] << (i*8);

    val = val & ((1 << Fq_bits) - 1);
  }

  return val;
}

void rnd_sys_mat(Fq *M, int M_r, int M_c, const uint8_t *seed, size_t seed_len)
{
  keccak_state shake;
  shake256_absorb_once(&shake, seed, seed_len);

  for (int r = 0; r < M_r; r++)
    for (int c = M_r; c < M_c; c++)
      pmod_mat_set_entry(M, M_r, M_c, r, c, rnd_GF(&shake));

  for (int r = 0; r < M_r; r++)
    for (int c = 0; c < M_r; c++)
      if (r == c)
        pmod_mat_set_entry(M, M_r, M_c, r, c, 1);
      else
        pmod_mat_set_entry(M, M_r, M_c, r, c, 0);
}

void rnd_inv_matrix(Fq *M, int M_r, int M_c, uint8_t *seed, size_t seed_len)
{
  if (M_r != M_c)
    return;

  keccak_state shake;
  shake256_absorb_once(&shake, seed, seed_len);

  while (1)
  {
    for (int r = 0; r < M_r; r++)
      for (int c = 0; c < M_c; c++)
        pmod_mat_set_entry(M, M_r, M_c, r, c, rnd_GF(&shake));

    if (pmod_mat_is_invertible_vartime(M, M_r))
      return;
  }
}

static unsigned bit_length_u32(uint32_t value)
{
  unsigned bits = 0;

  do
  {
    bits++;
    value >>= 1;
  }
  while (value != 0);

  return bits;
}

static int sample_masked_le(
    uint32_t *out,
    keccak_state *shake,
    unsigned bit_count)
{
  if (out == NULL || shake == NULL || bit_count == 0 || bit_count > 32)
    return -1;

  const size_t byte_count = MEDS_CEIL_DIV((size_t)bit_count, 8u);
  uint8_t buf[4] = {0};
  uint32_t value = 0;

  shake256_squeeze(buf, byte_count, shake);

  for (size_t i = 0; i < byte_count; i++)
    value |= (uint32_t)buf[i] << (8u * i);

  if (bit_count < 32)
    value &= ((uint32_t)1u << bit_count) - 1u;

  *out = value;
  return 0;
}

int trine_parse_hash(
    const uint8_t *digest,
    size_t digest_len,
    trine_challenge_t *out,
    size_t out_len)
{
  if (digest == NULL || out == NULL)
    return -1;

  if (digest_len != MEDS_digest_bytes || out_len != MEDS_r)
    return -1;

  trine_challenge_t parsed[MEDS_r];
  uint16_t raw[MEDS_r] = {0};
  keccak_state shake;
  const unsigned position_bits = bit_length_u32((uint32_t)(MEDS_r - 1));
  const unsigned value_bits = bit_length_u32((uint32_t)MEDS_X);

  shake256_init(&shake);
  shake256_absorb(&shake, digest, digest_len);
  shake256_finalize(&shake);

  for (int selected = 0; selected < MEDS_K; selected++)
  {
    uint32_t position = 0;
    uint32_t value = 0;

    do
    {
      if (sample_masked_le(&position, &shake, position_bits) != 0)
        return -1;
    }
    while (position >= (uint32_t)MEDS_r || raw[position] != 0);

    do
    {
      if (sample_masked_le(&value, &shake, value_bits) != 0)
        return -1;
    }
    while (value == 0 || value > (uint32_t)MEDS_X);

    raw[position] = (uint16_t)value;
  }

  for (int i = 0; i < MEDS_r; i++)
    parsed[i] = (trine_challenge_t)(MEDS_X - raw[i]);

  memcpy(out, parsed, sizeof(parsed));
  return 0;
}

int solve(Fq *A, Fq *B_inv, Fq *G0prime, Fq Amm)
{
  Fq P0prime0[MEDS_n*MEDS_n];
  Fq P0prime1[MEDS_n*MEDS_n];

  for (int i = 0; i < MEDS_n*MEDS_n; i++)
  {
    P0prime0[i] = G0prime[i];
    P0prime1[i] = G0prime[i + MEDS_n * MEDS_n];
  }

  Fq N[MEDS_n * MEDS_n];

  for (int i = 0; i < MEDS_n; i++)
    for (int j = 0; j < MEDS_n; j++)
      N[j*MEDS_n + i] = (MEDS_p - P0prime0[i*MEDS_n + j]) % MEDS_p;

  //LOG_MAT(N, MEDS_n, MEDS_n);


  Fq M[MEDS_n*(MEDS_n + MEDS_n + 2)] = {0};

  for (int i = 0; i < MEDS_n; i++)
    for (int j = 0; j < MEDS_n; j++)
      M[j*(MEDS_n + MEDS_n + 2) + i] = (MEDS_p - P0prime1[i*MEDS_n + j]) % MEDS_p;

  for (int i = 0; i < MEDS_n; i++)
    for (int j = 0; j < MEDS_n; j++)
      M[j*(MEDS_n + MEDS_n + 2) + i + MEDS_n] = P0prime0[i*MEDS_n + j];

  for (int j = 0; j < MEDS_n; j++)
    M[j*(MEDS_n + MEDS_n + 2) + MEDS_n + MEDS_n] = ((uint32_t)P0prime0[(MEDS_n-1)*MEDS_n + j] * (MEDS_p - (uint32_t)Amm)) % MEDS_p;

  for (int j = 0; j < MEDS_n; j++)
    M[j*(MEDS_n + MEDS_n + 2) + MEDS_n + MEDS_n + 1] = ((uint32_t)P0prime1[(MEDS_n-1)*MEDS_n + j] * (uint32_t)Amm) % MEDS_p;


  //LOG_MAT(M, MEDS_n, MEDS_n + MEDS_n + 2);

  if (pmod_mat_syst_ct(M, MEDS_n-1, MEDS_n + MEDS_n + 2) < 0)
    return -1;

  //LOG_MAT_FMT(M, MEDS_n, MEDS_n + MEDS_n + 2, "M part");

  // eliminate last row
  for (int r = 0; r < MEDS_n-1; r++)
  {
    uint64_t factor = pmod_mat_entry(M, MEDS_n, MEDS_n + MEDS_n + 2, MEDS_n-1, r);

    // ignore last column
    for (int c = MEDS_n-1; c < MEDS_n + MEDS_n + 1; c++)
    {
      uint64_t tmp0 = pmod_mat_entry(M, MEDS_n, MEDS_n + MEDS_n + 2, MEDS_n-1, c);
      uint64_t tmp1 = pmod_mat_entry(M, MEDS_n, MEDS_n + MEDS_n + 2, r, c);

      int64_t val = (tmp1 * factor) % MEDS_p;

      val = tmp0 - val;

      val += MEDS_p * (val < 0);

      pmod_mat_set_entry(M, MEDS_n, MEDS_n + MEDS_n + 2,  MEDS_n-1, c, val);
    }

    pmod_mat_set_entry(M, MEDS_n, MEDS_n + MEDS_n + 2, MEDS_n-1, r, 0);
  }

  // normalize last row
  {
    uint64_t val = pmod_mat_entry(M, MEDS_n, MEDS_n + MEDS_n + 2, MEDS_n-1, MEDS_n-1);

    if (val == 0)
      return -1;

    val = GF_inv(val);

    // ignore last column
    for (int c = MEDS_n; c < MEDS_n + MEDS_n + 1; c++)
    {
      uint64_t tmp = pmod_mat_entry(M, MEDS_n, MEDS_n + MEDS_n + 2, MEDS_n-1, c);

      tmp = (tmp * val) % MEDS_p;

      pmod_mat_set_entry(M, MEDS_n, MEDS_n + MEDS_n + 2, MEDS_n-1, c, tmp);
    }
  }

  pmod_mat_set_entry(M, MEDS_n, MEDS_n + MEDS_n + 2, MEDS_n-1, MEDS_n-1, 1);

  M[MEDS_n*(MEDS_n + MEDS_n + 2)-1] = 0;

  //LOG_MAT_FMT(M, MEDS_n, MEDS_n + MEDS_n + 2, "M red");

  // back substitute
  for (int r = 0; r < MEDS_n-1; r++)
  {
    uint64_t factor = pmod_mat_entry(M, MEDS_n, MEDS_n + MEDS_n + 2, r, MEDS_n-1);

    // ignore last column
    for (int c = MEDS_n; c < MEDS_n + MEDS_n + 1; c++)
    {
        uint64_t tmp0 = pmod_mat_entry(M, MEDS_n, MEDS_n + MEDS_n + 2, MEDS_n-1, c);
        uint64_t tmp1 = pmod_mat_entry(M, MEDS_n, MEDS_n + MEDS_n + 2, r, c);

        int64_t val = (tmp0 * factor) % MEDS_p;

        val = tmp1 - val;

        val += MEDS_p * (val < 0);

        pmod_mat_set_entry(M, MEDS_n, MEDS_n + MEDS_n + 2,  r, c, val);
    }

    pmod_mat_set_entry(M, M_r, MEDS_n + MEDS_n + 2, r, MEDS_n-1, 0);
  }


  //LOG_MAT_FMT(M, MEDS_n, MEDS_n + MEDS_n + 2, "M done");


  Fq sol[MEDS_n*MEDS_n + MEDS_n*MEDS_n] = {0};

  sol[MEDS_n*MEDS_n + MEDS_n * MEDS_n - 1] = Amm;

  for (int i = 0; i < MEDS_n - 1; i++)
    sol[MEDS_n*MEDS_n + MEDS_n*MEDS_n - MEDS_n + i] = M[(i+1) * (MEDS_n + MEDS_n + 2) - 1];

  for (int i = 0; i < MEDS_n; i++)
    sol[MEDS_n*MEDS_n + MEDS_n*MEDS_n - 2*MEDS_n + i] = M[(i+1) * (MEDS_n + MEDS_n + 2) - 2];

  for (int i = 0; i < MEDS_n; i++)
    sol[MEDS_n*MEDS_n - MEDS_n + i] = ((uint32_t)P0prime0[(MEDS_n-1)*MEDS_n + i] * (uint32_t)Amm) % MEDS_p;

  //LOG_VEC_FMT(sol, MEDS_n*MEDS_n + MEDS_n*MEDS_n, "initial sol");


  // incomplete blocks:

  for (int i = 0; i < MEDS_n; i++)
    for (int j = 0; j < MEDS_n-1; j++)
    {
      uint32_t tmp = (uint32_t)sol[MEDS_n*MEDS_n + MEDS_n*MEDS_n - 2*MEDS_n + i] + MEDS_p -
        ((uint32_t)M[i * (MEDS_n + MEDS_n + 2) + MEDS_n + MEDS_n-2 - j] * (uint32_t)sol[MEDS_n*MEDS_n + MEDS_n*MEDS_n - 2 - j]) % MEDS_p;
      sol[MEDS_n*MEDS_n + MEDS_n*MEDS_n - 2*MEDS_n + i] = tmp % MEDS_p;
    }

  for (int i = 0; i < MEDS_n; i++)
    for (int j = 0; j < MEDS_n-1; j++)
    {
      uint32_t tmp = (uint32_t)sol[MEDS_n*MEDS_n - MEDS_n + i] + MEDS_p - 
        ((uint32_t)N[i * (MEDS_n) + MEDS_n-2 - j] * (uint32_t)sol[MEDS_n*MEDS_n + MEDS_n*MEDS_n - 2 - j]) % MEDS_p;
      sol[MEDS_n*MEDS_n - MEDS_n + i] = tmp % MEDS_p;
    }

  //LOG_VEC_FMT(sol, MEDS_n*MEDS_n + MEDS_n*MEDS_n, "incomplete blocks");


  // complete blocks:

  for (int block = 3; block <= MEDS_n; block++)
    for (int i = 0; i < MEDS_n; i++)
      for (int j = 0; j < MEDS_n; j++)
      {
        uint32_t tmp = sol[MEDS_n*MEDS_n + MEDS_n*MEDS_n - block*MEDS_n + i] + MEDS_p -
          ((uint32_t)M[i * (MEDS_n + MEDS_n + 2) + MEDS_n + MEDS_n-1 - j] * (uint32_t)sol[MEDS_n*MEDS_n + MEDS_n*MEDS_n - 1 - (block-2)*MEDS_n - j]) % MEDS_p;
        sol[MEDS_n*MEDS_n + MEDS_n*MEDS_n - block*MEDS_n + i] = tmp % MEDS_p;
      }

  for (int block = 2; block <= MEDS_n; block++)
    for (int i = 0; i < MEDS_n; i++)
      for (int j = 0; j < MEDS_n; j++)
      {
        uint32_t tmp = sol[MEDS_n*MEDS_n - block*MEDS_n + i] + MEDS_p - 
          ((uint32_t)N[i * (MEDS_n) + MEDS_n-1 - j] * (uint32_t)sol[MEDS_n*MEDS_n + MEDS_n*MEDS_n - 1 - (block-1)*MEDS_n - j]) % MEDS_p;
        sol[MEDS_n*MEDS_n - block*MEDS_n + i] = tmp % MEDS_p;
      }

  //LOG_VEC_FMT(sol, MEDS_n*MEDS_n + MEDS_n*MEDS_n, "complete blocks");


  for (int i = 0; i < MEDS_n * MEDS_n; i++)
    A[i] = sol[i + MEDS_n * MEDS_n];

  for (int i = 0; i < MEDS_n * MEDS_n; i++)
    B_inv[i] = sol[i];

  //LOG_MAT(A, MEDS_n, MEDS_n);
  //LOG_MAT(B_inv, MEDS_n, MEDS_n);

  return 0;
}


void G_mat_init(Fq *G, Fq *Gsub[MEDS_n])
{
  for (int i = 0; i < MEDS_n; i++)
    Gsub[i] = G + i*MEDS_n*MEDS_n;
}

void pi(Fq *Gout, Fq *A, Fq *B, Fq *G)
{
  Fq *G0sub[MEDS_n];
  G_mat_init(G, G0sub);

  Fq *Gsub[MEDS_n];
  G_mat_init(Gout, Gsub);

  for (int i = 0; i < MEDS_n; i++)
  {
    pmod_mat_mul(Gsub[i], A, G0sub[i], MEDS_n);
    pmod_mat_mul(Gsub[i], Gsub[i], B, MEDS_n);
  }
}

#if MEDS_n != MEDS_n || MEDS_n != MEDS_n
#error "MEDS2endGen phi() assumes MEDS_n == MEDS_n == MEDS_n"
#endif

void phi(Fq *Gout, Fq *A, Fq *B, Fq *C, Fq *G)
{
  Fq tmp[MEDS_n * MEDS_n * MEDS_n];

  pi(tmp, A, B, G);

  pmod_mat_mul_rect(Gout,
      MEDS_n, MEDS_n * MEDS_n,
      C, MEDS_n, MEDS_n,
      tmp, MEDS_n, MEDS_n * MEDS_n);
}
