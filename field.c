#include "field.h"

GFq_t GF_inv(GFq_t val)
{
  if (MEDS_p == 8191)
  {
    uint64_t tmp_0  = val;
    uint64_t tmp_1  = (tmp_0 * tmp_0) % MEDS_p;
    uint64_t tmp_2  = (tmp_1 * tmp_0) % MEDS_p;
    uint64_t tmp_3  = (tmp_2 * tmp_1) % MEDS_p;
    uint64_t tmp_4  = (tmp_3 * tmp_3) % MEDS_p;
    uint64_t tmp_5  = (tmp_4 * tmp_3) % MEDS_p;
    uint64_t tmp_6  = (tmp_5 * tmp_5) % MEDS_p;
    uint64_t tmp_7  = (tmp_6 * tmp_6) % MEDS_p;
    uint64_t tmp_8  = (tmp_7 * tmp_7) % MEDS_p;
    uint64_t tmp_9  = (tmp_8 * tmp_8) % MEDS_p;
    uint64_t tmp_10 = (tmp_9 * tmp_5) % MEDS_p;
    uint64_t tmp_11 = (tmp_10 * tmp_10) % MEDS_p;
    uint64_t tmp_12 = (tmp_11 * tmp_11) % MEDS_p;
    uint64_t tmp_13 = (tmp_12 * tmp_2) % MEDS_p;
    uint64_t tmp_14 = (tmp_13 * tmp_13) % MEDS_p;
    uint64_t tmp_15 = (tmp_14 * tmp_14) % MEDS_p;
    uint64_t tmp_16 = (tmp_15 * tmp_15) % MEDS_p;
    uint64_t tmp_17 = (tmp_16 * tmp_3) % MEDS_p;

    return (GFq_t)tmp_17;
  }
  else
  {
    uint64_t exponent = MEDS_p - 2;
    uint64_t t = 1;

    while (exponent > 0)
    {
      if ((exponent & 1) != 0)
        t = (t * (uint64_t)val) % MEDS_p;

      val = (GFq_t)(((uint64_t)val * (uint64_t)val) % MEDS_p);

      exponent >>= 1;
    }

    return (GFq_t)t;
  }
}

int GF_inv_checked(GFq_t *out, GFq_t a)
{
  if (a == 0)
    return -1;

  *out = GF_inv(a);
  return 0;
}

int GF_batch_inv(
    GFq_t *out,
    const GFq_t *in,
    size_t count)
{
  if (count == 0)
    return 0;

  GFq_t prefix[count];

  if (in[0] == 0)
    return -1;

  prefix[0] = in[0];

  for (size_t i = 1; i < count; i++)
  {
    if (in[i] == 0)
      return -1;

    prefix[i] = GF_mul(prefix[i - 1], in[i]);
  }

  GFq_t acc;
  if (GF_inv_checked(&acc, prefix[count - 1]) != 0)
    return -1;

  for (size_t i = count; i-- > 1;)
  {
    out[i] = GF_mul(acc, prefix[i - 1]);
    acc = GF_mul(acc, in[i]);
  }

  out[0] = acc;
  return 0;
}
