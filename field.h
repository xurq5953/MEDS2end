#ifndef FIELD_H
#define FIELD_H

#include <stddef.h>
#include <stdint.h>

#include "params.h"

static inline GFq_t GF_add(GFq_t a, GFq_t b)
{
  uint32_t x = (uint32_t)a + b;
  x -= MEDS_p;
  x += MEDS_p & (uint32_t)-(int32_t)(x >> 31);
  return (GFq_t)x;
}

static inline GFq_t GF_sub(GFq_t a, GFq_t b)
{
  int32_t x = (int32_t)a - b;
  x += MEDS_p & -(x < 0);
  return (GFq_t)x;
}

static inline GFq_t GF_neg(GFq_t a)
{
  return (a == 0) ? 0 : (GFq_t)(MEDS_p - a);
}

static inline GFq_t GF_mul(GFq_t a, GFq_t b)
{
  return (GFq_t)(((uint32_t)a * b) % MEDS_p);
}

GFq_t GF_inv(GFq_t a);
int GF_inv_checked(GFq_t *out, GFq_t a);

int GF_batch_inv(
    GFq_t *out,
    const GFq_t *in,
    size_t count);

#endif
