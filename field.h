#ifndef FIELD_H
#define FIELD_H

#include <stddef.h>
#include <stdint.h>

#include "params.h"

static inline Fq GF_add(Fq a, Fq b)
{
  uint32_t x = (uint32_t)a + b;
  x -= MEDS_p;
  x += MEDS_p & (uint32_t)-(int32_t)(x >> 31);
  return (Fq)x;
}

static inline Fq GF_sub(Fq a, Fq b)
{
  int32_t x = (int32_t)a - b;
  x += MEDS_p & -(x < 0);
  return (Fq)x;
}

static inline Fq GF_neg(Fq a)
{
  return (a == 0) ? 0 : (Fq)(MEDS_p - a);
}

static inline Fq GF_mul(Fq a, Fq b)
{
  return (Fq)(((uint32_t)a * b) % MEDS_p);
}

Fq GF_inv(Fq a);
int GF_inv_checked(Fq *out, Fq a);

int GF_batch_inv(
    Fq *out,
    const Fq *in,
    size_t count);

#endif
