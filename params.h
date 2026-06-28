#ifndef PARAMS_H
#define PARAMS_H

#include <stdint.h>

#define MEDS_name "MEDS9923"

#define MEDS_digest_bytes 32
#define MEDS_pub_seed_bytes 32
#define MEDS_sec_seed_bytes 32
#define MEDS_round_seed_bytes 16
#define MEDS_salt_bytes 32

#define MEDS_p 4093
typedef uint16_t Fq;
typedef uint16_t trine_challenge_t;
#define Fq_bits 12
#define Fq_bytes 2

#define MEDS_n 22

#define MEDS_X 4
#define MEDS_r 144
#define MEDS_K 14

#define MEDS_FORM_COUNT (MEDS_X + 1)
#define MEDS_NONBASE_COUNT MEDS_X
#define MEDS_BASE_FORM_INDEX MEDS_X

#if MEDS_n < 5
#error "TRINE CF requires MEDS_n >= 5"
#endif

#if MEDS_K > MEDS_r
#error "MEDS_K must not exceed MEDS_r"
#endif

#if MEDS_X < 2 || MEDS_X > 256
#error "TRINE requires 2 <= MEDS_X <= 256"
#endif

#define MEDS_t_mask 0x000007FF
#define MEDS_t_bytes 2

#define MEDS_s_mask 0x00000003

#define MEDS_CEIL_DIV(x,y) (((x) + (y) - 1) / (y))
#define CEILING(x,y) MEDS_CEIL_DIV(x,y)

#define MEDS_MAT_BYTES MEDS_CEIL_DIV((MEDS_n * MEDS_n) * Fq_bits, 8)
#define MEDS_G_BYTES MEDS_CEIL_DIV((MEDS_n * MEDS_n * MEDS_n) * Fq_bits, 8)
#define MEDS_PK_BYTES (MEDS_pub_seed_bytes + (MEDS_X - 1) * MEDS_G_BYTES)
#define MEDS_SK_BYTES (MEDS_sec_seed_bytes + MEDS_pub_seed_bytes + 3 * (MEDS_X - 1) * MEDS_MAT_BYTES)
#define MEDS_RESPONSE_BYTES (MEDS_K * 3 * MEDS_MAT_BYTES)
#define MEDS_ZERO_SEED_COUNT (MEDS_r - MEDS_K)
#define MEDS_ZERO_SEED_BYTES (MEDS_ZERO_SEED_COUNT * MEDS_round_seed_bytes)
#define MEDS_RESPONSE_OFFSET 0
#define MEDS_ZERO_SEED_OFFSET (MEDS_RESPONSE_OFFSET + MEDS_RESPONSE_BYTES)
#define MEDS_DIGEST_OFFSET (MEDS_ZERO_SEED_OFFSET + MEDS_ZERO_SEED_BYTES)
#define MEDS_SALT_OFFSET (MEDS_DIGEST_OFFSET + MEDS_digest_bytes)
#define MEDS_SIG_BYTES (MEDS_SALT_OFFSET + MEDS_salt_bytes)

#define TRINE_VECTOR_BYTES MEDS_CEIL_DIV(MEDS_n * Fq_bits, 8)
#define TRINE_TRIFORM_BYTES MEDS_CEIL_DIV((MEDS_n * MEDS_n * MEDS_n) * Fq_bits, 8)
#define TRINE_PK_BYTES (MEDS_pub_seed_bytes + MEDS_X * TRINE_TRIFORM_BYTES)
#define TRINE_SK_BYTES MEDS_sec_seed_bytes
#define TRINE_RESPONSE_BYTES (MEDS_K * TRINE_VECTOR_BYTES)
#define TRINE_BASE_SEED_COUNT (MEDS_r - MEDS_K)
#define TRINE_BASE_SEED_BYTES (TRINE_BASE_SEED_COUNT * MEDS_round_seed_bytes)
#define TRINE_RESPONSE_OFFSET 0
#define TRINE_BASE_SEED_OFFSET (TRINE_RESPONSE_OFFSET + TRINE_RESPONSE_BYTES)
#define TRINE_DIGEST_OFFSET (TRINE_BASE_SEED_OFFSET + TRINE_BASE_SEED_BYTES)
#define TRINE_SALT_OFFSET (TRINE_DIGEST_OFFSET + MEDS_digest_bytes)
#define TRINE_SIG_BYTES (TRINE_SALT_OFFSET + MEDS_salt_bytes)

#endif
