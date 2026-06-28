#ifndef PARAMS_H
#define PARAMS_H

#include <stdint.h>

#define PARAMS 1

#if PARAMS == 1
#define TRINE_lambda 128
#define TRINE_n 22
#define TRINE_q 4093
#define TRINE_r 138
#define TRINE_K 52
#define TRINE_X 1
#elif PARAMS == 2
#define TRINE_lambda 256
#define TRINE_n 40
#define TRINE_q 4093
#define TRINE_r 271
#define TRINE_K 104
#define TRINE_X 1
#elif PARAMS == 3
#define TRINE_lambda 512
#define TRINE_n 81
#define TRINE_q 4093
#define TRINE_r 534
#define TRINE_K 211
#define TRINE_X 1
#elif PARAMS == 4
#define TRINE_lambda 128
#define TRINE_n 22
#define TRINE_q 4093
#define TRINE_r 61
#define TRINE_K 36
#define TRINE_X 4
#elif PARAMS == 5
#define TRINE_lambda 256
#define TRINE_n 40
#define TRINE_q 4093
#define TRINE_r 116
#define TRINE_K 76
#define TRINE_X 4
#elif PARAMS == 6
#define TRINE_lambda 512
#define TRINE_n 81
#define TRINE_q 4093
#define TRINE_r 232
#define TRINE_K 149
#define TRINE_X 4
#endif


#define TRINE_q_bits 12

typedef uint16_t Fq;
typedef uint16_t trine_challenge_t;

#define TRINE_FORM_COUNT (TRINE_X + 1)
#define TRINE_NONBASE_COUNT TRINE_X
#define TRINE_BASE_FORM_INDEX TRINE_X

#if TRINE_n < 5
#error "TRINE CF requires TRINE_n >= 5"
#endif

#if TRINE_K > TRINE_r
#error "TRINE_K must not exceed TRINE_r"
#endif

#if TRINE_X < 1 || TRINE_X > 256
#error "TRINE requires 1 <= TRINE_X <= 256"
#endif

#if TRINE_lambda % 8 != 0
#error "TRINE_lambda must be byte aligned"
#endif

_Static_assert(
    ((uint32_t)1u << (TRINE_q_bits - 1)) < TRINE_q,
    "TRINE_q_bits is too large");

_Static_assert(
    TRINE_q <= ((uint32_t)1u << TRINE_q_bits),
    "TRINE_q_bits is too small");

#define TRINE_CEIL_DIV(x,y) (((x) + (y) - 1) / (y))

#define TRINE_lambda_bytes (TRINE_lambda / 8)
#define TRINE_digest_bytes (2 * TRINE_lambda_bytes)
#define TRINE_public_seed_bytes (2 * TRINE_lambda_bytes)
#define TRINE_secret_seed_bytes (2 * TRINE_lambda_bytes)
#define TRINE_round_seed_bytes TRINE_lambda_bytes
#define TRINE_salt_bytes (2 * TRINE_lambda_bytes)

#define TRINE_VECTOR_BYTES TRINE_CEIL_DIV(TRINE_n * TRINE_q_bits, 8)
#define TRINE_MATRIX_BYTES TRINE_CEIL_DIV((TRINE_n * TRINE_n) * TRINE_q_bits, 8)
#define TRINE_TRIFORM_BYTES TRINE_CEIL_DIV((TRINE_n * TRINE_n * TRINE_n) * TRINE_q_bits, 8)
#define TRINE_PK_BYTES (TRINE_public_seed_bytes + TRINE_X * TRINE_TRIFORM_BYTES)
#define TRINE_SK_BYTES TRINE_secret_seed_bytes
#define TRINE_RESPONSE_BYTES TRINE_CEIL_DIV((TRINE_K * TRINE_n * TRINE_q_bits), 8)
#define TRINE_BASE_SEED_COUNT (TRINE_r - TRINE_K)
#define TRINE_BASE_SEED_BYTES (TRINE_BASE_SEED_COUNT * TRINE_round_seed_bytes)
#define TRINE_RESPONSE_OFFSET 0
#define TRINE_BASE_SEED_OFFSET (TRINE_RESPONSE_OFFSET + TRINE_RESPONSE_BYTES)
#define TRINE_DIGEST_OFFSET (TRINE_BASE_SEED_OFFSET + TRINE_BASE_SEED_BYTES)
#define TRINE_SALT_OFFSET (TRINE_DIGEST_OFFSET + TRINE_digest_bytes)
#define TRINE_SIG_BYTES (TRINE_SALT_OFFSET + TRINE_salt_bytes)

_Static_assert(TRINE_VECTOR_BYTES == 33, "unexpected vector size");
_Static_assert(TRINE_TRIFORM_BYTES == 15972, "unexpected triform size");
_Static_assert(TRINE_PK_BYTES == 16004, "unexpected public-key size");
_Static_assert(TRINE_SK_BYTES == 32, "unexpected secret-key size");
_Static_assert(TRINE_RESPONSE_BYTES == 1716, "unexpected response size");
_Static_assert(TRINE_BASE_SEED_COUNT == 86, "unexpected base-seed count");
_Static_assert(TRINE_BASE_SEED_BYTES == 1376, "unexpected base-seed segment size");
_Static_assert(TRINE_SIG_BYTES == 3156, "unexpected signature size");

#endif
