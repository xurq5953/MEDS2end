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
#define Fq_bits 12
#define Fq_bytes 2

#define MEDS_n 22

#define MEDS_X 4
#define MEDS_r 144
#define MEDS_K 14

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

#endif
