#ifndef PARAMS_H
#define PARAMS_H

#define MEDS_name "MEDS9923"

#define MEDS_digest_bytes 32
#define MEDS_pub_seed_bytes 32
#define MEDS_sec_seed_bytes 32
#define MEDS_st_seed_bytes 16
#define MEDS_st_salt_bytes 32

#define MEDS_p 4093
#define GFq_t uint16_t
#define GFq_bits 12
#define GFq_bytes 2

#define MEDS_m 14
#define MEDS_n 14
#define MEDS_k 14

#define MEDS_s 4
#define MEDS_t 1152
#define MEDS_w 14

#define MEDS_seed_tree_height 11
#define SEED_TREE_size 4095
#define MEDS_max_path_len 100

#define MEDS_t_mask 0x000007FF
#define MEDS_t_bytes 2

#define MEDS_s_mask 0x00000003

#define CEILING(x,y) (((x) + (y) - 1) / (y))

#define MEDS_G_BYTES CEILING((MEDS_k * (MEDS_m * MEDS_n - MEDS_k)) * GFq_bits, 8)
#define MEDS_PK_BYTES (MEDS_pub_seed_bytes + (MEDS_s - 1) * MEDS_G_BYTES)
#define MEDS_SK_BYTES 1828
#define MEDS_SIG_BYTES 9896

#endif
