#ifndef CPUCYCLES_H
#define CPUCYCLES_H
#include <stdint.h>


#ifdef __cplusplus
extern "C" {
#endif

uint64_t cpucycles(void);
uint64_t cpucycles_overhead(void);

#ifdef __cplusplus
}
#endif

#endif
