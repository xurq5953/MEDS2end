#ifndef API_H
#define API_H

#include "params.h"

#define CRYPTO_SECRETKEYBYTES TRINE_SK_BYTES
#define CRYPTO_PUBLICKEYBYTES TRINE_PK_BYTES
#define CRYPTO_BYTES TRINE_SIG_BYTES

#define CRYPTO_ALGNAME "MEDS9923"

int crypto_sign_keypair(
    unsigned char *pk,
    unsigned char *sk
  );

int crypto_sign(
    unsigned char *sm, unsigned long long *smlen,
    const unsigned char *m, unsigned long long mlen,
    const unsigned char *sk
  );

int crypto_sign_open(
    unsigned char *m, unsigned long long *mlen,
    const unsigned char *sm, unsigned long long smlen,
    const unsigned char *pk
  );

#endif
