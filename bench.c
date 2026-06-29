#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <sys/random.h>
#include <sys/time.h>

#include "params.h"
#include "api.h"
#include "meds.h"

#include <x86intrin.h>

double osfreq(void);

long long cpucycles(void)
{
  return __rdtsc();
}

int main(int argc, char *argv[])
{
  printf("name: %s\n", CRYPTO_ALGNAME);
  printf("parameter set: %s (PARAMS=%d)\n", TRINE_PARAMETER_SET_NAME, PARAMS);

  long long keygen_time = 0;
  long long sign_time = 0;
  long long verify_time = 0;

  int rounds = 1;

  if (argc > 1)
    rounds = atoi(argv[1]);

  if (rounds < 1)
    rounds = 1;

  char msg[16] = "TestTestTestTest";

  uint8_t *sk = calloc(CRYPTO_SECRETKEYBYTES, sizeof(*sk));
  uint8_t *pk = calloc(CRYPTO_PUBLICKEYBYTES, sizeof(*pk));
  uint8_t *sig = calloc((size_t)CRYPTO_BYTES + sizeof(msg), sizeof(*sig));
  uint8_t *msg_out = calloc(sizeof(msg), sizeof(*msg_out));
  unsigned long long sig_len = 0;
  int exit_code = EXIT_FAILURE;

  if (sk == NULL || pk == NULL || sig == NULL || msg_out == NULL)
  {
    fprintf(stderr, "failed to allocate benchmark buffers\n");
    goto cleanup;
  }


  printf(
      "parameters (n, q, r, K, X): (%d, %d, %d, %d, %d)\n",
      TRINE_n,
      TRINE_q,
      TRINE_r,
      TRINE_K,
      TRINE_X);
  printf("pk:   %i bytes\n", TRINE_PK_BYTES);
  printf("sk:   %i bytes\n", TRINE_SK_BYTES);
  printf("sig:  %i bytes\n", TRINE_SIG_BYTES);


  for (int round = 0; round < rounds; round++)
  {
    keygen_time = -cpucycles();
    if (crypto_sign_keypair(pk, sk) != 0)
    {
      fprintf(stderr, "crypto_sign_keypair failed\n");
      goto cleanup;
    }
    keygen_time += cpucycles();

    sign_time = -cpucycles();
    sig_len = 0;
    if (crypto_sign(sig, &sig_len, (const unsigned char *)msg, sizeof(msg), sk) != 0)
    {
      fprintf(stderr, "crypto_sign failed\n");
      goto cleanup;
    }
    if (sig_len != (unsigned long long)CRYPTO_BYTES + (unsigned long long)sizeof(msg))
    {
      fprintf(stderr, "crypto_sign produced an unexpected signed-message length\n");
      goto cleanup;
    }
    sign_time += cpucycles();

    unsigned long long msg_out_len = 0;

    verify_time = -cpucycles();
    int ret = crypto_sign_open(msg_out, &msg_out_len, sig, sig_len, pk);
    verify_time += cpucycles();

    if (ret != 0 ||
        msg_out_len != (unsigned long long)sizeof(msg) ||
        memcmp(msg_out, msg, sizeof(msg)) != 0)
    {
      fprintf(stderr, "crypto_sign_open failed or recovered the wrong message\n");
      goto cleanup;
    }

    double freq = osfreq() / 1000;
    printf("F: %f\n", freq);

    printf("%f (%llu cycles)  ", keygen_time / freq, keygen_time);
    printf("%f (%llu cycles)  ", sign_time / freq, sign_time);
    printf("%f (%llu cycles)  \n", verify_time / freq, verify_time);
  }

  exit_code = EXIT_SUCCESS;

cleanup:
  free(sk);
  free(pk);
  free(sig);
  free(msg_out);
  return exit_code;
}
