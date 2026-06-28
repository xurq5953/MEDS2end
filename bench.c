#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

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

  long long keygen_time = 0;
  long long sign_time = 0;
  long long verify_time = 0;

  int rounds = 128;

  if (argc > 1)
    rounds = atoi(argv[1]);

  char msg[16] = "TestTestTestTest";

  uint8_t sk[CRYPTO_SECRETKEYBYTES] = {0};
  uint8_t pk[CRYPTO_PUBLICKEYBYTES] = {0};

  uint8_t sig[CRYPTO_BYTES + sizeof(msg)] = {0};
  unsigned long long sig_len = sizeof(sig);


  printf("m:    %i\n", TRINE_n);
  printf("n:    %i\n", TRINE_n);
  printf("q:    %i\n", TRINE_q);
  printf("k:    %i\n", TRINE_n);
  printf("s:    %i\n", TRINE_X);
  printf("t:    %i\n", TRINE_r);
  printf("w:    %i\n", TRINE_K);

  printf("pk:   %i bytes\n", TRINE_PK_BYTES);

  printf("sig:  %i bytes\n", TRINE_SIG_BYTES);


  for (int round = 0; round < rounds; round++)
  {
    keygen_time = -cpucycles();
    crypto_sign_keypair(pk, sk);
    keygen_time += cpucycles();

    sign_time = -cpucycles();
    crypto_sign(sig, &sig_len, (const unsigned char *)msg, sizeof(msg), sk);
    sign_time += cpucycles();

    unsigned char msg_out[17];
    unsigned long long msg_out_len = sizeof(msg_out);

    verify_time = -cpucycles();
    int ret = crypto_sign_open(msg_out, &msg_out_len, sig, sizeof(sig), pk);
    verify_time += cpucycles();

    if (ret != 0)
    {
      fprintf(stderr, "\n\n   ERROR! Signature did not verify!\n\n\n");

      exit(-1);
    }

    double freq = osfreq() / 1000;
    printf("F: %f\n", freq);

    printf("%f (%llu cycles)  ", keygen_time / freq, keygen_time);
    printf("%f (%llu cycles)  ", sign_time / freq, sign_time);
    printf("%f (%llu cycles)  \n", verify_time / freq, verify_time);
  }

  return 0;
}

