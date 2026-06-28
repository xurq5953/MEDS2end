#include "cpucycles.h"

uint64_t cpucycles(void)
{
  uint64_t result;
  __asm__ volatile("rdtsc;shlq $32, %%rdx; orq %%rdx, %%rax":"=a"(result)::"%rdx");
  return result;
}

uint64_t cpucycles_overhead(void)
{
  uint64_t t0,t1,overhead = UINT64_MAX;
  int i;
  for(i=0;i<100000;i++)
    {
      t0 = cpucycles();
      __asm__ volatile("");
      t1=cpucycles();
      if(t1-t0<overhead)
	overhead = t1-t0;
    }
  return overhead;
}
