#ifndef __LWPF_SW_PCR
#define __LWPF_SW_PCR
#include "sw5_pcr.h"
typedef struct perf_config_t{
  long pcr0, pcr1, pcr2, pcrc;
} perf_config_t;

static inline void set_perf_mode(long pcr0, long pcr1, long pcr2, long pcrc){
  pcr0 <<= 59;
  pcr1 <<= 59;
  pcr2 <<= 59;

  asm volatile("wcsr %0, 5\n"
               "wcsr %1, 6\n"
               "wcsr %2, 7\n"
               "wcsr %3, 8\n": :
               "r"(pcr0), "r"(pcr1), "r"(pcr2), "r"(pcrc));
}
static inline long pcr0(){
  long pcr;
  long msk = (1LL << 59) - 1;
  asm volatile("rcsr %0, 5" : "=r"(pcr));
  return pcr & msk;
}
static inline long pcr1(){
  long pcr;
  long msk = (1LL << 59) - 1;
  asm volatile("rcsr %0, 6" : "=r"(pcr));
  return pcr & msk;
}
static inline long pcr2(){
  long pcr;
  long msk = (1LL << 59) - 1;
  asm volatile("rcsr %0, 7" : "=r"(pcr));
  return pcr & msk;
}
static inline long rpcc(){
  long rpc;
  asm volatile("rcsr %0, 4" : "=r"(rpc));
  return rpc;
}
#endif
#ifdef CPE
#ifndef __LWPF_CPE_H
#define __LWPF_CPE_H

#include <slave.h>

//define a enum of kernels
#define _K(x) x
#define K(x) , x
typedef enum {
  LWPF_KERNELS
} lwpf_kernel;
#undef _K
#undef K

//define an array of kernel names
#define _K(x) #x
#define K(x) , #x
#define U(x) lwpf_kernel_names_ ## x
char *LWPF_UNIT[] = {LWPF_KERNELS};
#undef _K
#undef K
#undef U
//define arrays of kernels
#define _K(x) 1
#define K(x) +1
#define U(x) lwpf_kernel_count_ ## x
const long LWPF_UNIT = LWPF_KERNELS;
__thread_local static long KCNT = LWPF_KERNELS;
#undef U
__thread_local static long __rpcc[LWPF_KERNELS];
__thread_local static long __pcr0[LWPF_KERNELS];
__thread_local static long __pcr1[LWPF_KERNELS];
__thread_local static long __pcr2[LWPF_KERNELS];
#undef _K
#undef K

/* static inline void lwpf_start(lwpf_kernel kernel){ */
/*   __rpcc[kernel] -= rpcc(); */
/*   __pcr0[kernel] -= pcr0(); */
/*   __pcr1[kernel] -= pcr1(); */
/*   __pcr2[kernel] -= pcr2(); */
/* } */

/* static inline void lwpf_stop(lwpf_kernel kernel){ */
/*   __rpcc[kernel] += rpcc(); */
/*   __pcr0[kernel] += pcr0(); */
/*   __pcr1[kernel] += pcr1(); */
/*   __pcr2[kernel] += pcr2(); */
/* } */
#define lwpf_start(kernel) {                                    \
    long rpc, pc0, pc1, pc2;                                    \
    asm volatile("rcsr %0, 4\n"                                 \
                 "rcsr %1, 5\n"                                 \
                 "rcsr %2, 6\n"                                 \
                 "rcsr %3, 7\n"                                 \
                 : "=r"(rpc), "=r"(pc0), "=r"(pc1), "=r"(pc2)); \
    __rpcc[kernel] -= rpc;                                      \
    __pcr0[kernel] -= pc0;                                      \
    __pcr1[kernel] -= pc1;                                      \
    __pcr2[kernel] -= pc2;                                      \
  }                                                             \

#define lwpf_stop(kernel) {                                     \
    long rpc, pc0, pc1, pc2;                                    \
    asm volatile("rcsr %0, 4\n"                                 \
                 "rcsr %1, 5\n"                                 \
                 "rcsr %2, 6\n"                                 \
                 "rcsr %3, 7\n"                                 \
                 : "=r"(rpc), "=r"(pc0), "=r"(pc1), "=r"(pc2)); \
    __rpcc[kernel] += rpc;                                      \
    __pcr0[kernel] += pc0;                                      \
    __pcr1[kernel] += pc1;                                      \
    __pcr2[kernel] += pc2;                                      \
  }

#define U(x) lwpf_init_ ## x
void LWPF_UNIT(perf_config_t *conf){
  int i;
  for (i = 0; i < KCNT; i ++){
    __rpcc[i] = 0;
    __pcr0[i] = 0;
    __pcr1[i] = 0;
    __pcr2[i] = 0;
  }  
  set_perf_mode(conf->pcr0, conf->pcr1, conf->pcr2, conf->pcrc);
}
#undef U
#define U(x) lwpf_finish_ ## x
void LWPF_UNIT(long *mem){
  volatile int reply = 0;
  long *mrpcc = mem;
  long *mpcr0 = mem + KCNT * 64;
  long *mpcr1 = mem + KCNT * 2 * 64;
  long *mpcr2 = mem + KCNT * 3 * 64;
  //if (_MYID == 0)
  //  printf("%d\n", KCNT);
  reply = 0;
  athread_put(PE_MODE, __rpcc, mrpcc + KCNT * _MYID, sizeof(long) * KCNT, (void *)&reply, 0, 0);
  athread_put(PE_MODE, __pcr0, mpcr0 + KCNT * _MYID, sizeof(long) * KCNT, (void *)&reply, 0, 0);
  athread_put(PE_MODE, __pcr1, mpcr1 + KCNT * _MYID, sizeof(long) * KCNT, (void *)&reply, 0, 0);
  athread_put(PE_MODE, __pcr2, mpcr2 + KCNT * _MYID, sizeof(long) * KCNT, (void *)&reply, 0, 0);
  while (reply != 4);
}
#undef U
#endif
#endif
#ifdef MPE
#ifndef __LWPF_MPE_H
#define __LWPF_MPE_H
#include <athread.h>

//extern kernel names
#define U(x) extern char *lwpf_kernel_names_ ## x[];
LWPF_UNITS
#undef U

//extern kernel counts
#define U(x) extern long lwpf_kernel_count_ ## x;
LWPF_UNITS
#undef U

//extern init functions
#define U(x) extern SLAVE_FUN(lwpf_init_ ## x)(void *);
LWPF_UNITS
#undef U

//extern finish functions
#define U(x) extern SLAVE_FUN(lwpf_finish_ ## x)();
LWPF_UNITS
#undef U

//create init function
#define U(x) athread_spawn(lwpf_init_ ## x, conf); athread_join();
static inline void lwpf_init(perf_config_t *conf){
  LWPF_UNITS
}
#undef U
inline int longlen(long in){
  if (in == 0)
    return 1;
  int p = in, ret = 0;
  while (p > 0){
    p /= 10;
    ret ++;
  }
  return ret;
}
inline void calc_unit(long max, char *unit, long *div){
  *unit = '1';
  *div = 1;
  if (max > 1000000) {
    *unit = 'K';
    *div = 1000;
  }
  if (max > 1000000000) {
    *unit = 'M';
    *div = 1000000;
  }
  if (max > 1000000000000) {
    *unit = 'G';
    *div = 1000000000;
  }
  if (max > 1000000000000000) {
    *unit = 'T';
    *div = 1000000000000;
  }
}
inline void lwpf_print(FILE *f, long *ctrs, int kcnt, const char *unit, char *kname[]){
  long *crpcc = ctrs;
  long *cpcr0 = ctrs + kcnt * 64;
  long *cpcr1 = ctrs + kcnt * 2 * 64;
  long *cpcr2 = ctrs + kcnt * 3 * 64;

  fprintf(f, "counter staticists for %s\n", unit);
  long rpcc_max = 0;
  long pcr0_max = 0;
  long pcr1_max = 0;
  long pcr2_max = 0;

  long *rpcc_maxk = malloc(sizeof(long) * kcnt);
  long *pcr0_maxk = malloc(sizeof(long) * kcnt);
  long *pcr1_maxk = malloc(sizeof(long) * kcnt);
  long *pcr2_maxk = malloc(sizeof(long) * kcnt);

  long *rpcc_mink = malloc(sizeof(long) * kcnt);
  long *pcr0_mink = malloc(sizeof(long) * kcnt);
  long *pcr1_mink = malloc(sizeof(long) * kcnt);
  long *pcr2_mink = malloc(sizeof(long) * kcnt);

  long *rpcc_avgk = malloc(sizeof(long) * kcnt);
  long *pcr0_avgk = malloc(sizeof(long) * kcnt);
  long *pcr1_avgk = malloc(sizeof(long) * kcnt);
  long *pcr2_avgk = malloc(sizeof(long) * kcnt);

  //printf("%d\n", kcnt);
  int i, j;
  for (i = 0; i < kcnt; i ++){
    rpcc_maxk[i] = 0; rpcc_mink[i] = 1L << 62; rpcc_avgk[i] = 0;
    pcr0_maxk[i] = 0; pcr0_mink[i] = 1L << 62; pcr0_avgk[i] = 0;
    pcr1_maxk[i] = 0; pcr1_mink[i] = 1L << 62; pcr1_avgk[i] = 0;
    pcr2_maxk[i] = 0; pcr2_mink[i] = 1L << 62; pcr2_avgk[i] = 0;

    for (j = 0; j < 64; j ++){
      int ic = j * kcnt + i;
      if (crpcc[ic] > rpcc_maxk[i]) rpcc_maxk[i] = crpcc[ic];
      if (cpcr0[ic] > pcr0_maxk[i]) pcr0_maxk[i] = cpcr0[ic];
      if (cpcr1[ic] > pcr1_maxk[i]) pcr1_maxk[i] = cpcr1[ic];
      if (cpcr2[ic] > pcr2_maxk[i]) pcr2_maxk[i] = cpcr2[ic];
      if (crpcc[ic] < rpcc_mink[i]) rpcc_mink[i] = crpcc[ic];
      if (cpcr0[ic] < pcr0_mink[i]) pcr0_mink[i] = cpcr0[ic];
      if (cpcr1[ic] < pcr1_mink[i]) pcr1_mink[i] = cpcr1[ic];
      if (cpcr2[ic] < pcr2_mink[i]) pcr2_mink[i] = cpcr2[ic];
      rpcc_avgk[i] += crpcc[ic];
      pcr0_avgk[i] += cpcr0[ic];
      pcr1_avgk[i] += cpcr1[ic];
      pcr2_avgk[i] += cpcr2[ic];
    }
    rpcc_avgk[i] = (rpcc_avgk[i] + 32) / 64;
    pcr0_avgk[i] = (pcr0_avgk[i] + 32) / 64;
    pcr1_avgk[i] = (pcr1_avgk[i] + 32) / 64;
    pcr2_avgk[i] = (pcr2_avgk[i] + 32) / 64;

    if (rpcc_maxk[i] > rpcc_max) rpcc_max = rpcc_maxk[i];
    if (pcr0_maxk[i] > pcr0_max) pcr0_max = pcr0_maxk[i];
    if (pcr1_maxk[i] > pcr1_max) pcr1_max = pcr1_maxk[i];
    if (pcr2_maxk[i] > pcr2_max) pcr2_max = pcr2_maxk[i];
  }

  char rpcc_unit, pcr0_unit, pcr1_unit, pcr2_unit;
  long rpcc_div, pcr0_div, pcr1_div, pcr2_div;
  calc_unit(rpcc_max, &rpcc_unit, &rpcc_div);
  calc_unit(pcr0_max, &pcr0_unit, &pcr0_div);
  calc_unit(pcr1_max, &pcr1_unit, &pcr1_div);
  calc_unit(pcr2_max, &pcr2_unit, &pcr2_div);

  char *fmth = "|%16s|%18s(%c)|%18s(%c)|%18s(%c)|%18s(%c)|\n";
  char *fmtl = "|%16s|%7d%7d%7d|%7d%7d%7d|%7d%7d%7d|%7d%7d%7d|\n";
  //puts(fmth);
  //puts(fmtl);
  fprintf(f, fmth, "kernel", "rpcc(avg/max/min)", rpcc_unit, "pcr0(avg/max/min)", pcr0_unit, "pcr1(avg/max/min)", pcr1_unit, "pcr2(avg/max/min)", pcr2_unit);
  char hz[128];
  memset(hz, '-', sizeof(hz));
  hz[21] = 0;
  fprintf(f, "+----------------+%s+%s+%s+%s+\n", hz, hz, hz, hz);
  for (i = 0; i < kcnt; i ++){
    fprintf(f, fmtl, kname[i]
            , rpcc_avgk[i] / rpcc_div, rpcc_maxk[i] / rpcc_div, rpcc_mink[i] / rpcc_div
            , pcr0_avgk[i] / pcr0_div, pcr0_maxk[i] / pcr0_div, pcr0_mink[i] / pcr0_div
            , pcr1_avgk[i] / pcr1_div, pcr1_maxk[i] / pcr1_div, pcr1_mink[i] / pcr1_div
            , pcr2_avgk[i] / pcr2_div, pcr2_maxk[i] / pcr2_div, pcr2_mink[i] / pcr2_div);
  }

  free(rpcc_maxk);free(rpcc_mink);free(rpcc_avgk);
  free(pcr0_maxk);free(pcr0_mink);free(pcr0_avgk);
  free(pcr1_maxk);free(pcr1_mink);free(pcr1_avgk);
  free(pcr2_maxk);free(pcr2_mink);free(pcr2_avgk);
}
#define U(x)                                                            \
  ctrs = malloc(sizeof(long) * lwpf_kernel_count_ ## x * 64 * 4);       \
  athread_spawn(lwpf_finish_ ##x, ctrs);                                \
  athread_join();                                                       \
  lwpf_print(f, ctrs, lwpf_kernel_count_ ## x, #x, lwpf_kernel_names_ ##x); \
  free(ctrs);                                                           \

static inline void lwpf_finish(FILE *f){
  long *ctrs;
  LWPF_UNITS;
  puts("finished");
}
#undef U
#endif
#endif
