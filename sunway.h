#ifndef SUNWAY_H_
#define SUNWAY_H_
#ifdef MPE
#include <athread.h>
#endif

#ifdef CPE
#include <slave.h>
#define pe_init() volatile int reply = 0; int pe_cnt = 0;
#define pe_get(mem, ldm, size) {athread_get(PE_MODE, (mem), (ldm), (size), (void*)&reply, 0, 0, 0); pe_cnt ++;}
#define pe_put(mem, ldm, size) {athread_put(PE_MODE, (ldm), (mem), (size), (void*)&reply, 0, 0); pe_cnt ++;}
#define pe_syn() {while (reply != pe_cnt); asm volatile("memb");}
#define dma_rpl(desc, mem, ldm, reply) asm("dma %0, %1, %2\n\t" : : "r"(desc), "r"(mem), "r"(ldm), "r"(&reply) : "memory");
#endif
#ifdef __cplusplus
extern "C"{
#endif
#ifndef NEIGHMASK
#define NEIGHMASK 0x3fffffff
#define SBBITS 30
#define sbmask(j) ((j >> SBBITS) & 3)
#endif
  typedef struct bin_pack_atom_t{
    double x[3];
    int type, id;
    int *firstneigh;
    //int numneigh;
  } bin_pack_atom_t;
  typedef struct atom_in_t{
    double x[3];
    int type, sbj;
  } atom_in_t;

  typedef struct neigh_param_t{
    //packed bins
    int *binpackhead, *binpacknn;
    /* int *binpack, *binpacktype, *binpackhead; */
    /* double (*binpackx)[3]; */
    bin_pack_atom_t *binpack;
    //neighpointers
    /* int **binpackfn, *binpacknn; */

    //params needed
    double *cutneighsq;
    
    //stencil params
    int nstencil, *stencil;


    //numbers
    int nlocal, nghost, ntotal, ntypes, mbins;
    int maxchunk;
  } neigh_param_t;
#ifdef CPE
#include <simd.h>
  static inline void reg_reduce_inplace_doublev4(doublev4 *arr, int len){
    int i, j;
    doublev4 tmp;
    for (i = 1; i < 8; i += i){
      if ((_ROW & i) == i){
        for (j = 0; j < len; j ++)
          asm("putc %0, %1": : "r"(arr[j]), "r"(_ROW ^ i));
      }
      if ((_ROW & i) == 0){
        for (j = 0; j < len; j ++){
          asm(
              "getc %0\n\t"
              "vaddd %0, %1, %1\n\t"
              : "=r"(tmp), "+r"(arr[j]));
          //arr[j] += tmp;
        }
      }
      athread_syn(COL_SCOPE, 0xff);
    }
    athread_syn(ARRAY_SCOPE, 0xffff);
    if (_ROW == 0){
      for (i = 1; i < 8; i += i){
        if ((_COL & i) == i){
          for (j = 0; j < len; j ++)
            asm("putr %0, %1": : "r"(arr[j]), "r"(_COL ^ i));
        }
        if ((_COL & i) == 0){
          for (j = 0; j < len; j ++){
            /* asm("getr %0\n" : "=r"(tmp)); */
            /* arr[j] += tmp; */
            asm(
                "getr %0\n\t" 
                "vaddd %0, %1, %1\n\t"
                : "=r"(tmp), "+r"(arr[j]));
          }
        }
      }
      athread_syn(ROW_SCOPE, 0xff);
    }
  }
#endif
#ifdef __cplusplus
}
#endif
#endif
