#ifdef MPE
#include <athread.h>
#endif

#ifdef CPE
#include <slave.h>
#define pe_init() volatile int reply = 0; int pe_cnt = 0;
#define pe_get(mem, ldm, size) {athread_get(PE_MODE, (mem), (ldm), (size), (void*)&reply, 0, 0, 0); pe_cnt ++;}
#define pe_put(mem, ldm, size) {athread_put(PE_MODE, (ldm), (mem), (size), (void*)&reply, 0, 0); pe_cnt ++;}
#define pe_syn() {while (reply != pe_cnt); asm volatile("memb");}
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
  typedef struct fix_nve_param_t{
    double **x, **v, **f;
    double *rmass, *mass;
    double dtv, dtf;
    int *type, *mask;
    int nlocal, groupbit;
  }fix_nve_param_t;

#ifdef __cplusplus
}
#endif
