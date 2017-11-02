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
  typedef struct atom_in_t{
    double x[3];
    int type, pad;
  } atom_in_t;
  typedef struct compute_param_t{
    int *ilist, *numneigh, **firstneigh;
    //vars to be put back
    double (*x)[3], (*f)[3], (*vatom)[6], *eatom;
    atom_in_t *atom_in;
    double *cutsq, *lj1, *lj2, *lj3, *lj4, *offset;
    int *type;
    int nlocal, nghost, ntotal, inum, ntypes, rank;
    double special_lj[4];
    //vars to be copy back
    double eng_vdwl, eng_coul, virial[6];
    double eflag, vflag, evflag;
    double eflag_global, vflag_global;
    double eflag_atom, vflag_atom;
    double eflag_either, vflag_either;
  } compute_param_t;
#define NEIGH_ISTEP 64
  typedef struct neigh_param_t{
    //packed bins
    int *binhead, *atom2bin;
    int *binpack, *binpacktype, *binpackhead;
    double (*binpackx)[3];

    //params needed
    double (*x)[3];
    int *type;
    double *cutneighsq;
    
    //stencil params
    int nstencil, *stencil;

    //neighpointers
    int **neighptr, (*neighoffset)[NEIGH_ISTEP + 32], *iilast;

    //numbers
    int nlocal, nghost, ntotal, istart, ntypes;
    int maxchunk;
  } neigh_param_t;
#ifdef __cplusplus
}
#endif
