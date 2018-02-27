#ifndef PAIR_REAXC_SW64_H
#define PAIR_REAXC_SW64_H
#include "reaxc_ctypes_sunway.h"
#include "reaxc_nsdef_sunway.h"
#ifdef __cplusplus
extern "C"{
#else
#endif
  typedef struct write_atoms_param_t{
    NSDEF(reax_system_c) *csys;
    int *map, *type;
    NSDEF(rc_tagint) *tag;
    double (*x)[3];
    double *q;
    int *num_bonds, *num_hbonds;
  } write_atoms_param_t;

  typedef struct write_lists_param_t{
    NSDEF(atom_pack_t) *patoms;
    NSDEF(reax_list) *far_nbrs;
    int maxfar;
    int *ilist;
    int *numneigh, **firstneigh;
    int numall;
    double nonb_cutsq;
  } write_lists_param_t;
  void write_reax_atoms_and_pack(write_atoms_param_t *pm);
  void write_reax_lists_c(write_lists_param_t *pm);
#ifdef __cplusplus
}
#endif
#endif

