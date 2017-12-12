#ifndef PAIR_LJ_CUT_SW64_H_
#define PAIR_LJ_CUT_SW64_H_
#ifdef __cplusplus
extern "C"{
#endif
  typedef struct atom_in_t{
    double x[3];
    int type, sbj;
  } atom_in_t;
  typedef struct pair_lj_cut_compute_param_t{
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
    int eflag, vflag, evflag;
    int eflag_global, vflag_global;
    int eflag_atom, vflag_atom;
    int eflag_either, vflag_either;
  } pair_lj_cut_compute_param_t;
#ifdef __cplusplus
}
#endif
#endif
