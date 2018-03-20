#ifndef PAIR_TERSOFF_SW64_H_
#define PAIR_TERSOFF_SW64_H_
#ifdef __cplusplus
extern "C"{
#endif
  typedef struct tersoff_param_t{
    double lam1,lam2,lam3;
    double c,d,h;
    double gamma,powerm;
    double powern,beta,half_powern_inv;
    double biga,bigb,bigd,bigr,bigdinv,c2divd2;
    double /* cut, */cutsq;
    double c1/* ,c2,c3 */,c4;
    /* int ielement,jelement,kelement; */
    int powermint;
    /* double Z_i,Z_j;              // added for TersoffZBL */
    /* double ZBLcut,ZBLexpscale; */
    /* double c5,ca1,ca4;           // added for TersoffMOD */
    /* double powern_del; */
    /* double c0;                   // added for TersoffMODC */
  } tersoff_param_t;
  typedef struct short_neigh_t{
    //double prefactor_fwd, prefactor_rev;
    double d[3], r2, rinv;
    long type, idx;
    double padding;
  } short_neigh_t;

  typedef struct pair_tersoff_compute_param_t{
    int *ilist, *numneigh, **firstneigh, *numshort;
    int *firstshort;
    double (*fend)[4], (*ftmp)[4];
    int *fdone;
    long fend_base, fdone_base, ftmp_base, atom_in_base;
    int maxshort, cutshortsq;
    short_neigh_t *shortlist;
    int *shortidx;
    //global params
    int *elem2param, *map;
    tersoff_param_t *params;
    //vars to be put back
    double (*x)[3], (*f)[3], (*vatom)[6], *eatom;
    atom_in_t *atom_in;
    //atom_in_t *atom_in;
    int *type;
    int nlocal, nghost, ntotal, inum, gnum, ntypes, rank, nelements, nparams;
    //vars to be copy back
    double eng_vdwl, eng_coul, virial[6];
    int eflag, vflag, evflag;
    int eflag_global, vflag_global;
    int eflag_atom, vflag_atom;
    int eflag_either, vflag_either;
  } pair_tersoff_compute_param_t;
  void pair_tersoff_compute_attractive(pair_tersoff_compute_param_t *pm);
  void pair_tersoff_compute_zeta(pair_tersoff_compute_param_t *pm);
#ifdef __cplusplus
}
#endif
#endif
