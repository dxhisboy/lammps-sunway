#include "sunway.h"
#include <stdlib.h>
#include <simd.h>

#ifdef MPE
void ev_tally_full(compute_param_t *pm, int i, double evdwl, double ecoul, double fpair,
                   double delx, double dely, double delz)
{
  double v[6];
  //printf("%f %f\n", pm->eatom[i], pm->eng_vdwl);
  if (pm->eflag_either) {
    if (pm->eflag_global) {
      pm->eng_vdwl += 0.5*evdwl;
      pm->eng_coul += 0.5*ecoul;
    }
    if (pm->eflag_atom) pm->eatom[i] += 0.5 * (evdwl + ecoul);
  }

  if (pm->vflag_either) {
    v[0] = 0.5*delx*delx*fpair;
    v[1] = 0.5*dely*dely*fpair;
    v[2] = 0.5*delz*delz*fpair;
    v[3] = 0.5*delx*dely*fpair;
    v[4] = 0.5*delx*delz*fpair;
    v[5] = 0.5*dely*delz*fpair;

    if (pm->vflag_global) {
      pm->virial[0] += v[0];
      pm->virial[1] += v[1];
      pm->virial[2] += v[2];
      pm->virial[3] += v[3];
      pm->virial[4] += v[4];
      pm->virial[5] += v[5];
    }

    if (pm->vflag_atom) {
      pm->vatom[i][0] += v[0];
      pm->vatom[i][1] += v[1];
      pm->vatom[i][2] += v[2];
      pm->vatom[i][3] += v[3];
      pm->vatom[i][4] += v[4];
      pm->vatom[i][5] += v[5];
    }
  }
}

extern SLAVE_FUN(pair_lj_cut_sunway_compute_para)(compute_param_t *);
void pair_lj_cut_sunway_compute(compute_param_t *pm){
  if (athread_idle() == 0)
    athread_init();
  athread_spawn(pair_lj_cut_sunway_compute_para, pm);
  //puts("spawned");
  athread_join();
  //puts("joined");
}
/* #define o(i, j) ((i) * (pm->ntypes + 1) + (j)) */
/* #define mcutsq(i, j) pm->cutsq[o(i, j)] */
/* #define mlj1(i, j) pm->lj1[o(i, j)] */
/* #define mlj2(i, j) pm->lj2[o(i, j)] */
/* #define mlj3(i, j) pm->lj3[o(i, j)] */
/* #define mlj4(i, j) pm->lj4[o(i, j)] */
/* #define moffset(i, j) pm->offset[o(i, j)] */

/* void pair_lj_cut_sunway_compute(compute_param_t *pm){ */
/*   int i,j,ii,jj,inum,jnum,itype,jtype; */
/*   double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair; */
/*   double rsq,r2inv,r6inv,forcelj,factor_lj; */
/*   int *ilist,*jlist,*numneigh,**firstneigh; */

/*   double (*x)[3] = pm->x; */
/*   double (*f)[3] = pm->f; */
/*   int *type = pm->type; */
/*   int nlocal = pm->nlocal; */
/*   double *special_lj = pm->special_lj; */

/*   inum = pm->inum; */
/*   ilist = pm->ilist; */
/*   numneigh = pm->numneigh; */
/*   firstneigh = pm->firstneigh; */

/*   // loop over neighbors of my atoms */

/*   for (ii = 0; ii < inum; ii++) { */
/*     i = ilist[ii]; */
/*     xtmp = x[i][0]; */
/*     ytmp = x[i][1]; */
/*     ztmp = x[i][2]; */
/*     itype = type[i]; */
/*     jlist = firstneigh[i]; */
/*     jnum = numneigh[i]; */

/*     for (jj = 0; jj < jnum; jj++) { */
/*       j = jlist[jj]; */
/*       factor_lj = special_lj[sbmask(j)]; */
/*       j &= NEIGHMASK; */

/*       delx = xtmp - x[j][0]; */
/*       dely = ytmp - x[j][1]; */
/*       delz = ztmp - x[j][2]; */
/*       rsq = delx*delx + dely*dely + delz*delz; */
/*       jtype = type[j]; */

/*       if (rsq < mcutsq(itype,jtype)) { */
/*         r2inv = 1.0/rsq; */
/*         r6inv = r2inv*r2inv*r2inv; */
/*         forcelj = r6inv * (mlj1(itype, jtype)*r6inv - mlj2(itype, jtype)); */
/*         fpair = factor_lj*forcelj*r2inv; */

/*         f[i][0] += delx*fpair; */
/*         f[i][1] += dely*fpair; */
/*         f[i][2] += delz*fpair; */

/*         if (pm->eflag) { */
/*           evdwl = r6inv*(mlj3(itype, jtype)*r6inv-mlj4(itype,jtype)) - */
/*             moffset(itype,jtype); */
/*           evdwl *= factor_lj; */
/*         } */
/*         //printf("%f %f %f\n", factor_lj, rsq, evdwl); */
/*         if (pm->evflag) ev_tally_full(pm, i,evdwl,0.0,fpair,delx,dely,delz); */
/*       } */
/*     } */
/*   } */
/* } */

#endif
#ifdef CPE
#define JPAGE_SIZE 64
#define JLIST_SIZE 256
#define IPAGE_SIZE 16
void ev_tally_full(compute_param_t *pm, int i, double evdwl, double ecoul, double fpair,
                   double delx, double dely, double delz,
                   double *eng_coul, double *eng_vdwl, double *virial,
                   double *eatom, double (*vatom)[6])
{
  double v[6];

  if (pm->eflag_either) {
    if (pm->eflag_global) {
      eng_vdwl[0] += 0.5*evdwl;
      eng_coul[0] += 0.5*ecoul;
    }
    if (pm->eflag_atom) eatom[i] += 0.5 * (evdwl + ecoul);
  }

  if (pm->vflag_either) {
    v[0] = 0.5*delx*delx*fpair;
    v[1] = 0.5*dely*dely*fpair;
    v[2] = 0.5*delz*delz*fpair;
    v[3] = 0.5*delx*dely*fpair;
    v[4] = 0.5*delx*delz*fpair;
    v[5] = 0.5*dely*delz*fpair;

    if (pm->vflag_global) {
      virial[0] += v[0];
      virial[1] += v[1];
      virial[2] += v[2];
      virial[3] += v[3];
      virial[4] += v[4];
      virial[5] += v[5];
    }

    if (pm->vflag_atom) {
      vatom[i][0] += v[0];
      vatom[i][1] += v[1];
      vatom[i][2] += v[2];
      vatom[i][3] += v[3];
      vatom[i][4] += v[4];
      vatom[i][5] += v[5];
    }
  }
}

/* inline void ev_tally_full(compute_param_t *pm, double evdwl, double ecoul, double fpair, */
/*                    double delx, double dely, double delz, */
/*                    double *eng_coul, double *eng_vdwl, double *virial, */
/*                    double *eatom, double *vatom) */
/* { */
/*   double v[6]; */

/*   if (pm->eflag_either) { */
/*     if (pm->eflag_global) { */
/*       eng_vdwl[0] += 0.5*evdwl; */
/*       eng_coul[0] += 0.5*ecoul; */
/*     } */
/*     if (pm->eflag_atom) *eatom += 0.5 * (evdwl + ecoul); */
/*   } */

/*   if (pm->vflag_either) { */
/*     v[0] = 0.5*delx*delx*fpair; */
/*     v[1] = 0.5*dely*dely*fpair; */
/*     v[2] = 0.5*delz*delz*fpair; */
/*     v[3] = 0.5*delx*dely*fpair; */
/*     v[4] = 0.5*delx*delz*fpair; */
/*     v[5] = 0.5*dely*delz*fpair; */

/*     if (pm->vflag_global) { */
/*       virial[0] += v[0]; */
/*       virial[1] += v[1]; */
/*       virial[2] += v[2]; */
/*       virial[3] += v[3]; */
/*       virial[4] += v[4]; */
/*       virial[5] += v[5]; */
/*     } */

/*     if (pm->vflag_atom) { */
/*       vatom[0] += v[0]; */
/*       vatom[1] += v[1]; */
/*       vatom[2] += v[2]; */
/*       vatom[3] += v[3]; */
/*       vatom[4] += v[4]; */
/*       vatom[5] += v[5]; */
/*     } */
/*   } */
/* } */

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
        asm("getc %0\n" : "=r"(tmp));
        arr[j] += tmp;
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
          asm("getr %0\n" : "=r"(tmp));
          arr[j] += tmp;
          //      "vaddd %0, %1, %0" : "=r"(arr[j]), "=r"(tmp));
        }
      }
    }
    athread_syn(ROW_SCOPE, 0xff);
  }
}

#define IPAGE_SIZE 16
#define o(i, j) ((i) * (l_pm.ntypes + 1) + (j))
#define mcutsq(i, j) cutsq[o(i, j)]
#define mlj1(i, j) lj1[o(i, j)]
#define mlj2(i, j) lj2[o(i, j)]
#define mlj3(i, j) lj3[o(i, j)]
#define mlj4(i, j) lj4[o(i, j)]
#define moffset(i, j) offset[o(i, j)]

#define MAXT2 16
#define JLIST_PAGESIZE 64
#define ILIST_PAGESIZE 16
#define JCACHE_LINESIZE 16
#define JCACHE_LINECNT 128
#define JCACHE_HBIT 11
#define JCACHE_MMASK 15
#define JCACHE_LMASK 127
#define JCACHE_SBIT 4
void pair_lj_cut_sunway_compute_para(compute_param_t *pm){
  pe_init();
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double (*x)[3] = pm->x;
  double (*f)[3] = pm->f;
  //double xi[IPAGE_SIZE][3], fi[IPAGE_SIZE][3];
  //double ti[IPAGE_SIZE];
  int *type = pm->type;
  int nlocal = pm->nlocal;
  double *special_lj;
  doublev4 eng_virial[2];
  double *eng_vdwl = (double*)(void*)eng_virial;
  double *eng_coul = eng_vdwl + 1;
  double *virial = eng_coul + 1;

  compute_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(compute_param_t));
  pe_syn();
  double cutsq[MAXT2], lj1[MAXT2], lj2[MAXT2], lj3[MAXT2], lj4[MAXT2], offset[MAXT2];
  pe_get(l_pm.cutsq, cutsq, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  pe_get(l_pm.lj1, lj1, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  pe_get(l_pm.lj2, lj2, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  pe_get(l_pm.lj3, lj3, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  pe_get(l_pm.lj4, lj4, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  pe_get(l_pm.offset, offset, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  special_lj = l_pm.special_lj;
  pe_syn();
  eng_virial[0] = 0;
  eng_virial[1] = eng_virial[0];
  inum = pm->inum;
  //ilist = pm->ilist;
  numneigh = pm->numneigh;
  firstneigh = pm->firstneigh;
  int jlist_buf[JLIST_PAGESIZE];
  double xj_cache[JCACHE_LINECNT][JCACHE_LINESIZE][3];
  int tj_cache[JCACHE_LINECNT][JCACHE_LINESIZE];
  int jcache_id[JCACHE_LINECNT];
  double *xjc, xj[3];
  int tj;
  double xi[ILIST_PAGESIZE][3], fi[ILIST_PAGESIZE][3];
  double ei[ILIST_PAGESIZE], vi[ILIST_PAGESIZE][6];
  int ti[ILIST_PAGESIZE], *fn[ILIST_PAGESIZE], nn[ILIST_PAGESIZE];
  int ipage_start;
  for (i = 0; i < JCACHE_LINECNT; i ++)
    jcache_id[i] = -1;
  for (ipage_start = ILIST_PAGESIZE * _MYID; ipage_start < inum; ipage_start += ILIST_PAGESIZE * 64){
    int ipage_end = ipage_start + ILIST_PAGESIZE;
    if (ipage_end > inum)
      ipage_end = inum;
    int ipage_size = ipage_end - ipage_start;
    pe_get(x[ipage_start], xi, ipage_size * sizeof(double) * 3);
    pe_get(type + ipage_start, ti, ipage_size * sizeof(int));
    pe_get(firstneigh + ipage_start, fn, ipage_size * sizeof(int*));
    pe_get(numneigh + ipage_start, nn, ipage_size * sizeof(int));
    /* if (l_pm.eflag_either && l_pm.eflag_atom) */
    /*   pe_get(l_pm.eatom + ipage_start, ei, sizeof(double) * ipage_size); */
    /* if (l_pm.vflag_either && l_pm.vflag_atom) */
    /*   pe_get(l_pm.vatom[ipage_start], vi[0], sizeof(double) * 6 * ipage_size); */

    pe_syn();
    for (ii = ipage_start; ii < ipage_end; ii ++) {
      i = ii;
      int ioff = i - ipage_start;
      fi[ioff][0] = fi[ioff][1] = fi[ioff][2] = 0;
      ei[ioff] = 0;
      vi[ioff][0] = vi[ioff][1] = vi[ioff][2] = 0;
      vi[ioff][3] = vi[ioff][4] = vi[ioff][5] = 0;
      itype = ti[ioff];
      jlist = fn[ioff];
      jnum = nn[ioff];
      int jpage_start;
      for (jpage_start = 0; jpage_start < jnum; jpage_start += JLIST_PAGESIZE){
        int jpage_size = JLIST_PAGESIZE;
        if (JLIST_PAGESIZE + jpage_start > jnum)
          jpage_size = jnum - jpage_start;
        pe_get(jlist + jpage_start, jlist_buf, jpage_size * sizeof(int));
        pe_syn();
        for (jj = 0; jj < jpage_size; jj++) {
          j = jlist_buf[jj];
          factor_lj = special_lj[sbmask(j)];
          j &= NEIGHMASK;
          /* pe_get(x[j], xj, sizeof(double) * 3); */
          /* pe_get(type + j, &jtype, sizeof(int)); */
          /* pe_syn(); */
          if (jcache_id[(j >> JCACHE_SBIT) & JCACHE_LMASK] != j >> JCACHE_SBIT){
            pe_get(x[j & ~JCACHE_MMASK], xj_cache[(j >> JCACHE_SBIT) & JCACHE_LMASK], sizeof(double) * 3 * JCACHE_LINESIZE);
            pe_get(type + (j & ~JCACHE_MMASK), tj_cache[(j >> JCACHE_SBIT) & JCACHE_LMASK], sizeof(int) * JCACHE_LINESIZE);
            pe_syn();
            jcache_id[(j >> JCACHE_SBIT) & JCACHE_LMASK] = j >> JCACHE_SBIT;
          }
          xjc = xj_cache[(j >> JCACHE_SBIT) & JCACHE_LMASK][j & JCACHE_MMASK];
          tj = tj_cache[(j >> JCACHE_SBIT) & JCACHE_LMASK][j & JCACHE_MMASK];
          
          /* pe_get(x[j], xj, sizeof(double) * 3); */
          /* pe_get(type + j, &jtype, sizeof(int)); */
          /* pe_syn(); */
          /* if (xjc[2] != xj[2]){ */
          /*   printf("%d %f %f\n", j, xjc[2], xj[2]); */
          /* } */
          xj[0] = xjc[0];
          xj[1] = xjc[1];
          xj[2] = xjc[2];
          jtype = tj;

          delx = xi[ioff][0] - xj[0];
          dely = xi[ioff][1] - xj[1];
          delz = xi[ioff][2] - xj[2];
          rsq = delx*delx + dely*dely + delz*delz;

          if (rsq < mcutsq(itype,jtype)) {
            r2inv = 1.0/rsq;
            r6inv = r2inv*r2inv*r2inv;
            forcelj = r6inv * (mlj1(itype, jtype)*r6inv - mlj2(itype, jtype));
            fpair = factor_lj*forcelj*r2inv;

            fi[ioff][0] += delx*fpair;
            fi[ioff][1] += dely*fpair;
            fi[ioff][2] += delz*fpair;

            if (l_pm.eflag) {
              evdwl = r6inv*(mlj3(itype, jtype)*r6inv-mlj4(itype,jtype)) -
                moffset(itype,jtype);
              evdwl *= factor_lj;
            }

            if (l_pm.evflag) ev_tally_full(&l_pm, ioff, evdwl,0.0,fpair,delx,dely,delz, eng_coul, eng_vdwl, virial, ei, vi);
          }
        }
      }
    }
    pe_put(f[ipage_start], fi[0], sizeof(double) * 3 * ipage_size);
    if (l_pm.eflag_either && l_pm.eflag_atom)
      pe_put(l_pm.eatom + ipage_start, ei, sizeof(double) * ipage_size);
    if (l_pm.vflag_either && l_pm.vflag_atom)
      pe_put(l_pm.vatom[ipage_start], vi[0], sizeof(double) * 6 * ipage_size);
    pe_syn();
  }
  reg_reduce_inplace_doublev4(eng_virial, 2);
  if (_MYID == 0){
    if (l_pm.eflag_either && l_pm.eflag_global) {
      pm->eng_vdwl = *eng_vdwl;
      pm->eng_coul = *eng_coul;
    }
    if (l_pm.vflag_either && l_pm.vflag_global) {
      pm->virial[0] = virial[0];
      pm->virial[1] = virial[1];
      pm->virial[2] = virial[2];
      pm->virial[3] = virial[3];
      pm->virial[4] = virial[4];
      pm->virial[5] = virial[5];
    }
  }
}

#endif
