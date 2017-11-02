#include "sunway.h"
#include <stdlib.h>
#include <simd.h>

#ifdef MPE
#define LWPF_UNITS U(LJCUT)
#include "lwpf.h"
int r = 0;
extern SLAVE_FUN(pair_lj_cut_sunway_compute_para)(compute_param_t *);
extern SLAVE_FUN(pair_lj_cut_sunway_compute_a2s)(compute_param_t *);
void pair_lj_cut_sunway_compute(compute_param_t *pm){
  if (athread_idle() == 0)
    athread_init();
  if (r == 0){
    perf_config_t conf;
    conf.pcr0 = PC0_CYC_CHNL_BLOCK;
    conf.pcr1 = PC1_CYCLE;
    conf.pcr2 = PC2_CNT_DMA_REQ;
    conf.pcrc = PCRC_ALL_PC;
    lwpf_init(&conf);
  }
  atom_in_t *atom_in = malloc(sizeof(atom_in_t) * (pm->nlocal + pm->nghost + 32));
  //printf("%d %d %d %d %d %d\n", pm->eflag_either, pm->vflag_either, pm->eflag_atom, pm->vflag_atom, pm->eflag_global, pm->vflag_global);
  pm->atom_in = ((long)atom_in | 255) + 1;
  athread_spawn(pair_lj_cut_sunway_compute_a2s, pm);
  athread_join();
  athread_spawn(pair_lj_cut_sunway_compute_para, pm);
  athread_join();
  if (r == 10 && pm->rank == 0){
    lwpf_finish(stdout);
  }
  r ++;
  free(atom_in);
}
#endif
#ifdef CPE
#define LWPF_KERNELS _K(COMP) K(CACHE) K(IREAD) K(IWRITE) K(ALL) K(GST) K(JLOOP)
#define LWPF_UNIT U(LJCUT)
#include "lwpf.h"

/* #define JPAGE_SIZE 64 */
/* #define JLIST_SIZE 256 */
/* #define IPAGE_SIZE 16 */
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

#define JLIST_PAGESIZE 128
#define ILIST_PAGESIZE 128

#define JCACHE_HBIT 10
#define JCACHE_SBIT 3
#define JCACHE_LINESIZE (1 << JCACHE_SBIT)
#define JCACHE_LINECNT (1 << (JCACHE_HBIT - JCACHE_SBIT))
#define JCACHE_MMASK (JCACHE_LINESIZE - 1)
#define JCACHE_LMASK (JCACHE_LINECNT - 1)

void pair_lj_cut_sunway_compute_a2s(compute_param_t *pm){
  pe_init();
  compute_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(compute_param_t));
  pe_syn();
  double x[ILIST_PAGESIZE][3];
  int type[ILIST_PAGESIZE];
  atom_in_t atom_in[ILIST_PAGESIZE];
  int ipage_start, i;
  int inum = l_pm.nlocal + l_pm.nghost;
  for (ipage_start = ILIST_PAGESIZE * _MYID; ipage_start < inum; ipage_start += ILIST_PAGESIZE * 64){
    int ipage_end = ipage_start + ILIST_PAGESIZE;
    if (ipage_end > inum)
      ipage_end = inum;
    int ipage_size = ipage_end - ipage_start;
    pe_get(l_pm.x[ipage_start], x, ipage_size * sizeof(double) * 3);
    pe_get(l_pm.type + ipage_start, type, ipage_size * sizeof(int));
    pe_syn();
    for (i = 0; i < ipage_size; i ++){
      atom_in[i].x[0] = x[i][0];
      atom_in[i].x[1] = x[i][1];
      atom_in[i].x[2] = x[i][2];
      atom_in[i].type = type[i];
    }
    pe_put(l_pm.atom_in + ipage_start, atom_in, sizeof(atom_in_t) * ipage_size);
    pe_syn();
  }
}

void pair_lj_cut_sunway_compute_para(compute_param_t *pm){
  lwpf_start(ALL);
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
  int ntp1 = l_pm.ntypes + 1;
  double cutsq[ntp1][ntp1], lj1[ntp1][ntp1], lj2[ntp1][ntp1], lj3[ntp1][ntp1], lj4[ntp1][ntp1], offset[ntp1][ntp1];
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
  numneigh = pm->numneigh;
  firstneigh = pm->firstneigh;
  int jlist_buf[JLIST_PAGESIZE];

  atom_in_t j_cache[JCACHE_LINECNT][JCACHE_LINESIZE];
  int jcache_tag[JCACHE_LINECNT];
  double *xjc, xj[3];
  int tj;
  double xi[ILIST_PAGESIZE][3], fi[ILIST_PAGESIZE][3];
  double ei[ILIST_PAGESIZE], vi[ILIST_PAGESIZE][6], v[6];
  int ti[ILIST_PAGESIZE], *fn[ILIST_PAGESIZE], nn[ILIST_PAGESIZE];
  int ipage_start;
  lwpf_start(COMP);
  for (i = 0; i < JCACHE_LINECNT; i ++)
    jcache_tag[i] = -1;
  for (ipage_start = ILIST_PAGESIZE * _MYID; ipage_start < inum; ipage_start += ILIST_PAGESIZE * 64){
    int ipage_end = ipage_start + ILIST_PAGESIZE;
    if (ipage_end > inum)
      ipage_end = inum;
    int ipage_size = ipage_end - ipage_start;
    lwpf_start(IREAD);
    pe_get(x[ipage_start], xi, ipage_size * sizeof(double) * 3);
    pe_get(type + ipage_start, ti, ipage_size * sizeof(int));
    pe_get(firstneigh + ipage_start, fn, ipage_size * sizeof(int*));
    pe_get(numneigh + ipage_start, nn, ipage_size * sizeof(int));

    pe_syn();
    lwpf_stop(IREAD);
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
        lwpf_start(JLOOP);
        for (jj = 0; jj < jpage_size; jj++) {
          j = jlist_buf[jj];
          factor_lj = special_lj[sbmask(j)];
          j &= NEIGHMASK;
          int line = j >> JCACHE_SBIT & JCACHE_LMASK;
          int tag = j >> JCACHE_HBIT;

          if (jcache_tag[line] != tag){
            int mem = j & ~JCACHE_MMASK;
            lwpf_start(CACHE);
            pe_get(l_pm.atom_in + mem, j_cache[line], sizeof(atom_in_t) * JCACHE_LINESIZE);
            pe_syn();
            jcache_tag[line] = tag;
            lwpf_stop(CACHE);
          }

          atom_in_t *jc = j_cache[line] + (j & JCACHE_MMASK);
          xj[0] = jc->x[0];
          xj[1] = jc->x[1];
          xj[2] = jc->x[2];
          jtype = jc->type;

          delx = xi[ioff][0] - xj[0];
          dely = xi[ioff][1] - xj[1];
          delz = xi[ioff][2] - xj[2];
          rsq = delx*delx + dely*dely + delz*delz;

          if (rsq < cutsq[itype][jtype]) {
            r2inv = 1.0/rsq;
            r6inv = r2inv*r2inv*r2inv;
            forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
            fpair = factor_lj*forcelj*r2inv;

            fi[ioff][0] += delx*fpair;
            fi[ioff][1] += dely*fpair;
            fi[ioff][2] += delz*fpair;

            if (l_pm.eflag) {
              evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
                offset[itype][jtype];
              evdwl *= factor_lj;
            }
            if (l_pm.evflag) ev_tally_full(&l_pm, ioff, evdwl,0.0,fpair,delx,dely,delz, eng_coul, eng_vdwl, virial, ei, vi);
          }
        }
        lwpf_stop(JLOOP);
      }
    }
    lwpf_start(IWRITE);
    pe_put(f[ipage_start], fi[0], sizeof(double) * 3 * ipage_size);
    if (l_pm.eflag_either && l_pm.eflag_atom)
      pe_put(l_pm.eatom + ipage_start, ei, sizeof(double) * ipage_size);
    if (l_pm.vflag_either && l_pm.vflag_atom)
      pe_put(l_pm.vatom[ipage_start], vi[0], sizeof(double) * 6 * ipage_size);
    pe_syn();
    lwpf_stop(IWRITE);
  }
  lwpf_stop(COMP);
  reg_reduce_inplace_doublev4(eng_virial, 2);
  lwpf_start(GST);
  if (_MYID == 0){
    pe_put(&(pm->eng_vdwl), eng_vdwl, sizeof(double) * 8);
    pe_syn();
    /* if (l_pm.eflag_either && l_pm.eflag_global) { */
    /*   pm->eng_vdwl = *eng_vdwl; */
    /*   pm->eng_coul = *eng_coul; */
    /* } */
    /* if (l_pm.vflag_either && l_pm.vflag_global) { */
    /*   pm->virial[0] = virial[0]; */
    /*   pm->virial[1] = virial[1]; */
    /*   pm->virial[2] = virial[2]; */
    /*   pm->virial[3] = virial[3]; */
    /*   pm->virial[4] = virial[4]; */
    /*   pm->virial[5] = virial[5]; */
    /* } */
  }
  lwpf_stop(GST);
  lwpf_stop(ALL);
}

#endif
