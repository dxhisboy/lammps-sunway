#include "sunway.h"
#include "pair_tersoff_sw64.h"
#include <stdio.h>
#include <stdlib.h>
#include "gptl.h"
#define ISTEP 64
#define JSTEP 32

#ifdef MPE
#define LWPF_UNITS U(TERSOFF)
#include "lwpf.h"
#include <simd.h>
extern SLAVE_FUN(pair_tersoff_compute_attractive_para)(void*);
extern SLAVE_FUN(pair_tersoff_a2s)(void*);
extern SLAVE_FUN(pair_tersoff_reduction_force)(void*);

int r = 0;
void waitint(int *ptr){
  while (1){
    int tmp;
    asm ("ldw_nc %0, 0(%1)": "=&r"(tmp): "r"(ptr));
    //volatile int i;
    if (tmp != -1)
      break;
  }
}


void pair_tersoff_compute_attractive(pair_tersoff_compute_param_t *pm){
  if (athread_idle() == 0)
    athread_init();

  long fend_base = (long)pm->fend_base;
  long ftmp_base = (long)pm->ftmp_base;
  long fdone_base = (long)pm->fdone_base;
  long atom_in_base = (long)pm->atom_in_base;
  pm->fend = (void*)((fend_base + 255) & (~255));
  pm->fdone = (void*)((fdone_base + 255) & (~255));
  pm->ftmp = (void*)((ftmp_base + 255) & (~255));
  pm->atom_in = (void*)((atom_in_base + 255) & (~255));
  GPTLstart("tersoff a2s");
  athread_spawn(pair_tersoff_a2s, pm);
  GPTLstart("tersoff fill");
  double (*ftmp)[4] = pm->ftmp;
  double (*fend)[4] = pm->fend;
  int *fdone = pm->fdone;
  
  int ist, ied, ii;
  int ntotal = pm->ntotal;
  int *firstshort = pm->firstshort;
  double (*f)[3] = pm->f;
  short_neigh_t *shortlist = pm->shortlist;
  int *shortidx = pm->shortidx;
  int nlocal = pm->nlocal;
  int maxshort = pm->maxshort;
  int *numshort = pm->numshort;
  intv8 v8_0 = 0;
  doublev4 v4_0 = 0.0;
  for (ii = 0; ii < nlocal; ii += 4){
    simd_store(v4_0, ftmp[ii]);
    simd_store(v4_0, ftmp[ii + 1]);
    simd_store(v4_0, ftmp[ii + 2]);
    simd_store(v4_0, ftmp[ii + 3]);
  }

  GPTLstop("tersoff fill");
  athread_join();
  GPTLstop("tersoff a2s");
  /* if (r == 0){ */
  /*   perf_config_t conf; */
  /*   conf.pcr0 = PC0_CNT_ADD_SUB_MUL; */
  /*   conf.pcr1 = PC1_CNT_DIV_SQRT; */
  /*   conf.pcr2 = PC2_CNT_DMA_REQ; */
  /*   conf.pcrc = PCRC_ALL; */
  /*   lwpf_init(&conf); */
  /* } */
  /* r++; */
  GPTLstart("attractive athread");
  athread_spawn(pair_tersoff_compute_attractive_para, pm);
  athread_join();
  GPTLstop("attractive athread");
  GPTLstart("attractive reduction bond");
  for (ist = 0; ist < ntotal; ist += ISTEP * 64){
    ied = ist + ISTEP * 64;
    if (ied > ntotal)
      ied = ntotal;
    waitint(numshort + ist);
    for (ii = ist; ii < ied; ii ++){
      int jj, jend;
      jend = numshort[ii] + ii * maxshort;
      for (jj = ii * maxshort; jj < jend; jj += 4){
        int j0 = shortidx[jj + 0];
        int j1 = shortidx[jj + 1];
        int j2 = shortidx[jj + 2];
        int j3 = shortidx[jj + 3];
        doublev4 fend0, fend1, fend2, fend3;
        simd_load(fend0, fend[jj + 0]);
        simd_load(fend1, fend[jj + 1]);
        simd_load(fend2, fend[jj + 2]);
        simd_load(fend3, fend[jj + 3]);
        if (j0 < nlocal){
          doublev4 ftmp_j;
          simd_load(ftmp_j, ftmp[j0]);
          simd_store(fend0 + ftmp_j, ftmp[j0]);
        }

        if (j1 < nlocal && jj + 1 < jend){
          doublev4 ftmp_j;
          simd_load(ftmp_j, ftmp[j1]);
          simd_store(fend1 + ftmp_j, ftmp[j1]);
        }
        if (j2 < nlocal && jj + 2 < jend){
          doublev4 ftmp_j;
          simd_load(ftmp_j, ftmp[j2]);
          simd_store(fend2 + ftmp_j, ftmp[j2]);
        }
        if (j3 < nlocal && jj + 3 < jend){
          doublev4 ftmp_j;
          simd_load(ftmp_j, ftmp[j3]);
          simd_store(fend3 + ftmp_j, ftmp[j3]);
        }
      }
    }
  }
  GPTLstop("attractive reduction bond");

  /* if (r == 10 && pm->rank == 0){ */
  /*   lwpf_finish(stdout); */
  /* } */


  GPTLstart("attractive reduction force");
  athread_spawn(pair_tersoff_reduction_force, pm);
  athread_join();
  GPTLstop("attractive reduction force");
}
#endif
#ifdef CPE

#include <math.h>
#include "poly_math.h"
__thread_local rank;
#define MY_PI2 1.57079632679489661923
#define MY_PI4 0.78539816339744830962

#define inv_sqrt(x, r) {                        \
    double y = x;                               \
    double xhalf = 0.5 * y;                     \
    long i = *(long*)(&y);                      \
    i = 0x5fe6ec85e7de30daLL - (i >> 1);        \
    y = *(double*)(&i);                         \
    y = y * (1.5 - xhalf * y * y );             \
    y = y * (1.5 - xhalf * y * y );             \
    y = y * (1.5 - xhalf * y * y );             \
    r = y;                                      \
  }

void print_4(double *in, char *name, int lineno){
  call_printf("%d: %s=%e %e %e %e\n", lineno, name, in[0], in[1], in[2], in[3]);
}
void print_f(double in, char *name, int lineno){
  call_printf("%d: %s=%e\n", lineno, name, in);
}
#define print_vec(x) print_4(&(x), #x, __LINE__)
#define print_sca(x) print_f(x, #x, __LINE__)

#define HALF 0.5
void ev_tally_global(int i, int j, int nlocal, double evdwl, double fpair,
                     double delx, double dely, double delz,
                     double *eng_vdwl, double *virial)
{
  double v[6];
  double factor = 0;
  if (i < nlocal) factor += HALF;
  if (j < nlocal) factor += HALF;
  *eng_vdwl += evdwl * factor;
  if (_MYID == 0){
    print_sca(factor);
    print_sca(evdwl);
    print_sca(fpair);
  }

  v[0] = delx*delx*fpair;
  v[1] = dely*dely*fpair;
  v[2] = delz*delz*fpair;
  v[3] = delx*dely*fpair;
  v[4] = delx*delz*fpair;
  v[5] = dely*delz*fpair;

  virial[0] += factor*v[0];
  virial[1] += factor*v[1];
  virial[2] += factor*v[2];
  virial[3] += factor*v[3];
  virial[4] += factor*v[4];
  virial[5] += factor*v[5];

}
void ev_tally_full_global(int i, double evdwl, double fpair,
                          double delx, double dely, double delz,
                          double *eng_vdwl, double *virial)
{
  double v[6];
  *eng_vdwl += evdwl * HALF;
  v[0] = delx*delx*fpair;
  v[1] = dely*dely*fpair;
  v[2] = delz*delz*fpair;
  v[3] = delx*dely*fpair;
  v[4] = delx*delz*fpair;
  v[5] = dely*delz*fpair;

  virial[0] += HALF*v[0];
  virial[1] += HALF*v[1];
  virial[2] += HALF*v[2];
  virial[3] += HALF*v[3];
  virial[4] += HALF*v[4];
  virial[5] += HALF*v[5];
}

#define THIRD 0.33333333333333333333333
//v_tally3rd(i, j, k, nlocal, vflag_global, vflag_atom, tfj, tfk, dij, dik, virial, vatom);
void v_tally3rd(int i, int j, int k, int nlocal, int vflag_global, int vflag_atom,
                double *fi, double *fj, double *drik, double *drjk, double *virial, double (*vatom)[6])
{
  double v[6];
  double factor = 0;
  if (i < nlocal) factor += THIRD;
  if (j < nlocal) factor += THIRD;
  if (k < nlocal) factor += THIRD;
  v[0] = factor * (drik[0]*fi[0] + drjk[0]*fj[0]);
  v[1] = factor * (drik[1]*fi[1] + drjk[1]*fj[1]);
  v[2] = factor * (drik[2]*fi[2] + drjk[2]*fj[2]);
  v[3] = factor * (drik[0]*fi[1] + drjk[0]*fj[1]);
  v[4] = factor * (drik[0]*fi[2] + drjk[0]*fj[2]);
  v[5] = factor * (drik[1]*fi[2] + drjk[1]*fj[2]);

  if (vflag_global){
    virial[0] += v[0];
    virial[1] += v[1];
    virial[2] += v[2];
    virial[3] += v[3];
    virial[4] += v[4];
    virial[5] += v[5];
  }
  if (vflag_atom){
    vatom[i][0] += v[0]; vatom[i][1] += v[1]; vatom[i][2] += v[2];
    vatom[i][3] += v[3]; vatom[i][4] += v[4]; vatom[i][5] += v[5];
  }
}

#define HBIT 10
#define SBIT 3
#define LINESIZE (1 << SBIT)
#define LINECNT (1 << HBIT - SBIT)
#define MMASK (LINESIZE - 1)
#define LMASK (LINECNT - 1)
#include <dma.h>
//biga, bigr, bigd, bigdinv, lam1, cutsq
typedef struct repulsive_param_t{
  double biga, bigr, bigd, bigdinv, lam1, cutsq, padding1, padding2;
} repulsive_param_t;
typedef struct tersoff_param3_v4_t{
  doublev4 lam3, bigr, bigd, bigb, bigdinv, c, d, c2divd2;
  doublev4 h, gamma;
  doublev4 power3, skip;
} tersoff_param3_v4_t;
typedef struct tersoff_param3_4_t{
  double lam3[4], bigr[4], bigd[4], bigb[4], bigdinv[4], c[4], d[4], c2divd2[4];
  double h[4], gamma[4], power3[4], skip[4];
} tersoff_param3_4_t;
typedef struct tersoff_thb_param_t{
  double lam3, bigr, bigd, bigb;
  double bigdinv, c2, d2, c2divd2;
  double h, gamma, power3, skip;
} tersoff_thb_param_t;
typedef union tersoff_param3_u{
  tersoff_param3_v4_t vec;
  tersoff_param3_4_t  sca;
} tersoff_param3_u;
typedef struct tersoff_param2_v4_t{
  doublev4 bigr, bigd, bigb, lam2, c1, c4, beta, powern, half_powern_inv, padding;
  doublev4 skip_f, skip_ev;
} tersoff_param2_v4_t;
typedef struct tersoff_param2_4_t{
  double bigr[4], bigd[4], bigb[4];
  double lam2[4], c1[4], c4[4], beta[4], powern[4], half_powern_inv[4], padding[4];
  double skip_f[4], skip_ev[4];
} tersoff_param2_4_t;
typedef union tersoff_param2_u{
  tersoff_param2_v4_t vec;
  tersoff_param2_4_t  sca;
} tersoff_param2_u;
#define simd_bcastf(x) simd_vshff((doublev4)(x), (doublev4)(x), 0)
#define align_ptr(x) (void*)(((long)(x)) + 31 & ~31)
#define transpose4x4_2x4(in0, in1, in2, in3, ot0, ot1) {                \
    doublev4 o0 = simd_vshff(in1,in0,0x44);                             \
    doublev4 o2 = simd_vshff(in3,in2,0x44);                             \
    ot0 = simd_vshff(o2,o0,0x88);                                       \
    ot1 = simd_vshff(o2,o0,0xDD);                                       \
  }

#define transpose4x4(in0, in1, in2, in3, ot0, ot1, ot2, ot3) {  \
    doublev4 o0 = simd_vshff(in1,in0,0x44);                     \
    doublev4 o1 = simd_vshff(in1,in0,0xEE);                     \
    doublev4 o2 = simd_vshff(in3,in2,0x44);                     \
    doublev4 o3 = simd_vshff(in3,in2,0xEE);                     \
    ot0 = simd_vshff(o2,o0,0x88);                               \
    ot1 = simd_vshff(o2,o0,0xDD);                               \
    ot2 = simd_vshff(o3,o1,0x88);                               \
    ot3 = simd_vshff(o3,o1,0xDD);                               \
  }
#define vshuffd_rc(a, b, c, d) (d | (c << 2) | (b << 4) | (a << 6))
#define simd_vsumd(x) {                                 \
    x += simd_vshff(x, x, vshuffd_rc(2, 3, 0, 1));      \
    x += simd_vshff(x, x, vshuffd_rc(1, 0, 3, 2));      \
  } 

//#define LWPF_KERNELS _K(all) K(repulsive) K(repulsive_load) K(repulsive_comp) K(repulsive_trans) K(attractive) K(zeta) K(force_zeta) K(parampack2) K(jboot) K(snwrite) K(jlist_load)
//#define LWPF_UNIT U(TERSOFF)
//#include "lwpf.h"
#define SQR(x) ((x) * (x))

void pair_tersoff_compute_attractive_para(pair_tersoff_compute_param_t *pm){
  pe_init();
  //lwpf_start(all);
  pair_tersoff_compute_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(pair_tersoff_compute_param_t));
  pe_syn();

  rank = l_pm.rank;
  int nlocal = l_pm.nlocal;
  int nghost = l_pm.nghost;
  int ntotal = l_pm.ntotal;
  int inum = l_pm.inum;
  int gnum = l_pm.gnum;
  int ntypes = l_pm.ntypes;
  int rank = l_pm.rank;
  int nelements = l_pm.nelements;
  int nparams = l_pm.nparams;
  int vflag_global = l_pm.vflag_global;
  int vflag_atom = l_pm.vflag_atom;
  int vflag_either = l_pm.vflag_either;
  int *ilist = l_pm.ilist;
  int *firstshort = l_pm.firstshort;
  int maxshort = l_pm.maxshort;
  short_neigh_t *shortlist = l_pm.shortlist;
  int map[ntypes + 1];
  int elem2param[nelements][nelements][nelements];
  tersoff_param_t params[nparams], params3[nelements][nelements][nelements];
  int nep3 = nelements * nelements * nelements;
  tersoff_thb_param_t thb_param_stor[nep3 + 1];

  tersoff_thb_param_t (*thb_param)[nelements][nelements] = align_ptr(&thb_param_stor);
  repulsive_param_t repul_param_stor[nelements * nelements + 1];
  repulsive_param_t (*repul_param)[nelements] = align_ptr(repul_param_stor);

  pe_get(l_pm.map, map, sizeof(int) * (ntypes + 1));
  pe_get(l_pm.elem2param, elem2param, sizeof(int) * nep3);
  pe_get(l_pm.params, params, sizeof(tersoff_param_t) * nparams);
  pe_syn();

  int t1, t2, t3;
  for (t1 = 0; t1 < nelements; t1 ++){
    for (t2 = 0; t2 < nelements; t2 ++){
      for (t3 = 0; t3 < nelements; t3 ++){
        params3[t1][t2][t3] = params[elem2param[t1][t2][t3]];
        thb_param[t1][t2][t3].lam3 = params3[t1][t2][t3].lam3;
        thb_param[t1][t2][t3].bigr = params3[t1][t2][t3].bigr;
        thb_param[t1][t2][t3].bigd = params3[t1][t2][t3].bigd;
        thb_param[t1][t2][t3].bigb = params3[t1][t2][t3].bigb;

        thb_param[t1][t2][t3].bigdinv = params3[t1][t2][t3].bigdinv;
        thb_param[t1][t2][t3].c2      = SQR(params3[t1][t2][t3].c);
        thb_param[t1][t2][t3].d2      = SQR(params3[t1][t2][t3].d);
        thb_param[t1][t2][t3].c2divd2 = params3[t1][t2][t3].c2divd2;

        thb_param[t1][t2][t3].h      = params3[t1][t2][t3].h             ;
        thb_param[t1][t2][t3].gamma  = params3[t1][t2][t3].gamma         ;
        thb_param[t1][t2][t3].power3 = params3[t1][t2][t3].powermint == 3;
        if (t2 == t3){
          repul_param[t1][t2].biga    = params3[t1][t2][t3].biga   ;
          repul_param[t1][t2].bigr    = params3[t1][t2][t3].bigr   ;
          repul_param[t1][t2].bigd    = params3[t1][t2][t3].bigd   ;
          repul_param[t1][t2].bigdinv = params3[t1][t2][t3].bigdinv;
          repul_param[t1][t2].lam1    = params3[t1][t2][t3].lam1   ;
          repul_param[t1][t2].cutsq   = params3[t1][t2][t3].cutsq  ;
        }
      }
    }
  }
  double (*x)[3] = l_pm.x;
  double (*f)[3] = l_pm.f;
  double (*vatom)[6] = l_pm.vatom;
  double *eatom = l_pm.eatom;

  int *type = l_pm.type;
  int ii;
  int ist, ied;
  int *firstneigh_local[ISTEP];
  int numneigh_local[ISTEP];
  int jlist_local[ISTEP];
  int shortidx_local[ISTEP];
  int numshort_local[ISTEP];
  double fi[ISTEP][3], xi[ISTEP][3];
  int ti[ISTEP];
  int *numshort = l_pm.numshort;
  short_neigh_t short_neigh_stor[maxshort + 1];
  short_neigh_t *short_neigh = align_ptr(short_neigh_stor);
  doublev4 eng_virial_v4[2];
  eng_virial_v4[0] = 0.0;
  eng_virial_v4[1] = 0.0;
  double *eng_vdwl = eng_virial_v4;
  double *eng_coul = eng_vdwl + 1;
  double *virial = eng_vdwl + 2;
  doublev4 eng_vdwl_v4 = simd_bcastf(0.0);
  doublev4 virial_v4[6];
  for (ii = 0; ii < 6; ii ++)
    virial_v4[ii] = simd_bcastf(0.0);
  atom_in_t j_cache[LINECNT][LINESIZE];
  int j_tag[LINECNT];
  dma_desc cache_get_desc = 0;
  volatile int cache_reply = 0;
  dma_set_mode(&cache_get_desc, PE_MODE);
  dma_set_size(&cache_get_desc, sizeof(atom_in_t) * LINESIZE);
  dma_set_op(&cache_get_desc, DMA_GET);
  dma_set_reply(&cache_get_desc, &cache_reply);
  doublev4 tf_mask_v4 = simd_set_doublev4(0, 1, 0, 0);
  doublev4 one_v4 = simd_bcastf(1.0);
  doublev4 half_v4 = simd_bcastf(0.5);
  doublev4 zero_v4 = simd_bcastf(0.0);
  doublev4 pi2_v4 = simd_bcastf(MY_PI2);
  doublev4 pi4_v4 = simd_bcastf(MY_PI4);
  for (ii = 0; ii < LINECNT; ii ++)
    j_tag[ii] = -1;
  for (ist = _MYID * ISTEP; ist < ntotal; ist += ISTEP * 64){
    ied = ist + ISTEP;
    if (ied > ntotal)
      ied = ntotal;
    int isz = ied - ist;
    int iwe = ied;
    if (iwe > inum)
      iwe = inum;
    int iwr = iwe - ist;
    int i;
    pe_get(x + ist, xi, sizeof(double) * 3 * isz);
    /* pe_get(numshort + ist, nn, sizeof(int) * isz); */
    pe_get(type + ist, ti, sizeof(int) * isz);
    pe_syn();
    pe_get(l_pm.firstneigh + ist, firstneigh_local, sizeof(int*) * isz);
    pe_get(l_pm.numneigh + ist, numneigh_local, sizeof(int) * isz);
    pe_syn();


    for (i = ist; i < ied; i ++){
      int ioff = i - ist;
      int itype = map[ti[ioff]];
      int *jlist = firstneigh_local[ioff];
      int jnum = numneigh_local[ioff];
      int *idx_short = l_pm.shortidx + i * maxshort;
      int nshort = 0;
      double fc_i_stor[maxshort + 3], fc_d_i_stor[maxshort + 3];
      double *fc_i = align_ptr(fc_i_stor);
      double *fc_d_i = align_ptr(fc_d_i_stor);
      double fxtmp = 0.0;
      double fytmp = 0.0;
      double fztmp = 0.0;
      doublev4 fxtmp_v4 = 0.0;
      doublev4 fytmp_v4 = 0.0;
      doublev4 fztmp_v4 = 0.0;
      doublev4 factor2_base_v4, factor3_base_v4;
      if (i < nlocal)
        factor2_base_v4 = simd_bcastf(0.5);
      else
        factor2_base_v4 = simd_bcastf(0.0);
      if (i < nlocal)
        factor3_base_v4 = simd_bcastf(THIRD);
      else
        factor3_base_v4 = simd_bcastf(0.0);
      int jj, jjj;
      int jst, jed;
      //lwpf_start(repulsive);
      for (jst = 0; jst < jnum; jst += JSTEP){
        //lwpf_start(jlist_load);
        jed = jst + JSTEP;
        if (jed > jnum)
          jed = jnum;
        int jsz = jed - jst;
        pe_get(jlist + jst, jlist_local, sizeof(int) * jsz);
        pe_syn();
        //lwpf_stop(jlist_load);
        doublev4 xj_v4[3], jatom_v4[4], padding_v4;
        int padding[8]; double skip_4[4];
        for (jjj = 0; jjj < jsz; jjj += 4){
          double j_4[4];
          atom_in_t *jatom;
          //lwpf_start(repulsive_load);
          for (jj = 0; jj < 4; jj ++){
            if (jj + jjj < jsz){
              int j = jlist_local[jj + jjj];
              j &= NEIGHMASK;
              j_4[jj] = j;
              int line = j >> SBIT & LMASK;
              int tag = j >> HBIT;
              if (j_tag[line] != tag){
                int mem = j & ~MMASK;
                cache_reply = 0;
                dma_rpl(cache_get_desc, l_pm.atom_in + mem, j_cache[line], cache_reply);
                while (cache_reply != 1);
                j_tag[line] = tag;
              }
              jatom = j_cache[line] + (j & MMASK);
              skip_4[jj] = 0;
            } else {
              jatom = j_cache[0];
              skip_4[jj] = 1;
            }
            simd_load(jatom_v4[jj], jatom);
          }
          //lwpf_stop(repulsive_load);
          //lwpf_start(repulsive_trans);
          transpose4x4(jatom_v4[0], jatom_v4[1], jatom_v4[2], jatom_v4[3],
                       xj_v4[0], xj_v4[1], xj_v4[2], padding_v4);
          simd_store(padding_v4, padding);
          doublev4 skip_v4;
          simd_load(skip_v4, skip_4);
          doublev4 xi_v4[3];
          xi_v4[0] = simd_bcastf(xi[ioff][0]);
          xi_v4[1] = simd_bcastf(xi[ioff][1]);
          xi_v4[2] = simd_bcastf(xi[ioff][2]);
          void *param_ij0 = repul_param[itype] + map[padding[0]];
          void *param_ij1 = repul_param[itype] + map[padding[2]];
          void *param_ij2 = repul_param[itype] + map[padding[4]];
          void *param_ij3 = repul_param[itype] + map[padding[6]];
          doublev4 param_lo0, param_lo1, param_lo2, param_lo3;
          doublev4 param_hi0, param_hi1, param_hi2, param_hi3;
          simd_load(param_lo0, param_ij0);
          simd_load(param_lo1, param_ij1);
          simd_load(param_lo2, param_ij2);
          simd_load(param_lo3, param_ij3);

          simd_load(param_hi0, param_ij0 + 32);
          simd_load(param_hi1, param_ij1 + 32);
          simd_load(param_hi2, param_ij2 + 32);
          simd_load(param_hi3, param_ij3 + 32);
          
          doublev4 biga_v4, bigr_v4, bigd_v4, bigdinv_v4, lam1_v4, cutsq_v4;
          transpose4x4(param_lo0, param_lo1, param_lo2, param_lo3,
                       biga_v4, bigr_v4, bigd_v4, bigdinv_v4);
          transpose4x4_2x4(param_hi0, param_hi1, param_hi2, param_hi3, lam1_v4, cutsq_v4);
          //lwpf_stop(repulsive_trans);
          //lwpf_start(repulsive_comp);
          doublev4 dij_v4[3];
          dij_v4[0] = xi_v4[0] - xj_v4[0];
          dij_v4[1] = xi_v4[1] - xj_v4[1];
          dij_v4[2] = xi_v4[2] - xj_v4[2];
          doublev4 r2ij_v4;
          r2ij_v4 = dij_v4[0] * dij_v4[0] + dij_v4[1] * dij_v4[1] + dij_v4[2] * dij_v4[2];
          doublev4 dcutsq_v4 = r2ij_v4 - cutsq_v4;
          doublev4 r_v4 = simd_vsqrtd(r2ij_v4);
          doublev4 rinv_v4 = one_v4 / r_v4;

          doublev4 angle_v4 = pi2_v4 * (r_v4 - bigr_v4) * bigdinv_v4;
          doublev4 fc_v4 = half_v4 * (one_v4 - simd_vsinnpi_pid(angle_v4));
          doublev4 fc_d_v4 = -(pi4_v4 * bigdinv_v4) * simd_vcosnpi_pid(angle_v4);
          doublev4 rhi_v4 = r_v4 - (bigr_v4 + bigd_v4);
          doublev4 rlo_v4 = r_v4 - (bigr_v4 - bigd_v4);
          fc_v4 = simd_vsellt(rlo_v4, one_v4, fc_v4);
          fc_d_v4 = simd_vsellt(rlo_v4, zero_v4, fc_d_v4);
          fc_v4 = simd_vsellt(rhi_v4, fc_v4, zero_v4);
          fc_d_v4 = simd_vsellt(rhi_v4, fc_d_v4, zero_v4);
          doublev4 tmp_exp_v4 = simd_vexpd(-lam1_v4 * r_v4);
          tmp_exp_v4 = simd_vsellt(dcutsq_v4, tmp_exp_v4, zero_v4);
          doublev4 fpair_v4 = -biga_v4 * tmp_exp_v4 * (fc_d_v4 - fc_v4 * lam1_v4) * rinv_v4;
          doublev4 evdwl_v4 = fc_v4 * biga_v4 * tmp_exp_v4;
          
          if (i < nlocal && l_pm.evflag){
            eng_vdwl_v4 += evdwl_v4 * half_v4;
            doublev4 v_v4[6];
            v_v4[0] = dij_v4[0] * dij_v4[0] * fpair_v4;
            v_v4[1] = dij_v4[1] * dij_v4[1] * fpair_v4;
            v_v4[2] = dij_v4[2] * dij_v4[2] * fpair_v4;
            v_v4[3] = dij_v4[0] * dij_v4[1] * fpair_v4;
            v_v4[4] = dij_v4[0] * dij_v4[2] * fpair_v4;
            v_v4[5] = dij_v4[1] * dij_v4[2] * fpair_v4;
            virial_v4[0] += half_v4 * v_v4[0];
            virial_v4[1] += half_v4 * v_v4[1];
            virial_v4[2] += half_v4 * v_v4[2];
            virial_v4[3] += half_v4 * v_v4[3];
            virial_v4[4] += half_v4 * v_v4[4];
            virial_v4[5] += half_v4 * v_v4[5];
          }
          fxtmp_v4 += fpair_v4 * dij_v4[0];
          fytmp_v4 += fpair_v4 * dij_v4[1];
          fztmp_v4 += fpair_v4 * dij_v4[2];
          //lwpf_start(snwrite);
          double fc_4[4];
          double fc_d_4[4];
          double r2_4[4];
          double rinv_4[4];
          double dij_4[3][4];
          simd_store(fc_v4, fc_4);
          simd_store(fc_d_v4, fc_d_4);
          simd_store(r2ij_v4, r2_4);
          simd_store(rinv_v4, rinv_4);
          simd_store(dij_v4[0], dij_4[0]);
          simd_store(dij_v4[1], dij_4[1]);
          simd_store(dij_v4[2], dij_4[2]);
          //lwpf_stop(repulsive_comp);
          for (jj = 0; jj < 4; jj ++){
            if (jjj + jj < jsz){
              double r2ij = r2_4[jj];
              if (r2ij < l_pm.cutshortsq){
                short_neigh[nshort].idx = j_4[jj];
                short_neigh[nshort].type = map[padding[jj << 1]];
                short_neigh[nshort].d[0] = -dij_4[0][jj];
                short_neigh[nshort].d[1] = -dij_4[1][jj];
                short_neigh[nshort].d[2] = -dij_4[2][jj];
                short_neigh[nshort].r2   = r2ij;
                short_neigh[nshort].padding = 0;
                short_neigh[nshort].rinv = rinv_4[jj];
                shortidx_local[nshort] = j_4[jj];
                fc_i[nshort] = fc_4[jj];
                fc_d_i[nshort] = fc_d_4[jj];
                nshort ++;
              }
            }
          }
          //lwpf_stop(snwrite);
        }
      }
      //lwpf_stop(repulsive);
      numshort_local[ioff] = nshort;
      pe_put(l_pm.shortidx + i * maxshort, shortidx_local, sizeof(int) * nshort);
      int nshortpad = nshort + 3 & ~3;
      for (jj = nshort; jj < nshortpad; jj ++){
        short_neigh[jj].idx = -1;
        short_neigh[jj].type = 0;
        short_neigh[jj].padding = 1;
      }

      doublev4 fend_stor[nshortpad + 1];
      doublev4 *fend_v4 = align_ptr(fend_stor);

      doublev4 v4_0 = 0.0;
      for (jj = 0; jj < nshort; jj ++){
        fend_v4[jj] = v4_0;
      }
      for (jjj = 0; jjj < nshortpad; jjj += 4){
        //lwpf_start(jboot);

        doublev4 dij_v4[3];
        doublev4 r2ij_v4, rij_v4, rijinv_v4;
        doublev4 jshort_lo0, jshort_lo1, jshort_lo2, jshort_lo3;
        doublev4 jshort_hi0, jshort_hi1, jshort_hi2, jshort_hi3;

        simd_load(jshort_lo0, (void*)(short_neigh + jjj + 0));
        simd_load(jshort_lo1, (void*)(short_neigh + jjj + 1));
        simd_load(jshort_lo2, (void*)(short_neigh + jjj + 2));
        simd_load(jshort_lo3, (void*)(short_neigh + jjj + 3));
        simd_load(jshort_hi0, (void*)(short_neigh + jjj + 0) + 32);
        simd_load(jshort_hi1, (void*)(short_neigh + jjj + 1) + 32);
        simd_load(jshort_hi2, (void*)(short_neigh + jjj + 2) + 32);
        simd_load(jshort_hi3, (void*)(short_neigh + jjj + 3) + 32);

        transpose4x4(jshort_lo0, jshort_lo1, jshort_lo2, jshort_lo3,
                     dij_v4[0], dij_v4[1], dij_v4[2], r2ij_v4);

        doublev4 idx_v4, type_v4, padding_v4;
        transpose4x4(jshort_hi0, jshort_hi1, jshort_hi2, jshort_hi3,
                     rijinv_v4, type_v4, idx_v4, padding_v4);
        long jtype_4[4], j_4[4];
        simd_store(type_v4, jtype_4);
        simd_store(idx_v4, j_4);
        rij_v4 = rijinv_v4 * r2ij_v4;

        doublev4 rij_hat_v4[3];
        rij_hat_v4[0] = dij_v4[0] * rijinv_v4;
        rij_hat_v4[1] = dij_v4[1] * rijinv_v4;
        rij_hat_v4[2] = dij_v4[2] * rijinv_v4;

        doublev4 zeta_ij_v4 = v4_0;
        //lwpf_start(parampack2);
        tersoff_param2_u param2_stor;
        tersoff_param2_u *param2 = (void*)(((long)(&param2_stor)) + 31 & ~31);
        for (jj = jjj; jj < jjj + 4; jj ++){
          int jjoff = jj - jjj;
          tersoff_param_t *param_ij = params3[itype][jtype_4[jjoff]] + jtype_4[jjoff];
          param2->sca.lam2           [jjoff] = param_ij->lam2           ;
          param2->sca.bigr           [jjoff] = param_ij->bigr           ;
          param2->sca.bigd           [jjoff] = param_ij->bigd           ;
          param2->sca.bigb           [jjoff] = param_ij->bigb           ;
          param2->sca.c1             [jjoff] = param_ij->c1             ;
          param2->sca.c4             [jjoff] = param_ij->c4             ;
          param2->sca.beta           [jjoff] = param_ij->beta           ;
          param2->sca.powern         [jjoff] = param_ij->powern         ;
          param2->sca.half_powern_inv[jjoff] = param_ij->half_powern_inv;
          param2->sca.skip_f         [jjoff] = j_4[jjoff] == -1         ;
          param2->sca.skip_ev        [jjoff] = j_4[jjoff] == -1 || j_4[jjoff] >= nlocal;
        }
        //lwpf_stop(parampack2);
        int kk;
        doublev4 ex_delr_j_stor[nshort + 1][4];
        doublev4 gijk_j_stor[nshort + 1][4];
        doublev4 gijk_d_j_stor[nshort + 1][4];
        doublev4 (*ex_delr_j_v4) = align_ptr(ex_delr_j_stor);
        doublev4 (*gijk_j_v4) = align_ptr(gijk_j_stor);
        doublev4 (*gijk_d_j_v4) = align_ptr(gijk_d_j_stor);
        tersoff_param3_u param3_stor[nshort + 1];
        tersoff_param3_u *param3 = (void*)(((long)(param3_stor)) + 31 & ~31);

        for (kk = 0; kk < nshort; kk ++){
          short_neigh_t *kshort = short_neigh + kk;
          int ktype = kshort->type;
          double dik[3];
          dik[0] = kshort->d[0];
          dik[1] = kshort->d[1];
          dik[2] = kshort->d[2];
          double r2ik = kshort->r2;
          double rikinv = kshort->rinv;
          double rik = rikinv * r2ik;
          void *param_ij0k = thb_param[itype][jtype_4[0]] + ktype;
          void *param_ij1k = thb_param[itype][jtype_4[1]] + ktype;
          void *param_ij2k = thb_param[itype][jtype_4[2]] + ktype;
          void *param_ij3k = thb_param[itype][jtype_4[3]] + ktype;
          doublev4 param_lo0, param_lo1, param_lo2, param_lo3;
          simd_load(param_lo0, param_ij0k);
          simd_load(param_lo1, param_ij1k);
          simd_load(param_lo2, param_ij2k);
          simd_load(param_lo3, param_ij3k);
          transpose4x4(param_lo0, param_lo1, param_lo2, param_lo3,
                       param3[kk].vec.lam3, param3[kk].vec.bigr,
                       param3[kk].vec.bigd, param3[kk].vec.bigb);
          doublev4 param_md0, param_md1, param_md2, param_md3;
          simd_load(param_md0, param_ij0k + 32);
          simd_load(param_md1, param_ij1k + 32);
          simd_load(param_md2, param_ij2k + 32);
          simd_load(param_md3, param_ij3k + 32);
          transpose4x4(param_md0, param_md1, param_md2, param_md3,
                       param3[kk].vec.bigdinv, param3[kk].vec.c,
                       param3[kk].vec.d, param3[kk].vec.c2divd2);
          doublev4 param_hi0, param_hi1, param_hi2, param_hi3;
          simd_load(param_hi0, param_ij0k + 64);
          simd_load(param_hi1, param_ij1k + 64);
          simd_load(param_hi2, param_ij2k + 64);
          simd_load(param_hi3, param_ij3k + 64);
          doublev4 skip_v4;
          transpose4x4(param_hi0, param_hi1, param_hi2, param_hi3,
                       param3[kk].vec.h, param3[kk].vec.gamma,
                       param3[kk].vec.power3, skip_v4);
          int kkoff = (kk - jjj) << 1;
          doublev4 skipjk_v4 = simd_vshff(tf_mask_v4, tf_mask_v4, (1 << kkoff) & 255);
          asm("vbisw %1, %2, %0": "=r"(param3[kk].vec.skip): "r"(skip_v4), "r"(skipjk_v4));
        }

        //lwpf_stop(jboot);
        //lwpf_start(zeta);
        for (kk = 0; kk < nshort; kk ++){
          short_neigh_t *kshort = short_neigh + kk;
          int ktype = kshort->type;
          double dik[3];
          dik[0] = kshort->d[0];
          dik[1] = kshort->d[1];
          dik[2] = kshort->d[2];
          double r2ik = kshort->r2;
          double rikinv = kshort->rinv;
          double rik = rikinv * r2ik;
          doublev4 dik_v4[3];
          dik_v4[0] = simd_bcastf(kshort->d[0]);
          dik_v4[1] = simd_bcastf(kshort->d[1]);
          dik_v4[2] = simd_bcastf(kshort->d[2]);
          doublev4 r2ik_v4 = simd_bcastf(kshort->r2);
          doublev4 rikinv_v4 = simd_bcastf(kshort->rinv);
          doublev4 rik_v4 = r2ik_v4 * rikinv_v4;
          doublev4 costheta_v4 = (dij_v4[0] * dik_v4[0] + dij_v4[1] * dik_v4[1] +
                                  dij_v4[2] * dik_v4[2]) * (rijinv_v4 * rikinv_v4);
          doublev4 arg_v4 = param3[kk].vec.lam3 * (rij_v4 - rik_v4);
          doublev4 argcb_v4 = arg_v4 * arg_v4 * arg_v4;
          arg_v4 = simd_vseleq(param3[kk].vec.power3, arg_v4, argcb_v4);
          ex_delr_j_v4[kk] = simd_vexpd(arg_v4);
          doublev4 ters_R_v4  = param3[kk].vec.bigr;
          doublev4 ters_D_v4  = param3[kk].vec.bigd;
          doublev4 ters_c_v4  = param3[kk].vec.c;
          doublev4 ters_d_v4  = param3[kk].vec.d;
          doublev4 hcth_v4    = param3[kk].vec.h - costheta_v4;
          doublev4 gamma_v4   = param3[kk].vec.gamma;
          doublev4 c2divd2_v4 = param3[kk].vec.c2divd2;
          doublev4 numerator_v4 = simd_bcastf(-2.0) * ters_c_v4 * hcth_v4;
          doublev4 denominator_v4 = simd_bcastf(1.0) / (ters_d_v4 + hcth_v4 * hcth_v4);
          gijk_j_v4[kk] = 
            gamma_v4*(simd_bcastf(1.0) + c2divd2_v4 - ters_c_v4*denominator_v4);
          gijk_d_j_v4[kk] = gamma_v4*numerator_v4*denominator_v4*denominator_v4;
          //simd_load(zeta_ij_v4, zeta_ij_4);
          doublev4 fc_v4 = simd_bcastf(fc_i[kk]);
          fc_v4 = simd_vseleq(param3[kk].vec.skip, fc_v4, simd_bcastf(0.0));
          zeta_ij_v4 = zeta_ij_v4 + fc_v4 * gijk_j_v4[kk] * ex_delr_j_v4[kk];
          //simd_store(zeta_ij_v4, zeta_ij_4);
          //simd_store(ex_delr_j_v4, ex_delr_j[kk]);
          //simd_store(gijk_j_v4, gijk_j[kk]);
          //simd_store(gijk_d_j_v4, gijk_d_j[kk]);
        }
        //lwpf_stop(zeta);
        //lwpf_start(force_zeta);
        double prefactor[4];
        doublev4 ters_R_v4 = param2->vec.bigr;
        doublev4 ters_D_v4 = param2->vec.bigd;
        doublev4 ters_B_v4 = param2->vec.bigb;
        doublev4 ters_lam2_v4 = param2->vec.lam2;
        doublev4 er_v4 = simd_vexpd(-ters_lam2_v4 * rij_v4);
        doublev4 fc_ij_v4;
        simd_load(fc_ij_v4, fc_i + jjj);
        doublev4 fc_d_ij_v4;
        simd_load(fc_d_ij_v4, fc_d_i + jjj);
        doublev4 fa_v4 = -ters_B_v4 * er_v4 * fc_ij_v4;
        doublev4 fa_d_v4 = ters_B_v4 * er_v4 * (ters_lam2_v4 * fc_ij_v4 - fc_d_ij_v4);
        doublev4 beta_v4 = param2->vec.beta;
        doublev4 tmp_v4 = beta_v4 * zeta_ij_v4;
        doublev4 powern_v4 = param2->vec.powern;
        doublev4 half_powern_inv_v4 = param2->vec.half_powern_inv;
        doublev4 c1_v4 = param2->vec.c1;
        doublev4 c4_v4 = param2->vec.c4;
        doublev4 tmp_n_v4 = simd_vpowd(tmp_v4, powern_v4);
        doublev4 tmp_n_p1_v4 = simd_bcastf(1.0) + tmp_n_v4;
        doublev4 bij_v4 = simd_vpowd(tmp_n_p1_v4, -half_powern_inv_v4);
        doublev4 bij_d_v4 = simd_bcastf(-0.5) * bij_v4 * tmp_n_v4 / (zeta_ij_v4 * tmp_n_p1_v4);
        doublev4 zero_v4 = simd_bcastf(0.0);
        bij_v4 = simd_vseleq(param2->vec.skip_f, bij_v4, zero_v4);
        doublev4 fpair_v4 = simd_bcastf(0.5) * bij_v4 * fa_d_v4 * rijinv_v4;
        doublev4 prefactor_v4 = simd_bcastf(-0.5) * fa_v4 * bij_d_v4;
        doublev4 evdwl_v4 = 0.5 * bij_v4 * fa_v4;

        doublev4 fij_v4[3];
        fij_v4[0] = dij_v4[0] * fpair_v4;
        fij_v4[1] = dij_v4[1] * fpair_v4;
        fij_v4[2] = dij_v4[2] * fpair_v4;
        doublev4 fend_ij_v4[4], pad_v4;
        transpose4x4(fij_v4[0], fij_v4[1], fij_v4[2], zero_v4,
                     fend_ij_v4[0], fend_ij_v4[1], fend_ij_v4[2], fend_ij_v4[3]);
        fend_v4[jjj + 0] -= fend_ij_v4[0];
        fend_v4[jjj + 1] -= fend_ij_v4[1];
        fend_v4[jjj + 2] -= fend_ij_v4[2];
        fend_v4[jjj + 3] -= fend_ij_v4[3];

        fxtmp_v4 += fij_v4[0];
        fytmp_v4 += fij_v4[1];
        fztmp_v4 += fij_v4[2];
        if (l_pm.evflag){
          doublev4 half_v4 = simd_bcastf(0.5);
          doublev4 factor_v4 = simd_vseleq(param2->vec.skip_ev, half_v4, zero_v4) + factor2_base_v4;
          eng_vdwl_v4 += evdwl_v4 * factor_v4;
          doublev4 v_v4[6];
          v_v4[0] = -dij_v4[0] * dij_v4[0] * fpair_v4;
          v_v4[1] = -dij_v4[1] * dij_v4[1] * fpair_v4;
          v_v4[2] = -dij_v4[2] * dij_v4[2] * fpair_v4;
          v_v4[3] = -dij_v4[0] * dij_v4[1] * fpair_v4;
          v_v4[4] = -dij_v4[0] * dij_v4[2] * fpair_v4;
          v_v4[5] = -dij_v4[1] * dij_v4[2] * fpair_v4;

          virial_v4[0] += factor_v4 * v_v4[0];
          virial_v4[1] += factor_v4 * v_v4[1];
          virial_v4[2] += factor_v4 * v_v4[2];
          virial_v4[3] += factor_v4 * v_v4[3];
          virial_v4[4] += factor_v4 * v_v4[4];
          virial_v4[5] += factor_v4 * v_v4[5];
        }
        //lwpf_stop(force_zeta);
        //lwpf_start(attractive);
        for (kk = 0; kk < nshort; kk ++){
          short_neigh_t *kshort = short_neigh + kk;
          int k = kshort->idx;
          doublev4 dik_v4[3];
          dik_v4[0] = simd_bcastf(kshort->d[0]);
          dik_v4[1] = simd_bcastf(kshort->d[1]);
          dik_v4[2] = simd_bcastf(kshort->d[2]);
          doublev4 r2ik_v4 = simd_bcastf(kshort->r2);
          doublev4 rikinv_v4 = simd_bcastf(kshort->rinv);
          doublev4 rik_v4 = r2ik_v4 * rikinv_v4;
          doublev4 rik_hat_v4[3];
          rik_hat_v4[0] = rikinv_v4 * dik_v4[0];
          rik_hat_v4[1] = rikinv_v4 * dik_v4[1];
          rik_hat_v4[2] = rikinv_v4 * dik_v4[2];
          doublev4 ters_R_v4 = param3[kk].vec.bigr;
          doublev4 ters_D_v4 = param3[kk].vec.bigd;
          doublev4 Dinv_v4   = param3[kk].vec.bigdinv;
          doublev4 ters_lam3_v4 = param3[kk].vec.lam3;
          doublev4 ex_delr_v4 = ex_delr_j_v4[kk];
          //simd_load(ex_delr_v4, ex_delr_j[kk]);
          doublev4 ex_delr_d_v4 = ters_lam3_v4 * ex_delr_v4;
          doublev4 ex_delr_d_cu_v4 = simd_bcastf(3.0) * 
            ters_lam3_v4 * ters_lam3_v4 * ters_lam3_v4 *
            (rij_v4 - rik_v4) * (rij_v4 - rik_v4) * ex_delr_v4;
          ex_delr_d_v4 = simd_vseleq(param3[kk].vec.power3, ex_delr_d_v4, ex_delr_d_cu_v4);
          doublev4 cos_theta_v4 = 
            rij_hat_v4[0] * rik_hat_v4[0] + 
            rij_hat_v4[1] * rik_hat_v4[1] +
            rij_hat_v4[2] * rik_hat_v4[2];
          doublev4 dcosdrj_v4[3], dcosdrk_v4[3];
          dcosdrj_v4[0] = (-cos_theta_v4 * rij_hat_v4[0] + rik_hat_v4[0]) * rijinv_v4;
          dcosdrj_v4[1] = (-cos_theta_v4 * rij_hat_v4[1] + rik_hat_v4[1]) * rijinv_v4;
          dcosdrj_v4[2] = (-cos_theta_v4 * rij_hat_v4[2] + rik_hat_v4[2]) * rijinv_v4;

          dcosdrk_v4[0] = (-cos_theta_v4 * rik_hat_v4[0] + rij_hat_v4[0]) * rikinv_v4;
          dcosdrk_v4[1] = (-cos_theta_v4 * rik_hat_v4[1] + rij_hat_v4[1]) * rikinv_v4;
          dcosdrk_v4[2] = (-cos_theta_v4 * rik_hat_v4[2] + rij_hat_v4[2]) * rikinv_v4;

          doublev4 gijk_d_v4 = gijk_d_j_v4[kk];
          doublev4 gijk_v4 = gijk_j_v4[kk];
          doublev4 fc_v4 = simd_bcastf(fc_i[kk]);
          doublev4 fc_d_v4 = simd_bcastf(fc_d_i[kk]);
          doublev4 tfj_v4[3], tfk_v4[3], tfi_v4[3];
          doublev4 lastfac_v4 = simd_vseleq(param3[kk].vec.skip, prefactor_v4, zero_v4);

          doublev4 tmp = gijk_d_v4 * ex_delr_v4;
          tfj_v4[0]  = tmp * dcosdrj_v4[0];
          tfj_v4[1]  = tmp * dcosdrj_v4[1];
          tfj_v4[2]  = tmp * dcosdrj_v4[2];

          tfj_v4[0] += gijk_v4 * ex_delr_d_v4 * rij_hat_v4[0];
          tfj_v4[1] += gijk_v4 * ex_delr_d_v4 * rij_hat_v4[1];
          tfj_v4[2] += gijk_v4 * ex_delr_d_v4 * rij_hat_v4[2];

          tfj_v4[0] *= fc_v4 * lastfac_v4;
          tfj_v4[1] *= fc_v4 * lastfac_v4;
          tfj_v4[2] *= fc_v4 * lastfac_v4;

          tfk_v4[0]  = gijk_v4*rik_hat_v4[0] * (fc_d_v4 * ex_delr_v4 - fc_v4 * ex_delr_d_v4);
          tfk_v4[1]  = gijk_v4*rik_hat_v4[1] * (fc_d_v4 * ex_delr_v4 - fc_v4 * ex_delr_d_v4);
          tfk_v4[2]  = gijk_v4*rik_hat_v4[2] * (fc_d_v4 * ex_delr_v4 - fc_v4 * ex_delr_d_v4);

          tfk_v4[0] += fc_v4 * gijk_d_v4 * ex_delr_v4 * dcosdrk_v4[0];
          tfk_v4[1] += fc_v4 * gijk_d_v4 * ex_delr_v4 * dcosdrk_v4[1];
          tfk_v4[2] += fc_v4 * gijk_d_v4 * ex_delr_v4 * dcosdrk_v4[2];

          tfk_v4[0] *= lastfac_v4;
          tfk_v4[1] *= lastfac_v4;
          tfk_v4[2] *= lastfac_v4;
          
          tfi_v4[0] = -tfj_v4[0] - tfk_v4[0];
          tfi_v4[1] = -tfj_v4[1] - tfk_v4[1];
          tfi_v4[2] = -tfj_v4[2] - tfk_v4[2];

          doublev4 fend_ij_v4[4], pad_v4;
          transpose4x4(tfj_v4[0], tfj_v4[1], tfj_v4[2], zero_v4,
                       fend_ij_v4[0], fend_ij_v4[1], fend_ij_v4[2], fend_ij_v4[3]);
          doublev4 fend_jj_v4[4];
          fend_v4[jjj + 0] += fend_ij_v4[0];
          fend_v4[jjj + 1] += fend_ij_v4[1];
          fend_v4[jjj + 2] += fend_ij_v4[2];
          fend_v4[jjj + 3] += fend_ij_v4[3];

          doublev4 fend_ik_v4[4];
          transpose4x4(tfk_v4[0], tfk_v4[1], tfk_v4[2], zero_v4,
                       fend_ik_v4[0], fend_ik_v4[1], fend_ik_v4[2], fend_ik_v4[3]);
          doublev4 fend_kk_v4;
          fend_v4[kk] += fend_ik_v4[0] + fend_ik_v4[1] + fend_ik_v4[2] + fend_ik_v4[3];

          fxtmp_v4 += tfi_v4[0];
          fytmp_v4 += tfi_v4[1];
          fztmp_v4 += tfi_v4[2];
          
          doublev4 factor3_v4;
          doublev4 zero_v4 = simd_bcastf(0.0);
          doublev4 third_v4 = simd_bcastf(THIRD);
          if (k < nlocal)
            factor3_v4 = factor3_base_v4 + simd_bcastf(THIRD);
          else
            factor3_v4 = factor3_base_v4;
          factor3_v4 += simd_vseleq(param2->vec.skip_ev, third_v4, zero_v4);
          if (l_pm.evflag){
            doublev4 v_v4[6];
            v_v4[0] = dij_v4[0] * tfj_v4[0] + dik_v4[0] * tfk_v4[0];
            v_v4[1] = dij_v4[1] * tfj_v4[1] + dik_v4[1] * tfk_v4[1];
            v_v4[2] = dij_v4[2] * tfj_v4[2] + dik_v4[2] * tfk_v4[2];
            v_v4[3] = dij_v4[0] * tfj_v4[1] + dik_v4[0] * tfk_v4[1];
            v_v4[4] = dij_v4[0] * tfj_v4[2] + dik_v4[0] * tfk_v4[2];
            v_v4[5] = dij_v4[1] * tfj_v4[2] + dik_v4[1] * tfk_v4[2];

            virial_v4[0] += factor3_v4 * v_v4[0];
            virial_v4[1] += factor3_v4 * v_v4[1];
            virial_v4[2] += factor3_v4 * v_v4[2];
            virial_v4[3] += factor3_v4 * v_v4[3];
            virial_v4[4] += factor3_v4 * v_v4[4];
            virial_v4[5] += factor3_v4 * v_v4[5];
          }
        }
        //lwpf_stop(attractive);
      }
      pe_put(l_pm.fend[i * maxshort], fend_v4, sizeof(doublev4) * nshort);
      simd_vsumd(fxtmp_v4);
      simd_vsumd(fytmp_v4);
      simd_vsumd(fztmp_v4);
      fi[ioff][0] = fxtmp + fxtmp_v4;
      fi[ioff][1] = fytmp + fytmp_v4;
      fi[ioff][2] = fztmp + fztmp_v4;
      pe_syn();
    }
    if (iwr > 0){
      pe_put(f + ist, fi, sizeof(double) * iwr * 3);
    }
    pe_put(l_pm.numshort + ist, numshort_local, sizeof(int) * isz);
    pe_syn();
  }

  simd_vsumd(eng_vdwl_v4);
  *eng_vdwl += eng_vdwl_v4;
  for (ii = 0; ii < 6; ii ++){
    simd_vsumd(virial_v4[ii]);
    virial[ii] += virial_v4[ii];
  }
  reg_reduce_inplace_doublev4(eng_virial_v4, 2);

  if (_MYID == 0){
    pe_put(&(pm->eng_vdwl), eng_virial_v4, sizeof(double) * 8);
    pe_syn();
  }
  //lwpf_stop(all);
}

void pair_tersoff_a2s(pair_tersoff_compute_param_t *pm){
  pe_init();
  pair_tersoff_compute_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(pair_tersoff_compute_param_t));
  pe_syn();
  double x[ISTEP][3];
  int type[ISTEP];
  atom_in_t atom_in[ISTEP];
  int n1s[ISTEP];
  doublev4 zeros[ISTEP];
  int ipage_start, i;
  int inum = l_pm.nlocal + l_pm.nghost;
  intv8 v8_n1 = -1;
  doublev4 v4_0 = 0;
  for (i = 0; i < ISTEP; i += 8){
    simd_store(v8_n1, n1s + i);
  }
  for (ipage_start = ISTEP * _MYID; ipage_start < inum; ipage_start += ISTEP * 64){
    int ipage_end = ipage_start + ISTEP;
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
    pe_put(l_pm.numshort + ipage_start, n1s, sizeof(int) * ipage_size);
    //pe_put(l_pm.ftmp[ipage_start], zeros, sizeof(doublev4) * ipage_size);
    pe_syn();
  }
}

void pair_tersoff_reduction_force(pair_tersoff_compute_param_t *pm){
  pe_init();
  pair_tersoff_compute_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(pair_tersoff_compute_param_t));
  pe_syn();
  double f[ISTEP][3];
  double ftmp[ISTEP][4];
  int ist, ied;
  for (ist = _MYID * ISTEP; ist < l_pm.nlocal; ist += ISTEP * 64){
    ied = ist + ISTEP;
    if (ied > l_pm.nlocal)
      ied = l_pm.nlocal;
    int isz = ied - ist;
    pe_get(l_pm.f[ist], f, sizeof(double) * 3 * isz);
    pe_get(l_pm.ftmp[ist], ftmp, sizeof(double) * 4 * isz);
    pe_syn();
    int ii;
    for (ii = 0; ii < isz; ii ++){
      f[ii][0] += ftmp[ii][0];
      f[ii][1] += ftmp[ii][1];
      f[ii][2] += ftmp[ii][2];
    }
    pe_put(l_pm.f[ist], f, sizeof(double) * 3 * isz);
    pe_syn();
  }
}

#endif
