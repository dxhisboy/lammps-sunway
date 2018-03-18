#include "sunway.h"
#include "pair_tersoff_sw64.h"
#include <stdio.h>
#include <stdlib.h>
#include "gptl.h"
#define ISTEP 128
#define JSTEP 64

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

  GPTLstop("tersoff fill");
  athread_join();
  GPTLstop("tersoff a2s");
  GPTLstart("attractive athread");
  athread_spawn(pair_tersoff_compute_attractive_para, pm);
  GPTLstart("attractive reduction bond");
  for (ii = 0; ii < nlocal; ii += 4){
    simd_store(v4_0, ftmp[ii]);
    simd_store(v4_0, ftmp[ii + 1]);
    simd_store(v4_0, ftmp[ii + 2]);
    simd_store(v4_0, ftmp[ii + 3]);
  }
  for (ist = 0; ist < ntotal; ist += ISTEP){
    ied = ist + ISTEP;
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

  athread_join();
  GPTLstop("attractive athread");

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
typedef struct tersoff_param3_v4_t{
  doublev4 lam3, bigr, bigd, bigb, bigdinv, c, d, c2divd2;
  doublev4 h, gamma;
  doublev4 power3, skip;
} tersoff_param3_v4_t;
typedef struct tersoff_param3_4_t{
  double lam3[4], bigr[4], bigd[4], bigb[4], bigdinv[4], c[4], d[4], c2divd2[4];
  double h[4], gamma[4], power3[4], skip[4];
} tersoff_param3_4_t;
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

void pair_tersoff_compute_attractive_para(pair_tersoff_compute_param_t *pm){
  pe_init();
  //lwpf_start(ALL);
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
  int mask_swap_hi_lo = 0x10325476;
  int map[ntypes + 1];
  int elem2param[nelements][nelements][nelements];
  tersoff_param_t params[nparams], params3[nelements][nelements][nelements];
  int nep3 = nelements * nelements * nelements;
  pe_get(l_pm.map, map, sizeof(int) * (ntypes + 1));
  pe_get(l_pm.elem2param, elem2param, sizeof(int) * nep3);
  pe_get(l_pm.params, params, sizeof(tersoff_param_t) * nparams);
  pe_syn();

  int t1, t2, t3;
  for (t1 = 0; t1 < nelements; t1 ++){
    for (t2 = 0; t2 < nelements; t2 ++){
      for (t3 = 0; t3 < nelements; t3 ++){
        params3[t1][t2][t3] = params[elem2param[t1][t2][t3]];
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
  int nn[ISTEP];
  //short_neigh_t js[maxshort];
  short_neigh_t short_neigh[maxshort];
  double fend_stor[maxshort][4];
  double (*fend)[4] = (void*)(((long)fend_stor) + 31 & ~31);
  int fdone[ISTEP];
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

  for (ii = 0; ii < LINECNT; ii ++)
    j_tag[ii] = -1;
  for (ii = 0; ii < ISTEP; ii ++)
    fdone[ii] = 1;
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
    pe_syn();
    pe_get(numshort + ist, nn, sizeof(int) * isz);
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
      int jj;
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
      doublev4 factor_base_v4;
      if (i < nlocal)
        factor_base_v4 = simd_bcastf(0.5);
      else
        factor_base_v4 = simd_bcastf(0.0);
      int jst, jed;
      for (jst = 0; jst < jnum; jst += JSTEP){
        jed = jst + JSTEP;
        if (jed > jnum)
          jed = jnum;
        int jsz = jed - jst;
        pe_get(jlist + jst, jlist_local, sizeof(int) * jsz);
        pe_syn();
        for (jj = 0; jj < jsz; jj ++){
          int j = jlist_local[jj];
          j &= NEIGHMASK;
          int line = j >> SBIT & LMASK;
          int tag = j >> HBIT;
          if (j_tag[line] != tag){
            int mem = j & ~MMASK;
            cache_reply = 0;
            dma_rpl(cache_get_desc, l_pm.atom_in + mem, j_cache[line], cache_reply);
            while (cache_reply != 1);
            j_tag[line] = tag;
          }
          atom_in_t *jatom = j_cache[line] + (j & MMASK);
          int jtype = map[jatom->type];
          double dij[3];
          dij[0] = xi[ioff][0] - jatom->x[0];
          dij[1] = xi[ioff][1] - jatom->x[1];
          dij[2] = xi[ioff][2] - jatom->x[2];
          double r2ij = dij[0] * dij[0] + dij[1] * dij[1] + dij[2] * dij[2];
          
          tersoff_param_t *param_ij = params3[itype][jtype] + jtype;
          if (r2ij > param_ij->cutsq) continue;
          double fpair, evdwl;
          double tmp_exp;
          double ters_R = param_ij->bigr;
          double ters_D = param_ij->bigd;
          double Dinv = param_ij->bigdinv;
          double r, rinv;
          inv_sqrt(r2ij, rinv);
          r = rinv * r2ij;

          if (r < ters_R-ters_D)
            {
              fc_i[nshort] = 1.0;
              fc_d_i[nshort] = 0.0;
            }
          else if (r > ters_R+ters_D)
            {
              fc_i[nshort] =  0.0;
              fc_d_i[nshort] =  0.0;
            }
          else
            {
              fc_i[nshort] = 0.5*(1.0 - p_sind(MY_PI2*(r - ters_R) * Dinv));
              fc_d_i[nshort] = -(MY_PI4 * Dinv) * p_cosd(MY_PI2*(r - ters_R) * Dinv);
            }

          tmp_exp = p_expd(-param_ij->lam1 * r);
          fpair = -param_ij->biga * tmp_exp * 
            (fc_d_i[nshort] - fc_i[nshort]*param_ij->lam1) * rinv;
          evdwl = fc_i[nshort] * param_ij->biga * tmp_exp;
          
          if (r2ij < l_pm.cutshortsq){
            short_neigh[nshort].idx = j;
            short_neigh[nshort].type = jtype;
            short_neigh[nshort].d[0] = -dij[0];
            short_neigh[nshort].d[1] = -dij[1];
            short_neigh[nshort].d[2] = -dij[2];
            short_neigh[nshort].r2   = r2ij;
            inv_sqrt(r2ij, short_neigh[nshort].rinv);
            shortidx_local[nshort] = j;
            nshort ++;
          }

          if (i < nlocal){
            fxtmp += dij[0] * fpair;
            fytmp += dij[1] * fpair;
            fztmp += dij[2] * fpair;
            if (l_pm.evflag) ev_tally_full_global(i, evdwl, fpair, dij[0], dij[1], dij[2], eng_vdwl, virial);
          }
        }
      }
      pe_put(l_pm.shortidx + i * maxshort, shortidx_local, sizeof(int) * nshort);
      numshort_local[ioff] = nshort;
      int nshortpad = nshort + 3 & ~3;
      for (jj = nshort; jj < nshortpad; jj ++){
        short_neigh[jj].idx = -1;
        short_neigh[jj].type = 0;
      }
      doublev4 v4_0 = 0.0;
      for (jj = 0; jj < nshort; jj ++){
        simd_store(v4_0, fend[jj]);
      }
      int jjj;
      for (jjj = 0; jjj < nshortpad; jjj += 4){
        long jtype_4[4];
        jtype_4[0] = short_neigh[jjj + 0].type;
        jtype_4[1] = short_neigh[jjj + 1].type;
        jtype_4[2] = short_neigh[jjj + 2].type;
        jtype_4[3] = short_neigh[jjj + 3].type;
        long j_4[4];
        j_4[0] = short_neigh[jjj + 0].idx;
        j_4[1] = short_neigh[jjj + 1].idx;
        j_4[2] = short_neigh[jjj + 2].idx;
        j_4[3] = short_neigh[jjj + 3].idx;
        double dij_4[3][4];
        dij_4[0][0] = short_neigh[jjj + 0].d[0];
        dij_4[0][1] = short_neigh[jjj + 1].d[0];
        dij_4[0][2] = short_neigh[jjj + 2].d[0];
        dij_4[0][3] = short_neigh[jjj + 3].d[0];

        dij_4[1][0] = short_neigh[jjj + 0].d[1];
        dij_4[1][1] = short_neigh[jjj + 1].d[1];
        dij_4[1][2] = short_neigh[jjj + 2].d[1];
        dij_4[1][3] = short_neigh[jjj + 3].d[1];

        dij_4[2][0] = short_neigh[jjj + 0].d[2];
        dij_4[2][1] = short_neigh[jjj + 1].d[2];
        dij_4[2][2] = short_neigh[jjj + 2].d[2];
        dij_4[2][3] = short_neigh[jjj + 3].d[2];

        double r2_4[4];
        r2_4[0] = short_neigh[jjj + 0].r2;
        r2_4[1] = short_neigh[jjj + 1].r2;
        r2_4[2] = short_neigh[jjj + 2].r2;
        r2_4[3] = short_neigh[jjj + 3].r2;
        
        double rinv_4[4];
        rinv_4[0] = short_neigh[jjj + 0].rinv;
        rinv_4[1] = short_neigh[jjj + 1].rinv;
        rinv_4[2] = short_neigh[jjj + 2].rinv;
        rinv_4[3] = short_neigh[jjj + 3].rinv;
        double rij_4[4];
        rij_4[0] = r2_4[0] * rinv_4[0];
        rij_4[1] = r2_4[1] * rinv_4[1];
        rij_4[2] = r2_4[2] * rinv_4[2];
        rij_4[3] = r2_4[3] * rinv_4[3];
        double zeta_ij_4[4];
        doublev4 zeta_ij_v4;
        zeta_ij_4[0] = 0.0;
        zeta_ij_4[1] = 0.0;
        zeta_ij_4[2] = 0.0;
        zeta_ij_4[3] = 0.0;

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
        int kk;
        double ex_delr_j_stor[nshort + 1][4];
        double gijk_j_stor[nshort + 1][4];
        double gijk_d_j_stor[nshort + 1][4];
        double (*ex_delr_j)[4] = align_ptr(ex_delr_j_stor);
        double (*gijk_j)[4] = align_ptr(gijk_j_stor);
        double (*gijk_d_j)[4] = align_ptr(gijk_d_j_stor);
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
          for (jj = jjj; jj < jjj + 4; jj ++){
            int jjoff = jj - jjj;
            tersoff_param_t *param_ijk = params3[itype][jtype_4[jjoff]] + ktype;
            param3[kk].sca.lam3           [jjoff] = param_ijk->lam3           ;
            param3[kk].sca.bigr           [jjoff] = param_ijk->bigr           ;
            param3[kk].sca.bigd           [jjoff] = param_ijk->bigd           ;
            param3[kk].sca.bigdinv        [jjoff] = param_ijk->bigdinv        ;
            param3[kk].sca.c              [jjoff] = param_ijk->c*param_ijk->c ;
            param3[kk].sca.d              [jjoff] = param_ijk->d*param_ijk->d ;
            param3[kk].sca.c2divd2        [jjoff] = param_ijk->c2divd2        ;
            param3[kk].sca.h              [jjoff] = param_ijk->h              ;
            param3[kk].sca.gamma          [jjoff] = param_ijk->gamma          ;
            param3[kk].sca.power3         [jjoff] = param_ijk->powermint == 3 ;
            param3[kk].sca.skip           [jjoff] = jj == kk || j_4[jjoff] == -1;
          }
        }
        doublev4 dij_v4[3];
        simd_load(dij_v4[0], dij_4[0]);
        simd_load(dij_v4[1], dij_4[1]);
        simd_load(dij_v4[2], dij_4[2]);
        doublev4 r2ij_v4, rij_v4, rijinv_v4;
        simd_load(r2ij_v4, r2_4);
        simd_load(rij_v4, rij_4);
        simd_load(rijinv_v4, rinv_4);

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
          doublev4 ex_delr_j_v4 = simd_vexpd(arg_v4);
          doublev4 ters_R_v4  = param3[kk].vec.bigr;
          doublev4 ters_D_v4  = param3[kk].vec.bigd;
          doublev4 ters_c_v4  = param3[kk].vec.c;
          doublev4 ters_d_v4  = param3[kk].vec.d;
          doublev4 hcth_v4    = param3[kk].vec.h - costheta_v4;
          doublev4 gamma_v4   = param3[kk].vec.gamma;
          doublev4 c2divd2_v4 = param3[kk].vec.c2divd2;
          doublev4 numerator_v4 = simd_bcastf(-2.0) * ters_c_v4 * hcth_v4;
          doublev4 denominator_v4 = simd_bcastf(1.0) / (ters_d_v4 + hcth_v4 * hcth_v4);
          doublev4 
            gijk_j_v4 = gamma_v4*(simd_bcastf(1.0) + c2divd2_v4 - ters_c_v4*denominator_v4);
          doublev4 gijk_d_j_v4 = gamma_v4*numerator_v4*denominator_v4*denominator_v4;
          simd_load(zeta_ij_v4, zeta_ij_4);
          doublev4 fc_v4 = simd_bcastf(fc_i[kk]);
          fc_v4 = simd_vseleq(param3[kk].vec.skip, fc_v4, simd_bcastf(0.0));
          zeta_ij_v4 += fc_v4 * gijk_j_v4 * ex_delr_j_v4;
          simd_store(zeta_ij_v4, zeta_ij_4);
          simd_store(ex_delr_j_v4, ex_delr_j[kk]);
          simd_store(gijk_j_v4, gijk_j[kk]);
          simd_store(gijk_d_j_v4, gijk_d_j[kk]);
        }
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
        double fa_ij_4[4], fa_d_ij_4[4];
        simd_store(fa_v4, fa_ij_4);
        simd_store(fa_d_v4, fa_d_ij_4);
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
        double bij_4[4], bij_d_4[4];
        simd_store(bij_v4, bij_4);
        simd_store(bij_d_v4, bij_d_4);
        doublev4 zero_v4 = simd_bcastf(0.0);
        bij_v4 = simd_vseleq(param2->vec.skip_f, bij_v4, zero_v4);
        doublev4 fpair_v4 = simd_bcastf(0.5) * bij_v4 * fa_d_v4 * rijinv_v4;
        doublev4 prefactor_v4 = simd_bcastf(-0.5) * fa_v4 * bij_d_v4;
        doublev4 evdwl_v4 = 0.5 * bij_v4 * fa_v4;
        simd_store(prefactor_v4, prefactor);
        doublev4 fij_v4[3];
        fij_v4[0] = dij_v4[0] * fpair_v4;
        fij_v4[1] = dij_v4[1] * fpair_v4;
        fij_v4[2] = dij_v4[2] * fpair_v4;
        doublev4 fend_ij_v4[4], pad_v4;
        transpose4x4(fij_v4[0], fij_v4[1], fij_v4[2], zero_v4,
                     fend_ij_v4[0], fend_ij_v4[1], fend_ij_v4[2], fend_ij_v4[3]);
        doublev4 fend_v4[4];
        simd_load(fend_v4[0], fend[jjj + 0]);
        simd_load(fend_v4[1], fend[jjj + 1]);
        simd_load(fend_v4[2], fend[jjj + 2]);
        simd_load(fend_v4[3], fend[jjj + 3]);
        fend_v4[0] -= fend_ij_v4[0];
        fend_v4[1] -= fend_ij_v4[1];
        fend_v4[2] -= fend_ij_v4[2];
        fend_v4[3] -= fend_ij_v4[3];
        simd_store(fend_v4[0], fend[jjj + 0]);
        simd_store(fend_v4[1], fend[jjj + 1]);
        simd_store(fend_v4[2], fend[jjj + 2]);
        simd_store(fend_v4[3], fend[jjj + 3]);

        fxtmp_v4 += fij_v4[0];
        fytmp_v4 += fij_v4[1];
        fztmp_v4 += fij_v4[2];
        if (l_pm.evflag){
          doublev4 half_v4 = simd_bcastf(0.5);
          doublev4 factor_v4 = simd_vseleq(param2->vec.skip_ev, half_v4, zero_v4) + factor_base_v4;
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
        for (jj = jjj; jj < jjj + 4; jj ++){
          int jjoff = jj - jjj;
          int j = j_4[jjoff];
          int jtype = jtype_4[jjoff];

          double r2ij = r2_4[jjoff];
          double rij = rij_4[jjoff];
          double rijinv = rinv_4[jjoff];
          double dij[3];
          dij[0] = dij_4[0][jjoff];
          dij[1] = dij_4[1][jjoff];
          dij[2] = dij_4[2][jjoff];
          double evdwl, fpair;
          double fa, fa_d, bij, bij_d;

          double ters_R = param2->sca.bigr[jjoff];
          double ters_D = param2->sca.bigd[jjoff];
          double ters_B = param2->sca.bigb[jjoff];
          double ters_lam2 = param2->sca.lam2[jjoff];

          /* if (rij > ters_R + ters_D) { */
          /*   fa = 0.0; fa_d = 0.0; */
          /* } else { */
          /*   double er = p_expd(-ters_lam2 * rij); */
          /*   fa   = -ters_B * er * fc_i[jj]; */
          /*   fa_d = ters_B * er * (ters_lam2 * fc_i[jj] - fc_d_i[jj]); */
          /* } */
          fa = fa_ij_4[jjoff];
          fa_d = fa_d_ij_4[jjoff];
          double beta = param2->sca.beta[jjoff];
          double tmp = beta * zeta_ij_4[jjoff];
          double powern = param2->sca.powern[jjoff];
          double half_powern_inv = param2->sca.half_powern_inv[jjoff];
          double c1 = param2->sca.c1[jjoff];
          double c4 = param2->sca.c4[jjoff];
          /* if (tmp > c1) { */
          /*   inv_sqrt(tmp, bij); */
          /*   bij_d = beta * -0.5 * bij * bij * bij; */
          /* } */
          /* else if (tmp < c4) { */
          /*   bij = 1.0; bij_d = 0.0; */
          /* } */
          /* else { */
          /*   double tmp_n = p_powd(tmp,powern); */
          /*   bij   = p_powd(1.0 + tmp_n, -half_powern_inv); */
          /*   bij_d = -0.5 * bij * tmp_n / (zeta_ij_4[jjoff] * (1 + tmp_n)); */
          /* } */
          bij = bij_4[jjoff];
          bij_d = bij_d_4[jjoff];
          fpair = 0.5*bij*fa_d * rijinv;
          //prefactor[jjoff] = -0.5*fa * bij_d;
          evdwl = 0.5*bij*fa;

          /* fend[jj][0] -= dij[0] * fpair; */
          /* fend[jj][1] -= dij[1] * fpair; */
          /* fend[jj][2] -= dij[2] * fpair; */
          /* if (j >= 0){ */
          /*   if (i < nlocal){ */
          /*     fxtmp += dij[0] * fpair; */
          /*     fytmp += dij[1] * fpair; */
          /*     fztmp += dij[2] * fpair; */
          /*   } */
          /*   /\* if (l_pm.evflag) *\/ */
          /*   /\*   ev_tally_global(i, j, nlocal, evdwl, -fpair, *\/ */
          /*   /\*                   -dij[0], -dij[1], -dij[2], eng_vdwl, virial); *\/ */
          /* } */
        }
        for (kk = 0; kk < nshort; kk ++){
          for (jj = jjj; jj < jjj + 4; jj ++){
            if (jj == kk) continue;
            int jjoff = jj - jjj;
            int j = j_4[jjoff];
            short_neigh_t *kshort = short_neigh + kk;
            int ktype = kshort->type;
            int k = kshort->idx;
            //tersoff_param_t *param_ijk = params3[itype][jtype_4[jjoff]] + ktype;
            double dij[3];
            dij[0] = dij_4[0][jjoff];
            dij[1] = dij_4[1][jjoff];
            dij[2] = dij_4[2][jjoff];
            double dik[3];
            dik[0] = kshort->d[0];
            dik[1] = kshort->d[1];
            dik[2] = kshort->d[2];
            double r2ij = r2_4[jjoff], rijinv = rinv_4[jjoff];
            double rij = r2ij * rijinv;
            double r2ik = kshort->r2;

            double tfi[3], tfj[3], tfk[3];
            double rij_hat[3], rik_hat[3];
            double rikinv = kshort->rinv, rik;
            rik = rikinv * r2ik;

            rij_hat[0] = rijinv * dij[0];
            rij_hat[1] = rijinv * dij[1];
            rij_hat[2] = rijinv * dij[2];

            rik_hat[0] = rikinv * dik[0];
            rik_hat[1] = rikinv * dik[1];
            rik_hat[2] = rikinv * dik[2];

            double cos_theta, ex_delr_d;
            double dcosdri[3],dcosdrj[3],dcosdrk[3];

            double ters_R = param3[kk].sca.bigr[jjoff];
            double ters_D = param3[kk].sca.bigd[jjoff];
            double Dinv   = param3[kk].sca.bigdinv[jjoff];
            double ters_lam3 = param3[kk].sca.lam3[jjoff];
            if (param3[kk].sca.power3[jjoff] != 0)
              ex_delr_d = 3.0 *
                ters_lam3 * ters_lam3 * ters_lam3 *
                (rij - rik) * (rij - rik) * ex_delr_j[kk][jjoff];
            else
              ex_delr_d = ters_lam3 * ex_delr_j[kk][jjoff];
            cos_theta = rij_hat[0] * rik_hat[0] + rij_hat[1] * rik_hat[1] + rij_hat[2] * rik_hat[2];

            dcosdrj[0] = (-cos_theta * rij_hat[0] + rik_hat[0]) * rijinv;
            dcosdrj[1] = (-cos_theta * rij_hat[1] + rik_hat[1]) * rijinv;
            dcosdrj[2] = (-cos_theta * rij_hat[2] + rik_hat[2]) * rijinv;

            dcosdrk[0] = (-cos_theta * rik_hat[0] + rij_hat[0]) * rikinv;
            dcosdrk[1] = (-cos_theta * rik_hat[1] + rij_hat[1]) * rikinv;
            dcosdrk[2] = (-cos_theta * rik_hat[2] + rij_hat[2]) * rikinv;

            tfj[0]=(gijk_d_j[kk][jjoff]*ex_delr_j[kk][jjoff]*dcosdrj[0] + gijk_j[kk][jjoff]*ex_delr_d*rij_hat[0]) * fc_i[kk] * prefactor[jjoff];
            tfj[1]=(gijk_d_j[kk][jjoff]*ex_delr_j[kk][jjoff]*dcosdrj[1] + gijk_j[kk][jjoff]*ex_delr_d*rij_hat[1]) * fc_i[kk] * prefactor[jjoff];
            tfj[2]=(gijk_d_j[kk][jjoff]*ex_delr_j[kk][jjoff]*dcosdrj[2] + gijk_j[kk][jjoff]*ex_delr_d*rij_hat[2]) * fc_i[kk] * prefactor[jjoff];

            tfk[0] = (gijk_j[kk][jjoff]*(ex_delr_j[kk][jjoff]*fc_d_i[kk] - fc_i[kk]*ex_delr_d)*rik_hat[0] + fc_i[kk]*gijk_d_j[kk][jjoff]*ex_delr_j[kk][jjoff]*dcosdrk[0]) * prefactor[jjoff];
            tfk[1] = (gijk_j[kk][jjoff]*(ex_delr_j[kk][jjoff]*fc_d_i[kk] - fc_i[kk]*ex_delr_d)*rik_hat[1] + fc_i[kk]*gijk_d_j[kk][jjoff]*ex_delr_j[kk][jjoff]*dcosdrk[1]) * prefactor[jjoff];
            tfk[2] = (gijk_j[kk][jjoff]*(ex_delr_j[kk][jjoff]*fc_d_i[kk] - fc_i[kk]*ex_delr_d)*rik_hat[2] + fc_i[kk]*gijk_d_j[kk][jjoff]*ex_delr_j[kk][jjoff]*dcosdrk[2]) * prefactor[jjoff];

            tfi[0] = -tfj[0] - tfk[0];
            tfi[1] = -tfj[1] - tfk[1];
            tfi[2] = -tfj[2] - tfk[2];

            if (j >= 0){
              fend[jj][0] += tfj[0];
              fend[jj][1] += tfj[1];
              fend[jj][2] += tfj[2];
              fend[kk][0] += tfk[0];
              fend[kk][1] += tfk[1];
              fend[kk][2] += tfk[2];
              fxtmp += tfi[0];
              fytmp += tfi[1];
              fztmp += tfi[2];
              if (vflag_either) v_tally3rd(i, j, k, nlocal, vflag_global, vflag_atom, tfj, tfk, dij, dik, virial, vatom);
            }
          }
        }
      }
      pe_put(l_pm.fend[i * maxshort], fend, sizeof(double) * nshort * 4);
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
  int ipage_start, i;
  int inum = l_pm.nlocal + l_pm.nghost;
  intv8 v8_n1 = -1;
  for (i = 0; i < ISTEP; i += 8)
    simd_store(v8_n1, n1s + i);
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
