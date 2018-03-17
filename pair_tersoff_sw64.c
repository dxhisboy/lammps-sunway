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
      double fc_i[maxshort], fc_d_i[maxshort];
      double fxtmp = 0.0;
      double fytmp = 0.0;
      double fztmp = 0.0;
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
      for (jj = 0; jj < nshortpad; jj ++){
        short_neigh_t *jshort = short_neigh + jj;
        int jtype = jshort->type;
        int j = jshort->idx;
        double dij[3];
        dij[0] = jshort->d[0];
        dij[1] = jshort->d[1];
        dij[2] = jshort->d[2];
        double r2ij = jshort->r2;
        double rijinv = jshort->rinv, rij;
        rij = rijinv * r2ij;
        int iparam_ij = elem2param[itype][jtype][jtype];
        tersoff_param_t *param_ij_base = params3[itype][jtype];
        tersoff_param_t *param_ij = param_ij_base + jtype;
        if (r2ij >= param_ij->cutsq) continue;

        double zeta_ij = 0.0;
        int kk;
        double ex_delr_j[maxshort], gijk_j[maxshort], gijk_d_j[maxshort];
        for (kk = 0; kk < nshort; kk ++){
          if (jj == kk) continue;
          short_neigh_t *kshort = short_neigh + kk;
          int ktype = kshort->type;//map[type[k]];
          tersoff_param_t *param_ijk = param_ij_base + ktype;
          double dik[3];
          dik[0] = kshort->d[0];
          dik[1] = kshort->d[1];
          dik[2] = kshort->d[2];
          double r2ik = kshort->r2;
          if (r2ik >= param_ijk->cutsq) continue;
          double rik,costheta,arg;
          double rikinv = kshort->rinv;
          rik = rikinv * r2ik;
          costheta = (dij[0]*dik[0] + dij[1]*dik[1] +
                      dij[2]*dik[2]) * (rijinv*rikinv);

          arg = param_ijk->lam3 * (rij - rik);
          if (param_ijk->powermint == 3) arg = arg * arg * arg;

          if (arg > 69.0776) ex_delr_j[kk] = 1.e30;
          else if (arg < -69.0776) ex_delr_j[kk] = 0.0;
          else ex_delr_j[kk] = p_expd(arg);
  
          double ters_R = param_ijk->bigr;
          double ters_D = param_ijk->bigd;
 
          const double ters_c = param_ijk->c * param_ijk->c;
          const double ters_d = param_ijk->d * param_ijk->d;
          const double hcth = param_ijk->h - costheta;
          double numerator = -2.0 * ters_c * hcth;
          double denominator = 1.0/(ters_d + hcth*hcth);
          gijk_j[kk] = param_ijk->gamma * (1.0 + param_ijk->c2divd2 - ters_c * denominator);
          gijk_d_j[kk] = param_ijk->gamma * numerator * denominator * denominator;

          zeta_ij += fc_i[kk] * gijk_j[kk] * ex_delr_j[kk];
        }
        double evdwl, fpair;
        double prefactor;
        double fa, fa_d, bij, bij_d;

        double ters_R = param_ij->bigr;
        double ters_D = param_ij->bigd;
        double ters_B = param_ij->bigb;
        double ters_lam2 = param_ij->lam2;

        if (rij > ters_R + ters_D) {
          fa = 0.0; fa_d = 0.0;
        } else {
          double er = p_expd(-param_ij->lam2 * rij);
          fa   = -ters_B * er * fc_i[jj];
          fa_d = ters_B * er * (param_ij->lam2 * fc_i[jj] - fc_d_i[jj]);
        }
  
        double tmp = param_ij->beta * zeta_ij;
        double powern = param_ij->powern;
        double half_powern_inv = param_ij->half_powern_inv;
    
        if (tmp > param_ij->c1) {
          inv_sqrt(tmp, bij);
          bij_d = param_ij->beta * -0.5 * bij * bij * bij;
        }
        else if (tmp < param_ij->c4) {
          bij = 1.0; bij_d = 0.0;
        }
        else {
          double tmp_n = p_powd(tmp,powern);
          bij   = p_powd(1.0 + tmp_n, -half_powern_inv);
          bij_d = -0.5 * bij * tmp_n / (zeta_ij * (1 + tmp_n));
        }
  
        fpair = 0.5*bij*fa_d * rijinv;
        prefactor = -0.5*fa * bij_d;
        evdwl = 0.5*bij*fa;

        fend[jj][0] -= dij[0] * fpair;
        fend[jj][1] -= dij[1] * fpair;
        fend[jj][2] -= dij[2] * fpair;
        if (j >= 0){
          if (i < nlocal){
            fxtmp += dij[0] * fpair;
            fytmp += dij[1] * fpair;
            fztmp += dij[2] * fpair;
          }
          if (l_pm.evflag)
            ev_tally_global(i, j, nlocal, evdwl, -fpair,
                            -dij[0], -dij[1], -dij[2], eng_vdwl, virial);
        }

        for (kk = 0; kk < nshort; kk ++){
          if (jj == kk) continue;
          short_neigh_t *kshort = short_neigh + kk;
          int ktype = kshort->type;
          int k = kshort->idx;
          tersoff_param_t *param_ijk = param_ij_base + ktype;
          double dik[3];
          dik[0] = kshort->d[0];
          dik[1] = kshort->d[1];
          dik[2] = kshort->d[2];
          double r2ik = kshort->r2;
          if (r2ik >= param_ijk->cutsq) continue;
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

          double ters_R = param_ijk->bigr;
          double ters_D = param_ijk->bigd;
          double Dinv = param_ijk->bigdinv;

          if (param_ijk->powermint == 3)
            ex_delr_d = 3.0 *
              param_ijk->lam3 * param_ijk->lam3 * param_ijk->lam3 *
              (rij - rik) * (rij - rik) * ex_delr_j[kk];
          else ex_delr_d = param_ijk->lam3 * ex_delr_j[kk];

          cos_theta = rij_hat[0] * rik_hat[0] + rij_hat[1] * rik_hat[1] + rij_hat[2] * rik_hat[2];

          dcosdrj[0] = (-cos_theta * rij_hat[0] + rik_hat[0]) * rijinv;
          dcosdrj[1] = (-cos_theta * rij_hat[1] + rik_hat[1]) * rijinv;
          dcosdrj[2] = (-cos_theta * rij_hat[2] + rik_hat[2]) * rijinv;

          dcosdrk[0] = (-cos_theta * rik_hat[0] + rij_hat[0]) * rikinv;
          dcosdrk[1] = (-cos_theta * rik_hat[1] + rij_hat[1]) * rikinv;
          dcosdrk[2] = (-cos_theta * rik_hat[2] + rij_hat[2]) * rikinv;

          tfj[0]=(gijk_d_j[kk]*ex_delr_j[kk]*dcosdrj[0] + gijk_j[kk]*ex_delr_d*rij_hat[0]) * fc_i[kk] * prefactor;
          tfj[1]=(gijk_d_j[kk]*ex_delr_j[kk]*dcosdrj[1] + gijk_j[kk]*ex_delr_d*rij_hat[1]) * fc_i[kk] * prefactor;
          tfj[2]=(gijk_d_j[kk]*ex_delr_j[kk]*dcosdrj[2] + gijk_j[kk]*ex_delr_d*rij_hat[2]) * fc_i[kk] * prefactor;

          tfk[0] = (gijk_j[kk]*(ex_delr_j[kk]*fc_d_i[kk] - fc_i[kk]*ex_delr_d)*rik_hat[0] + fc_i[kk]*gijk_d_j[kk]*ex_delr_j[kk]*dcosdrk[0]) * prefactor;
          tfk[1] = (gijk_j[kk]*(ex_delr_j[kk]*fc_d_i[kk] - fc_i[kk]*ex_delr_d)*rik_hat[1] + fc_i[kk]*gijk_d_j[kk]*ex_delr_j[kk]*dcosdrk[1]) * prefactor;
          tfk[2] = (gijk_j[kk]*(ex_delr_j[kk]*fc_d_i[kk] - fc_i[kk]*ex_delr_d)*rik_hat[2] + fc_i[kk]*gijk_d_j[kk]*ex_delr_j[kk]*dcosdrk[2]) * prefactor;

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
      pe_put(l_pm.fend[i * maxshort], fend, sizeof(double) * nshort * 4);
      fi[ioff][0] = fxtmp;
      fi[ioff][1] = fytmp;
      fi[ioff][2] = fztmp;
      pe_syn();
    }
    if (iwr > 0){
      pe_put(f + ist, fi, sizeof(double) * iwr * 3);
    }
    pe_put(l_pm.numshort + ist, numshort_local, sizeof(int) * isz);
    pe_syn();
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
