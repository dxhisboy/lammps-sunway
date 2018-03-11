#include "sunway.h"
#include "pair_tersoff_sw64.h"
#include <stdio.h>
#include <stdlib.h>
#include "gptl.h"
#define ISTEP 128


#ifdef MPE
#define LWPF_UNITS U(TERSOFF)
#include "lwpf.h"
#include <simd.h>
extern SLAVE_FUN(pair_tersoff_compute_attractive_para)(void*);

int r = 0;
void waitint(int *ptr){
  while (1){
    int tmp;
    asm ("ldw_nc %0, 0(%1)": "=&r"(tmp): "r"(ptr));
    if (tmp)
      break;
  }
}


void pair_tersoff_compute_attractive(pair_tersoff_compute_param_t *pm){
  if (athread_idle() == 0)
    athread_init();
  long fend_base = (long)malloc(sizeof(double) * pm->firstshort[pm->ntotal] * 4 + 32);
  long ftmp_base = (long)calloc(pm->nlocal * 4 + 4, sizeof(double));
  /* long fend_base = (long)pm->fend; */
  /* long ftmp_base = (long)pm->ftmp; */
  pm->fend = (void*)((fend_base + 31) & (~31));
  double (*ftmp)[4] = (void*)((ftmp_base + 31) & (~31));
  pm->fdone = calloc(pm->ntotal, sizeof(int));
  GPTLstart("attractive athread");
  athread_spawn(pair_tersoff_compute_attractive_para, pm);
  int ist, ied, ii;
  int ntotal = pm->ntotal;
  int stotal = pm->firstshort[ntotal];
  int *firstshort = pm->firstshort;
  double (*f)[3] = pm->f;
  short_neigh_t *shortlist = pm->shortlist;
  int *shortidx = pm->shortidx;
  int nlocal = pm->nlocal;
  double (*fend)[4] = pm->fend;
  int *fdone = pm->fdone;
  GPTLstart("attractive reduction bond");
  for (ist = 0; ist < ntotal; ist += ISTEP){
    ied = ist + ISTEP;
    if (ied > ntotal)
      ied = ntotal;
    waitint(fdone + ist);
    for (ii = ist; ii < ied; ii ++){
      int jj, jend;
      jend = firstshort[ii + 1];
      for (jj = firstshort[ii]; jj < jend; jj += 4){
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
          //simd_print_doublev4(fend3);
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
  for (ii = 0; ii < pm->nlocal; ii ++){
    pm->f[ii][0] += ftmp[ii][0];
    pm->f[ii][1] += ftmp[ii][1];
    pm->f[ii][2] += ftmp[ii][2];
  }
  GPTLstop("attractive reduction force");
  free((void*)fend_base);
  free((void*)ftmp_base);
  free(ftmp);
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

inline double vec3_dot(const double *x, const double *y){
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

inline void vec3_add(const double *x, const double *y,
                     double * const z) {
  z[0] = x[0]+y[0];  z[1] = x[1]+y[1];  z[2] = x[2]+y[2];
}

inline void vec3_scale(const double k, const double *x,
                       double *y) {
  y[0] = k*x[0];  y[1] = k*x[1];  y[2] = k*x[2];
}

inline void vec3_scaleadd(const double k, const double *x,
                          const double *y, double * const z) {
  z[0] = k*x[0]+y[0];
  z[1] = k*x[1]+y[1];
  z[2] = k*x[2]+y[2];
}

inline void ters_attractive_unroll(double prefactor,
                            double rsqij, double rsqik,
                            double *delrij, double *delrik,
                            double *dri, double *drj, double *drk,
                            tersoff_param_t *param)
{
  double rij_hat[3], rik_hat[3];
  double rijinv, rikinv, rij, rik;

  inv_sqrt(rsqij, rijinv);
  inv_sqrt(rsqik, rikinv);

  rij = rijinv * rsqij;
  rik = rikinv * rsqik;

  vec3_scale(rijinv, delrij, rij_hat);
  vec3_scale(rikinv, delrik, rik_hat);
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  double ters_R = param->bigr;
  double ters_D = param->bigd;
  double Dinv = param->bigdinv;
  if (rik < ters_R - ters_D){
    fc = 1.0;
    dfc = 0.0;
  } else if (rik > ters_R + ters_D){
    fc = 0.0;
    dfc = 0.0;
  } else {
    fc = 0.5*(1.0 - p_sinnpi_pid(MY_PI2*(rik - ters_R)*Dinv));
    dfc = -(MY_PI4*Dinv) * p_cosnpi_pid(MY_PI2*(rik - ters_R)*Dinv);
  }

  tmp = param->lam3 * (rij-rik);
  if (param->powermint == 3) tmp = tmp * tmp * tmp;

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = p_expd(tmp);

  if (param->powermint == 3)
    ex_delr_d = 3.0 *
      param->lam3 * param->lam3 * param->lam3 *
      (rij - rik) * (rij - rik) * ex_delr;
  else ex_delr_d = param->lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  double ters_c = param->c * param->c;
  double ters_d = param->d * param->d;
  double hcth = param->h - cos_theta;
  double numerator = -2.0 * ters_c * hcth;
  double denominator = 1.0/(ters_d + hcth*hcth);
  gijk = param->gamma * (1.0 + param->c2divd2 - ters_c * denominator);
  gijk_d = param->gamma * numerator * denominator * denominator;

  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,dcosdrj);
  vec3_scale(rijinv,dcosdrj,dcosdrj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,dcosdrk);
  vec3_scale(rikinv,dcosdrk,dcosdrk);

  vec3_scale(fc*gijk_d*ex_delr,dcosdrj,drj);
  vec3_scaleadd(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  vec3_scale(prefactor,drj,drj);

  vec3_scale(dfc*gijk*ex_delr,rik_hat,drk);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdrk,drk,drk);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  vec3_scale(prefactor,drk,drk);

  vec3_add(drj, drk, dri);
  vec3_scale(-1, dri, dri);
}


double zeta_unroll(tersoff_param_t *param, double rsqij, double rsqik,
                   double *delrij, double *delrik)
{
  double rij,rik,costheta,arg,ex_delr;
  double rijinv, rikinv;
  inv_sqrt(rsqij, rijinv);
  inv_sqrt(rsqik, rikinv);
  rij = rijinv * rsqij;
  rik = rikinv * rsqik;
  costheta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] +
              delrij[2]*delrik[2]) * (rijinv*rikinv);

  arg = param->lam3 * (rij - rik);
  if (param->powermint == 3) arg = arg * arg * arg;

  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = p_expd(arg);
  
  double ters_R = param->bigr;
  double ters_D = param->bigd;
  double fc;
  if (rik < ters_R-ters_D) fc = 1.0;
  else if (rik > ters_R+ters_D) fc = 0.0;
  else fc =  0.5*(1.0 - p_sinnpi_pid(MY_PI2*(rik - ters_R) * param->bigdinv));
 
  double gijk;
  const double ters_c = param->c * param->c;
  const double ters_d = param->d * param->d;
  const double hcth = param->h - costheta;
  gijk =  param->gamma*(1.0 + param->c2divd2/* ters_c/ters_d */ - ters_c / (ters_d + hcth*hcth));

  return fc * gijk * ex_delr;
}

inline void force_zeta_unroll(tersoff_param_t *param, double rsq, double zeta_ij,
                       double *fforce, double *prefactor, double *eng)
{
  double r, rinv, fa, fa_d, bij, bij_d, fc, fc_d;
  inv_sqrt(rsq, rinv);
  r = rsq * rinv;

  double ters_R = param->bigr;
  double ters_D = param->bigd;
  double ters_B = param->bigb;
  double ters_lam2 = param->lam2;

  if (r < ters_R-ters_D) 
  {
    fc = 1.0; fc_d = 0.0;
  }
  else if (r > ters_R+ters_D) 
  {
    fc = 0.0; fc_d = 0.0;
  }
  else
  {
    fc = 0.5*(1.0 - p_sinnpi_pid(MY_PI2*(r - ters_R)* param->bigdinv));
    fc_d = -(MY_PI4* param->bigdinv) * p_cosnpi_pid(MY_PI2*(r - ters_R)* param->bigdinv);
  }

  if (r > ters_R + ters_D)
  {
    fa = 0.0; fa_d = 0.0;
  }
  else 
  {
    double er = p_expd(-param->lam2 * r);
    fa   = -ters_B * er * fc;
    fa_d = ters_B * er * (param->lam2 * fc - fc_d);
  }
  
  double tmp = param->beta * zeta_ij;
  double powern = param->powern;
  double half_powern_inv = param->half_powern_inv;
  /* double n1 = -1.0/(2.0*powern); */
  /* double n2 = -1.0-(1.0/(2.0*powern)); */
    
  if (tmp > param->c1) 
  {
    inv_sqrt(tmp, bij);
    bij_d = param->beta * -0.5 * bij * bij * bij;
  }
  else if (tmp > param->c2)
  { 
    double tmpinvsqrt;
    double tmp_nn = p_powd(tmp, -powern);
    inv_sqrt(tmp, tmpinvsqrt);
    bij   = (1.0 - tmp_nn * half_powern_inv) * tmpinvsqrt;
    bij_d = param->beta * (-0.5 * tmpinvsqrt * tmpinvsqrt * tmpinvsqrt *
            (1.0 - (1.0 + half_powern_inv) * tmp_nn));
  }
  else if (tmp < param->c4) 
  {
    bij = 1.0; bij_d = 0.0;
  }
  else if (tmp < param->c3)
  {
    double tmp_n = p_powd(tmp, powern);
    bij   = 1.0 - tmp_n * half_powern_inv;
    bij_d = -0.5*param->beta * tmp_n;
  }
  else
  { 
    double tmp_n = p_powd(tmp,powern);
    bij   = p_powd(1.0 + tmp_n, -half_powern_inv);
    //double tmp_n1inv = 1 / (1 + tmp_n);
    //bij_d = -0.5 * p_powd(1.0+tmp_n, -1.0-(half_powern_inv))*tmp_n / zeta_ij;
    bij_d = -0.5 * bij * tmp_n / (zeta_ij * (1 + tmp_n));
  }
  
  *fforce = 0.5*bij*fa_d * rinv;
  *prefactor = -0.5*fa * bij_d;
  *eng = 0.5*bij*fa;
}


#define HALF 0.5
void ev_tally_global(int i, int j, int nlocal, double evdwl, double fpair,
                     double delx, double dely, double delz,
                     double *eng_vdwl, double *virial, int eflag_global, int vflag_global)
{
  double v[6];
  double factor = 0;
  if (i < nlocal) factor += HALF;
  if (j < nlocal) factor += HALF;
  if (eflag_global){
    *eng_vdwl += evdwl * factor;
  }
  if (vflag_global){
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

#define SNSTEP 64
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
  short_neigh_t *shortlist = l_pm.shortlist;

  int map[ntypes + 1];
  int elem2param[nelements][nelements][nelements];
  tersoff_param_t params[nparams];
  int nep3 = nelements * nelements * nelements;
  pe_get(l_pm.map, map, sizeof(int) * (ntypes + 1));
  pe_get(l_pm.elem2param, elem2param, sizeof(int) * nep3);
  pe_get(l_pm.params, params, sizeof(tersoff_param_t) * nparams);
  pe_syn();

  double (*x)[3] = l_pm.x;
  double (*f)[3] = l_pm.f;
  double (*vatom)[6] = l_pm.vatom;
  double *eatom = l_pm.eatom;

  int *type = l_pm.type;
  int ii;
  int ist, ied;
  double fi[ISTEP][3];
  int ti[ISTEP], fs[ISTEP + 1];
  short_neigh_t js[SNSTEP], ks[SNSTEP];
  double fend[SNSTEP][4];
  int fdone[ISTEP];
  doublev4 eng_virial_v4[2];
  eng_virial_v4[0] = 0.0;
  eng_virial_v4[1] = 0.0;
  double *eng_vdwl = eng_virial_v4;
  double *eng_coul = eng_vdwl + 1;
  double *virial = eng_vdwl + 2;
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
    //pe_get(x + ist, xi, sizeof(double) * 3 * isz);
    pe_get(f + ist, fi, sizeof(double) * 3 * isz);
    pe_get(firstshort + ist, fs, sizeof(int) * (isz + 1));
    pe_get(type + ist, ti, sizeof(int) * isz);
    pe_syn();

    for (i = ist; i < ied; i ++){
      int ioff = i - ist;
      int itype = map[ti[ioff]];
      short_neigh_t *jlist_short = shortlist + fs[ioff];
      int jnum = fs[ioff + 1] - fs[ioff];
      double fxtmp = 0;
      double fytmp = 0;
      double fztmp = 0;
      int jj;
      pe_get(jlist_short, js, sizeof(short_neigh_t) * jnum);
      for (jj = 0; jj < jnum; jj ++){
        fend[jj][0] = 0;
        fend[jj][1] = 0;
        fend[jj][2] = 0;
      }
      pe_syn();
      for (jj = 0; jj < jnum; jj ++){
        short_neigh_t *jshort = js + jj;
        int jtype = jshort->type;
        int j = jshort->idx;
        double dij[3];
        dij[0] = jshort->d[0];
        dij[1] = jshort->d[1];
        dij[2] = jshort->d[2];
        double r2ij = jshort->r2;
        int iparam_ij = elem2param[itype][jtype][jtype];
        if (r2ij >= params[iparam_ij].cutsq) continue;

        short_neigh_t *klist_short = js;
        int knum = jnum;
        double zeta_ij = 0.0;

        int kk;
        for (kk = 0; kk < knum; kk ++){
          if (jj == kk) continue;
          short_neigh_t *kshort = klist_short + kk;
          int ktype = kshort->type;//map[type[k]];
          int iparam_ijk = elem2param[itype][jtype][ktype];
          double dik[3];
          dik[0] = kshort->d[0];
          dik[1] = kshort->d[1];
          dik[2] = kshort->d[2];
          double r2ik = kshort->r2;//dik[0] * dik[0] + dik[1] * dik[1] + dik[2] * dik[2];
          if (r2ik >= params[iparam_ijk].cutsq) continue;
          zeta_ij += zeta_unroll(params + iparam_ijk, r2ij, r2ik, dij, dik);
        }
        double evdwl, fpair;
        double prefactor;
        force_zeta_unroll(params + iparam_ij, r2ij, zeta_ij, &fpair, &prefactor, &evdwl);
        fend[jj][0] -= dij[0] * fpair;
        fend[jj][1] -= dij[1] * fpair;
        fend[jj][2] -= dij[2] * fpair;
        if (i < nlocal){
          fxtmp += dij[0] * fpair;
          fytmp += dij[1] * fpair;
          fztmp += dij[2] * fpair;
        }
        if (l_pm.evflag)
          ev_tally_global(i, j, nlocal, evdwl, -fpair,
                          -dij[0], -dij[1], -dij[2], eng_vdwl, virial,
                          l_pm.eflag_global, l_pm.vflag_global);


        for (kk = 0; kk < knum; kk ++){
          if (jj == kk) continue;
          short_neigh_t *kshort = klist_short + kk;
          int ktype = kshort->type;
          int k = kshort->idx;
          int iparam_ijk = elem2param[itype][jtype][ktype];
          double dik[3];
          dik[0] = kshort->d[0];
          dik[1] = kshort->d[1];
          dik[2] = kshort->d[2];
          double r2ik = kshort->r2;
          if (r2ik >= params[iparam_ijk].cutsq) continue;
          double tfi[3], tfj[3], tfk[3];
          ters_attractive_unroll(prefactor, r2ij, r2ik, dij, dik, tfi, tfj, tfk, params + iparam_ijk);
          //lwpf_stop(ATTRACTIVE);
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
      pe_put(l_pm.fend[fs[ioff]], fend, sizeof(double) * jnum * 4);
      fi[ioff][0] += fxtmp;
      fi[ioff][1] += fytmp;
      fi[ioff][2] += fztmp;
      pe_syn();
    }
    if (iwr > 0){
      pe_put(f + ist, fi, sizeof(double) * iwr * 3);
    }
    pe_put(l_pm.fdone + ist, fdone, sizeof(int) * isz);
    pe_syn();

  }

  reg_reduce_inplace_doublev4(eng_virial_v4, 2);

  if (_MYID == 0){
    pe_put(&(pm->eng_vdwl), eng_virial_v4, sizeof(double) * 8);
    pe_syn();
  }
  //lwpf_stop(ALL);
}

#endif
