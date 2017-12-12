#include <math.h>

#include "sunway.h"
#include "pair_tersoff_sw64.h"
#define MY_PI2 1.57079632679489661923
#define MY_PI4 0.78539816339744830962
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

inline double ters_gijk(const double costheta,
                        const tersoff_param_t * const param) {
  const double ters_c = param->c * param->c;
  const double ters_d = param->d * param->d;
  const double hcth = param->h - costheta;

  return param->gamma*(1.0 + ters_c/ters_d - ters_c / (ters_d + hcth*hcth));
}

inline double ters_gijk_d(const double costheta,
                          const tersoff_param_t * const param) {
  const double ters_c = param->c * param->c;
  const double ters_d = param->d * param->d;
  const double hcth = param->h - costheta;
  const double numerator = -2.0 * ters_c * hcth;
  const double denominator = 1.0/(ters_d + hcth*hcth);
  return param->gamma*numerator*denominator*denominator;
}
inline void costheta_d(double *rij_hat, double rij,
                             double *rik_hat, double rik,
                             double *dri, double *drj, double *drk)
{
  // first element is devative wrt Ri, second wrt Rj, third wrt Rk

  double cos_theta = vec3_dot(rij_hat,rik_hat);

  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,drj);
  vec3_scale(1.0/rij,drj,drj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,drk);
  vec3_scale(1.0/rik,drk,drk);
  vec3_add(drj,drk,dri);
  vec3_scale(-1.0,dri,dri);
}

inline double ters_fc(double r, tersoff_param_t *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r < ters_R-ters_D) return 1.0;
  if (r > ters_R+ters_D) return 0.0;
  return 0.5*(1.0 - sin(MY_PI2*(r - ters_R)/ters_D));
}

/* ---------------------------------------------------------------------- */

inline double ters_fc_d(double r, tersoff_param_t *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r < ters_R-ters_D) return 0.0;
  if (r > ters_R+ters_D) return 0.0;
  return -(MY_PI4/ters_D) * cos(MY_PI2*(r - ters_R)/ters_D);
}

inline void ters_zetaterm_d(double prefactor,
                            double *rij_hat, double rij,
                            double *rik_hat, double rik,
                            double *dri, double *drj, double *drk,
                            tersoff_param_t *param)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc(rik,param);
  dfc = ters_fc_d(rik,param);
  if (param->powermint == 3) tmp = pow(param->lam3 * (rij-rik),3.0);
  else tmp = param->lam3 * (rij-rik);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (param->powermint == 3)
    ex_delr_d = 3.0*pow(param->lam3,3.0) * pow(rij-rik,2.0)*ex_delr;
  else ex_delr_d = param->lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(cos_theta,param);
  gijk_d = ters_gijk_d(cos_theta,param);
  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Ri
  // dri = -dfc*gijk*ex_delr*rik_hat;
  // dri += fc*gijk_d*ex_delr*dcosdri;
  // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);

  vec3_scale(-dfc*gijk*ex_delr,rik_hat,dri);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdri,dri,dri);
  vec3_scaleadd(fc*gijk*ex_delr_d,rik_hat,dri,dri);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rij_hat,dri,dri);
  vec3_scale(prefactor,dri,dri);

  // compute the derivative wrt Rj
  // drj = fc*gijk_d*ex_delr*dcosdrj;
  // drj += fc*gijk*ex_delr_d*rij_hat;

  vec3_scale(fc*gijk_d*ex_delr,dcosdrj,drj);
  vec3_scaleadd(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  vec3_scale(prefactor,drj,drj);

  // compute the derivative wrt Rk
  // drk = dfc*gijk*ex_delr*rik_hat;
  // drk += fc*gijk_d*ex_delr*dcosdrk;
  // drk += -fc*gijk*ex_delr_d*rik_hat;

  vec3_scale(dfc*gijk*ex_delr,rik_hat,drk);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdrk,drk,drk);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  vec3_scale(prefactor,drk,drk);
}

void attractive(tersoff_param_t *param, double prefactor,
                             double rsqij, double rsqik,
                             double *delrij, double *delrik,
                             double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  rij = sqrt(rsqij);
  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);

  rik = sqrt(rsqik);
  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  ters_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,param);
}

#ifdef MPE
#define THIRD 0.3333333333333333333
void v_tally3rd(int i, int vflag_global, int vflag_atom,
                double *fi, double *fj, double *drik, double *drjk,
                double *virial, double (*vatom)[6])
{
  double v[6];

  v[0] = THIRD * (drik[0]*fi[0] + drjk[0]*fj[0]);
  v[1] = THIRD * (drik[1]*fi[1] + drjk[1]*fj[1]);
  v[2] = THIRD * (drik[2]*fi[2] + drjk[2]*fj[2]);
  v[3] = THIRD * (drik[0]*fi[1] + drjk[0]*fj[1]);
  v[4] = THIRD * (drik[0]*fi[2] + drjk[0]*fj[2]);
  v[5] = THIRD * (drik[1]*fi[2] + drjk[1]*fj[2]);

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

void pair_tersoff_compute_attractive(pair_tersoff_compute_param_t *pm){
  int nlocal = pm->nlocal;
  int nghost = pm->nghost;
  int ntotal = pm->ntotal;
  int inum = pm->inum;
  int gnum = pm->gnum;
  int ntypes = pm->ntypes;
  int rank = pm->rank;
  int nelements = pm->nelements;
  int vflag_global = pm->vflag_global;
  int vflag_atom = pm->vflag_atom;
  int vflag_either = pm->vflag_either;
  int *ilist = pm->ilist;
  int *firstshort = pm->firstshort;
  short_neigh_t *shortlist = pm->shortlist;
  int *map = pm->map;
  int (*elem2param)[nelements][nelements] = pm->elem2param;
  tersoff_param_t *params = pm->params;
  double (*x)[3] = pm->x;
  double (*f)[3] = pm->f;
  double (*vatom)[6] = pm->vatom;
  double *eatom = pm->eatom;
  double *virial = pm->virial;
  int *type = pm->type;
  int ii;

  for (ii = 0; ii < pm->inum; ii ++){
    int i = pm->ilist[ii];
    int itype = pm->map[type[i]];
    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];
    short_neigh_t *jlist_short = shortlist + firstshort[i];
    int jnum = firstshort[i + 1] - firstshort[i];
    double fxtmp = 0;
    double fytmp = 0;
    double fztmp = 0;
    int jj;
    for (jj = 0; jj < jnum; jj ++){
      short_neigh_t *jshort = jlist_short + jj;
      int jtype = jshort->type;//map[type[j]];
      int j = jshort->idx;
      int iparam_ij = elem2param[itype][jtype][jtype];
      double dij[3];
      dij[0] = jshort->x[0] - xtmp;
      dij[1] = jshort->x[1] - ytmp;
      dij[2] = jshort->x[2] - ztmp;
      double dji[3];
      dji[0] = -dij[0];
      dji[1] = -dij[1];
      dji[2] = -dij[2];

      double r2ij = jshort->r2;
      double prefactor_ij = jshort->prefactor_fwd;
      double prefactor_ji = jshort->prefactor_rev;

      double zeta_ij = 0;
      double zeta_ji = 0.0;

      short_neigh_t *iklist_short = shortlist + firstshort[i];
      int iknum = firstshort[i + 1] - firstshort[i];

      int kk;
      for (kk = 0; kk < iknum; kk ++){
        if (jj == kk) continue;
        short_neigh_t *kshort = iklist_short + kk;
        int ktype = kshort->type;
        int iparam_ijk = elem2param[itype][jtype][ktype];
        double dik[3];
        dik[0] = kshort->x[0] - xtmp;
        dik[1] = kshort->x[1] - ytmp;
        dik[2] = kshort->x[2] - ztmp;
        double r2ik = kshort->r2;
        double fi[3], fj[3], fk[3];
        attractive(params + iparam_ijk, prefactor_ij, r2ij, r2ik, dij, dik, fi, fj, fk);
        fxtmp += fi[0];
        fytmp += fi[1];
        fztmp += fi[2];
        if (vflag_either) v_tally3rd(i, vflag_global, vflag_atom, fj, fk, dij, dik, virial, vatom);
      }
      short_neigh_t *jklist_short = shortlist + firstshort[j];
      int jknum = firstshort[j + 1] - firstshort[j];
      for (kk = 0; kk < jknum; kk ++){
        short_neigh_t *kshort = jklist_short + kk;
        if (kshort->idx == i) continue;
        int ktype = kshort->type;//map[type[k]];
        int iparam_jik = elem2param[jtype][itype][ktype];
        int iparam_jki = elem2param[jtype][ktype][itype];
        double djk[3];
        djk[0] = kshort->x[0] - x[j][0];
        djk[1] = kshort->x[1] - x[j][1];
        djk[2] = kshort->x[2] - x[j][2];
        double r2jk = kshort->r2;
        double fi[3], fj[3], fk[3];
        attractive(params + iparam_jik, prefactor_ji, r2ij, r2jk, dji, djk, fj, fi, fk);
        fxtmp += fi[0];
        fytmp += fi[1];
        fztmp += fi[2];
        if (vflag_either) v_tally3rd(i, vflag_global, vflag_atom, fi, fk, dji, djk, virial, vatom);
        double prefactor_jk = kshort->prefactor_fwd;//[firstshort[j] + kk][0];

        attractive(params + iparam_jki, prefactor_jk, r2jk, r2ij, djk, dji, fj, fk, fi);
        fxtmp += fi[0];
        fytmp += fi[1];
        fztmp += fi[2];
        if (vflag_either) v_tally3rd(i, vflag_global, vflag_atom, fk, fi, djk, dji, virial, vatom);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}
#endif
