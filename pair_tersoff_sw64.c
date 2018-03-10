#include "sunway.h"
#include "pair_tersoff_sw64.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gptl.h"
#define ISTEP 128
#ifdef CPE
__thread_local rank;
#define MY_PI2 1.57079632679489661923
#define MY_PI4 0.78539816339744830962

#define EXP_PPC0 1.00000000000969824220931059244e+00
#define EXP_PPC1 9.99999996155887638238368708699e-01
#define EXP_PPC2 4.99999998889176233696218787372e-01
#define EXP_PPC3 1.66666783941218449305310400632e-01
#define EXP_PPC4 4.16666866967536769772451066274e-02
#define EXP_PPC5 8.33238238257071571479794869219e-03
#define EXP_PPC6 1.38876391271098567278818869397e-03
#define EXP_PPC7 2.01234277126192251773997843323e-04
#define EXP_PPC8 2.51169087281361337897402780106e-05
#define EXP_PPC9 0.007812500000000000000000000000000


inline double exp_4_tersoff(double x){
  double x_128th = x * EXP_PPC9;

  double expx = x_128th * EXP_PPC8 + EXP_PPC7;
  expx = expx * x_128th + EXP_PPC6;
  expx = expx * x_128th + EXP_PPC5;
  expx = expx * x_128th + EXP_PPC4;
  expx = expx * x_128th + EXP_PPC3;
  expx = expx * x_128th + EXP_PPC2;
  expx = expx * x_128th + EXP_PPC1;
  expx = expx * x_128th + EXP_PPC0;
  expx *= expx;
  expx *= expx;
  expx *= expx;
  expx *= expx;
  expx *= expx;
  expx *= expx;
  expx *= expx;
  return expx;
}

#define SIN_PPC0  9.99999999919516047164336214337e-01
#define SIN_PPC1 -1.66666665687026344100374330992e-01
#define SIN_PPC2  8.33332995332082202000201220926e-03
#define SIN_PPC3 -1.98407729348946363246569690730e-04
#define SIN_PPC4  2.75219407932635322134055574916e-06
#define SIN_PPC5 -2.38436847087956781783966877901e-08

#define COS_PPC0  9.99999999999460209565427248890e-01
#define COS_PPC1 -4.99999999977592313182839234287e-01
#define COS_PPC2  4.16666665120830198576484804107e-02
#define COS_PPC3 -1.38888849148028747232386237442e-03
#define COS_PPC4  2.48011030436940119387519837479e-05
#define COS_PPC5 -2.75271172247895381065456790748e-07
#define COS_PPC6  1.99428286541087023330572854704e-09

inline double sin_4_tersoff(double x){
  double x2 = x * x;
  double sinx = x2 * SIN_PPC5 + SIN_PPC4;
  sinx = sinx * x2 + SIN_PPC3;
  sinx = sinx * x2 + SIN_PPC2;
  sinx = sinx * x2 + SIN_PPC1;
  sinx = sinx * x2 + SIN_PPC0;
  sinx *= x;
  return sinx;
}
inline double cos_4_tersoff(double x){
  double x2 = x * x;
  double cosx = x2 * COS_PPC6 + COS_PPC5;
  cosx = cosx * x2 + COS_PPC4;
  cosx = cosx * x2 + COS_PPC3;
  cosx = cosx * x2 + COS_PPC2;
  cosx = cosx * x2 + COS_PPC1;
  cosx = cosx * x2 + COS_PPC0;
  return cosx;
}
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

/* #define inv_sqrt(x, r) {r = 1.0 / sqrt(x);} */
/* #define cos_4_tersoff cos */
/* #define sin_4_tersoff sin */
/* #define exp_4_tersoff exp */
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


/* inline void attractive(tersoff_param_t *param, double prefactor, */
/*                        double rsqij, double rsqik, */
/*                        double *delrij, double *delrik, */
/*                        double *fi, double *fj, double *fk) */
/* { */
/*   double rij_hat[3],rik_hat[3]; */
/*   double rij,rijinv,rik,rikinv; */

/*   inv_sqrt(rsqij, rijinv); */
/*   rij = rsqij * rijinv; */
/*   vec3_scale(rijinv,delrij,rij_hat); */

/*   inv_sqrt(rsqik, rikinv); */
/*   rik = rsqik * rikinv; */

/*   vec3_scale(rikinv,delrik,rik_hat); */

/*   //ters_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,param); */
/*   double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp; */
/*   double dcosdri[3],dcosdrj[3],dcosdrk[3]; */
 
/*   double ters_R = param->bigr; */
/*   double ters_D = param->bigd; */
/*   double Dinv = param->bigdinv; */
/*   if (rik < ters_R - ters_D){ */
/*     fc = 1.0; */
/*     dfc = 0.0; */
/*   } else if (rik > ters_R + ters_D){ */
/*     fc = 0.0; */
/*     dfc = 0.0; */
/*   } else { */
/*     fc = 0.5 * (1.0 - sin_4_tersoff(MY_PI2 * (rik - ters_R) * Dinv)); */
/*     dfc = -(MY_PI4 * Dinv) * cos_4_tersoff(MY_PI2 * (rik - ters_R) * Dinv); */
/*   } */
/*   /\* fc = ters_fc(rik,param); *\/ */
/*   /\* dfc = ters_fc_d(rik,param); *\/ */
/*   tmp = param->lam3 * (rij-rik); */
/*   if (param->powermint == 3) tmp = tmp * tmp * tmp; */

/*   if (tmp > 69.0776) ex_delr = 1.e30; */
/*   else if (tmp < -69.0776) ex_delr = 0.0; */
/*   else ex_delr = exp_4_tersoff(tmp); */

/*   if (param->powermint == 3) */
/*     ex_delr_d = 3.0*param->lam3*param->lam3*param->lam3 *(rij-rik)*(rij-rik)*ex_delr; */
/*   else ex_delr_d = param->lam3 * ex_delr; */

/*   cos_theta = vec3_dot(rij_hat,rik_hat); */
/*   //ters_gijk(d) below */
/*   double ters_c = param->c * param->c; */
/*   double ters_d = param->d * param->d; */
/*   double hcth = param->h - cos_theta; */
/*   double numerator = -2.0 * ters_c * hcth; */
/*   double denominator = 1.0 / (ters_d + hcth * hcth); */
/*   gijk = param->gamma * (1.0 + param->c2divd2 - ters_c * denominator); */
/*   gijk_d = param->gamma * numerator * denominator * denominator; */
/*   /\* gijk = ters_gijk(cos_theta,param); *\/ */
/*   /\* gijk_d = ters_gijk_d(cos_theta,param); *\/ */
  

/*   //costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk); */

/*   vec3_scaleadd(-cos_theta,rij_hat,rik_hat,dcosdrj); */
/*   vec3_scale(rijinv,dcosdrj,dcosdrj); */
/*   vec3_scaleadd(-cos_theta,rik_hat,rij_hat,dcosdrk); */
/*   vec3_scale(rikinv,dcosdrk,dcosdrk); */
/*   /\* vec3_add(dcosdrj,dcosdrk,dcosdri); *\/ */
/*   /\* vec3_scale(-1.0,dcosdri,dcosdri); *\/ */

/*   // compute the derivative wrt Ri */
/*   // dri = -dfc*gijk*ex_delr*rik_hat; */
/*   // dri += fc*gijk_d*ex_delr*dcosdri; */
/*   // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat); */

/*   vec3_scale(-dfc*gijk*ex_delr,rik_hat,fi); */
/*   vec3_scaleadd(fc*gijk_d*ex_delr,dcosdri,fi,fi); */
/*   vec3_scaleadd(fc*gijk*ex_delr_d,rik_hat,fi,fi); */
/*   vec3_scaleadd(-fc*gijk*ex_delr_d,rij_hat,fi,fi); */
/*   vec3_scale(prefactor,fi,fi); */

/*   // compute the derivative wrt Rj */
/*   // drj = fc*gijk_d*ex_delr*dcosdrj; */
/*   // drj += fc*gijk*ex_delr_d*rij_hat; */
/*   //Used below: fc, gijk_d, ex_delr, dcosdrj, gijk, ex_delr_d, dfc, rik_hat, dcosdrk */
/*   /\* vec3_scale(fc*gijk_d*ex_delr,dcosdrj,fj); *\/ */
/*   /\* vec3_scaleadd(fc*gijk*ex_delr_d,rij_hat,fj,fj); *\/ */
/*   /\* vec3_scale(prefactor,fj,fj); *\/ */

/*   /\* // compute the derivative wrt Rk *\/ */
/*   /\* // drk = dfc*gijk*ex_delr*rik_hat; *\/ */
/*   /\* // drk += fc*gijk_d*ex_delr*dcosdrk; *\/ */
/*   /\* // drk += -fc*gijk*ex_delr_d*rik_hat; *\/ */

/*   /\* vec3_scale(dfc*gijk*ex_delr,rik_hat,fk); *\/ */
/*   /\* vec3_scaleadd(fc*gijk_d*ex_delr,dcosdrk,fk,fk); *\/ */
/*   /\* vec3_scaleadd(-fc*gijk*ex_delr_d,rik_hat,fk,fk); *\/ */
/*   /\* vec3_scale(prefactor,fk,fk); *\/ */

/*   /\* vec3_add(fj, fk, fi); *\/ */
/*   /\* vec3_scale(-1.0, fi, fi); *\/ */
/* } */

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

double ters_fc(double r, tersoff_param_t *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r < ters_R-ters_D) return 1.0;
  if (r > ters_R+ters_D) return 0.0;
  return 0.5*(1.0 - sin_4_tersoff(MY_PI2*(r - ters_R)/ters_D));
}

double ters_fc_d(double r, tersoff_param_t *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r < ters_R-ters_D) return 0.0;
  if (r > ters_R+ters_D) return 0.0;
  return -(MY_PI4/ters_D) * cos_4_tersoff(MY_PI2*(r - ters_R)/ters_D);
}

/* ---------------------------------------------------------------------- */
void costheta_d(double *rij_hat, double rij,
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

double ters_fa(double r, tersoff_param_t *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return -param->bigb * exp_4_tersoff(-param->lam2 * r) * ters_fc(r,param);
}

/* ---------------------------------------------------------------------- */

double ters_fa_d(double r, tersoff_param_t *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return param->bigb * exp_4_tersoff(-param->lam2 * r) *
    (param->lam2 * ters_fc(r,param) - ters_fc_d(r,param));
}

/* ---------------------------------------------------------------------- */

double ters_bij(double zeta, tersoff_param_t *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->c1) return 1.0/sqrt(tmp);
  if (tmp > param->c2)
    return (1.0 - pow(tmp,-param->powern) / (2.0*param->powern))/sqrt(tmp);
  if (tmp < param->c4) return 1.0;
  if (tmp < param->c3)
    return 1.0 - pow(tmp,param->powern)/(2.0*param->powern);
  return pow(1.0 + pow(tmp,param->powern), -1.0/(2.0*param->powern));
}

/* ---------------------------------------------------------------------- */

double ters_bij_d(double zeta, tersoff_param_t *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->c1) return param->beta * -0.5*pow(tmp,-1.5);
  if (tmp > param->c2)
    return param->beta * (-0.5*pow(tmp,-1.5) *
			  // error in negligible 2nd term fixed 9/30/2015
			  // (1.0 - 0.5*(1.0 +  1.0/(2.0*param->powern)) *
                          (1.0 - (1.0 +  1.0/(2.0*param->powern)) *
                           pow(tmp,-param->powern)));
  if (tmp < param->c4) return 0.0;
  if (tmp < param->c3)
    return -0.5*param->beta * pow(tmp,param->powern-1.0);

  double tmp_n = pow(tmp,param->powern);
  return -0.5 * pow(1.0+tmp_n, -1.0-(1.0/(2.0*param->powern)))*tmp_n / zeta;
}

/* ---------------------------------------------------------------------- */

void ters_zetaterm_d(double prefactor,
                     double *rij_hat, double rij,
                     double *rik_hat, double rik,
                     double *dri, double *drj, double *drk,
                     tersoff_param_t *param)
{
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
    fc = 0.5*(1.0 - sin_4_tersoff(MY_PI2*(rik - ters_R)*Dinv));
    dfc = -(MY_PI4*Dinv) * cos_4_tersoff(MY_PI2*(rik - ters_R)*Dinv);
  }

  tmp = param->lam3 * (rij-rik);
  if (param->powermint == 3) tmp = tmp * tmp * tmp;

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp_4_tersoff(tmp);

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
  gijk = param->gamma * (1.0 + ters_c / ters_d - ters_c * denominator);
  gijk_d = param->gamma * numerator * denominator * denominator;

  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,dcosdrj);
  vec3_scale(1.0/rij,dcosdrj,dcosdrj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,dcosdrk);
  vec3_scale(1.0/rik,dcosdrk,dcosdrk);

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

/* ---------------------------------------------------------------------- */

inline void ters_attractive_unroll(double prefactor,
                            /* double *rij_hat, double rij, */
                            /* double *rik_hat, double rik, */
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
    fc = 0.5*(1.0 - sin_4_tersoff(MY_PI2*(rik - ters_R)*Dinv));
    dfc = -(MY_PI4*Dinv) * cos_4_tersoff(MY_PI2*(rik - ters_R)*Dinv);
  }

  tmp = param->lam3 * (rij-rik);
  if (param->powermint == 3) tmp = tmp * tmp * tmp;

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp_4_tersoff(tmp);

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
  gijk = param->gamma * (1.0 + ters_c / ters_d - ters_c * denominator);
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


inline void ters_attractive_unroll_pair(double prefactor_ij, double prefactor_ik,
                                        double rsqij, double rsqik,
                                        double *delrij, double *delrik,
                                        double *dri, double *drj, double *drk,
                                        tersoff_param_t *param_ijk,
                                        tersoff_param_t *param_ikj)
{
  double rij_hat[3], rik_hat[3];
  double rijinv, rikinv, rij, rik;

  inv_sqrt(rsqij, rijinv);
  inv_sqrt(rsqik, rikinv);

  rij = rijinv * rsqij;
  rik = rikinv * rsqik;

  vec3_scale(rijinv, delrij, rij_hat);
  vec3_scale(rikinv, delrik, rik_hat);

  double cos_theta = vec3_dot(rij_hat,rik_hat);
  double dcosdri[3],dcosdrj[3],dcosdrk[3];
  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,dcosdrj);
  vec3_scale(rijinv,dcosdrj,dcosdrj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,dcosdrk);
  vec3_scale(rikinv,dcosdrk,dcosdrk);

  double ters_R_ij = param_ijk->bigr;
  double ters_D_ij = param_ijk->bigd;
  double Dinv_ij = param_ijk->bigdinv;
  double fc_ij, dfc_ij;
  if (rik < ters_R_ij - ters_D_ij){
    fc_ij = 1.0;
    dfc_ij = 0.0;
  } else if (rik > ters_R_ij + ters_D_ij){
    fc_ij = 0.0;
    dfc_ij = 0.0;
  } else {
    fc_ij = 0.5*(1.0 - sin_4_tersoff(MY_PI2*(rik - ters_R_ij)*Dinv_ij));
    dfc_ij = -(MY_PI4*Dinv_ij)*cos_4_tersoff(MY_PI2*(rik-ters_R_ij)*Dinv_ij);
  }

  double ters_R_ik = param_ikj->bigr;
  double ters_D_ik = param_ikj->bigd;
  double Dinv_ik = param_ikj->bigdinv;
  double fc_ik, dfc_ik;
  if (rij < ters_R_ik - ters_D_ik){
    fc_ik = 1.0;
    dfc_ik = 0.0;
  } else if (rij > ters_R_ik + ters_D_ik){
    fc_ik = 0.0;
    dfc_ik = 0.0;
  } else {
    fc_ik = 0.5*(1.0 - sin_4_tersoff(MY_PI2*(rij - ters_R_ik)*Dinv_ik));
    dfc_ik = -(MY_PI4*Dinv_ik)*cos_4_tersoff(MY_PI2*(rij-ters_R_ik)*Dinv_ik);
  }

  double tmp_ij = param_ijk->lam3 * (rij-rik);
  if (param_ijk->powermint == 3) tmp_ij = tmp_ij * tmp_ij * tmp_ij;
  /* double tmp_ik = param_ikj->lam3 * (rik-rij); */
  /* if (param_ikj->powermint == 3) tmp_ik = tmp_ik * tmp_ik * tmp_ik; */
  /* double tmp_ik = -tmp_ij; */
  double ex_delr_ij, ex_delr_ik;
  if (tmp_ij > 69.0776) {
    ex_delr_ij = 1.e30;
    ex_delr_ik = 0.0;
  } else if (tmp_ij < -69.0776) {
    ex_delr_ij = 0.0;
    ex_delr_ik = 1.e30;
  } else {
    ex_delr_ij = exp_4_tersoff(tmp_ij);
    ex_delr_ik = exp_4_tersoff(-tmp_ij);
  }

  double lam3_ij = param_ijk->lam3;
  double ex_delr_d_ij;
  if (param_ijk->powermint == 3)
    ex_delr_d_ij = 3.0 * lam3_ij * lam3_ij * lam3_ij * (rij - rik) * (rij - rik) * ex_delr_ij;
  else ex_delr_d_ij = param_ijk->lam3 * ex_delr_ij;

  double lam3_ik = param_ikj->lam3;
  double ex_delr_d_ik;
  if (param_ikj->powermint == 3)
    ex_delr_d_ik = 3.0 * lam3_ik * lam3_ik * lam3_ik * (rik - rij) * (rik - rij) * ex_delr_ik;
  else ex_delr_d_ik = param_ikj->lam3 * ex_delr_ik;

  double ters_c_ij = param_ijk->c * param_ijk->c;
  double ters_c_ik = param_ikj->c * param_ikj->c;
  double ters_d_ij = param_ijk->d * param_ijk->d;
  double ters_d_ik = param_ikj->d * param_ikj->d;
  double hcth_ij = param_ijk->h - cos_theta;
  double hcth_ik = param_ikj->h - cos_theta;
  double numerator_ij = -2.0 * ters_c_ij * hcth_ij;
  double numerator_ik = -2.0 * ters_c_ik * hcth_ik;
  double denominator_ij = 1.0/(ters_d_ij + hcth_ij*hcth_ij);
  double denominator_ik = 1.0/(ters_d_ik + hcth_ik*hcth_ik);
  double gijk_ij = param_ijk->gamma * (1.0 + ters_c_ij/ters_d_ij - ters_c_ij*denominator_ij);
  double gijk_ik = param_ikj->gamma * (1.0 + ters_c_ik/ters_d_ik - ters_c_ik*denominator_ik);
  double gijk_d_ij = param_ijk->gamma * numerator_ij * denominator_ij * denominator_ij;
  double gijk_d_ik = param_ikj->gamma * numerator_ik * denominator_ik * denominator_ik;

  double fdrj[3], fdrk[3];
  vec3_scale(fc_ij*gijk_d_ij*ex_delr_ij,dcosdrj,fdrj);
  vec3_scaleadd(fc_ij*gijk_ij*ex_delr_d_ij,rij_hat,fdrj,fdrj);
  vec3_scale(prefactor_ij,fdrj,drj);

  vec3_scale(dfc_ij*gijk_ij*ex_delr_ij,rik_hat,fdrk);
  vec3_scaleadd(fc_ij*gijk_d_ij*ex_delr_ij,dcosdrk,fdrk,fdrk);
  vec3_scaleadd(-fc_ij*gijk_ij*ex_delr_d_ij,rik_hat,fdrk,fdrk);
  vec3_scale(prefactor_ij,fdrk,drk);

  double rdrj[3], rdrk[3];
  vec3_scale(fc_ik*gijk_d_ik*ex_delr_ik,dcosdrk,rdrk);
  vec3_scaleadd(fc_ik*gijk_ik*ex_delr_d_ik,rik_hat,rdrk,rdrk);
  vec3_scaleadd(prefactor_ik,rdrk,drk,drk);

  vec3_scale(dfc_ik*gijk_ik*ex_delr_ik,rij_hat,rdrj);
  vec3_scaleadd(fc_ik*gijk_d_ik*ex_delr_ik,dcosdrj,rdrj,rdrj);
  vec3_scaleadd(-fc_ik*gijk_ik*ex_delr_d_ik,rij_hat,rdrj,rdrj);
  vec3_scaleadd(prefactor_ik,rdrj,drj,drj);

  /* vec3_add(fdrj, rdrj, drj); */
  /* vec3_add(fdrk, rdrk, drk); */
  vec3_add(drj, drk, dri);
  vec3_scale(-1, dri, dri);
}

void attractive(tersoff_param_t *param, double prefactor,
                double rsqij, double rsqik,
                double *delrij, double *delrik,
                double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  inv_sqrt(rsqij, rijinv);
  rij = rijinv * rsqij;
  vec3_scale(rijinv,delrij,rij_hat);

  inv_sqrt(rsqik, rikinv);
  rik = rikinv * rsqik;
  vec3_scale(rikinv,delrik,rik_hat);

  ters_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,param);
}

void attractive_p(tersoff_param_t *param, double *prefactorp,
                double *rsqijp, double *rsqikp,
                double *delrij, double *delrik,
                double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  double prefactor = *prefactorp;
  double rsqij = *rsqijp;
  double rsqik = *rsqikp;

  rij = sqrt(rsqij);
  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);

  rik = sqrt(rsqik);
  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  ters_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,param);
}


void attractive_i(tersoff_param_t *param, double *prefactorp,
                  double *rsqijp, double *rsqikp,
                  double *delrij, double *delrik,
                  double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  double prefactor = *prefactorp;
  double rsqij = *rsqijp;
  double rsqik = *rsqikp;

  attractive(param, prefactor, rsqij, rsqik, delrij, delrik, fi, fj, fk);
}


void attractive_jk_dt(tersoff_param_t *param_ij, 
                      tersoff_param_t *param_ik,
                      double *prefactor_ijp, double *prefactor_ikp,
                      double *rsqijp, double *rsqikp,
                      double *delrij, double *delrik,
                      double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  double prefactor_ij = *prefactor_ijp;
  double prefactor_ik = *prefactor_ikp;
  double rsqij = *rsqijp;
  double rsqik = *rsqikp;
  double ffi[3], ffj[3], ffk[3];
  double rfi[3], rfj[3], rfk[3];
  attractive(param_ij, prefactor_ij, rsqij, rsqik, delrij, delrik, ffi, ffj, ffk);
  attractive(param_ik, prefactor_ik, rsqik, rsqij, delrik, delrij, rfi, rfk, rfj);
  vec3_add(ffi, rfi, fi);
  vec3_add(ffj, rfj, fj);
  vec3_add(ffk, rfk, fk);
}

#endif
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


#ifdef MPE
#define LWPF_UNITS U(TERSOFF)
#include "lwpf.h"
#include <simd.h>
extern SLAVE_FUN(pair_tersoff_compute_attractive_para)(void*);

                                                      //#include <unistd.h>
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
  long fend_base = malloc(sizeof(double) * pm->firstshort[pm->ntotal] * 4 + 32);
  long ftmp_base = calloc(pm->nlocal * 4 + 4, sizeof(double));
  pm->fend = (void*)((fend_base + 31) & (~31));
  pm->fdone = calloc(pm->ntotal, sizeof(int));
  double (*ftmp)[4] = (void*)((ftmp_base + 31) & (~31));
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
          simd_load(ftmp_j, ftmp + j0);
          simd_store(fend0 + ftmp_j, ftmp + j0);
        }

        if (j1 < nlocal && jj + 1 < jend){
          doublev4 ftmp_j;
          simd_load(ftmp_j, ftmp + j1);
          simd_store(fend1 + ftmp_j, ftmp + j1);
        }
        if (j2 < nlocal && jj + 2 < jend){
          doublev4 ftmp_j;
          simd_load(ftmp_j, ftmp + j2);
          simd_store(fend2 + ftmp_j, ftmp + j2);
        }
        if (j3 < nlocal && jj + 3 < jend){
          doublev4 ftmp_j;
          simd_load(ftmp_j, ftmp + j3);
          simd_store(fend3 + ftmp_j, ftmp + j3);
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
  free(fend_base);
  free(ftmp_base);
  free(ftmp);
  /* if (++r == 10 && pm->rank == 0) */
  /*   lwpf_finish(stdout); */
}
#endif
#ifdef CPE
/* #define LWPF_KERNELS _K(ALL) K(JLOOP) K(JKLOOP) K(JKLOAD) K(ATTRACTIVE) K(ATTRACTIVE_JK) */
/* #define LWPF_UNIT U(TERSOFF) */
/* #include "lwpf.h" */
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
  double fend[SNSTEP][4], prefactor[SNSTEP];
  int fdone[ISTEP];
  doublev4 virial_v4[2];
  virial_v4[0] = 0.0;
  virial_v4[1] = 0.0;
  double *virial = virial_v4;
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
      /* for (jj = 0; jj < jnum; jj ++){ */
      /*   short_neigh_t *jshort = js + jj; */
      /*   int jtype = jshort->type; */
      /*   int j = jshort->idx; */
      /*   double dij[3]; */
      /*   dij[0] = jshort->d[0]; */
      /*   dij[1] = jshort->d[1]; */
      /*   dij[2] = jshort->d[2]; */
      /*   double r2ij = jshort->r2; */
      /*   int iparam_ij = elem2param[itype][jtype][jtype]; */
      /*   if (r2ij >= params[iparam_ij].cutsq) continue; */
      /*   double zeta_ij = 0.0; */
      /*   short_neigh_t *klist_short = js; */
      /*   int knum = jnum; */

      /*   int kk; */
      /*   for (kk = 0; kk < knum; kk ++){ */
      /*     if (jj == kk) continue; */
      /*     kshort = klist_short + kk; */
      /*     ktype = kshort->type;//map[type[k]]; */
      /*     iparam_ijk = elem2param[itype][jtype][ktype]; */

      /*     dik[0] = kshort->d[0]; */
      /*     dik[1] = kshort->d[1]; */
      /*     dik[2] = kshort->d[2]; */
      /*     r2ik = kshort->r2;//dik[0] * dik[0] + dik[1] * dik[1] + dik[2] * dik[2]; */
      /*     if (r2ik >= params[iparam_ijk].cutsq) continue; */
      /*     zeta_ij += zeta(params + iparam_ijk, r2ij, r2ik, dij, dik); */
      /*   } */
      /*   if (i < nlocal){ */
      /*     fxtmp += dij[0] * fpair; */
      /*     fytmp += dij[1] * fpair; */
      /*     fztmp += dij[2] * fpair; */
       
      /*     if (evflag) ev_tally_full(i, evdwl, 0.0, -fpair, -dij[0], -dij[1], -dij[2]); */
      /*   } */

      /* } */
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
        double prefactor_ij = jshort->prefactor_fwd;

        short_neigh_t *klist_short = js;
        int knum = jnum;

        int kk;
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
          ters_attractive_unroll(prefactor_ij, r2ij, r2ik, dij, dik, tfi, tfj, tfk, params + iparam_ijk);
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
      /* if (i < 64 && rank == 2) */
      /*   printf("%d %f %f %f\n", i, fxtmp, fytmp, fztmp); */

    }
    if (iwr > 0){
      pe_put(f + ist, fi, sizeof(double) * iwr * 3);
    }
    pe_put(l_pm.fdone + ist, fdone, sizeof(int) * isz);
    pe_syn();

  }
  /* int i; */
  /* if (rank == 3){ */
  /*   for (i = 0; i < 2; i ++){ */
  /*     if (_MYID == i) */
  /*       printf("%d %f %f %f\n", i, virial[0], virial[1], virial[2]); */
  /*     athread_syn(ARRAY_SCOPE, 0xffff); */
  /*   } */
  /* } */
  reg_reduce_inplace_doublev4(virial_v4, 2);

  /* if (rank == 3){ */
  /*   for (i = 0; i < 2; i ++){ */
  /*     if (_MYID == i) */
  /*       printf("%d %f %f %f\n", i, virial[0], virial[1], virial[2]); */
  /*     athread_syn(ARRAY_SCOPE, 0xffff); */
  /*   } */
  /* } */
  if (_MYID == 0){
    pe_put(pm->virial, virial, sizeof(double) * 6);
    pe_syn();
  }
  //lwpf_stop(ALL);
}

#endif
/*
inline void attractive_jk(tersoff_param_t *param,
                          double prefactor_ij, double prefactor_ik,
                          double rsqij, double rsqik,
                          double *delrij, double *delrik,
                          double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  inv_sqrt(rsqij, rijinv);
  rij = rsqij * rijinv;
  vec3_scale(rijinv,delrij,rij_hat);

  inv_sqrt(rsqik, rikinv);
  rik = rsqik * rikinv;

  vec3_scale(rikinv,delrik,rik_hat);

  double gijk,gijk_d,cos_theta;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];
 
  double ters_R = param->bigr;
  double ters_D = param->bigd;
  double Dinv = param->bigdinv;
  double fc_ik,dfc_ik, fc_ij, dfc_ij;
  if (rik < ters_R - ters_D){
    fc_ik = 1.0;
    dfc_ik = 0.0;
  } else if (rik > ters_R + ters_D){
    fc_ik = 0.0;
    dfc_ik = 0.0;
  } else {
    fc_ik = 0.5 * (1.0 - sin_4_tersoff(MY_PI2 * (rik - ters_R) * Dinv));
    dfc_ik = -(MY_PI4 * Dinv) * cos_4_tersoff(MY_PI2 * (rik - ters_R) * Dinv);
  }
  if (rij < ters_R - ters_D){
    fc_ij = 1.0;
    dfc_ij = 0.0;
  } else if (rij > ters_R + ters_D){
    fc_ij = 0.0;
    dfc_ij = 0.0;
  } else {
    fc_ij = 0.5 * (1.0 - sin_4_tersoff(MY_PI2 * (rij - ters_R) * Dinv));
    dfc_ij = -(MY_PI4 * Dinv) * cos_4_tersoff(MY_PI2 * (rij - ters_R) * Dinv);
  }
  double tmp = param->lam3 * (rij-rik);
  if (param->powermint == 3) tmp = tmp * tmp * tmp;
  double ex_delr_ij, ex_delr_ik, ex_delr_d_ij, ex_delr_d_ik;
  if (tmp > 69.0776) {
    ex_delr_ik = 1.e30;
    ex_delr_ij = 0;
  } else if (tmp < -69.0776) {
    ex_delr_ik = 0.0;
    ex_delr_ij = 1.e30;
  } else {
    ex_delr_ik = exp_4_tersoff(tmp);
    ex_delr_ij = exp_4_tersoff(-tmp);
  }
  if (param->powermint == 3){
    ex_delr_d_ij = 3.0*param->lam3*param->lam3*param->lam3 *(rij-rik)*(rij-rik)*ex_delr_ij;
    ex_delr_d_ik = 3.0*param->lam3*param->lam3*param->lam3 *(rij-rik)*(rij-rik)*ex_delr_ik;
  }
  else {
    ex_delr_d_ij = param->lam3 * ex_delr_ij;
    ex_delr_d_ik = param->lam3 * ex_delr_ik;
  }
  cos_theta = vec3_dot(rij_hat,rik_hat);

  double ters_c = param->c * param->c;
  double ters_d = param->d * param->d;
  double hcth = param->h - cos_theta;
  double numerator = -2.0 * ters_c * hcth;
  double denominator = 1.0 / (ters_d + hcth * hcth);
  gijk = param->gamma * (1.0 + param->c2divd2 - ters_c * denominator);
  gijk_d = param->gamma * numerator * denominator * denominator;

  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,dcosdrj);
  vec3_scale(rijinv,dcosdrj,dcosdrj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,dcosdrk);
  vec3_scale(rikinv,dcosdrk,dcosdrk);

  double fjtmp[3], fktmp[3];
  vec3_scale(fc_ik*gijk_d*ex_delr_ik,dcosdrj,fjtmp);
  vec3_scaleadd(fc_ik*gijk*ex_delr_d_ik,rij_hat,fjtmp,fjtmp);
  vec3_scale(prefactor_ij,fjtmp,fj);

  vec3_scale(dfc_ik*gijk*ex_delr_ik,rik_hat,fktmp);
  vec3_scaleadd(fc_ik*gijk_d*ex_delr_ik,dcosdrk,fktmp,fktmp);
  vec3_scaleadd(-fc_ik*gijk*ex_delr_d_ik,rik_hat,fktmp,fktmp);
  vec3_scale(prefactor_ij,fktmp,fk);

  vec3_scale(fc_ij*gijk_d*ex_delr_ij,dcosdrk,fktmp);
  vec3_scaleadd(fc_ij*gijk*ex_delr_d_ij,rik_hat,fktmp,fktmp);
  vec3_scaleadd(prefactor_ik,fktmp,fk,fk);

  vec3_scale(dfc_ij*gijk*ex_delr_ij,rik_hat,fjtmp);
  vec3_scaleadd(fc_ij*gijk_d*ex_delr_ij,dcosdrj,fjtmp,fjtmp);
  vec3_scaleadd(-fc_ij*gijk*ex_delr_d_ij,rij_hat,fjtmp,fjtmp);
  vec3_scaleadd(prefactor_ik,fjtmp,fj,fj);

}

 */
