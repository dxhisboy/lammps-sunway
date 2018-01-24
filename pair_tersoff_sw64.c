#include "sunway.h"
#include "pair_tersoff_sw64.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef CPE
__thread_local rank;
#endif
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

inline void sincos_4_tersoff(double x, double *sin, double *cos){
  double x2 = x * x;
  double cosx = x2 * COS_PPC6 + COS_PPC5;
  double sinx = x2 * SIN_PPC5 + SIN_PPC4;
  cosx = cosx * x2 + COS_PPC4;
  sinx = sinx * x2 + SIN_PPC3;
  cosx = cosx * x2 + COS_PPC3;
  sinx = sinx * x2 + SIN_PPC2;
  cosx = cosx * x2 + COS_PPC2;
  sinx = sinx * x2 + SIN_PPC1;
  cosx = cosx * x2 + COS_PPC1;
  sinx = sinx * x2 + SIN_PPC0;
  *cos = cosx * x2 + COS_PPC0;
  *sin = sinx * x;
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

double ters_fa(double r, tersoff_param_t *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return -param->bigb * exp(-param->lam2 * r) * ters_fc(r,param);
}

/* ---------------------------------------------------------------------- */

double ters_fa_d(double r, tersoff_param_t *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return param->bigb * exp(-param->lam2 * r) *
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
                          (1.0 - (1.0 +  1.0/(2.0*param->powern)) *
                           pow(tmp,-param->powern)));
  if (tmp < param->c4) return 0.0;
  if (tmp < param->c3)
    return -0.5*param->beta * pow(tmp,param->powern-1.0);

  double tmp_n = pow(tmp,param->powern);
  return -0.5 * pow(1.0+tmp_n, -1.0-(1.0/(2.0*param->powern)))*tmp_n / zeta;
}

double zeta(tersoff_param_t *param, double rsqij, double rsqik,
            double *delrij, double *delrik)
{
  double rij,rik,costheta,arg,ex_delr;
  double rijinv, rikinv;
  inv_sqrt(rsqij, rijinv);
  inv_sqrt(rsqik, rikinv);
  rij = rijinv * rsqij;
  rik = rikinv * rsqik;
  /* rij = sqrt(rsqij); */
  /* rik = sqrt(rsqik); */
  costheta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] +
              delrij[2]*delrik[2]) * (rijinv*rikinv);

  arg = param->lam3 * (rij - rik);
  if (param->powermint == 3) arg = arg * arg * arg;
  //else arg = param->lam3 * (rij-rik);

  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp_4_tersoff(arg);

  return ters_fc(rik,param) * ters_gijk(costheta,param) * ex_delr;
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
  /* rij = sqrt(rsqij); */
  /* rik = sqrt(rsqik); */
  costheta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] +
              delrij[2]*delrik[2]) * (rijinv*rikinv);

  arg = param->lam3 * (rij - rik);
  if (param->powermint == 3) arg = arg * arg * arg;
  //else arg = param->lam3 * (rij-rik);

  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp_4_tersoff(arg);
  
  double ters_R = param->bigr;
  double ters_D = param->bigd;
  double fc;
  if (rik < ters_R-ters_D) fc = 1.0;
  else if (rik > ters_R+ters_D) fc = 0.0;
  else fc =  0.5*(1.0 - sin_4_tersoff(MY_PI2*(rik - ters_R)/ters_D));
 
  double gijk;
  const double ters_c = param->c * param->c;
  const double ters_d = param->d * param->d;
  const double hcth = param->h - costheta;
  gijk =  param->gamma*(1.0 + ters_c/ters_d - ters_c / (ters_d + hcth*hcth));

  return fc * gijk * ex_delr;
  //return ters_fc(rik,param) * ters_gijk(costheta,param) * ex_delr;
}


/* ---------------------------------------------------------------------- */

void force_zeta(tersoff_param_t *param, double rsq, double zeta_ij,
                double *fforce, double *prefactor,
                int eflag, double *eng)
{
  double r,fa,fa_d,bij;

  r = sqrt(rsq);
  fa = ters_fa(r,param);
  fa_d = ters_fa_d(r,param);
  bij = ters_bij(zeta_ij,param);
  *fforce = 0.5*bij*fa_d / r;
  *prefactor = -0.5*fa * ters_bij_d(zeta_ij,param);
  if (eflag) *eng = 0.5*bij*fa;
}

void force_zeta_unroll(tersoff_param_t *param, double rsq, double zeta_ij,
                       double *fforce, double *prefactor,
                       int eflag, double *eng)
{
  double r,fa,fa_d,bij, bij_d, fc, fc_d;
  double rinv;
  inv_sqrt(rsq, rinv);
  r = rinv * rsq;
  /* fa = ters_fa(r,param); */
  /* fa_d = ters_fa_d(r,param); */

  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r > ters_R + ters_D){
    fa = 0.0;
    fa_d = 0.0;
  } else {
    double er = exp(-param->lam2 * r);
    /* double Dinv = param->bigdinv; */
    /* double fc, dfc; */
    /* if (r < ters_R - ters_D){ */
    /*   fc = 1.0; */
    /*   dfc = 0.0; */
    /* } else if (r > ters_R + ters_D){ */
    /*   fc = dfc = 0.0; */
    /* } else { */
    /*   double cosr, sinr; */
    /*   /\* sincos_4_tersoff(MY_PI2 * (r - ters_R) / ters_D, &cosr, &sinr); *\/ */
    /*   /\* fc = 0.5 * (1.0 - sinr); *\/ */
    /*   /\* dfc = -(MY_PI4 / ters_D) * cosr; *\/ */
    /*   fc = 0.5*(1.0 - sin_4_tersoff(MY_PI2*(r-ters_R)/ters_D)); */
    /*   dfc = -(MY_PI4/ters_D)*cos_4_tersoff(MY_PI2*(r-ters_R)/ters_D); */
    /* } */

    /* if (fc - ters_fc(r, param) != 0) */
    /*   puts("error"); */
    fa = -param->bigb * er * ters_fc(r, param);
    fa_d = param->bigb * er *
      //(param->lam2 * fc - dfc);//ters_fc_d(r, param));
      (param->lam2 * ters_fc(r, param) - ters_fc_d(r, param));
    /* fa = -param->bigb * er * fc; */
    /* double lam2fc = param->lam2 * fc; */
    /* fa_d = param->bigb * er * (lam2fc - dfc); */
  }
  /* double tmp = param->beta * zeta_ij; */
  /* double tmp_n = pow(tmp, param->powern); */
  /* double tmp_p = pow(1.0 + tmp_n, -1.0-(1.0/(2.0*param->powern))); */
  /* bij = pow(1.0 + tmp_n, -(1.0/(2.0*param->powern))); */
  /* //bij = tmp_pp1; */
  /* bij_d = -0.5 * tmp_p * tmp_n / zeta_ij; */
  bij = ters_bij(zeta_ij, param);
  bij_d = ters_bij_d(zeta_ij, param);
  *fforce = 0.5*bij*fa_d * rinv;
  *prefactor = -0.5*fa * bij_d;
  if (eflag) *eng = 0.5*bij*fa;
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
  double cdivd = param->c2divd2;
  double hcth = param->h - cos_theta;
  double numerator = -2.0 * ters_c * hcth;
  double denominator = 1.0/(ters_d + hcth*hcth);
  gijk = param->gamma * (1.0 + cdivd - ters_c * denominator);
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
  double cdivd_ij = param_ijk->c2divd2;
  double cdivd_ik = param_ikj->c2divd2;
  double gijk_ij = param_ijk->gamma * (1.0 + cdivd_ij - ters_c_ij*denominator_ij);
  double gijk_ik = param_ikj->gamma * (1.0 + cdivd_ik - ters_c_ik*denominator_ik);
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

//#endif
#define THIRD 0.33333333333333333333333
void v_tally3rd(int i, int vflag_global, int vflag_atom,
                double *fi, double *fj, double *drik, double *drjk, double *virial, double (*vatom)[6])
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

static void ev_tally_full(pair_tersoff_compute_param_t *pm,
                          int i, double evdwl, double ecoul, double fpair,
                          double delx, double dely, double delz,
                          double *eng_coul, double *eng_vdwl, double *virial,
                          double *eatom, double (*vatom)[6])
{
  double v[6];
  if (pm->eflag_either) 
  {
    if (pm->eflag_global) 
    {
      eng_vdwl[0] += 0.5*evdwl;
      eng_coul[0] += 0.5*ecoul;
    }
    if (pm->eflag_atom) eatom[i] += 0.5 * (evdwl + ecoul);
  }//if-eflag-either
  if (pm->vflag_either) 
  {
    v[0] = 0.5*delx*delx*fpair;
    v[1] = 0.5*dely*dely*fpair;
    v[2] = 0.5*delz*delz*fpair;
    v[3] = 0.5*delx*dely*fpair;
    v[4] = 0.5*delx*delz*fpair;
    v[5] = 0.5*dely*delz*fpair;
    if (pm->vflag_global) 
    {
      virial[0] += v[0];
      virial[1] += v[1];
      virial[2] += v[2];
      virial[3] += v[3];
      virial[4] += v[4];
      virial[5] += v[5];
    }
    if (pm->vflag_atom) 
    {
      vatom[i][0] += v[0];
      vatom[i][1] += v[1];
      vatom[i][2] += v[2];
      vatom[i][3] += v[3];
      vatom[i][4] += v[4];
      vatom[i][5] += v[5];
    }
  }//if-vflag-eigher
}


#ifdef MPE
#define LWPF_UNITS U(TERSOFF)
#include "lwpf.h"
#include <athread.h>
extern SLAVE_FUN(pair_tersoff_compute_attractive_para)(void*);
extern SLAVE_FUN(pair_tersoff_compute_zeta_para)(void*);


void pair_tersoff_compute_attractive(pair_tersoff_compute_param_t *pm)
{
  if (athread_idle() == 0)
    athread_init();
  athread_spawn(pair_tersoff_compute_attractive_para, pm);
  athread_join();
}
void pair_tersoff_compute_zeta(pair_tersoff_compute_param_t *pm)
{
  if (athread_idle() == 0)
    athread_init();
  athread_spawn(pair_tersoff_compute_zeta_para, pm);
  athread_join();
}
#endif
#ifdef CPE
#define LWPF_KERNELS _K(ALL) K(JLOOP) K(JKLOOP) K(JKLOAD) K(ATTRACTIVE) K(ATTRACTIVE_JK)
#define LWPF_UNIT U(TERSOFF)
#include "lwpf.h"
#include "dma.h"
#include <simd.h>
#include<math.h>
#define ISTEP 64
#define SNSTEP 128

/************version3***unroll-evtally-zeta************/
void pair_tersoff_compute_zeta_para(pair_tersoff_compute_param_t *pm)
{
  pe_init();
  //lwpf_start(ALL);
  pair_tersoff_compute_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(pair_tersoff_compute_param_t));
  pe_syn();
  rank = l_pm.rank;
  int inum      = l_pm.inum;
  int nlocal    = l_pm.nlocal;
  int nghost    = l_pm.nghost;
  int allnum    = nlocal + nghost;
  int nelements = l_pm.nelements;
  int nep3      = nelements * nelements * nelements;
  int ntypes    = l_pm.ntypes;
  int nparams   = l_pm.nparams;
  int map[ntypes + 1];
  int elem2param[nelements][nelements][nelements];
  tersoff_param_t params[nparams];
  pe_get(l_pm.map, map, sizeof(int) * (ntypes + 1));
  pe_get(l_pm.elem2param, elem2param, sizeof(int) * nep3);
  pe_get(l_pm.params, params, sizeof(tersoff_param_t) * nparams);
  pe_syn();

  int *type           = l_pm.type;
  int *firstshort     = l_pm.firstshort;
  double (*vatom)[6]  = l_pm.vatom;
  double (*f)[3]      = l_pm.f;
  double *eatom       = l_pm.eatom;

  short_neigh_t *shortlist = l_pm.shortlist;

  int eflag        = l_pm.eflag;
  int vflag        = l_pm.vflag;
  int evflag       = l_pm.evflag;
  int eflag_global = l_pm.eflag_global;
  int vflag_global = l_pm.vflag_global;
  int eflag_atom   = l_pm.eflag_atom;
  int vflag_atom   = l_pm.vflag_atom;
  int eflag_either = l_pm.eflag_either;
  int vflag_either = l_pm.vflag_either;
  
  doublev4 eng_virial[2];
  eng_virial[0] = 0;
  eng_virial[1] = eng_virial[0];
  double *eng_vdwl = (double*)(void*)(eng_virial);
  double *eng_coul = eng_vdwl + 1;
  double *virial = eng_coul + 1;

  int ist, ied;
  double fi[ISTEP][3], ei[ISTEP], vi[ISTEP][6];
  int ti[ISTEP], fs[ISTEP + 1];
  short_neigh_t js[SNSTEP], ks[SNSTEP];
  for(ist = _MYID * ISTEP; ist < allnum; ist += ISTEP * 64)
  {
    ied = ist + ISTEP;
    if(ied > allnum)
      ied = allnum;
    int isz = ied - ist;
    int i;
    pe_get(f + ist, fi, sizeof(double) * 3 * isz);
    pe_get(firstshort + ist, fs, sizeof(int) * (isz + 1));//fs<=>firstshort;
    pe_get(type + ist, ti, sizeof(int) * isz);
    pe_syn();

    for (i = ist; i < ied; i ++)
    {
      int ioff = i - ist;
      int itype = map[ti[ioff]];
      double fxtmp = 0;
      double fytmp = 0;
      double fztmp = 0;
      short_neigh_t *jlist_short = shortlist + fs[ioff];
      int jnum = fs[ioff + 1] - fs[ioff];
      pe_get(jlist_short, js, sizeof(short_neigh_t) * jnum);
      pe_syn();
      
      ei[ioff] = 0;
      vi[ioff][0] = vi[ioff][1] = vi[ioff][2] = 0;
      vi[ioff][3] = vi[ioff][4] = vi[ioff][5] = 0;

      int jj, kk;
      for (jj = 0; jj < jnum; jj ++)
      {
        short_neigh_t *jshort = js + jj;
        int j = jshort->idx;
        int jtype = jshort->type;
        int iparam_ij = elem2param[itype][jtype][jtype];
        
        double dij[3], dji[3];
        dij[0] = jshort->d[0];
        dij[1] = jshort->d[1];
        dij[2] = jshort->d[2];
        dji[0] = -dij[0];
        dji[1] = -dij[1];
        dji[2] = -dij[2];
        
        double r2ij = jshort->r2;
        if (r2ij >= params[iparam_ij].cutsq) continue;
        
        double zeta_ij, zeta_ji;
        zeta_ij = zeta_ji = 0.0;
        
        short_neigh_t *klist_short = shortlist + fs[ioff];
        int iknum = fs[ioff + 1] - fs[ioff];
        pe_get(klist_short, ks, sizeof(short_neigh_t) * iknum);
        pe_syn();

        //compute zeta_ij
        for (kk = 0; kk < iknum; kk ++)
        {
          if (jj == kk) continue;
          short_neigh_t *kshort = ks + kk;
          int ktype = kshort->type;
          int iparam_ijk = elem2param[itype][jtype][ktype];

          double dik[3];
          dik[0] = kshort->d[0];
          dik[1] = kshort->d[1];
          dik[2] = kshort->d[2];
          double r2ik = kshort->r2;
          if (r2ik >= params[iparam_ijk].cutsq) continue;
          //lwpf_start(ZETA);
          zeta_ij += zeta_unroll(params + iparam_ijk, r2ij, r2ik, dij, dik);
          //lwpf_stop(ZETA);
        }//for-kk

        //lwpf_start(JKLOAD);
        int fsj[2];
        pe_get(firstshort + j, fsj, sizeof(int)*2);
        pe_syn();
        short_neigh_t *jklist_short = shortlist + fsj[0];
        int jknum = fsj[1] - fsj[0];
        pe_get(jklist_short, ks, sizeof(short_neigh_t) * jknum);
        pe_syn();

        //lwpf_stop(JKLOAD);
        //lwpf_start(JKLOOP);

        for (kk = 0; kk < jknum; kk ++)
        {
          short_neigh_t *kshort = ks + kk;
          if (kshort->idx == i) continue;
          int ktype = kshort->type;
          int iparam_jik = elem2param[jtype][itype][ktype];
          int iparam_jki = elem2param[jtype][ktype][itype];

          double djk[3];
          djk[0] = kshort->d[0];
          djk[1] = kshort->d[1];
          djk[2] = kshort->d[2];
          double r2jk = kshort->r2;

          if (r2jk >= params[iparam_jik].cutsq) continue;
          //lwpf_start(ZETA_JK);
          zeta_ji += zeta_unroll(params + iparam_jik, r2ij, r2jk, dji, djk);
          //lwpf_stop(ZETA_JK);
        }//for-kk

        double fpair, evdwl, prefactor_ij, prefactor_ji;
        force_zeta_unroll(params + iparam_ij, r2ij, zeta_ij,
                          &fpair, &prefactor_ij, eflag, &evdwl);
        //fpair = prefactor_ij = evdwl = 0;
        if (i < nlocal)
        {
          fxtmp += dij[0] * fpair;
          fytmp += dij[1] * fpair;
          fztmp += dij[2] * fpair;
          if(evflag)
          {
            //ev_tally_full(&l_pm, ioff, evdwl, 0.0, -fpair, -dij[0],-dij[1],-dij[2], 
            //              eng_coul, eng_vdwl, virial, ei, vi);
            //ev_tally_full(pm, i, evdwl, 0.0, -fpair, -dij[0], -dij[1], -dij[2], 
            //              eng_coul, eng_vdwl, virial, eatom, vatom);
            double v[6];
            double delx = -dij[0], dely = -dij[1], delz = -dij[2];
            if (l_pm.eflag_either) 
            {
              if (l_pm.eflag_global) 
              {
                eng_vdwl[0] += 0.5*evdwl;
              }
              if (l_pm.eflag_atom) ei[ioff] += 0.5 * (evdwl);
            }//if-eflag-either
            if (l_pm.vflag_either) 
            {
              v[0] = 0.5*delx*delx*(-fpair);
              v[1] = 0.5*dely*dely*(-fpair);
              v[2] = 0.5*delz*delz*(-fpair);
              v[3] = 0.5*delx*dely*(-fpair);
              v[4] = 0.5*delx*delz*(-fpair);
              v[5] = 0.5*dely*delz*(-fpair);
              if (l_pm.vflag_global) 
              {
                virial[0] += v[0];
                virial[1] += v[1];
                virial[2] += v[2];
                virial[3] += v[3];
                virial[4] += v[4];
                virial[5] += v[5];
              }
              if (l_pm.vflag_atom) 
              {
                vi[ioff][0] += v[0];
                vi[ioff][1] += v[1];
                vi[ioff][2] += v[2];
                vi[ioff][3] += v[3];
                vi[ioff][4] += v[4];
                vi[ioff][5] += v[5];
              }
            }//if-vflag-eigher
            
          }//if-evflag
        }//if-nlocal
      
        force_zeta_unroll(params + iparam_ij, r2ij, zeta_ji, 
                          &fpair, &prefactor_ji, eflag, &evdwl);
        //fpair = prefactor_ji = evdwl = 0;
        if (i < nlocal)
        {
          fxtmp += dij[0] * fpair;
          fytmp += dij[1] * fpair;
          fztmp += dij[2] * fpair;
          if (evflag)
          {
            //ev_tally_full(&l_pm, ioff, evdwl, 0.0, -fpair, -dij[0],-dij[1],-dij[2], 
            //              eng_coul, eng_vdwl, virial, ei, vi);
            //ev_tally_full(pm, i, evdwl, 0.0, -fpair, -dij[0], -dij[1], -dij[2], 
            //              eng_coul, eng_vdwl, virial, eatom, vatom);            
            double v[6];
            double delx = -dij[0], dely = -dij[1], delz = -dij[2];
            if (l_pm.eflag_either) 
            {
              if (l_pm.eflag_global) 
              {
                eng_vdwl[0] += 0.5*evdwl;
              }
              if (l_pm.eflag_atom) ei[ioff] += 0.5 * (evdwl);
            }//if-eflag-either
            if (l_pm.vflag_either) 
            {
              v[0] = 0.5*delx*delx*(-fpair);
              v[1] = 0.5*dely*dely*(-fpair);
              v[2] = 0.5*delz*delz*(-fpair);
              v[3] = 0.5*delx*dely*(-fpair);
              v[4] = 0.5*delx*delz*(-fpair);
              v[5] = 0.5*dely*delz*(-fpair);
              if (l_pm.vflag_global) 
              {
                virial[0] += v[0];
                virial[1] += v[1];
                virial[2] += v[2];
                virial[3] += v[3];
                virial[4] += v[4];
                virial[5] += v[5];
              }
              if (l_pm.vflag_atom) 
              {
                vi[ioff][0] += v[0];
                vi[ioff][1] += v[1];
                vi[ioff][2] += v[2];
                vi[ioff][3] += v[3];
                vi[ioff][4] += v[4];
                vi[ioff][5] += v[5];
              }
            }//if-vflag-eigher
  
          }//if-evflag
        }//if-nlocal

        jshort->prefactor_fwd = prefactor_ij;
        jshort->prefactor_rev = prefactor_ji;
        //lwpf_stop(JKLOOP);

      }//for-jj
      pe_put(jlist_short, js, sizeof(short_neigh_t) * jnum);
      pe_syn();

      if (i < nlocal)
      {
        fi[ioff][0] += fxtmp;
        fi[ioff][1] += fytmp;
        fi[ioff][2] += fztmp;
      }
    }//for-i
    pe_put(f + ist, fi,  sizeof(double) * 3 * isz);
    if (eflag_either && eflag_atom)
    {
      pe_put(eatom + ist, ei,  sizeof(double) * isz);
    }
    if(vflag_either && vflag_atom)
    {
      pe_put(vatom + ist, vi,  sizeof(double) * 6 * isz);
    }
    pe_syn();
  }//for-ist
  reg_reduce_inplace_doublev4(eng_virial, 2);
  if (_MYID == 0)
  {
    pe_put(&(pm->eng_vdwl), eng_vdwl, sizeof(double) * 8);
    pe_syn();
  }
  //lwpf_stop(ALL);
}


/************version2***athread_get**1.45 or put:1.4148*************/
/*
void pair_tersoff_compute_zeta_para(pair_tersoff_compute_param_t *pm)
{
  pe_init();
  //lwpf_start(ALL);
  pair_tersoff_compute_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(pair_tersoff_compute_param_t));
  pe_syn();

  int inum      = l_pm.inum;
  int nlocal    = l_pm.nlocal;
  int nghost    = l_pm.nghost;
  int allnum    = nlocal + nghost;
  int nelements = l_pm.nelements;
  int nep3      = nelements * nelements * nelements;
  int ntypes    = l_pm.ntypes;
  int nparams   = l_pm.nparams;
  int map[ntypes + 1];
  int elem2param[nelements][nelements][nelements];
  tersoff_param_t params[nparams];
  pe_get(l_pm.map, map, sizeof(int) * (ntypes + 1));
  pe_get(l_pm.elem2param, elem2param, sizeof(int) * nep3);
  pe_get(l_pm.params, params, sizeof(tersoff_param_t) * nparams);
  pe_syn();

  int *type           = l_pm.type;
  int *firstshort     = l_pm.firstshort;
  double (*vatom)[6]  = l_pm.vatom;
  double (*f)[3]      = l_pm.f;
  double *eatom       = l_pm.eatom;

  short_neigh_t *shortlist = l_pm.shortlist;

  int eflag        = l_pm.eflag;
  int vflag        = l_pm.vflag;
  int evflag       = l_pm.evflag;
  int eflag_global = l_pm.eflag_global;
  int vflag_global = l_pm.vflag_global;
  int eflag_atom   = l_pm.eflag_atom;
  int vflag_atom   = l_pm.vflag_atom;
  int eflag_either = l_pm.eflag_either;
  int vflag_either = l_pm.vflag_either;
  
  doublev4 eng_virial[2];
  eng_virial[0] = 0;
  eng_virial[1] = eng_virial[0];
  double *eng_vdwl = (double*)(void*)(eng_virial);
  double *eng_coul = eng_vdwl + 1;
  double *virial = eng_coul + 1;

  int ist, ied;
  double fi[ISTEP][3], ei[ISTEP], vi[ISTEP][6];
  int ti[ISTEP], fs[ISTEP + 1];
  short_neigh_t js[SNSTEP], ks[SNSTEP];
  for(ist = _MYID * ISTEP; ist < allnum; ist += ISTEP * 64)
  {
    ied = ist + ISTEP;
    if(ied > allnum)
      ied = allnum;
    int isz = ied - ist;
    int i;
    pe_get(f + ist, fi, sizeof(double) * 3 * isz);
    pe_get(firstshort + ist, fs, sizeof(int) * (isz + 1));//fs<=>firstshort;
    pe_get(type + ist, ti, sizeof(int) * isz);
    pe_syn();

    for (i = ist; i < ied; i ++)
    {
      int ioff = i - ist;
      int itype = map[ti[ioff]];
      double fxtmp = 0;
      double fytmp = 0;
      double fztmp = 0;
      short_neigh_t *jlist_short = shortlist + fs[ioff];
      int jnum = fs[ioff + 1] - fs[ioff];
      pe_get(jlist_short, js, sizeof(short_neigh_t) * jnum);
      pe_syn();
      
      ei[ioff] = 0;
      vi[ioff][0] = vi[ioff][1] = vi[ioff][2] = 0;
      vi[ioff][3] = vi[ioff][4] = vi[ioff][5] = 0;

      int jj, kk;
      for (jj = 0; jj < jnum; jj ++)
      {
        short_neigh_t *jshort = js + jj;
        int j = jshort->idx;
        int jtype = jshort->type;
        int iparam_ij = elem2param[itype][jtype][jtype];
        
        double dij[3], dji[3];
        dij[0] = jshort->d[0];
        dij[1] = jshort->d[1];
        dij[2] = jshort->d[2];
        dji[0] = -dij[0];
        dji[1] = -dij[1];
        dji[2] = -dij[2];
        
        double r2ij = jshort->r2;
        if (r2ij >= params[iparam_ij].cutsq) continue;
        
        double zeta_ij, zeta_ji;
        zeta_ij = zeta_ji = 0.0;
        
        short_neigh_t *klist_short = shortlist + fs[ioff];
        int iknum = fs[ioff + 1] - fs[ioff];
        pe_get(klist_short, ks, sizeof(short_neigh_t) * iknum);
        pe_syn();

        //compute zeta_ij
        for (kk = 0; kk < iknum; kk ++)
        {
          if (jj == kk) continue;
          short_neigh_t *kshort = ks + kk;
          int ktype = kshort->type;
          int iparam_ijk = elem2param[itype][jtype][ktype];

          double dik[3];
          dik[0] = kshort->d[0];
          dik[1] = kshort->d[1];
          dik[2] = kshort->d[2];
          double r2ik = kshort->r2;
          if (r2ik >= params[iparam_ijk].cutsq) continue;
          //lwpf_start(ZETA);
          zeta_ij += zeta(params + iparam_ijk, r2ij, r2ik, dij, dik);
          //lwpf_stop(ZETA);
        }//for-kk

        //lwpf_start(JKLOAD);
        int fsj[2];
        pe_get(firstshort + j, fsj, sizeof(int)*2);
        pe_syn();
        short_neigh_t *jklist_short = shortlist + fsj[0];
        int jknum = fsj[1] - fsj[0];
        pe_get(jklist_short, ks, sizeof(short_neigh_t) * jknum);
        pe_syn();

        //lwpf_stop(JKLOAD);
        //lwpf_start(JKLOOP);

        for (kk = 0; kk < jknum; kk ++)
        {
          short_neigh_t *kshort = ks + kk;
          if (kshort->idx == i) continue;
          int ktype = kshort->type;
          int iparam_jik = elem2param[jtype][itype][ktype];
          int iparam_jki = elem2param[jtype][ktype][itype];

          double djk[3];
          djk[0] = kshort->d[0];
          djk[1] = kshort->d[1];
          djk[2] = kshort->d[2];
          double r2jk = kshort->r2;

          if (r2jk >= params[iparam_jik].cutsq) continue;
          //lwpf_start(ZETA_JK);
          zeta_ji += zeta(params + iparam_jik, r2ij, r2jk, dji, djk);
          //lwpf_stop(ZETA_JK);
        }//for-kk

        double fpair, evdwl, prefactor_ij, prefactor_ji;
        force_zeta(params + iparam_ij, r2ij, zeta_ij,
                          &fpair, &prefactor_ij, eflag, &evdwl);

        if (i < nlocal)
        {
          fxtmp += dij[0] * fpair;
          fytmp += dij[1] * fpair;
          fztmp += dij[2] * fpair;
          if(evflag)
          {
            ev_tally_full(&l_pm, ioff, evdwl, 0.0, -fpair, -dij[0],-dij[1],-dij[2], 
                          eng_coul, eng_vdwl, virial, ei, vi);
            //ev_tally_full(pm, i, evdwl, 0.0, -fpair, -dij[0], -dij[1], -dij[2], 
            //              eng_coul, eng_vdwl, virial, eatom, vatom);
          }//if-evflag
        }//if-nlocal
      
        force_zeta(params + iparam_ij, r2ij, zeta_ji, 
                          &fpair, &prefactor_ji, eflag, &evdwl);

        if (i < nlocal)
        {
          fxtmp += dij[0] * fpair;
          fytmp += dij[1] * fpair;
          fztmp += dij[2] * fpair;
          if (evflag)
          {
            ev_tally_full(&l_pm, ioff, evdwl, 0.0, -fpair, -dij[0],-dij[1],-dij[2], 
                          eng_coul, eng_vdwl, virial, ei, vi);
            //ev_tally_full(pm, i, evdwl, 0.0, -fpair, -dij[0], -dij[1], -dij[2], 
            //              eng_coul, eng_vdwl, virial, eatom, vatom);            
          }//if-evflag
        }//if-nlocal

        jshort->prefactor_fwd = prefactor_ij;
        jshort->prefactor_rev = prefactor_ji;
        //lwpf_stop(JKLOOP);

      }//for-jj
      pe_put(jlist_short, js, sizeof(short_neigh_t) * jnum);
      pe_syn();

      if (i < nlocal)
      {
        fi[ioff][0] += fxtmp;
        fi[ioff][1] += fytmp;
        fi[ioff][2] += fztmp;
      }
    }//for-i
    pe_put(f + ist, fi,  sizeof(double) * 3 * isz);
    if (eflag_either && eflag_atom)
    {
      pe_put(eatom + ist, ei,  sizeof(double) * isz);
    }
    if(vflag_either && vflag_atom)
    {
      pe_put(vatom + ist, vi,  sizeof(double) * 6 * isz);
    }
    pe_syn();
  }//for-ist
  reg_reduce_inplace_doublev4(eng_virial, 2);
  if (_MYID == 0)
  {
    pe_put(&(pm->eng_vdwl), eng_vdwl, sizeof(double) * 8);
    pe_syn();
  }
  //lwpf_stop(ALL);
}
*/

/**********zeta_para---version1---Partition**6.0777********/
/*
void pair_tersoff_compute_zeta_para(pair_tersoff_compute_param_t *pm)
{
  pe_init();
  int inum      = pm->inum;
  int nlocal    = pm->nlocal;
  int nghost    = pm->nghost;
  int allnum    = nlocal + nghost;
  int nelements = pm->nelements;
  int nep3 = nelements * nelements * nelements;
  int ntypes    = pm->ntypes;
  int nparams   = pm->nparams;
  int *map      = pm->map;
  int (*elem2param)[nelements][nelements] = pm->elem2param;
  tersoff_param_t *params = pm->params;
  
  int *type           = pm->type;
  int *firstshort     = pm->firstshort;
  double (*vatom)[6]  = pm->vatom;
  double (*f)[3]      = pm->f;
  double *eatom       = pm->eatom;
  short_neigh_t *shortlist = pm->shortlist;

  int eflag        = pm->eflag;
  int vflag        = pm->vflag;
  int evflag       = pm->evflag;
  int eflag_global = pm->eflag_global;
  int vflag_global = pm->vflag_global;
  int eflag_atom   = pm->eflag_atom;
  int vflag_atom   = pm->vflag_atom;
  int eflag_either = pm->eflag_either;
  int vflag_either = pm->vflag_either;
  
  doublev4 eng_virial[2];
  eng_virial[0] = 0;
  eng_virial[1] = eng_virial[0];
  double *eng_vdwl = (double*)(void*)(eng_virial);
  double *eng_coul = eng_vdwl + 1;
  double *virial = eng_coul + 1;

  int ist, ied;
  short_neigh_t js[SNSTEP];
  for(ist = _MYID * ISTEP; ist < allnum; ist += ISTEP * 64)
  {
    ied = ist + ISTEP;
    if(ied > allnum)
      ied = allnum;
    int isz = ied - ist;
    int i;
    for (i = ist; i < ied; i ++)
    {
      int itype = map[type[i]];
      double fxtmp, fytmp, fztmp;
      fxtmp = fytmp = fztmp = 0.0;
      int jnum = firstshort[i + 1] - firstshort[i];
      short_neigh_t *jlist_short = shortlist + firstshort[i];

      int jj;
      for (jj = 0; jj < jnum; jj ++)
      {
        short_neigh_t *jshort = jlist_short + jj;
        int j = jshort->idx;
        int jtype = jshort->type;
        int iparam_ij = elem2param[itype][jtype][jtype];

        double dij[3], dji[3];
        dij[0] = jshort->d[0];
        dij[1] = jshort->d[1];
        dij[2] = jshort->d[2];
        dji[0] = -dij[0];
        dji[1] = -dij[1];
        dji[2] = -dij[2];

        double r2ij = jshort->r2;
        if (r2ij >= params[iparam_ij].cutsq) continue;

        double zeta_ij, zeta_ji;
        zeta_ij = zeta_ji = 0.0;

        short_neigh_t *klist_short = shortlist + firstshort[i];
        int knum = firstshort[i + 1] - firstshort[i];

        int kk;
        //compute zeta_ij
        for (kk = 0; kk < knum; kk ++)
        {
          if (jj == kk) continue;
          short_neigh_t *kshort = klist_short + kk;
          int ktype = kshort->type;
          int iparam_ijk = elem2param[itype][jtype][ktype];

          double dik[3];
          dik[0] = kshort->d[0];
          dik[1] = kshort->d[1];
          dik[2] = kshort->d[2];
          double r2ik = kshort->r2;
          if (r2ik >= params[iparam_ijk].cutsq) continue;
          zeta_ij += zeta(params + iparam_ijk, r2ij, r2ik, dij, dik);
        }//for-kk
        klist_short = shortlist + firstshort[j];
        knum = firstshort[j + 1] - firstshort[j];

        for (kk = 0; kk < knum; kk ++)
        {
          short_neigh_t *kshort = klist_short + kk;
          if (kshort->idx == i) continue;
          int ktype = kshort->type;
          int iparam_jik = elem2param[jtype][itype][ktype];
          int iparam_jki = elem2param[jtype][ktype][itype];

          double djk[3];
          djk[0] = kshort->d[0];
          djk[1] = kshort->d[1];
          djk[2] = kshort->d[2];
          double r2jk = kshort->r2;

          if (r2jk >= params[iparam_jik].cutsq) continue;
          zeta_ji += zeta(params + iparam_jik, r2ij, r2jk, dji, djk);
        }//for-kk
        double fpair, evdwl, prefactor_ij, prefactor_ji;
        force_zeta(params+iparam_ij, r2ij, zeta_ij,&fpair, &prefactor_ij,eflag,&evdwl);

        if (i < nlocal)
        {
          fxtmp += dij[0] * fpair;
          fytmp += dij[1] * fpair;
          fztmp += dij[2] * fpair;
          if(evflag) 
            ev_tally_full(pm, i, evdwl, 0.0, -fpair, -dij[0], -dij[1], -dij[2], 
                          eng_coul, eng_vdwl, virial, eatom, vatom);
        }//if
        force_zeta(params + iparam_ij, r2ij,zeta_ji,&fpair,&prefactor_ji,eflag,&evdwl);

        if (i < nlocal)
        {
          fxtmp += dij[0] * fpair;
          fytmp += dij[1] * fpair;
          fztmp += dij[2] * fpair;
          if (evflag)
          {
            ev_tally_full(pm, i, evdwl, 0.0, -fpair, -dij[0], -dij[1], -dij[2], 
                          eng_coul, eng_vdwl, virial, eatom, vatom);
          }
        }//if
        jshort->prefactor_fwd = prefactor_ij;
        jshort->prefactor_rev = prefactor_ji;
      }//for-jj

      if (i < nlocal)
      {
        f[i][0] += fxtmp;
        f[i][1] += fytmp;
        f[i][2] += fztmp;
      }
    }//for-ii
  }//for-ist
  reg_reduce_inplace_doublev4(eng_virial, 2);
  if (_MYID == 0)
  {
    pe_put(&(pm->eng_vdwl), eng_vdwl, sizeof(double) * 8);
    pe_syn();
  }
}
*/


void pair_tersoff_compute_attractive_para(pair_tersoff_compute_param_t *pm)
{
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
  doublev4 virial_v4[2];
  virial_v4[0] = 0.0;
  virial_v4[1] = 0.0;
  double *virial = virial_v4;

  for (ist = _MYID * ISTEP; ist < inum; ist += ISTEP * 64){
    ied = ist + ISTEP;
    if (ied > inum)
      ied = inum;
    int isz = ied - ist;
    int i;
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
      pe_get(jlist_short, js, sizeof(short_neigh_t) * jnum);
      pe_syn();
      int jj;
      for (jj = 0; jj < jnum; jj ++){
        short_neigh_t *jshort = js + jj;
        int jtype = jshort->type;
        int j = jshort->idx;
        double dij[3];
        dij[0] = jshort->d[0];
        dij[1] = jshort->d[1];
        dij[2] = jshort->d[2];
        double dji[3];
        dji[0] = -dij[0];
        dji[1] = -dij[1];
        dji[2] = -dij[2];
        double r2ij = jshort->r2;
        int iparam_ij = elem2param[itype][jtype][jtype];
        if (r2ij >= params[iparam_ij].cutsq) continue;
        double prefactor_ij = jshort->prefactor_fwd;
        double prefactor_ji = jshort->prefactor_rev;

        short_neigh_t *iklist_short = js;
        int iknum = jnum;

        int kk;
        for (kk = 0; kk < iknum; kk ++){
          if (jj == kk) continue;
          short_neigh_t *kshort = iklist_short + kk;
          int ktype = kshort->type;
          int iparam_ijk = elem2param[itype][jtype][ktype];
          double dik[3];
          dik[0] = kshort->d[0];
          dik[1] = kshort->d[1];
          dik[2] = kshort->d[2];
          double r2ik = kshort->r2;
          if (r2ik >= params[iparam_ijk].cutsq) continue;
          double tfi[3], tfj[3], tfk[3];
          //lwpf_start(ATTRACTIVE);
          ters_attractive_unroll(prefactor_ij, r2ij, r2ik, dij, dik, tfi, tfj, tfk, params + iparam_ijk);
          //lwpf_stop(ATTRACTIVE);
          fxtmp += tfi[0];
          fytmp += tfi[1];
          fztmp += tfi[2];
          if (vflag_either) v_tally3rd(i, vflag_global, vflag_atom, tfj, tfk, dij, dik, virial, vatom);
        }
        //lwpf_start(JKLOAD);
        int fsj[2];
        pe_get(firstshort + j, fsj, sizeof(int) * 2);
        pe_syn();
        short_neigh_t *jklist_short = shortlist + fsj[0];
        int jknum = fsj[1] - fsj[0];
        pe_get(jklist_short, ks, sizeof(short_neigh_t) * jknum);
        pe_syn();
        //lwpf_stop(JKLOAD);
        //lwpf_start(JKLOOP);
        for (kk = 0; kk < jknum; kk ++){
          short_neigh_t *kshort = ks + kk;
          if (kshort->idx == i) continue;
          int ktype = kshort->type;
          int iparam_jik = elem2param[jtype][itype][ktype];
          int iparam_jki = elem2param[jtype][ktype][itype];
          double djk[3];
          djk[0] = kshort->d[0];
          djk[1] = kshort->d[1];
          djk[2] = kshort->d[2];
          double r2jk = kshort->r2;
          double ffi[3], ffj[3], ffk[3];
          double rfi[3], rfj[3], rfk[3];
          double tfi[3], tfj[3], tfk[3];
          double prefactor_jk = kshort->prefactor_fwd;

          if (r2jk >= params[iparam_jik].cutsq) continue;

          //lwpf_start(ATTRACTIVE_JK);
          ters_attractive_unroll_pair(prefactor_ji, prefactor_jk, r2ij, r2jk, dji, djk, tfj, tfi, tfk, params + iparam_jik, params + iparam_jki);
          //lwpf_stop(ATTRACTIVE_JK);
          fxtmp += tfi[0];
          fytmp += tfi[1];
          fztmp += tfi[2];
          if (vflag_either) v_tally3rd(i, vflag_global, vflag_atom, tfk, tfi, djk, dji, virial, vatom);
        }
        //lwpf_stop(JKLOOP);
      }
      fi[ioff][0] += fxtmp;
      fi[ioff][1] += fytmp;
      fi[ioff][2] += fztmp;
    }
    pe_put(f + ist, fi, sizeof(double) * isz * 3);
    pe_syn();
  }
  reg_reduce_inplace_doublev4(virial_v4, 2);

  if (_MYID == 0){
    pe_put(pm->virial, virial, sizeof(double) * 6);
    pe_syn();
  }
  //lwpf_stop(ALL);
}

#endif
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
