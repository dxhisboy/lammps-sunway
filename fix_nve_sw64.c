#include "sunway.h"
#include "fix_nve_sw64.h"
#ifdef MPE
extern SLAVE_FUN(fix_nve_initial_integrate_para)(fix_nve_param_t *);
void fix_nve_initial_integrate(fix_nve_param_t *pm){
  if (athread_idle() == 0)
    athread_init();
  athread_spawn(fix_nve_initial_integrate_para, pm);
  athread_join();
}
extern SLAVE_FUN(fix_nve_final_integrate_para)(fix_nve_param_t *);
void fix_nve_final_integrate(fix_nve_param_t *pm){
  if (athread_idle() == 0)
    athread_init();
  athread_spawn(fix_nve_final_integrate_para, pm);
  athread_join();
}

#endif
#ifdef CPE
#define ISTEP 512
#define NCOMP 64
void fix_nve_initial_integrate_para(fix_nve_param_t *pm){
  pe_init();
  fix_nve_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(fix_nve_param_t));
  pe_syn();
  double (*x)[3] = l_pm.x, (*v)[3] = l_pm.v, (*f)[3] = l_pm.f;
  double *rmass = l_pm.rmass;
  double dtv = l_pm.dtv, dtf = l_pm.dtf;
  int nlocal = l_pm.nlocal, groupbit = l_pm.groupbit;
  int *type = l_pm.type, *mask = l_pm.mask;
  int i, ist;
  int mi[ISTEP];
  int ti[ISTEP];
  double vi[ISTEP][3];
  double xi[ISTEP][3];
  double fi[ISTEP][3];
  if (_MYID >= NCOMP) return;
  if (rmass) {
    for (ist = _MYID * ISTEP; ist < nlocal; ist += NCOMP * ISTEP){
      int ied = ist + ISTEP;
      if (ied > nlocal) ied = nlocal;
      int isz = ied - ist;
      for(i = ist; i < ied; i ++){
        if (mask[i] & groupbit) {
          double dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
        }
      }
    }
  } else {
    double mass[l_pm.ntypes + 1];
    pe_get(l_pm.mass, mass, (l_pm.ntypes + 1) * sizeof(double));
    pe_syn();
    double dtfmt[l_pm.ntypes + 1];
    for (i = 0; i <= l_pm.ntypes; i ++)
      dtfmt[i] = dtf / mass[i];
    for (ist = _MYID * ISTEP; ist < nlocal; ist += NCOMP * ISTEP){
      int ied = ist + ISTEP;
      if (ied > nlocal) ied = nlocal;
      int isz = ied - ist;
      pe_get(x[ist], xi[0], isz * 3 * sizeof(double));
      pe_syn();
      pe_get(v[ist], vi[0], isz * 3 * sizeof(double));
      pe_syn();
      pe_get(f[ist], fi[0], isz * 3 * sizeof(double));
      pe_syn();
      pe_get(mask + ist, mi, isz * sizeof(int));
      pe_syn();
      pe_get(type + ist, ti, isz * sizeof(int));
      pe_syn();

      for(i = 0; i < isz; i ++){
        if (mi[i] & groupbit) {
          double dtfm = dtfmt[ti[i]];
          vi[i][0] += dtfm * fi[i][0];
          vi[i][1] += dtfm * fi[i][1];
          vi[i][2] += dtfm * fi[i][2];
          xi[i][0] += dtv * vi[i][0];
          xi[i][1] += dtv * vi[i][1];
          xi[i][2] += dtv * vi[i][2];
        }
      }
      pe_put(x[ist], xi[0], isz * 3 * sizeof(double));
      pe_put(v[ist], vi[0], isz * 3 * sizeof(double));
      pe_syn();
    }
  }
}
#undef NCOMP
#define NCOMP 64
void fix_nve_final_integrate_para(fix_nve_param_t *pm){
  pe_init();
  fix_nve_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(fix_nve_param_t));
  pe_syn();
  double (*x)[3] = l_pm.x, (*v)[3] = l_pm.v, (*f)[3] = l_pm.f;
  double *rmass = l_pm.rmass;
  double dtv = l_pm.dtv, dtf = l_pm.dtf;
  int nlocal = l_pm.nlocal, groupbit = l_pm.groupbit;
  int *type = l_pm.type, *mask = l_pm.mask;
  int i, ist;
  int mi[ISTEP];
  int ti[ISTEP];
  double vi[ISTEP][3];
  double fi[ISTEP][3];
  if (_MYID >= NCOMP) return;
  if (rmass) {
    for (ist = _MYID * ISTEP; ist < nlocal; ist += NCOMP * ISTEP){
      int ied = ist + ISTEP;
      if (ied > nlocal) ied = nlocal;
      int isz = ied - ist;
      for(i = ist; i < ied; i ++){
        if (mask[i] & groupbit) {
          double dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
        }
      }
    }
  } else {
    double mass[l_pm.ntypes + 1];
    pe_get(l_pm.mass, mass, (l_pm.ntypes + 1) * sizeof(double));
    pe_syn();
    double dtfmt[l_pm.ntypes + 1];
    for (i = 0; i <= l_pm.ntypes; i ++)
      dtfmt[i] = dtf / mass[i];
    for (ist = _MYID * ISTEP; ist < nlocal; ist += NCOMP * ISTEP){
      int ied = ist + ISTEP;
      if (ied > nlocal) ied = nlocal;
      int isz = ied - ist;
      pe_get(v[ist], vi[0], isz * 3 * sizeof(double));
      pe_syn();
      pe_get(f[ist], fi[0], isz * 3 * sizeof(double));
      pe_syn();
      pe_get(mask + ist, mi, isz * sizeof(int));
      pe_syn();
      pe_get(type + ist, ti, isz * sizeof(int));
      pe_syn();

      for(i = 0; i < isz; i ++){
        if (mi[i] & groupbit) {
          double dtfm = dtfmt[ti[i]];
          vi[i][0] += dtfm * fi[i][0];
          vi[i][1] += dtfm * fi[i][1];
          vi[i][2] += dtfm * fi[i][2];
        }
      }
      pe_put(v[ist], vi[0], isz * 3 * sizeof(double));
      pe_syn();
    }
  }
}

#endif
