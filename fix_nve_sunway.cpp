/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include "fix_nve_sunway.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "sunway.h"
using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVESunway::FixNVESunway(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"nve/sphere") != 0 && narg < 3)
    error->all(FLERR,"Illegal fix nve command");

  dynamic_group_allow = 1;
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixNVESunway::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVESunway::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  if (strstr(update->integrate_style,"respa"))
    step_respa = ((Respa *) update->integrate)->step;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */
extern "C"{
  void fix_nve_initial_integrate_sunway_compute(fix_nve_param_t *pm);
  void fix_nve_final_integrate_sunway_compute(fix_nve_param_t *pm);
}



void FixNVESunway::initial_integrate(int vflag)
{
  double dtfm;

  // update v and x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;


  fix_nve_param_t pm;


  pm.x = atom->x;
  pm.v = atom->v;
  pm.f = atom->f;
  pm.rmass = atom->rmass;
  pm.mass = atom->mass;
  pm.type = atom->type;
  pm.mask = atom->mask;
  pm.nlocal = nlocal;
  pm.groupbit = groupbit;
  pm.dtv = dtv;
  pm.dtf = dtf;

  double st = -MPI_Wtime();
  fix_nve_initial_integrate_sunway_compute(&pm);
  st += MPI_Wtime();
//  printf("the time %.10f\n",st);
//
//  int max = 0,min = 0xffff;
//  for(int i = 0;i < nlocal;i++)
//  {
//    if(max < type[i])
//      max = type[i];
//    if(min > type[i])
//      min = type[i];
//  }
//
//  printf("ntype = %d\n",ntype);
//  printf("max = %d ,min = %d\n",max,min);
//
//  st = -MPI_Wtime();
//  if (rmass) {
//
//    for (int i = 0; i < nlocal; i++)
//      if (mask[i] & groupbit) {
//        dtfm = dtf / rmass[i];
//
//        v[i][0] += dtfm * f[i][0];
//        v[i][1] += dtfm * f[i][1];
//        v[i][2] += dtfm * f[i][2];
//        x[i][0] += dtv * v[i][0];
//        x[i][1] += dtv * v[i][1];
//        x[i][2] += dtv * v[i][2];
//      }
//
//  } else {
//    for (int i = 0; i < nlocal; i++)
//      if (mask[i] & groupbit) {
//
//        dtfm = dtf / mass[type[i]];
//        v[i][0] += dtfm * f[i][0];
//        v[i][1] += dtfm * f[i][1];
//        v[i][2] += dtfm * f[i][2];
//        x[i][0] += dtv * v[i][0];
//        x[i][1] += dtv * v[i][1];
//        x[i][2] += dtv * v[i][2];
//      }
//  }
//  
//  st += MPI_Wtime();
//  printf("the original time %.10f\n",st);
}

/* ---------------------------------------------------------------------- */

void FixNVESunway::final_integrate()
{
  double dtfm;

  // update v of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  fix_nve_param_t pm;


  pm.x = atom->x;
  pm.v = atom->v;
  pm.f = atom->f;
  pm.rmass = atom->rmass;
  pm.mass = atom->mass;
  pm.type = atom->type;
  pm.mask = atom->mask;
  pm.nlocal = nlocal;
  pm.groupbit = groupbit;
  pm.dtv = dtv;
  pm.dtf = dtf;

  fix_nve_final_integrate_sunway_compute(&pm);

//
//  if (rmass) {
//    for (int i = 0; i < nlocal; i++)
//      if (mask[i] & groupbit) {
//        dtfm = dtf / rmass[i];
//
//        v[i][0] += dtfm * f[i][0];
//        v[i][1] += dtfm * f[i][1];
//        v[i][2] += dtfm * f[i][2];
//      }
//
//  } else {
//    for (int i = 0; i < nlocal; i++)
//      if (mask[i] & groupbit) {
//        dtfm = dtf / mass[type[i]];
//        v[i][0] += dtfm * f[i][0];
//        v[i][1] += dtfm * f[i][1];
//        v[i][2] += dtfm * f[i][2];
//      }
//  }
}

/* ---------------------------------------------------------------------- */

void FixNVESunway::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVESunway::final_integrate_respa(int ilevel, int iloop)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVESunway::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
