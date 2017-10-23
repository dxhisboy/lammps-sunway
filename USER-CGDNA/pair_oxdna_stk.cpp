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
/* ----------------------------------------------------------------------
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_oxdna_stk.h"
#include "mf_oxdna.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "atom_vec_ellipsoid.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MFOxdna;

/* ---------------------------------------------------------------------- */

PairOxdnaStk::PairOxdnaStk(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairOxdnaStk::~PairOxdnaStk()
{
  if (allocated) {

    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(epsilon_st);
    memory->destroy(a_st);
    memory->destroy(cut_st_0);
    memory->destroy(cut_st_c);
    memory->destroy(cut_st_lo);
    memory->destroy(cut_st_hi);
    memory->destroy(cut_st_lc);
    memory->destroy(cut_st_hc);
    memory->destroy(b_st_lo);
    memory->destroy(b_st_hi);
    memory->destroy(shift_st);
    memory->destroy(cutsq_st_hc);

    memory->destroy(a_st4);
    memory->destroy(theta_st4_0);
    memory->destroy(dtheta_st4_ast);
    memory->destroy(b_st4);
    memory->destroy(dtheta_st4_c);

    memory->destroy(a_st5);
    memory->destroy(theta_st5_0);
    memory->destroy(dtheta_st5_ast);
    memory->destroy(b_st5);
    memory->destroy(dtheta_st5_c);

    memory->destroy(a_st6);
    memory->destroy(theta_st6_0);
    memory->destroy(dtheta_st6_ast);
    memory->destroy(b_st6);
    memory->destroy(dtheta_st6_c);

    memory->destroy(a_st1);
    memory->destroy(cosphi_st1_ast);
    memory->destroy(b_st1);
    memory->destroy(cosphi_st1_c);
    memory->destroy(a_st2);
    memory->destroy(cosphi_st2_ast);
    memory->destroy(b_st2);
    memory->destroy(cosphi_st2_c);

  }
}

/* ----------------------------------------------------------------------
   compute function for oxDNA pair interactions
   s=sugar-phosphate backbone site, b=base site, st=stacking site
------------------------------------------------------------------------- */

void PairOxdnaStk::compute(int eflag, int vflag)
{

  double delf[3],delta[3],deltb[3]; // force, torque increment;
  double evdwl,fpair,finc,tpair;
  double delr_ss[3],delr_ss_norm[3],rsq_ss,r_ss,rinv_ss;
  double delr_st[3],delr_st_norm[3],rsq_st,r_st,rinv_st;
  double theta4,t4dir[3],cost4;
  double theta5p,t5pdir[3],cost5p;
  double theta6p,t6pdir[3],cost6p;
  double cosphi1,cosphi2,cosphi1dir[3],cosphi2dir[3];

  // distances COM-backbone site, COM-stacking site
  double d_cs=-0.4, d_cst=+0.34;
  // vectors COM-backbone site, COM-stacking site in lab frame
  double ra_cs[3],ra_cst[3];
  double rb_cs[3],rb_cst[3];

  // quaternions and Cartesian unit vectors in lab frame
  double *qa,ax[3],ay[3],az[3];
  double *qb,bx[3],by[3],bz[3];

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  int *type = atom->type;

  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;

  int a,b,in,atype,btype;

  double f1,f4t4,f4t5,f4t6,f5c1,f5c2;
  double df1,df4t4,df4t5,df4t6,df5c1,df5c2;
  double tptofp;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  // loop over stacking interaction neighours using bond topology

  for (in = 0; in < nbondlist; in++) {

    a = bondlist[in][1];
    b = bondlist[in][0];

    qa=bonus[a].quat;
    MathExtra::q_to_exyz(qa,ax,ay,az);
    qb=bonus[b].quat;
    MathExtra::q_to_exyz(qb,bx,by,bz);

    // vector COM a - stacking site a
    ra_cst[0] = d_cst*ax[0];
    ra_cst[1] = d_cst*ax[1];
    ra_cst[2] = d_cst*ax[2];

    // vector COM b - stacking site b
    rb_cst[0] = d_cst*bx[0];
    rb_cst[1] = d_cst*bx[1];
    rb_cst[2] = d_cst*bx[2];

    // vector stacking site b to a
    delr_st[0] = x[a][0] + ra_cst[0] - x[b][0] - rb_cst[0];
    delr_st[1] = x[a][1] + ra_cst[1] - x[b][1] - rb_cst[1];
    delr_st[2] = x[a][2] + ra_cst[2] - x[b][2] - rb_cst[2];

    // test for directionality of vector b to a
    tptofp = MFOxdna::is_3pto5p(delr_st,bz);

    // if b to a is 5' to 3' we need to swap roles of a and b
    if (tptofp == -1) {

      std::swap(a,b);
      std::swap(ax,bx);
      std::swap(ay,by);
      std::swap(az,bz);
      std::swap(ra_cst,rb_cst);

      delr_st[0] *= -1;
      delr_st[1] *= -1;
      delr_st[2] *= -1;

    }

    atype = type[a];
    btype = type[b];

    rsq_st = delr_st[0]*delr_st[0] + delr_st[1]*delr_st[1] + delr_st[2]*delr_st[2];
    r_st = sqrt(rsq_st);
    rinv_st = 1.0/r_st;

    delr_st_norm[0] = delr_st[0] * rinv_st;
    delr_st_norm[1] = delr_st[1] * rinv_st;
    delr_st_norm[2] = delr_st[2] * rinv_st;

    // vector COM a - backbone site a
    ra_cs[0] = d_cs*ax[0];
    ra_cs[1] = d_cs*ax[1];
    ra_cs[2] = d_cs*ax[2];

    // vector COM b - backbone site b
    rb_cs[0] = d_cs*bx[0];
    rb_cs[1] = d_cs*bx[1];
    rb_cs[2] = d_cs*bx[2];

    // vector backbone site b to a
    delr_ss[0] = (x[a][0] + ra_cs[0] - x[b][0] - rb_cs[0]);
    delr_ss[1] = (x[a][1] + ra_cs[1] - x[b][1] - rb_cs[1]);
    delr_ss[2] = (x[a][2] + ra_cs[2] - x[b][2] - rb_cs[2]);

    rsq_ss = delr_ss[0]*delr_ss[0] + delr_ss[1]*delr_ss[1] + delr_ss[2]*delr_ss[2];
    r_ss = sqrt(rsq_ss);
    rinv_ss = 1.0/r_ss;

    delr_ss_norm[0] = delr_ss[0] * rinv_ss;
    delr_ss_norm[1] = delr_ss[1] * rinv_ss;
    delr_ss_norm[2] = delr_ss[2] * rinv_ss;

    f1 = F1(r_st, epsilon_st[atype][btype], a_st[atype][btype], cut_st_0[atype][btype],
        cut_st_lc[atype][btype], cut_st_hc[atype][btype], cut_st_lo[atype][btype], cut_st_hi[atype][btype],
        b_st_lo[atype][btype], b_st_hi[atype][btype], shift_st[atype][btype]);

    // early rejection criterium
    if (f1) {

    // theta4 angle and correction
    cost4 = MathExtra::dot3(az,bz);
    if (cost4 >  1.0) cost4 =  1.0;
    if (cost4 < -1.0) cost4 = -1.0;
    theta4 = acos(cost4);

    f4t4 = F4(theta4, a_st4[atype][btype], theta_st4_0[atype][btype], dtheta_st4_ast[atype][btype],
        b_st4[atype][btype], dtheta_st4_c[atype][btype]);

    // early rejection criterium
    if (f4t4) {

    // theta5 angle and correction
    cost5p  = MathExtra::dot3(delr_st_norm,az);
    if (cost5p >  1.0) cost5p =  1.0;
    if (cost5p < -1.0) cost5p = -1.0;
    theta5p = acos(cost5p);

    f4t5 = F4(theta5p, a_st5[atype][btype], theta_st5_0[atype][btype], dtheta_st5_ast[atype][btype],
        b_st5[atype][btype], dtheta_st5_c[atype][btype]);

    // early rejection criterium
    if (f4t5) {

    cost6p = MathExtra::dot3(delr_st_norm,bz);
    if (cost6p >  1.0) cost6p =  1.0;
    if (cost6p < -1.0) cost6p = -1.0;
    theta6p = acos(cost6p);

    cosphi1 = MathExtra::dot3(delr_ss_norm,ay);
    if (cosphi1 >  1.0) cosphi1 =  1.0;
    if (cosphi1 < -1.0) cosphi1 = -1.0;

    cosphi2 = MathExtra::dot3(delr_ss_norm,by);
    if (cosphi2 >  1.0) cosphi2 =  1.0;
    if (cosphi2 < -1.0) cosphi2 = -1.0;

    f4t6 = F4(theta6p, a_st6[atype][btype], theta_st6_0[atype][btype], dtheta_st6_ast[atype][btype],
        b_st6[atype][btype], dtheta_st6_c[atype][btype]);

    f5c1 = F5(-cosphi1, a_st1[atype][btype], -cosphi_st1_ast[atype][btype], b_st1[atype][btype],
        cosphi_st1_c[atype][btype]);

    f5c2 = F5(-cosphi2, a_st2[atype][btype], -cosphi_st2_ast[atype][btype], b_st2[atype][btype],
        cosphi_st2_c[atype][btype]);


    evdwl = f1 * f4t4 * f4t5 * f4t6 * f5c1 * f5c2;

    // early rejection criterium
    if (evdwl) {

    // increment energy
    if (evflag) ev_tally(a,b,nlocal,newton_bond,evdwl,0.0,0.0,0.0,0.0,0.0);

    df1 = DF1(r_st, epsilon_st[atype][btype], a_st[atype][btype], cut_st_0[atype][btype],
        cut_st_lc[atype][btype], cut_st_hc[atype][btype], cut_st_lo[atype][btype], cut_st_hi[atype][btype],
        b_st_lo[atype][btype], b_st_hi[atype][btype]);

    df4t4 = DF4(theta4, a_st4[atype][btype], theta_st4_0[atype][btype], dtheta_st4_ast[atype][btype],
        b_st4[atype][btype], dtheta_st4_c[atype][btype])/sin(theta4);

    df4t5 = DF4(theta5p, a_st5[atype][btype], theta_st5_0[atype][btype], dtheta_st5_ast[atype][btype],
        b_st5[atype][btype], dtheta_st5_c[atype][btype])/sin(theta5p);

    df4t6 = DF4(theta6p, a_st6[atype][btype], theta_st6_0[atype][btype], dtheta_st6_ast[atype][btype],
        b_st6[atype][btype], dtheta_st6_c[atype][btype])/sin(theta6p);

    df5c1 = DF5(-cosphi1, a_st1[atype][btype], -cosphi_st1_ast[atype][btype], b_st1[atype][btype],
        cosphi_st1_c[atype][btype]);

    df5c2 = DF5(-cosphi2, a_st2[atype][btype], -cosphi_st2_ast[atype][btype], b_st2[atype][btype],
        cosphi_st2_c[atype][btype]);


    // force, torque and virial contribution for forces between stacking sites

    fpair = 0.0;

    delf[0] = 0.0;
    delf[1] = 0.0;
    delf[2] = 0.0;

    delta[0] = 0.0;
    delta[1] = 0.0;
    delta[2] = 0.0;

    deltb[0] = 0.0;
    deltb[1] = 0.0;
    deltb[2] = 0.0;

    // radial force
    finc  = -df1 * f4t4 * f4t5 * f4t6 * f5c1 * f5c2;
    fpair += finc;

    delf[0] += delr_st[0] * finc;
    delf[1] += delr_st[1] * finc;
    delf[2] += delr_st[2] * finc;

    // theta5p force
    if (theta5p) {

      finc   = -f1 * f4t4 * df4t5 * f4t6 * f5c1 * f5c2 * rinv_st;
      fpair += finc;

      delf[0] += (delr_st_norm[0]*cost5p - az[0]) * finc;
      delf[1] += (delr_st_norm[1]*cost5p - az[1]) * finc;
      delf[2] += (delr_st_norm[2]*cost5p - az[2]) * finc;

    }

    // theta6p force
    if (theta6p) {

      finc   = -f1 * f4t4 * f4t5 * df4t6 * f5c1 * f5c2 * rinv_st;
      fpair += finc;

      delf[0] += (delr_st_norm[0]*cost6p - bz[0]) * finc;
      delf[1] += (delr_st_norm[1]*cost6p - bz[1]) * finc;
      delf[2] += (delr_st_norm[2]*cost6p - bz[2]) * finc;

    }

    // increment forces, torques and virial

    if (newton_bond || a < nlocal) {

      f[a][0] += delf[0];
      f[a][1] += delf[1];
      f[a][2] += delf[2];

      MathExtra::cross3(ra_cst,delf,delta);

    }
    if (newton_bond || b < nlocal) {

      f[b][0] -= delf[0];
      f[b][1] -= delf[1];
      f[b][2] -= delf[2];

      MathExtra::cross3(rb_cst,delf,deltb);

    }

    if (newton_bond || a < nlocal) {

      torque[a][0] += delta[0];
      torque[a][1] += delta[1];
      torque[a][2] += delta[2];

    }
    if (newton_bond || b < nlocal) {

      torque[b][0] -= deltb[0];
      torque[b][1] -= deltb[1];
      torque[b][2] -= deltb[2];

    }

    if (evflag) ev_tally(a,b,nlocal,newton_bond,0.0,0.0,fpair,delr_st[0],delr_st[1],delr_st[2]);

    // force, torque and virial contribution for forces between backbone sites

    fpair = 0.0;

    delf[0] = 0.0;
    delf[1] = 0.0;
    delf[2] = 0.0;

    delta[0] = 0.0;
    delta[1] = 0.0;
    delta[2] = 0.0;

    deltb[0] = 0.0;
    deltb[1] = 0.0;
    deltb[2] = 0.0;

    // cosphi1 force
    if (cosphi1) {

      finc   = -f1 * f4t4 * f4t5 * f4t6 * df5c1 * f5c2 * rinv_ss;
      fpair += finc;

      delf[0] += (delr_ss_norm[0]*cosphi1 - ay[0]) * finc;
      delf[1] += (delr_ss_norm[1]*cosphi1 - ay[1]) * finc;
      delf[2] += (delr_ss_norm[2]*cosphi1 - ay[2]) * finc;

    }

    // cosphi2 force
    if (cosphi2) {

      finc   = -f1 * f4t4 * f4t5 * f4t6 * f5c1 * df5c2 * rinv_ss;
      fpair += finc;

      delf[0] += (delr_ss_norm[0]*cosphi2 - by[0]) * finc;
      delf[1] += (delr_ss_norm[1]*cosphi2 - by[1]) * finc;
      delf[2] += (delr_ss_norm[2]*cosphi2 - by[2]) * finc;

    }

    // increment forces, torques and virial

    if (newton_bond || a < nlocal) {

      f[a][0] += delf[0];
      f[a][1] += delf[1];
      f[a][2] += delf[2];

      MathExtra::cross3(ra_cs,delf,delta);

    }
    if (newton_bond || b < nlocal) {

      f[b][0] -= delf[0];
      f[b][1] -= delf[1];
      f[b][2] -= delf[2];

      MathExtra::cross3(rb_cs,delf,deltb);

    }

    if (newton_bond || a < nlocal) {

      torque[a][0] += delta[0];
      torque[a][1] += delta[1];
      torque[a][2] += delta[2];

    }
    if (newton_bond || b < nlocal) {

      torque[b][0] -= deltb[0];
      torque[b][1] -= deltb[1];
      torque[b][2] -= deltb[2];

    }

    if (evflag) ev_tally(a,b,nlocal,newton_bond,0.0,0.0,fpair,delr_ss[0],delr_ss[1],delr_ss[2]);


    // pure torques not expressible as r x f

    delta[0] = 0.0;
    delta[1] = 0.0;
    delta[2] = 0.0;
    deltb[0] = 0.0;
    deltb[1] = 0.0;
    deltb[2] = 0.0;

    // theta4 torque
    if (theta4) {

      tpair = -f1 * df4t4 * f4t5 * f4t6 * f5c1 * f5c2;
      MathExtra::cross3(bz,az,t4dir);

      delta[0] += t4dir[0]*tpair;
      delta[1] += t4dir[1]*tpair;
      delta[2] += t4dir[2]*tpair;

      deltb[0] += t4dir[0]*tpair;
      deltb[1] += t4dir[1]*tpair;
      deltb[2] += t4dir[2]*tpair;

    }

    // theta5p torque
    if (theta5p) {

      tpair = -f1 * f4t4 * df4t5 * f4t6 * f5c1 * f5c2;
      MathExtra::cross3(delr_st_norm,az,t5pdir);

      delta[0] += t5pdir[0] * tpair;
      delta[1] += t5pdir[1] * tpair;
      delta[2] += t5pdir[2] * tpair;

    }

    // theta6p torque
    if (theta6p) {

      tpair = -f1 * f4t4 * f4t5 * df4t6 * f5c1 * f5c2;
      MathExtra::cross3(delr_st_norm,bz,t6pdir);

      deltb[0] -= t6pdir[0] * tpair;
      deltb[1] -= t6pdir[1] * tpair;
      deltb[2] -= t6pdir[2] * tpair;

    }

    // cosphi1 torque
    if (cosphi1) {

      tpair   = -f1 * f4t4 * f4t5 * f4t6 * df5c1 * f5c2;
      MathExtra::cross3(delr_ss_norm,ay,cosphi1dir);

      delta[0] += cosphi1dir[0] * tpair;
      delta[1] += cosphi1dir[1] * tpair;
      delta[2] += cosphi1dir[2] * tpair;

    }

    // cosphi2 torque
    if (cosphi2) {

      tpair   = -f1 * f4t4 * f4t5 * f4t6 * f5c1 * df5c2;
      MathExtra::cross3(delr_ss_norm,by,cosphi2dir);

      deltb[0] -= cosphi2dir[0] * tpair;
      deltb[1] -= cosphi2dir[1] * tpair;
      deltb[2] -= cosphi2dir[2] * tpair;

    }

    // increment torques
    if (newton_bond || a < nlocal) {

      torque[a][0] += delta[0];
      torque[a][1] += delta[1];
      torque[a][2] += delta[2];

    }
    if (newton_bond || b < nlocal) {

      torque[b][0] -= deltb[0];
      torque[b][1] -= deltb[1];
      torque[b][2] -= deltb[2];

    }

    }
    }
    }
    }
    // end early rejection criteria

  }
  // end stacking interaction

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairOxdnaStk::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(epsilon_st,n+1,n+1,"pair:epsilon_st");
  memory->create(a_st,n+1,n+1,"pair:a_st");
  memory->create(cut_st_0,n+1,n+1,"pair:cut_st_0");
  memory->create(cut_st_c,n+1,n+1,"pair:cut_st_c");
  memory->create(cut_st_lo,n+1,n+1,"pair:cut_st_lo");
  memory->create(cut_st_hi,n+1,n+1,"pair:cut_st_hi");
  memory->create(cut_st_lc,n+1,n+1,"pair:cut_st_lc");
  memory->create(cut_st_hc,n+1,n+1,"pair:cut_st_hc");
  memory->create(b_st_lo,n+1,n+1,"pair:b_st_lo");
  memory->create(b_st_hi,n+1,n+1,"pair:b_st_hi");
  memory->create(shift_st,n+1,n+1,"pair:shift_st");
  memory->create(cutsq_st_hc,n+1,n+1,"pair:cutsq_st_hc");

  memory->create(a_st4,n+1,n+1,"pair:a_st4");
  memory->create(theta_st4_0,n+1,n+1,"pair:theta_st4_0");
  memory->create(dtheta_st4_ast,n+1,n+1,"pair:dtheta_st4_ast");
  memory->create(b_st4,n+1,n+1,"pair:b_st4");
  memory->create(dtheta_st4_c,n+1,n+1,"pair:dtheta_st4_c");

  memory->create(a_st5,n+1,n+1,"pair:a_st5");
  memory->create(theta_st5_0,n+1,n+1,"pair:theta_st5_0");
  memory->create(dtheta_st5_ast,n+1,n+1,"pair:dtheta_st5_ast");
  memory->create(b_st5,n+1,n+1,"pair:b_st5");
  memory->create(dtheta_st5_c,n+1,n+1,"pair:dtheta_st5_c");

  memory->create(a_st6,n+1,n+1,"pair:a_st6");
  memory->create(theta_st6_0,n+1,n+1,"pair:theta_st6_0");
  memory->create(dtheta_st6_ast,n+1,n+1,"pair:dtheta_st6_ast");
  memory->create(b_st6,n+1,n+1,"pair:b_st6");
  memory->create(dtheta_st6_c,n+1,n+1,"pair:dtheta_st6_c");

  memory->create(a_st1,n+1,n+1,"pair:a_st1");
  memory->create(cosphi_st1_ast,n+1,n+1,"pair:cosphi_st1_ast");
  memory->create(b_st1,n+1,n+1,"pair:b_st1");
  memory->create(cosphi_st1_c,n+1,n+1,"pair:cosphi_st1_c");
  memory->create(a_st2,n+1,n+1,"pair:a_st2");
  memory->create(cosphi_st2_ast,n+1,n+1,"pair:cosphi_st2_ast");
  memory->create(b_st2,n+1,n+1,"pair:b_st2");
  memory->create(cosphi_st2_c,n+1,n+1,"pair:cosphi_st2_c");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairOxdnaStk::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

}

/* ----------------------------------------------------------------------
   return temperature dependent oxDNA stacking strength
------------------------------------------------------------------------- */

double PairOxdnaStk::stacking_strength(double T)
{
  double eps;

  eps = 1.3448 + 2.6568 * T;

  return eps;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairOxdnaStk::coeff(int narg, char **arg)
{
  int count;

  if (narg != 21) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/stk");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  // stacking interaction
  count = 0;

  double T, epsilon_st_one, a_st_one, b_st_lo_one, b_st_hi_one;
  double cut_st_0_one, cut_st_c_one, cut_st_lo_one, cut_st_hi_one;
  double cut_st_lc_one, cut_st_hc_one, tmp, shift_st_one;

  double a_st4_one, theta_st4_0_one, dtheta_st4_ast_one;
  double b_st4_one, dtheta_st4_c_one;

  double a_st5_one, theta_st5_0_one, dtheta_st5_ast_one;
  double b_st5_one, dtheta_st5_c_one;

  double a_st6_one, theta_st6_0_one, dtheta_st6_ast_one;
  double b_st6_one, dtheta_st6_c_one;

  double a_st1_one, cosphi_st1_ast_one, b_st1_one, cosphi_st1_c_one;
  double a_st2_one, cosphi_st2_ast_one, b_st2_one, cosphi_st2_c_one;

  T = force->numeric(FLERR,arg[2]);
  epsilon_st_one = stacking_strength(T);

  a_st_one = force->numeric(FLERR,arg[3]);
  cut_st_0_one = force->numeric(FLERR,arg[4]);
  cut_st_c_one = force->numeric(FLERR,arg[5]);
  cut_st_lo_one = force->numeric(FLERR,arg[6]);
  cut_st_hi_one = force->numeric(FLERR,arg[7]);

  a_st4_one = force->numeric(FLERR,arg[8]);
  theta_st4_0_one = force->numeric(FLERR,arg[9]);
  dtheta_st4_ast_one = force->numeric(FLERR,arg[10]);
  a_st5_one = force->numeric(FLERR,arg[11]);
  theta_st5_0_one = force->numeric(FLERR,arg[12]);
  dtheta_st5_ast_one = force->numeric(FLERR,arg[13]);
  a_st6_one = force->numeric(FLERR,arg[14]);
  theta_st6_0_one = force->numeric(FLERR,arg[15]);
  dtheta_st6_ast_one = force->numeric(FLERR,arg[16]);
  a_st1_one = force->numeric(FLERR,arg[17]);
  cosphi_st1_ast_one = force->numeric(FLERR,arg[18]);
  a_st2_one = force->numeric(FLERR,arg[19]);
  cosphi_st2_ast_one = force->numeric(FLERR,arg[20]);

  b_st_lo_one = 2*a_st_one*exp(-a_st_one*(cut_st_lo_one-cut_st_0_one))*
        2*a_st_one*exp(-a_st_one*(cut_st_lo_one-cut_st_0_one))*
        (1-exp(-a_st_one*(cut_st_lo_one-cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_lo_one-cut_st_0_one)))/
        (4*((1-exp(-a_st_one*(cut_st_lo_one -cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_lo_one-cut_st_0_one)))-
        (1-exp(-a_st_one*(cut_st_c_one -cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_c_one-cut_st_0_one)))));

  cut_st_lc_one = cut_st_lo_one - a_st_one*exp(-a_st_one*(cut_st_lo_one-cut_st_0_one))*
        (1-exp(-a_st_one*(cut_st_lo_one-cut_st_0_one)))/b_st_lo_one;

  b_st_hi_one = 2*a_st_one*exp(-a_st_one*(cut_st_hi_one-cut_st_0_one))*
        2*a_st_one*exp(-a_st_one*(cut_st_hi_one-cut_st_0_one))*
        (1-exp(-a_st_one*(cut_st_hi_one-cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_hi_one-cut_st_0_one)))/
        (4*((1-exp(-a_st_one*(cut_st_hi_one -cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_hi_one-cut_st_0_one)))-
        (1-exp(-a_st_one*(cut_st_c_one -cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_c_one-cut_st_0_one)))));

  cut_st_hc_one = cut_st_hi_one - a_st_one*exp(-a_st_one*(cut_st_hi_one-cut_st_0_one))*
        (1-exp(-a_st_one*(cut_st_hi_one-cut_st_0_one)))/b_st_hi_one;

  tmp = 1 - exp(-(cut_st_c_one-cut_st_0_one) * a_st_one);
  shift_st_one = epsilon_st_one * tmp * tmp;

  b_st4_one = a_st4_one*a_st4_one*dtheta_st4_ast_one*dtheta_st4_ast_one/(1-a_st4_one*dtheta_st4_ast_one*dtheta_st4_ast_one);
  dtheta_st4_c_one = 1/(a_st4_one*dtheta_st4_ast_one);

  b_st5_one = a_st5_one*a_st5_one*dtheta_st5_ast_one*dtheta_st5_ast_one/(1-a_st5_one*dtheta_st5_ast_one*dtheta_st5_ast_one);
  dtheta_st5_c_one = 1/(a_st5_one*dtheta_st5_ast_one);

  b_st6_one = a_st6_one*a_st6_one*dtheta_st6_ast_one*dtheta_st6_ast_one/(1-a_st6_one*dtheta_st6_ast_one*dtheta_st6_ast_one);
  dtheta_st6_c_one = 1/(a_st6_one*dtheta_st6_ast_one);

  b_st1_one = a_st1_one*a_st1_one*cosphi_st1_ast_one*cosphi_st1_ast_one/(1-a_st1_one*cosphi_st1_ast_one*cosphi_st1_ast_one);
  cosphi_st1_c_one=1/(a_st1_one*cosphi_st1_ast_one);

  b_st2_one = a_st2_one*a_st2_one*cosphi_st2_ast_one*cosphi_st2_ast_one/(1-a_st2_one*cosphi_st2_ast_one*cosphi_st2_ast_one);
  cosphi_st2_c_one=1/(a_st2_one*cosphi_st2_ast_one);


  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {

      epsilon_st[i][j] = epsilon_st_one;
      a_st[i][j] = a_st_one;
      cut_st_0[i][j] = cut_st_0_one;
      cut_st_c[i][j] = cut_st_c_one;
      cut_st_lo[i][j] = cut_st_lo_one;
      cut_st_hi[i][j] = cut_st_hi_one;
      cut_st_lc[i][j] = cut_st_lc_one;
      cut_st_hc[i][j] = cut_st_hc_one;
      b_st_lo[i][j] = b_st_lo_one;
      b_st_hi[i][j] = b_st_hi_one;
      shift_st[i][j] = shift_st_one;

      a_st4[i][j] = a_st4_one;
      theta_st4_0[i][j] = theta_st4_0_one;
      dtheta_st4_ast[i][j] = dtheta_st4_ast_one;
      b_st4[i][j] = b_st4_one;
      dtheta_st4_c[i][j] = dtheta_st4_c_one;

      a_st5[i][j] = a_st5_one;
      theta_st5_0[i][j] = theta_st5_0_one;
      dtheta_st5_ast[i][j] = dtheta_st5_ast_one;
      b_st5[i][j] = b_st5_one;
      dtheta_st5_c[i][j] = dtheta_st5_c_one;

      a_st6[i][j] = a_st6_one;
      theta_st6_0[i][j] = theta_st6_0_one;
      dtheta_st6_ast[i][j] = dtheta_st6_ast_one;
      b_st6[i][j] = b_st6_one;
      dtheta_st6_c[i][j] = dtheta_st6_c_one;

      a_st1[i][j] = a_st1_one;
      cosphi_st1_ast[i][j] = cosphi_st1_ast_one;
      b_st1[i][j] = b_st1_one;
      cosphi_st1_c[i][j] = cosphi_st1_c_one;

      a_st2[i][j] = a_st2_one;
      cosphi_st2_ast[i][j] = cosphi_st2_ast_one;
      b_st2[i][j] = b_st2_one;
      cosphi_st2_c[i][j] = cosphi_st2_c_one;

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/stk");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairOxdnaStk::init_style()
{
  int irequest;

  // request regular neighbor lists

  irequest = neighbor->request(this,instance_me);

}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use regular
------------------------------------------------------------------------- */

void PairOxdnaStk::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  if (id  > 0) error->all(FLERR,"Respa not supported");

}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairOxdnaStk::init_one(int i, int j)
{

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Coefficient mixing not defined in oxDNA");
  }
  if (offset_flag) {
    error->all(FLERR,"Offset not supported in oxDNA");
  }

  epsilon_st[j][i] = epsilon_st[i][j];
  a_st[j][i] = a_st[i][j];
  b_st_lo[j][i] = b_st_lo[i][j];
  b_st_hi[j][i] = b_st_hi[i][j];
  cut_st_0[j][i] = cut_st_0[i][j];
  cut_st_c[j][i] = cut_st_c[i][j];
  cut_st_lo[j][i] = cut_st_lo[i][j];
  cut_st_hi[j][i] = cut_st_hi[i][j];
  cut_st_lc[j][i] = cut_st_lc[i][j];
  cut_st_hc[j][i] = cut_st_hc[i][j];
  shift_st[j][i] = shift_st[i][j];

  a_st4[j][i] = a_st4[i][j];
  theta_st4_0[j][i] = theta_st4_0[i][j];
  dtheta_st4_ast[j][i] = dtheta_st4_ast[i][j];
  b_st4[j][i] = b_st4[i][j];
  dtheta_st4_c[j][i] = dtheta_st4_c[i][j];

  a_st5[j][i] = a_st5[i][j];
  theta_st5_0[j][i] = theta_st5_0[i][j];
  dtheta_st5_ast[j][i] = dtheta_st5_ast[i][j];
  b_st5[j][i] = b_st5[i][j];
  dtheta_st5_c[j][i] = dtheta_st5_c[i][j];

  a_st6[j][i] = a_st6[i][j];
  theta_st6_0[j][i] = theta_st6_0[i][j];
  dtheta_st6_ast[j][i] = dtheta_st6_ast[i][j];
  b_st6[j][i] = b_st6[i][j];
  dtheta_st6_c[j][i] = dtheta_st6_c[i][j];

  a_st1[j][i] = a_st1[i][j];
  cosphi_st1_ast[j][i] = cosphi_st1_ast[i][j];
  b_st1[j][i] = b_st1[i][j];
  cosphi_st1_c[j][i] = cosphi_st1_c[i][j];

  a_st2[j][i] = a_st2[i][j];
  cosphi_st2_ast[j][i] = cosphi_st2_ast[i][j];
  b_st2[j][i] = b_st2[i][j];
  cosphi_st2_c[j][i] = cosphi_st2_c[i][j];

  cutsq_st_hc[i][j] = cut_st_hc[i][j]*cut_st_hc[i][j];
  cutsq_st_hc[j][i] = cutsq_st_hc[i][j];

  // set the master list distance cutoff
  return cut_st_hc[i][j];

}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdnaStk::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {

        fwrite(&epsilon_st[i][j],sizeof(double),1,fp);
        fwrite(&a_st[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_0[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_c[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_lo[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_hi[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_lc[i][j],sizeof(double),1,fp);
        fwrite(&cut_st_hc[i][j],sizeof(double),1,fp);
        fwrite(&b_st_lo[i][j],sizeof(double),1,fp);
        fwrite(&b_st_hi[i][j],sizeof(double),1,fp);
        fwrite(&shift_st[i][j],sizeof(double),1,fp);

        fwrite(&a_st4[i][j],sizeof(double),1,fp);
        fwrite(&theta_st4_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st4_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_st4[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st4_c[i][j],sizeof(double),1,fp);

        fwrite(&a_st5[i][j],sizeof(double),1,fp);
        fwrite(&theta_st5_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st5_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_st5[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st5_c[i][j],sizeof(double),1,fp);

        fwrite(&a_st6[i][j],sizeof(double),1,fp);
        fwrite(&theta_st6_0[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st6_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_st6[i][j],sizeof(double),1,fp);
        fwrite(&dtheta_st6_c[i][j],sizeof(double),1,fp);

        fwrite(&a_st1[i][j],sizeof(double),1,fp);
        fwrite(&cosphi_st1_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_st1[i][j],sizeof(double),1,fp);
        fwrite(&cosphi_st1_c[i][j],sizeof(double),1,fp);
        fwrite(&a_st2[i][j],sizeof(double),1,fp);
        fwrite(&cosphi_st2_ast[i][j],sizeof(double),1,fp);
        fwrite(&b_st2[i][j],sizeof(double),1,fp);
        fwrite(&cosphi_st2_c[i][j],sizeof(double),1,fp);

    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdnaStk::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {

          fread(&epsilon_st[i][j],sizeof(double),1,fp);
          fread(&a_st[i][j],sizeof(double),1,fp);
          fread(&cut_st_0[i][j],sizeof(double),1,fp);
          fread(&cut_st_c[i][j],sizeof(double),1,fp);
          fread(&cut_st_lo[i][j],sizeof(double),1,fp);
          fread(&cut_st_hi[i][j],sizeof(double),1,fp);
          fread(&cut_st_lc[i][j],sizeof(double),1,fp);
          fread(&cut_st_hc[i][j],sizeof(double),1,fp);
          fread(&b_st_lo[i][j],sizeof(double),1,fp);
          fread(&b_st_hi[i][j],sizeof(double),1,fp);
          fread(&shift_st[i][j],sizeof(double),1,fp);

          fread(&a_st4[i][j],sizeof(double),1,fp);
          fread(&theta_st4_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_st4_ast[i][j],sizeof(double),1,fp);
          fread(&b_st4[i][j],sizeof(double),1,fp);
          fread(&dtheta_st4_c[i][j],sizeof(double),1,fp);

          fread(&a_st5[i][j],sizeof(double),1,fp);
          fread(&theta_st5_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_st5_ast[i][j],sizeof(double),1,fp);
          fread(&b_st5[i][j],sizeof(double),1,fp);
          fread(&dtheta_st5_c[i][j],sizeof(double),1,fp);

          fread(&a_st6[i][j],sizeof(double),1,fp);
          fread(&theta_st6_0[i][j],sizeof(double),1,fp);
          fread(&dtheta_st6_ast[i][j],sizeof(double),1,fp);
          fread(&b_st6[i][j],sizeof(double),1,fp);
          fread(&dtheta_st6_c[i][j],sizeof(double),1,fp);

          fread(&a_st1[i][j],sizeof(double),1,fp);
          fread(&cosphi_st1_ast[i][j],sizeof(double),1,fp);
          fread(&b_st1[i][j],sizeof(double),1,fp);
          fread(&cosphi_st1_c[i][j],sizeof(double),1,fp);
          fread(&a_st2[i][j],sizeof(double),1,fp);
          fread(&cosphi_st2_ast[i][j],sizeof(double),1,fp);
          fread(&b_st2[i][j],sizeof(double),1,fp);
          fread(&cosphi_st2_c[i][j],sizeof(double),1,fp);

        }

        MPI_Bcast(&epsilon_st[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a_st[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_lo[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_hi[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_lc[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_st_hc[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st_lo[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st_hi[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&shift_st[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_st4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_st4_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st4_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st4_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_st5[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_st5_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st5_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st5[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st5_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_st6[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta_st6_0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st6_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st6[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&dtheta_st6_c[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&a_st1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cosphi_st1_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cosphi_st1_c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a_st2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cosphi_st2_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_st2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cosphi_st2_c[i][j],1,MPI_DOUBLE,0,world);

      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdnaStk::write_restart_settings(FILE *fp)
{
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdnaStk::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairOxdnaStk::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d\
         %g %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g\
         %g %g %g %g\
         \n",i,
        epsilon_st[i][i],a_st[i][i],cut_st_0[i][i],cut_st_c[i][i],cut_st_lo[i][i],cut_st_hi[i][i],
        cut_st_lc[i][i],cut_st_hc[i][i],b_st_lo[i][i],b_st_hi[i][i],shift_st[i][i],
        a_st4[i][i],theta_st4_0[i][i],dtheta_st4_ast[i][i],b_st4[i][i],dtheta_st4_c[i][i],
        a_st5[i][i],theta_st5_0[i][i],dtheta_st5_ast[i][i],b_st5[i][i],dtheta_st5_c[i][i],
        a_st6[i][i],theta_st6_0[i][i],dtheta_st6_ast[i][i],b_st6[i][i],dtheta_st6_c[i][i],
        a_st1[i][i],cosphi_st1_ast[i][i],b_st1[i][i], cosphi_st1_c[i][i],
        a_st2[i][i],cosphi_st2_ast[i][i],b_st2[i][i], cosphi_st2_c[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairOxdnaStk::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d\
         %g %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g %g\
         %g %g %g %g\
         %g %g %g %g\
         \n",i,j,
        epsilon_st[i][j],a_st[i][j],cut_st_0[i][j],cut_st_c[i][j],cut_st_lo[i][j],cut_st_hi[i][j],
        cut_st_lc[i][j],cut_st_hc[i][j],b_st_lo[i][j],b_st_hi[i][j],shift_st[i][j],
        a_st4[i][j],theta_st4_0[i][j],dtheta_st4_ast[i][j],b_st4[i][j],dtheta_st4_c[i][j],
        a_st5[i][j],theta_st5_0[i][j],dtheta_st5_ast[i][j],b_st5[i][j],dtheta_st5_c[i][j],
        a_st6[i][j],theta_st6_0[i][j],dtheta_st6_ast[i][j],b_st6[i][j],dtheta_st6_c[i][j],
        a_st1[i][j],cosphi_st1_ast[i][j],b_st1[i][j],cosphi_st1_c[i][j],
        a_st2[i][j],cosphi_st2_ast[i][j],b_st2[i][j],cosphi_st2_c[i][j]);
}

/* ---------------------------------------------------------------------- */

void *PairOxdnaStk::extract(const char *str, int &dim)
{
  dim = 2;

  if (strcmp(str,"epsilon_st") == 0) return (void *) epsilon_st;
  if (strcmp(str,"a_st") == 0) return (void *) a_st;
  if (strcmp(str,"cut_st_0") == 0) return (void *) cut_st_0;
  if (strcmp(str,"cut_st_c") == 0) return (void *) cut_st_c;
  if (strcmp(str,"cut_st_lo") == 0) return (void *) cut_st_lo;
  if (strcmp(str,"cut_st_hi") == 0) return (void *) cut_st_hi;
  if (strcmp(str,"cut_st_lc") == 0) return (void *) cut_st_lc;
  if (strcmp(str,"cut_st_hc") == 0) return (void *) cut_st_hc;
  if (strcmp(str,"b_st_lo") == 0) return (void *) b_st_lo;
  if (strcmp(str,"b_st_hi") == 0) return (void *) b_st_hi;
  if (strcmp(str,"shift_st") == 0) return (void *) shift_st;

  if (strcmp(str,"a_st4") == 0) return (void *) a_st4;
  if (strcmp(str,"theta_st4_0") == 0) return (void *) theta_st4_0;
  if (strcmp(str,"dtheta_st4_ast") == 0) return (void *) dtheta_st4_ast;
  if (strcmp(str,"b_st4") == 0) return (void *) b_st4;
  if (strcmp(str,"dtheta_st4_c") == 0) return (void *) dtheta_st4_c;

  if (strcmp(str,"a_st5") == 0) return (void *) a_st5;
  if (strcmp(str,"theta_st5_0") == 0) return (void *) theta_st5_0;
  if (strcmp(str,"dtheta_st5_ast") == 0) return (void *) dtheta_st5_ast;
  if (strcmp(str,"b_st5") == 0) return (void *) b_st5;
  if (strcmp(str,"dtheta_st5_c") == 0) return (void *) dtheta_st5_c;

  if (strcmp(str,"a_st6") == 0) return (void *) a_st6;
  if (strcmp(str,"theta_st6_0") == 0) return (void *) theta_st6_0;
  if (strcmp(str,"dtheta_st6_ast") == 0) return (void *) dtheta_st6_ast;
  if (strcmp(str,"b_st6") == 0) return (void *) b_st6;
  if (strcmp(str,"dtheta_st6_c") == 0) return (void *) dtheta_st6_c;

  if (strcmp(str,"a_st1") == 0) return (void *) a_st1;
  if (strcmp(str,"cosphi_st1_ast") == 0) return (void *) cosphi_st1_ast;
  if (strcmp(str,"b_st1") == 0) return (void *) b_st1;
  if (strcmp(str,"cosphi_st1_c") == 0) return (void *) cosphi_st1_c;

  if (strcmp(str,"a_st2") == 0) return (void *) a_st2;
  if (strcmp(str,"cosphi_st2_ast") == 0) return (void *) cosphi_st2_ast;
  if (strcmp(str,"b_st2") == 0) return (void *) b_st2;
  if (strcmp(str,"cosphi_st2_c") == 0) return (void *) cosphi_st2_c;

  return NULL;
}
