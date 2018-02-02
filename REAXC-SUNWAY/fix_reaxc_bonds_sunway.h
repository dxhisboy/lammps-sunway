/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(reax/c/bonds/sunway,FixReaxCBondsSunway)

#else

#ifndef LMP_FIX_REAXC_BONDS_SUNWAY_H
#define LMP_FIX_REAXC_BONDS_SUNWAY_H

#include <stdio.h>
#include "fix.h"
#include "pointers.h"

namespace LAMMPS_NS {

class FixReaxCBondsSunway : public Fix {
 public:
  FixReaxCBondsSunway(class LAMMPS *, int, char **);
  virtual ~FixReaxCBondsSunway();
  int setmask();
  virtual void init();
  void setup(int);
  void end_of_step();

 protected:
  int me, nprocs, nmax, ntypes, maxsize;
  int *numneigh;
  tagint **neighid;
  double **abo;
  FILE *fp;

  void allocate();
  void destroy();
  virtual void Output_ReaxC_Bonds(bigint, FILE *);
  void FindBond(struct _reax_list*, int &);
  void PassBuffer(double *, int &);
  void RecvBuffer(double *, int, int, int, int);
  int nint(const double &);
  virtual double memory_usage();

  bigint nvalid, nextvalid();
  struct _reax_list *lists;
  class PairReaxCSunway *reaxc;
  class NeighList *list;
};
}

#endif
#endif
